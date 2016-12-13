/*
taxator-tk predicts the taxon for DNA sequences based on sequence alignment.

Copyright (C) 2010 Johannes Dröge

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef taxonpredictionmodelsequence_hh_
#define taxonpredictionmodelsequence_hh_

#include <seqan/align.h>
#include <seqan/graph_msa.h>
#include <seqan/file.h>
#include <seqan/basic.h>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <boost/format.hpp>
#include <assert.h>
#include <limits>
#include <set>
#include <ostream>
#include <cmath>
#include "taxonpredictionmodel.hh"
#include "sequencestorage.hh"
#include "profiling.hh"

// alignment object
template <typename StringType>
struct alignment{
                int score; 
                int matches;
                int mmatches;
                int gaps;
                seqan::Align<StringType, seqan::ArrayGaps> alignment;
                };

// helper class
class BandFactor {
public:
    BandFactor(TaxonomyInterface& taxinter, uint reserve = 0) :
        bandfactor_(-1),
        taxinter_(taxinter) {
        if(reserve) data_.reserve(reserve);
    }

    void addSequence(const int score, const TaxonNode* node) {
        data_.push_back(boost::make_tuple(score, node));
    }

    float getFactor() {
        if(bandfactor_ < 0) {
            sort();
            setBandFactor();
        }
        return sqrt(bandfactor_);
    }

private:
    void setBandFactor(const float min_bandfactor = 1., const float max_bandfactor = std::numeric_limits<  int >::max()) { //data_ must be sorted TODO: optimize
        float bandfactor = min_bandfactor;
        int score;
        const TaxonNode* anchor;
        const TaxonNode* node;
        std::map< small_unsigned_int, int > worstscore_per_rank;
        boost::tie(score, anchor) = data_[0];
        small_unsigned_int last_rank = anchor->data->root_pathlength;
        worstscore_per_rank[ last_rank ] = score; //TODO: change name

        for (uint i = 1; i < data_.size(); ++i) {
            boost::tie(score, node) = data_[i];
            node = taxinter_.getLCA(node, anchor);
            small_unsigned_int rank = node->data->root_pathlength;
            if (rank == last_rank) {
            } else if (rank < last_rank) {
                worstscore_per_rank[rank] = score;
                last_rank = rank;
            } else { //disorder
                int refscore;
                small_unsigned_int r = rank - 1;
                do {
                    std::map< small_unsigned_int, int >::iterator it = worstscore_per_rank.find(r);
                    if (it != worstscore_per_rank.end()) {
                        refscore = it->second;
                        if (refscore) bandfactor = std::max(bandfactor, float(score)/float(refscore));
                    }
                } while (r--);
            }
        }
        bandfactor_ = std::min(bandfactor, max_bandfactor);
    }

    void sort() { //sort by increasing order
        data_type_::iterator start_it = data_.begin();
        std::sort(++start_it, data_.end(), comparator_);
    }

    compareTupleFirstLT< boost::tuple< int, const TaxonNode* >, 0 > comparator_;
    typedef std::vector< boost::tuple< int, const TaxonNode* > > data_type_;
    float bandfactor_;
    data_type_ data_;
    TaxonomyInterface taxinter_;

};

// TODO: make timers thread-safe
template< typename ContainerT, typename QStorType, typename DBStorType ,typename StringType>
class RPAPredictionModel : public TaxonPredictionModel< ContainerT > {
public:
    RPAPredictionModel(const Taxonomy* tax, QStorType& q_storage, const DBStorType& db_storage, float exclude_factor ,float reeval_bandwidth = .1) :
        TaxonPredictionModel< ContainerT >(tax),
        query_sequences_(q_storage),
        db_sequences_(db_storage),
        exclude_alignments_factor_(exclude_factor),
        reeval_bandwidth_factor_(1. - reeval_bandwidth),
        measure_sequence_retrieval_("sequence retrieval using index"),
        measure_pass_0_alignment_("best reference re-evaluation alignments (pass 0)"),
        measure_pass_1_alignment_("best reference anchor alignments (pass 1)"),
        measure_pass_2_alignment_("distant anchor alignments (pass 2)")
    {};

    void predict(ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink) {
        this->initPredictionRecord(recordset, prec);  // set query name and length
        const std::string& qid = prec.getQueryIdentifier();
        boost::format seqname_fmt("%d:%d@%s");  // local variable because not thread-safe
        StopWatchCPUTime stopwatch_init("initializing this record");  // log overall time for this predict phase
        stopwatch_init.start();
	// prepare alignment/ score type
	//seqan::Score<int> myBlosum = seqan::Blosum80(); 
	//align_method = &myBlosum;

        // push records into active_records  TODO: remove intermediate active_records?
        active_list_type_ active_records;
        {
            for(typename ContainerT::iterator rec_it = recordset.begin(); rec_it != recordset.end(); ++rec_it) {
                if(!(*rec_it)->isFiltered()) active_records.push_back(*rec_it);
            }
        }
        const uint n = active_records.size();
        
        // with no unmasked alignment, set to unclassified and return
        if(n==0) {  //TODO: record should not be reported at all in GFF3
            const std::string qrseqname = boost::str(seqname_fmt % -1 % -1 % qid);
            logsink << "ID" << tab << qrseqname << std::endl;
            logsink << "  NUMREF" << tab << n << std::endl << std::endl;
            logsink << "    RANGE\t" << this->taxinter_.getRoot()->data->annotation->name << tab << this->taxinter_.getRoot()->data->annotation->name << tab << this->taxinter_.getRoot()->data->annotation->name << std::endl << std::endl;
            logsink << "STATS" << tab << qrseqname << tab << n << "\t0\t0\t0\t0\t0\t0\t0\t.0" << std::endl << std::endl;
            
            TaxonPredictionModel< ContainerT >::setUnclassified(prec);
            return;
        }
        
        // with one alignment, don't align and return
        if(n==1) {
            typename ContainerT::value_type rec = active_records.front();
            large_unsigned_int qrstart = rec->getQueryStart();
            large_unsigned_int qrstop = rec->getQueryStop();
            const std::string qrseqname = boost::str(seqname_fmt % qrstart % qrstop % qid);
            
            logsink << "ID" << tab << qrseqname << std::endl;
            logsink << "  NUMREF" << tab << n << std::endl;
            logsink << "  RANGE\t" << rec->getReferenceNode()->data->annotation->name << tab << rec->getReferenceNode()->data->annotation->name << tab << this->taxinter_.getRoot()->data->annotation->name << std::endl << std::endl;
            logsink << "STATS" << tab << qrseqname << tab << n << "\t0\t0\t0\t0\t0\t0\t0\t.0" << std::endl << std::endl;
            
            prec.setQueryFeatureBegin(qrstart);
            prec.setQueryFeatureEnd(qrstop);
            prec.setInterpolationValue(1.);
            prec.setNodeRange(rec->getReferenceNode(), this->taxinter_.getRoot(), rec->getIdentities());
            prec.setBestReferenceTaxon(rec->getReferenceNode());
            return;
        }
        
        // n>1 -> screen alignments and determine query range
        large_unsigned_int qrstart;
        large_unsigned_int qrstop;
        {
            assert(n>1);
            typename active_list_type_::iterator rec_it = active_records.begin();
            qrstart = (*rec_it)->getQueryStart();
            qrstop = (*rec_it)->getQueryStop();
            ++rec_it;
            do {  // because n >= 2
                qrstart = std::min((*rec_it)->getQueryStart(), qrstart);
                qrstop = std::max((*rec_it)->getQueryStop(), qrstop);
            } while(++rec_it != active_records.end());
        }
        const large_unsigned_int qrlength = qrstop - qrstart + 1;
        
        // logging
        const std::string qrseqname = boost::str(seqname_fmt % qrstart % qrstop % qid);
        logsink << "ID" << tab << qrseqname << std::endl;
        logsink << "  NUMREF" << tab << n << std::endl;
        
        // sort the list by score
        sort_.filter(active_records);
        
        // data storage  TODO: maybe use Boost ptr containers
        const StringType qrseq = query_sequences_.getSequence(qid, qrstart, qrstop);
        
        std::vector< typename ContainerT::value_type > records(n);  //TODO: move below next section and do not create records if q==r_best
        {
            typename active_list_type_::iterator rec_it = active_records.begin();
            uint i = 0;
            while(rec_it != active_records.end()) {
                records[i] = *rec_it;
                ++rec_it;
                ++i;
            }
        }
        const float qmaxscore = records[0]->getScore();
        
        // n>1 and query is identical to reference, we will use local alignment scores only
        // TODO: add option to recalculate (github issue #24)
        if(records[0]->getAlignmentLength() == qrlength && records[0]->getIdentities() == qrlength) {
            float score_best = records[0]->getScore();
            const TaxonNode* lnode = records[0]->getReferenceNode();
            const TaxonNode* unode = nullptr;
            uint i = 1;

            while (true) {  // 2 breaks
                if(i == n) {  // 1 break
                    unode = this->taxinter_.getRoot();
                    break;
                }
                
                float score = records[i]->getScore();
                
                if(score == score_best) {
                    const TaxonNode* cnode = records[i]->getReferenceNode();
                    lnode = this->taxinter_.getLCA(lnode, cnode);
                    logsink << "    current ref/lower node: " << "("<< score <<") "<< lnode->data->annotation->name << " (+ " << cnode->data->annotation->name << " )" << std::endl;
                }
                else {  // 1 break
                    float uscore = score;
                    unode = lnode;
                    do {
                        const TaxonNode* cnode = records[i]->getReferenceNode();
                        unode = this->taxinter_.getLCA(unode, cnode);
                        logsink << "    current upper node: " << "("<< uscore <<") "<< unode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, lnode)->data->root_pathlength) << " )" << std::endl;
                    } while (++i < n && records[i]->getScore() == uscore);
                    break;
                }
                ++i;
            }
                        
            logsink << "  RANGE\t" << lnode->data->annotation->name << tab << lnode->data->annotation->name << tab << unode->data->annotation->name << std::endl << std::endl;
            logsink << "STATS" << tab << qrseqname << tab << n << "\t0\t0\t0\t0\t" << stopwatch_init.read() << "\t0\t0\t.0" << std::endl << std::endl;
            
            prec.setQueryFeatureBegin(qrstart);
            prec.setQueryFeatureEnd(qrstop);
            prec.setInterpolationValue(.0);
            prec.setNodeRange(lnode, unode, qrlength);
            prec.setBestReferenceTaxon(lnode);
            return;
        }
        
        std::vector< StringType > segments(n);    //TODO: don't call element constructors
        std::vector< int > queryscores(n, std::numeric_limits< int >::max());
        std::vector< large_int > querymatches(n, 0);   //TODO: value is not really relevant
        stopwatch_init.stop();
        
        // count number of alignment calculations in each of the three passes
        uint gcounter = 0;
        uint pass_0_counter = 0;
        uint pass_0_counter_naive = 0;
        uint pass_1_counter = 0;
        uint pass_1_counter_naive = 0;
        uint pass_2_counter = 0;
        uint pass_2_counter_naive = 0;
        
        StopWatchCPUTime stopwatch_seqret("retrieving sequences for this record");  // log overall time for this predict phase
        StopWatchCPUTime stopwatch_process("processing this record");  // log overall time for this predict phase
        stopwatch_process.start();

        std::set<uint> qgroup;
        large_unsigned_int anchors_support = 0;
        const TaxonNode* rtax = NULL;  // taxon of closest evolutionary neighbor(s)
        const TaxonNode* lca_allnodes = records.front()->getReferenceNode();  // used for optimization
        
        {   // pass 0 (re-alignment to most similar reference segments)
            logsink << std::endl << "  PASS\t0" << std::endl;
            float dbalignment_score_threshold = reeval_bandwidth_factor_*qmaxscore;
            uint index_best = 0;

            for (uint i = 0; i < n; ++i) { //calculate scores for best-scoring references
                int score;
                large_int matches;
                alignment<StringType> queryalignment;
                //int matches;
                const float qlscore = records[i]->getScore();
                const large_unsigned_int qlmatch = records[i]->getIdentities();
                const double qlpid = static_cast<double>(qlmatch)/qrlength;

                if(records[i]->getAlignmentLength() == qrlength && records[i]->getIdentities() == qrlength) {

                    qgroup.insert(i);
                    score = 0;
                    matches = records[i]->getIdentities();
                    logsink << std::setprecision(2) << "    *ALN " << i << " <=> query" << tab  << "qlscore=" << qlscore << "; qlmatch=" << qlmatch << "; score=" << score << "; match=" << matches << "; qpid=1.0" << std::endl;
                    ++pass_0_counter_naive;
                } else if (records[i]->getScore() >= dbalignment_score_threshold) {

                    qgroup.insert(i);
                    
                    
                    stopwatch_seqret.start();
                    if(seqan::empty(segments[i])) segments[i] = getSequence(records[i]->getReferenceIdentifier(),  records[i]->getReferenceStart(), records[i]->getReferenceStop(), records[i]->getQueryStart() - qrstart, qrstop - records[i]->getQueryStop());
                    stopwatch_seqret.stop();                   
                    //score = -seqan::globalAlignmentScore(segments[i], qrseq, seqan::MyersBitVector());
                    //score = -seqan::globalAlignmentScore(segments[i], qrseq, seqan::Blosum30());
                    queryalignment = getAlignment<seqan::Blosum80>(segments[i],qrseq);
                    score = queryalignment.score;
                    ++pass_0_counter;
                    ++pass_0_counter_naive;
                   
                    //matches = std::max(static_cast<large_int>(std::max(seqan::length(segments[i]), seqan::length(qrseq)) - score), static_cast<large_int>(records[i]->getIdentities()));
                    matches = std::max(static_cast<large_int>(queryalignment.matches), static_cast<large_int>(records[i]->getIdentities()));
                    //matches = std::max((std::max(seqan::length(segments[i]), seqan::length(qrseq)) - score), static_cast<int>records[i]->getIdentities());
                    double qpid = static_cast<double>(matches)/qrlength;
                    logsink << std::setprecision(2) << "    +ALN " << i << " <=> query" << tab  << "qlscore=" << qlscore << "; qlmatch=" << qlmatch << "; qlpid=" << qlpid << "; score=" << score << "; match=" << matches << "; qpid=" << qpid << std::endl;
                    //logsink << queryalignment.alignment << std::endl;
                    
                } else {  // not similar -> fill in some dummy values
                    score = std::numeric_limits< int >::max();
                    matches = records[i]->getIdentities();
                }
                queryscores[i] = score;
                querymatches[i] = matches;
                if (score < queryscores[index_best]) index_best = i;
                else if (score == queryscores[index_best]) {
                    if (matches > querymatches[index_best]) index_best = i;
                    else if (matches == querymatches[index_best] && qlscore > records[index_best]->getScore()) index_best = i;
                }
                //anchors_support = std::max(anchors_support, matches);
                anchors_support = std::max(static_cast<large_int>(anchors_support), matches);  //TODO: move to previous if-statement?
                lca_allnodes = this->taxinter_.getLCA(lca_allnodes, records[i]->getReferenceNode());
            }
            
            // TODO: resort alignments by local score (already done), and by secondary vector (matches)
            // only keep and use the best-scoring reference sequences
            rtax = records[index_best]->getReferenceNode();
            for (std::set< uint >::iterator it = qgroup.begin(); it != qgroup.end();) {
                if (queryscores[*it] != queryscores[index_best] || querymatches[*it] != querymatches[index_best] || records[*it]->getScore() != records[index_best]->getScore()) qgroup.erase(it++);
                else {
                    const TaxonNode* cnode = records[*it]->getReferenceNode();
                    rtax = this->taxinter_.getLCA(rtax, cnode);
                    logsink << "      current ref node: " << "("<< queryscores[*it] <<") "<< rtax->data->annotation->name << " (+ " << cnode->data->annotation->name << " )" << std::endl;
                    ++it;
                }
            }
            assert(! qgroup.empty());  // TODO: only in debug mode
            
            logsink << "    NUMALN\t" << pass_0_counter << tab << pass_0_counter_naive - pass_0_counter << std::endl << std::endl;
        }


        float anchors_taxsig = 1.;  // a measure of tree-like scores  
        float ival_global = 0.;
        const TaxonNode* lnode_global = rtax;
        const TaxonNode* unode_global = rtax;
        std::set<uint> outgroup;
        float bandfactor_max = 1.;
		
        {   // pass 1 (best reference alignment)
            logsink << "  PASS\t1" << std::endl;

            small_unsigned_int lca_root_dist_min = std::numeric_limits<small_unsigned_int>::max();
            do {  // for each most similar reference segment
                BandFactor bandfactor1(this->taxinter_, n);
                const uint index_anchor = *qgroup.begin();
                qgroup.erase(qgroup.begin());
                const int qscore = queryscores[index_anchor];
                const TaxonNode* rnode = records[index_anchor]->getReferenceNode();
                bandfactor1.addSequence(0, rnode);
                const TaxonNode* lnode = rtax;
                const TaxonNode* unode = NULL;
                int lscore = 0;
                int uscore = std::numeric_limits<int>::max();

                std::list< boost::tuple< uint, int > > outgroup_tmp;

                // align all others <=> anchor TODO: adaptive cut-off
                logsink << "      query: (" << qscore << ") unknown" << std::endl;
                pass_1_counter_naive += n - 1;
                
                // TODO: implement heuristic cut-off
                double qpid_upper = 0.;
                double qpid_thresh_guarantee = 0.;
                double qpid_thresh_heuristic = 0.;
                int qlscore_thresh_heuristic = 0.;
                
                for(uint i = 0; lnode != this->taxinter_.getRoot() && i < n && records[i]->getScore() >= qlscore_thresh_heuristic; ++i) {  //TODO: break loop when qlscore < qlscore_thresh_heuristic
                    const TaxonNode* cnode = records[i]->getReferenceNode();
                    const large_unsigned_int qlmatch = records[i]->getIdentities();
                    const double qlpid = static_cast<double>(qlmatch)/qrlength;  // must be used for cutoff for stability reasons
                    const double qpid = static_cast<double>(querymatches[i])/qrlength;
                    const float qlscore = records[i]->getScore();
                    double qpid_thresh = std::max(qpid_thresh_guarantee, qpid_thresh_heuristic);
                    
                    if(qpid >= qpid_thresh) {  //TODO: implement command line option
                        int score;
                        large_unsigned_int matches;
                        alignment<StringType> segmentalignment;
                        
                        if (i == index_anchor) score = 0;
                        else {
                            // use triangle relation to avoid alignment
                            if (queryscores[i] == 0) { // && queryscores[index_anchor] == 0 ) { //&& querymatches[i]) { // TODO: correct?
                                score = queryscores[index_anchor];
                                matches = querymatches[index_anchor];
                            }
                            else {
                                stopwatch_seqret.start();
                                if(seqan::empty(segments[index_anchor])) segments[index_anchor] = getSequence(records[index_anchor]->getReferenceIdentifier(),  records[index_anchor]->getReferenceStart(), records[index_anchor]->getReferenceStop(), records[index_anchor]->getQueryStart() - qrstart, qrstop - records[index_anchor]->getQueryStop());
                                if(seqan::empty(segments[i])) segments[i] = getSequence(records[i]->getReferenceIdentifier(),  records[i]->getReferenceStart(), records[i]->getReferenceStop(), records[i]->getQueryStart() - qrstart, qrstop - records[i]->getQueryStop());
                                stopwatch_seqret.stop();
                                
                                //score = -seqan::globalAlignmentScore(segments[i], segments[index_anchor], seqan::MyersBitVector());
                                //score = getAlignment(segments[i],segments[index_anchor]);
                                segmentalignment = getAlignment<seqan::Blosum80>(segments[i],segments[index_anchor]);
                                score = segmentalignment.score;
                                
                                ++pass_1_counter;
                                //matches = std::max(seqan::length(segments[i]), seqan::length(segments[index_anchor])) - score;
                                large_int matches = segmentalignment.matches;
                                
                                logsink << std::setprecision(2) << "    +ALN " << i << " <=> " << index_anchor << tab << "qlscore=" << qlscore << "; qlmatch=" << qlmatch << "; qlpid=" << qlpid << "; score=" << score << "; match=" << matches << "; qpid=" << qpid << "; qlscore_cut=" << qlscore_thresh_heuristic << "; qpid_cutg=" << qpid_thresh_guarantee << "; qpid_cut_h=" << qpid_thresh_heuristic << std::endl;
                                //logsink << segmentalignment.alignment << std::endl;
                            }
                        }
                        
                        bandfactor1.addSequence(score, cnode);

                        // place sequence
                        if (score == 0) qgroup.erase(i);  // remove this from list of qnodes because it is sequence-identical
                        else {
                            if(score <= qscore) {
                                lnode = this->taxinter_.getLCA(lnode, cnode);
                                if(score > lscore) lscore = score;
                                logsink << "      current lower node: " << "("<< score <<") "<<lnode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
                            }
                            else {
                                if(score < uscore) {  // true if we find a segment with a lower score than query
                                    uscore = score;
                                    if(qlpid > qpid_upper) {
                                        qpid_upper = qlpid;
                                        qpid_thresh_guarantee = qlpid*2. - 1.;  // hardcoded inequation: qpid+1.-qpid_up < qpid_up
                                        qpid_thresh_heuristic = qlpid*exclude_alignments_factor_;
                                    }
                                    if(!qlscore_thresh_heuristic) qlscore_thresh_heuristic = records[i]->getScore()*exclude_alignments_factor_;
                                }
                                outgroup_tmp.push_back(boost::make_tuple(i,score));
                            }
                        }
                    }
                }
                
                float bandfactor = bandfactor1.getFactor();  //TODO: limit and check
                bandfactor_max = std::max(bandfactor_max, bandfactor);
                int qscore_ex = qscore * bandfactor;
                int min_upper_score = std::numeric_limits< int >::max();
                
                logsink << std::endl << "    EXT\tqueryscore = " << qscore << "; threshold = " << qscore_ex << "; bandfactor = " << bandfactor << std::endl;
                for(std::list< boost::tuple<uint,int> >::iterator it = outgroup_tmp.begin(); it != outgroup_tmp.end();) {
                    int score = it->get<1>();

                    if(score > qscore_ex) {
                        if (score > min_upper_score) it = outgroup_tmp.erase(it);
                        else {
                            if(score < min_upper_score) min_upper_score = score;
                            ++it;
                        }
                    } else {
                        if(min_upper_score > qscore_ex) min_upper_score = score;
                        else min_upper_score = std::max(min_upper_score, score);
                        ++it;
                    }
                }
                
                // push elements from temporary to outgroup set
                if(min_upper_score != std::numeric_limits< int >::max()) unode = lnode;
                for(std::list< boost::tuple<uint,int> >::iterator it = outgroup_tmp.begin(); it != outgroup_tmp.end(); ++it) {
                    uint i;
                    int score;
                    boost::tie(i, score) = *it;
                    const TaxonNode* cnode = records[i]->getReferenceNode();
                    
                    if(score > min_upper_score) continue;
                    
                    // add to upper node if(score <= min_upper_score)
                    unode = this->taxinter_.getLCA(cnode, unode);
                    logsink << "      current upper node: " << "("<< score <<") "<< unode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;

                    // curate minimal outgroup TODO: only keep score == min_upper_score in outgroup?
                    const small_unsigned_int lca_root_dist = this->taxinter_.getLCA(cnode, rtax)->data->root_pathlength;
                    if(lca_root_dist > lca_root_dist_min) continue;
                    else if(lca_root_dist < lca_root_dist_min) {
                        lca_root_dist_min = lca_root_dist;
                        outgroup.clear();
                    }
                    outgroup.insert(i);
                }

                // adjust interpolation value and upper node
                float ival = 0.;                                 //  TODO: initialize
                if(!unode) {
                    unode = this->taxinter_.getRoot();
                    uscore = -1;
                    ival = 1.;
                } else if(unode != lnode && lscore < qscore) ival = (qscore - lscore)/static_cast<float>(uscore - lscore);
                
                logsink << std::endl << "    SCORE\tlscore = " << lscore << "; uscore = " << uscore << "; queryscore = " << qscore << "; queryscore_ex = " << qscore_ex << "; ival = " << ival  << std::endl << std::endl;
                const float taxsig = .0;  // TODO: placer.getTaxSignal(qscore);

                ival_global = std::max(ival, ival_global);  // combine interpolation values conservatively
                anchors_taxsig = std::min(taxsig, anchors_taxsig);  // combine taxonomic signal values conservatively
                unode_global = this->taxinter_.getLCA(unode_global, unode);
                lnode_global = this->taxinter_.getLCA(lnode_global, lnode);
                
            } while (! qgroup.empty() && lnode_global != this->taxinter_.getRoot());

            logsink << "    NUMALN\t" << pass_1_counter << tab << pass_1_counter_naive - pass_1_counter << std::endl;
            logsink << "    NUMOUTGRP\t" << outgroup.size() << std::endl;
        }

        logsink << "    RANGE\t" << rtax->data->annotation->name << tab << lnode_global->data->annotation->name << tab << unode_global->data->annotation->name << std::endl << std::endl;
        
        {   // pass 2 (stable upper node estimation alignment)
            logsink << "  PASS\t2" << std::endl;
            while (! outgroup.empty()) {
                const uint index_anchor = *outgroup.begin();
                outgroup.erase(outgroup.begin());
                
                if( unode_global == lca_allnodes ) {
                    if( queryscores[index_anchor] == std::numeric_limits<int>::max() ) pass_2_counter_naive += n;
                    else pass_2_counter_naive += n - 1;
                    continue;
                }

                // align all others <=> anchor TODO: heuristic
                const double qpid_anchor = static_cast<double>(querymatches[index_anchor])/qrlength;
                const double qpid_thresh_guarantee = qpid_anchor*2. - 1.;  // hardcoded inequation: qpid+1.-qpid_up < qpid_up;
                const double qpid_thresh_heuristic = qpid_anchor*exclude_alignments_factor_;
                const double qpid_thresh = std::max(qpid_thresh_guarantee, qpid_thresh_heuristic);
                const float qlscore_thresh_heuristic = records[index_anchor]->getScore()*exclude_alignments_factor_;
                ++pass_2_counter_naive; // query angainst reference alignment

                for (uint i = 0; i < n && records[i]->getScore() >= qlscore_thresh_heuristic; ++i) {
                    const double qpid = static_cast<double>(querymatches[i])/qrlength;
                    if(qpid >= qpid_thresh) {

                        const TaxonNode* cnode = records[i]->getReferenceNode();
                        const float qlscore = records[i]->getScore();
                        const large_unsigned_int qlmatch = records[i]->getIdentities();
                        int score;
                        alignment<StringType> segmentalignment;

                        if (i == index_anchor) score = 0;
                        else {
                            ++pass_2_counter_naive;
                            if( this->taxinter_.isParentOf(unode_global, cnode) || cnode == unode_global ) continue;
                            else {
                                stopwatch_seqret.start();
                                if(seqan::empty(segments[index_anchor])) segments[index_anchor] = getSequence(records[index_anchor]->getReferenceIdentifier(),  records[index_anchor]->getReferenceStart(), records[index_anchor]->getReferenceStop(), records[index_anchor]->getQueryStart() - qrstart, qrstop - records[index_anchor]->getQueryStop());
                                if(seqan::empty(segments[i])) segments[i] = getSequence(records[i]->getReferenceIdentifier(),  records[i]->getReferenceStart(), records[i]->getReferenceStop(), records[i]->getQueryStart() - qrstart, qrstop - records[i]->getQueryStop());
                                stopwatch_seqret.stop();
                                
                                //score = -seqan::globalAlignmentScore(segments[i], segments[index_anchor], seqan::MyersBitVector());
                                //score = getAlignment(segments[i], segments[index_anchor]);
                                segmentalignment = getAlignment<seqan::Blosum80>(segments[i], segments[index_anchor]);
                                score = segmentalignment.score;
                                
                                
                                logsink << std::setprecision(2) << "    +ALN " << i << " <=> " << index_anchor << tab << "qlscore=" << qlscore << "; qlmatch=" << qlmatch << "; score=" << score << "; qpid=" << qpid << std::endl;
                                //logsink << segmentalignment.alignment << std::endl;
                                ++pass_2_counter;
                                queryscores[i] = score;
                            }
                        }

                        if (score == 0) outgroup.erase(i);
                        else {
                            int qscore_ex;
                            if (queryscores[index_anchor] == std::numeric_limits<int>::max()) { //need to align query <=> anchor
                                stopwatch_seqret.start();
                                if(seqan::empty(segments[index_anchor])) segments[index_anchor] = getSequence(records[index_anchor]->getReferenceIdentifier(),  records[index_anchor]->getReferenceStart(), records[index_anchor]->getReferenceStop(), records[index_anchor]->getQueryStart() - qrstart, qrstop - records[index_anchor]->getQueryStop());
                                stopwatch_seqret.stop();
                                
                                //int score = -seqan::globalAlignmentScore(segments[index_anchor], qrseq, seqan::MyersBitVector());
                                //int score = getAlignment(segments[index_anchor],qrseq);
                                segmentalignment =  getAlignment<seqan::Blosum80>(segments[index_anchor], qrseq);
                                int score = segmentalignment.score;
                                
                                //large_int matches = std::max(static_cast<large_int>(std::max(seqan::length(segments[index_anchor]), seqan::length(qrseq)) - score), querymatches[index_anchor]);
                                //TODO
                                large_int matches = std::max(static_cast<large_int>(segmentalignment.matches), querymatches[index_anchor]);
                                
                                double qpid = static_cast<double>(matches)/qrlength;
                                logsink << std::setprecision(2) << "    +ALN query <=> " << index_anchor << tab << "qlscore=" << records[index_anchor]->getScore() << "; qlmatch=" << qlmatch << "; score=" << score << "; match=" << matches << "; qpid=" << qpid << std::endl;
                                queryscores[index_anchor] = score;
                                querymatches[index_anchor] = matches;
                                qscore_ex = score*bandfactor_max;
                                logsink << "      query: (" << qscore_ex << ") unknown" << std::endl;
                                //logsink << segmentalignment.alignment << std::endl;
                                ++pass_2_counter;
                            } else qscore_ex = queryscores[index_anchor]*bandfactor_max;

                            if(score <= qscore_ex) {
                                const TaxonNode* rnode = records[index_anchor]->getReferenceNode();
                                unode_global = this->taxinter_.getLCA(unode_global, cnode);
                                logsink << "      current upper node: " << "("<< score <<") "<< unode_global->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
                            }
                        }
                    }
                }
                logsink << std::endl;
            }
            logsink << "    NUMALN\t" << pass_2_counter << tab << pass_2_counter_naive - pass_2_counter << std::endl;
        }

        if(unode_global == lnode_global) ival_global = 1.;
        logsink << "    RANGE\t" << rtax->data->annotation->name << tab << lnode_global->data->annotation->name << tab << unode_global->data->annotation->name << std::endl << std::endl;

        prec.setSignalStrength(anchors_taxsig);
        prec.setQueryFeatureBegin(qrstart);
        prec.setQueryFeatureEnd(qrstop);
        prec.setInterpolationValue(ival_global);
        prec.setNodeRange(lnode_global, unode_global, anchors_support);
        prec.setBestReferenceTaxon(rtax);
        gcounter = pass_0_counter + pass_1_counter + pass_2_counter;
        float normalised_rt = (float)gcounter/(float)n;
        stopwatch_process.stop();
        logsink << "STATS" << tab << qrseqname << tab << n << tab << pass_0_counter << tab << pass_1_counter << tab << pass_2_counter << tab << gcounter << tab << stopwatch_init.read() << tab << stopwatch_seqret.read() << tab << stopwatch_process.read() << tab << std::setprecision(2) << std::fixed << normalised_rt << std::endl << std::endl;
    }
    
    const StringType getSequence(const std::string& id, const large_unsigned_int start, const large_unsigned_int stop, const large_unsigned_int left_ext = 0, const large_unsigned_int right_ext = 0 ) {
        if(typeid(StringType) == typeid(seqan::String<seqan::Dna5>)){
            if(start <= stop) {
            large_unsigned_int newstart = left_ext < start ? start - left_ext : 1;
            large_unsigned_int newstop = stop + right_ext;
            return db_sequences_.getSequence(id, newstart, newstop); //TODO: can we avoid copying
            }
            large_unsigned_int newstart = right_ext < stop ? stop - right_ext : 1;
            large_unsigned_int newstop = start + left_ext;
            return db_sequences_.getSequenceReverseComplement(id, newstart, newstop); 
        }else if(typeid(StringType) == typeid(seqan::String<seqan::AminoAcid>)){
        assert(start <= stop);
        large_unsigned_int newstart = left_ext < start ? start - left_ext : 1;
        large_unsigned_int newstop = stop + right_ext;
        return db_sequences_.getSequence(id, newstart, newstop); //TODO: can we avoid copying
        }
        else{
            assert(false);
        }
    }

protected:
    typedef std::list<typename ContainerT::value_type> active_list_type_;
    QStorType& query_sequences_;
    const DBStorType& db_sequences_;
    SortFilter< active_list_type_ > sort_;
    compareTupleFirstLT< boost::tuple< int, uint >, 0 > tuple_1_cmp_le_;

private:
    const float exclude_alignments_factor_;
    const float reeval_bandwidth_factor_;
    StopWatchCPUTime measure_sequence_retrieval_;
    StopWatchCPUTime measure_pass_0_alignment_;
    StopWatchCPUTime measure_pass_1_alignment_;
    StopWatchCPUTime measure_pass_2_alignment_;
    //const seqan::Score<int> *align_method;
    
//    int getAlignment(StringType A, StringType B){
//        int selfcomp = seqan::globalAlignmentScore(B,B,seqan::Blosum80());
//        int alignscore = seqan::globalAlignmentScore(A,B,seqan::Blosum80());
//        int returnscore = selfcomp - alignscore;
//        //std::cerr << A << std::endl;
//        //std::cerr << B << std::endl;
//        //std::cerr << selfcomp <<" "<< alignscore <<" "<< returnscore <<" "<< std::endl;
//        assert(selfcomp >= alignscore);
//        assert(returnscore >= 0);
//        return returnscore;
//        //return -seqan::globalAlignmentScore(A,B,seqan::Blosum30());
//    }
template<typename AlignMethod>
alignment<StringType> getAlignment(StringType A, StringType B){
	
	alignment<StringType> returnalignment;
            
    typedef typename seqan::Align<StringType, seqan::ArrayGaps> TAlign;
    typedef typename seqan::Row<TAlign>::Type TRow;
    typedef typename seqan::Iterator<TRow>::Type TRowIterator;	

 	TAlign selfalignA;
    TAlign selfalignB; 
	TAlign diffalign;
        
	int selfcomp;
	int alignscore;

	if(typeid(StringType)==typeid(seqan::String<seqan::AminoAcid>)){ 
        
        //align sequence to itself
       
        //int selfcomp = seqan::globalAlignment(selfalignA,AlignMethod());
        //selfcomp += seqan::globalAlignment(selfalignB,AlignMethod());
        
        //int alignscore = seqan::globalAlignment(diffalign,AlignMethod()); 
	

        resize(rows(selfalignA), 2);
        assignSource(row(selfalignA, 0), A);
        assignSource(row(selfalignA, 1), A);

	selfcomp = seqan::globalAlignment(selfalignA,seqan::Blosum80());

	resize(rows(selfalignB), 2);
        assignSource(row(selfalignB, 0), B);
        assignSource(row(selfalignB, 1), B);

        selfcomp += seqan::globalAlignment(selfalignB,seqan::Blosum80());
        

        resize(rows(diffalign), 2);
        assignSource(row(diffalign, 0), A);
        assignSource(row(diffalign, 1), B);

        alignscore = seqan::globalAlignment(diffalign,seqan::Blosum80()); 
	
	}
	else if(typeid(StringType)==typeid(seqan::String<seqan::Dna5>)){
        //align sequence to itself
		
        resize(rows(selfalignA), 2);
        assignSource(row(selfalignA, 0), A);
        assignSource(row(selfalignA, 1), A);
        
        selfcomp = seqan::globalAlignment(selfalignA,seqan::MyersHirschberg());

		resize(rows(selfalignB), 2);
        assignSource(row(selfalignB, 0), B);
        assignSource(row(selfalignB, 1), B);

        selfcomp += seqan::globalAlignment(selfalignB,seqan::MyersHirschberg());
        
	// score without alignment

        resize(rows(diffalign), 2);
        assignSource(row(diffalign, 0), A);
        assignSource(row(diffalign, 1), B);
        alignscore = seqan::globalAlignment(diffalign,seqan::MyersHirschberg()); 	
	}
	else{
		std::cerr << "Something went wrong";
		assert(false);
	}
        
        //logsink << diffalign << std::endl;
        
        TRow & row1 = row(diffalign, 0);
        TRow & row2 = row(diffalign, 1);
        
        TRowIterator itRow1 = begin(row1);
        TRowIterator itEndRow1 = end(row1);
        TRowIterator itRow2 = begin(row2);
        
        int gapCount = 0;
        int matchCount = 0;
        int missmatchCount = 0;
        
        
        
        for (; itRow1 != itEndRow1; ++itRow1, ++itRow2){
        if(seqan::isGap(itRow1) || seqan::isGap(itRow2)){  
            gapCount ++;
        }
        else if(*itRow1==*itRow2){
                matchCount ++;
            }
            else{
                missmatchCount ++;
            }
        
        }
        
        returnalignment.score = selfcomp-alignscore;
        returnalignment.matches = matchCount;
        returnalignment.mmatches = missmatchCount;
        returnalignment.gaps = gapCount;
        returnalignment.alignment = diffalign;
        //std::cerr << A << std::endl;
        //std::cerr << B << std::endl;
        //std::cerr << selfcomp <<" "<< alignscore << std::endl;
	//std::cerr << diffalign << std::endl;
	//std::cerr << "score: " << returnalignment.score << std::endl;
	//std::cerr << "match: " << matchCount << std::endl;
	//std::cerr << "mm: " << missmatchCount << std::endl;
	//std::cerr << "gaps :" << gapCount << std::endl; 
	//assert(selfcomp >= alignscore);
        //assert(returnalignment.score >= 0);
	
        return returnalignment;
        //return returnscore;
        //return -seqan::globalAlignmentScore(A,B,seqan::Blosum30());
    }
};

#endif // taxonpredictionmodelsequence_hh_
