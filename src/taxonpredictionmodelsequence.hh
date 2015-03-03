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
template< typename ContainerT, typename QStorType, typename DBStorType >
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


        // determine best local alignment score and number of unmasked records
        
        
//         for(typename ContainerT::iterator rec_it = recordset.begin(); rec_it != recordset.end(); ++rec_it) {
//             if(!(*rec_it)->isFiltered()) {
//                 qmaxscore = std::max(qmaxscore, (*rec_it)->getScore());
//                 ++n_pre;
//             }
//         }

        // go through records and determine maximum local alignment score
        active_list_type_ active_records;
        {
//             const float threshold = qmaxscore*exclude_alignments_factor_;
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
        sort_.filter(active_records);  //TODO: optimize sorting algo
        
        // data storage  TODO: use Boost ptr containers
        const seqan::Dna5String qrseq = query_sequences_.getSequence(qid, qrstart, qrstop);
        
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
        if(records[0]->getAlignmentLength() == qrlength && records[0]->getIdentities() == qrlength) {
            float score_best = records[0]->getScore();
            const TaxonNode* lnode = records[0]->getReferenceNode();
            float uscore = score_best;
            uint i = 1;
            while(i < n) {
                float score = records[i]->getScore();
                if(score == score_best) lnode = this->taxinter_.getLCA(lnode, records[i++]->getReferenceNode());
                else {
                    uscore = score;
                    ++i;
                    break;
                }
            }
            
            const TaxonNode* unode;
            if(i < n) {
                unode = lnode;
                do {
                    if(records[i]->getScore() == uscore) unode = this->taxinter_.getLCA(unode, records[i]->getReferenceNode());
                    else break;
                } while(++i<n);
            } else {
                unode = this->taxinter_.getRoot();
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
        
//         records.reserve(n);
        
//         std::vector< seqan::Dna5String > segments; //TODO: boost ptr_container/seqan::StringSet/set to detect equal sequences
//         segments.reserve(n);
        
        std::vector<seqan::Dna5String> segments(n);    //TODO: don't call element constructors
//         std::vector< std::string > segments_names;
        std::vector< int > queryscores(n, std::numeric_limits< int >::max());
//         queryscores.reserve(n);
        std::vector< large_unsigned_int > querymatches(n, 0);   //TODO: value is not really relevant
//         querymatches.reserve(n);
//         std::cerr << "initialization completed." << std::endl;
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
/*        {   // retrieve segment sequences
            std::cerr << "Retrieving " << n << " sequences for record " << qrseqname << std::endl;
            stopwatch_seqret.start();

            for(typename active_list_type_::iterator rec_it = active_records.begin(); rec_it != active_records.end(); ++rec_it) {
                records.push_back(*rec_it);
                AlignmentRecordTaxonomy& rec = **rec_it;
                large_unsigned_int left_ext = rec.getQueryStart() - qrstart;  //TODO: check for error
                large_unsigned_int right_ext = qrstop - rec.getQueryStop();
                segments.push_back(getSequence(rec.getReferenceIdentifier(),  rec.getReferenceStart(), rec.getReferenceStop(), left_ext, right_ext));

                // cut out reference region
//                 if(rstart <= rstop) {
//                     large_unsigned_int start = left_ext < rstart ? rstart - left_ext : 1;
//                     large_unsigned_int stop = rstop + right_ext;
//                     segments.push_back(db_sequences_.getSequence(rid, start, stop)); //TODO: avoid copying
// //                         segments_names.push_back(boost::str(seqname_fmt % start % stop % rid));
//                 } else {
//                     large_unsigned_int start = right_ext < rstop ? rstop - right_ext : 1;
//                     large_unsigned_int stop = rstart + left_ext;
//                     segments.push_back(db_sequences_.getSequenceReverseComplement(rid, start, stop)); //TODO: avoid copying
// //                         segments_names.push_back(boost::str(seqname_fmt % stop % start % rid));
//                 }
            }
            stopwatch_seqret.stop();
        }*/
        
        StopWatchCPUTime stopwatch_process("processing this record");  // log overall time for this predict phase
        stopwatch_process.start();

        std::set<uint> qgroup;
        large_unsigned_int anchors_support = 0;
        const TaxonNode* rtax = NULL;  // taxon of closest evolutionary neighbor(s)
        const TaxonNode* lca_allnodes = records.front()->getReferenceNode();  // used for optimization
        
        {   // pass 0 (re-alignment to most similar reference segments)
            logsink << std::endl << "  PASS\t0" << std::endl;
//             std::cerr << "Working on record " << qrseqname << std::endl;
            float dbalignment_score_threshold = reeval_bandwidth_factor_*qmaxscore;
            uint index_best = 0;

            for (uint i = 0; i < n; ++i) { //calculate scores for best-scoring references
                int score;
                large_unsigned_int matches;
                const float qlscore = records[i]->getScore();
                
                if(records[i]->getAlignmentLength() == qrlength && records[i]->getIdentities() == qrlength) {
                    qgroup.insert(i);
                    score = 0;
                    matches = records[i]->getIdentities();
                    double qpid = static_cast<double>(matches)/qrlength; 
                    logsink << std::setprecision(2) << "    *ALN " << i << " <=> query" << tab  << "qscore_loc = " << qlscore << "; qpid = " << qpid << "; score = " << score << "; matches = " << matches << std::endl;
                    ++pass_0_counter_naive;
                } else if (records[i]->getScore() >= dbalignment_score_threshold) {
                    qgroup.insert(i);
                    
                    
                    stopwatch_seqret.start();
                    if(seqan::empty(segments[i])) segments[i] = getSequence(records[i]->getReferenceIdentifier(),  records[i]->getReferenceStart(), records[i]->getReferenceStop(), records[i]->getQueryStart() - qrstart, qrstop - records[i]->getQueryStop());
                    stopwatch_seqret.stop();                   
                    score = -seqan::globalAlignmentScore(segments[i], qrseq, seqan::MyersBitVector());
                    
                    ++pass_0_counter;
                    ++pass_0_counter_naive;
                    matches = std::max(static_cast<large_unsigned_int>(std::max(seqan::length(segments[i]), seqan::length(qrseq)) - score), records[i]->getIdentities());
                    double qpid = static_cast<double>(matches)/qrlength;
                    logsink << std::setprecision(2) << "    +ALN " << i << " <=> query" << tab  << "qscore_loc = " << qlscore << "; qpid = " << qpid << "; score = " << score << "; matches = " << matches << std::endl;
                } else {  // not similar -> fill in some dummy values
                    score = std::numeric_limits< int >::max();
                    matches = records[i]->getIdentities();
                }
                queryscores[i] = score;
                querymatches[i] = matches;
                
                if (score < queryscores[index_best]) index_best = i; // || (score == queryscores[index_best] && matches > querymatches[index_best])) {
                else if(score == queryscores[index_best]) {
                    if(matches > querymatches[index_best] || (matches == querymatches[index_best] && qlscore < records[index_best]->getScore())) index_best = i;
                }
                
                anchors_support = std::max(anchors_support, matches);  //TODO: move to previous if-statement?
                
                lca_allnodes = this->taxinter_.getLCA(lca_allnodes, records[i]->getReferenceNode());
            }

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
            
            assert(! qgroup.empty());  //debug code
            
            logsink << "    NUMALN\t" << pass_0_counter << tab << pass_0_counter_naive - pass_0_counter << std::endl << std::endl;
        }

        float anchors_taxsig = 1.;                          // a measure of tree-like scores  
        float ival_global = 0.;
        const TaxonNode* lnode_global = rtax;
        const TaxonNode* unode_global = rtax;
        std::set<uint> outgroup;                 // outgroup sequences
        float bandfactor_max = 1.;

        {   // pass 1 (best reference alignment)
            logsink << "  PASS\t1" << std::endl;
//             const std::size_t num_qnodes = qgroup.size();
//             anchors_lnode.reserve(num_qnodes);  //TODO: performance check
//             anchors_unode.reserve(num_qnodes + 1);  // minimum

//             uint alignments_counter = 0;
//             uint alignments_counter_naive = 0;
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
                    const double qpid = static_cast<double>(querymatches[i])/qrlength;
                    const float qlscore = records[i]->getScore();
//                     bool skip_alignment = qpid < qpid_thresh_guarantee;
//                     float qpid_est = qpid + 1.0 - qpid_upper;
//                     bool skip_alignment = qpid_est < qpid_upper;
                    double qpid_thresh = std::max(qpid_thresh_guarantee, qpid_thresh_heuristic);
                    
                    if(qpid >= qpid_thresh) {  //TODO: implement algorithm parameter and command line option
//                         std::cerr << "Working on record " << qrseqname << std::endl;
                        int score;
                        large_unsigned_int matches;
                        
                        if (i == index_anchor) score = 0;
                        else {
                            // use triangle relation to avoid alignment
                            if (queryscores[i] == 0) { // && queryscores[index_anchor] == 0 ) { //&& querymatches[i]) {
                                score = queryscores[index_anchor];
                                matches = querymatches[index_anchor];
                            }
                            else {
                                stopwatch_seqret.start();
                                if(seqan::empty(segments[index_anchor])) segments[index_anchor] = getSequence(records[index_anchor]->getReferenceIdentifier(),  records[index_anchor]->getReferenceStart(), records[index_anchor]->getReferenceStop(), records[index_anchor]->getQueryStart() - qrstart, qrstop - records[index_anchor]->getQueryStop());
                                if(seqan::empty(segments[i])) segments[i] = getSequence(records[i]->getReferenceIdentifier(),  records[i]->getReferenceStart(), records[i]->getReferenceStop(), records[i]->getQueryStart() - qrstart, qrstop - records[i]->getQueryStop());
                                stopwatch_seqret.stop();
                                
                                score = -seqan::globalAlignmentScore(segments[i], segments[index_anchor], seqan::MyersBitVector());
                                ++pass_1_counter;
//                                 ++alignments_counter;
                                matches = std::max(seqan::length(segments[i]), seqan::length(segments[index_anchor])) - score;
                                logsink << std::setprecision(2) << "    +ALN " << i << " <=> " << index_anchor << tab << "qscore_loc = " << qlscore << "; qpid = " << qpid << "; score = " << score << "; matches = " << matches << std::endl;
                                
                                // update query alignment scores using triangle relation
//                                 if(score==0) {
//                                     queryscores[i] = score;
//                                     querymatches[i] = matches;
//                                 }
//                                 if (queryscores[index_anchor] == 0 && querymatches[i]) {  // TODO: why querymatches[i]? 
//                                     queryscores[i] = score;
//                                     querymatches[i] = matches;
//                                 }
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
    //                             if(lnode == this->taxinter_.getRoot()) {
    //                                 unode = this->taxinter_.getRoot();
    //                                 break;
    //                             }
//                                 if(skip_alignment) logsink << "    ERROR" << tab << "score <= qscore but excluded!" << std::endl;
                            }
                            else {
                                if(score < uscore) {  // true if we find a segment with a lower score than query
                                    uscore = score;
                                    if(qpid > qpid_upper) {
                                        qpid_upper = qpid;
                                        qpid_thresh_guarantee = qpid*2. - 1.;  // hardcoded inequation: qpid+1.-qpid_up < qpid_up
                                        qpid_thresh_heuristic = qpid*exclude_alignments_factor_;
                                    }
                                    if(!qlscore_thresh_heuristic) qlscore_thresh_heuristic = records[i]->getScore()*exclude_alignments_factor_;
                                    
//                                     qpid_upper = std::max(qpid_upper, querymatches[i]/float(qrlength));  //TODO: check for errors
                                }
                                outgroup_tmp.push_back(boost::make_tuple(i,score));
                            }

//                             logsink << std::setprecision(2) << "      EXCL_ALN" << tab << "qpid=" << qpid << "; qpid_thresh=" << qpid_thresh << "; excl=" << skip_alignment << std::endl;
//                             logsink << "      EXCL_ALN" << tab << "qpid=" << qpid << "; qpid_est=" << qpid_est << "; qpid_upper=" << qpid_upper << "; excl=" << skip_alignment << std::endl;
                            
    //                             // TODO: change logic in outgroup; keep highest taxon with worst score
    //                             if(score == min_upper_score) {
    //                                 unode = this->taxinter_.getLCA(cnode, this->taxinter_.getLCA(lnode, unode));
    //                                 logsink << "      current upper node: " << "("<< score <<") "<< unode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
    //                             }
    //                             else if (score < min_upper_score) {
    //                                 uscore = score;
    //                                 min_upper_score = score;
    //                                 unode = this->taxinter_.getLCA(cnode, lnode);
    //                                 outgroup_tmp.clear();
    //                                 logsink << "      current upper node: " << "("<< score <<") "<< unode->data->annotation->name << " (* " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
    //                             }
    //                         }
                        }
                    }
//                     } else logsink << std::setprecision(2) << "    -ALN " << i << " <=> " << index_anchor << tab << "qscore_loc = " << qlscore << "; qpid = " << qpid << "; qpid_thresh = " << qpid_thresh << "; qlscore_thresh = " << qlscore_thresh_heuristic << std::endl;
                }
                
                float bandfactor = bandfactor1.getFactor();  //TODO: limit and check
                bandfactor_max = std::max(bandfactor_max, bandfactor);
                int qscore_ex = qscore * bandfactor;
                int min_upper_score = std::numeric_limits< int >::max();
                
                logsink << std::endl << "    EXT\tqueryscore = " << qscore << "; threshold = " << qscore_ex << "; bandfactor = " << bandfactor << std::endl;
                for(std::list< boost::tuple<uint,int> >::iterator it = outgroup_tmp.begin(); it != outgroup_tmp.end();) {
//                     uint i;
                    int score = it->get<1>();
//                     boost::tie(i, score) = *it;
//                     const TaxonNode* cnode = records[i]->getReferenceNode();
//                     if(score <= qscore_ex) {
//                         unode = this->taxinter_.getLCA(unode, cnode);
//                         min_upper_score = std::min(min_upper_score, score);
//                         if(score > lscore) lscore = score;
//                         logsink << "      current upper node: " << "("<< score <<") "<<unode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
//                         ++it;
//                     }
//                     else {
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
//                     else {
//                         if(score == min_upper_score) {
//                             unode = this->taxinter_.getLCA(cnode, unode);
//                             logsink << "      current upper node: " << "("<< score <<") "<< unode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
//                         }
//                         else {  // if (score > min_upper_score)
//                             min_upper_score = uscore = score;
//                             unode = this->taxinter_.getLCA(cnode, lnode);
//                             logsink << "      current upper node: " << "("<< score <<") "<< unode->data->annotation->name << " (* " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
//                         }
//                         ++it;
//                     }
//                     }
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
//                 anchors_lnode.push_back(lnode);
//                 anchors_unode.push_back(unode);
                
            } while (! qgroup.empty() && lnode_global != this->taxinter_.getRoot());

            logsink << "    NUMALN\t" << pass_1_counter << tab << pass_1_counter_naive - pass_1_counter << std::endl;
            logsink << "    NUMOUTGRP\t" << outgroup.size() << std::endl;
        
//             { // retain only taxa in outgroup which are closest to root
//                 std::size_t tmp_size = outgroup.size();
//                 for (std::map<uint, small_unsigned_int>::iterator it = outgroup.begin(); it !=  outgroup.end();) {
//                     if(it->second > lca_root_dist_min) outgroup.erase(it++);  // TODO: check if this works
//                     else ++it;
//                 }
//                 logsink << "    NUMOUTGRP\t" << outgroup.size() << tab << tmp_size - outgroup.size() << std::endl;
//             }
            
//         anchors_unode.push_back(unode);
        }

        // avoid unneccessary computation if root node is reached TODO: where?
//         unode = this->taxinter_.getLCA(anchors_unode);
//         anchors_unode.clear();  // TODO: remove; code simplification
        logsink << "    RANGE\t" << rtax->data->annotation->name << tab << lnode_global->data->annotation->name << tab << unode_global->data->annotation->name << std::endl << std::endl;
        
        {   // pass 2 (stable upper node estimation alignment)
            logsink << "  PASS\t2" << std::endl;
//             std::cerr << "Working on record " << qrseqname << std::endl;
//             uint alignments_counter = 0;
//             uint alignments_counter_naive = 0;
            while (! outgroup.empty()) {
//                 std::list<uint> reevaluate;
//                 BandFactor bandfactor2(this->taxinter_, n);
                const uint index_anchor = *outgroup.begin();
                outgroup.erase(outgroup.begin());
                
//                 std::cerr << "after erase: " << outgroup.size() << std::endl;
//                 unode = records[index_anchor]->getReferenceNode();
//                 lnode = NULL;
//                 bandfactor2.addSequence(0, rnode);
                
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
                        int score;

                        if (i == index_anchor) score = 0;
                        else {
                            ++pass_2_counter_naive;
                            if( this->taxinter_.isParentOf(unode_global, cnode) || cnode == unode_global ) continue;
                            else {
                                stopwatch_seqret.start();
                                if(seqan::empty(segments[index_anchor])) segments[index_anchor] = getSequence(records[index_anchor]->getReferenceIdentifier(),  records[index_anchor]->getReferenceStart(), records[index_anchor]->getReferenceStop(), records[index_anchor]->getQueryStart() - qrstart, qrstop - records[index_anchor]->getQueryStop());
                                if(seqan::empty(segments[i])) segments[i] = getSequence(records[i]->getReferenceIdentifier(),  records[i]->getReferenceStart(), records[i]->getReferenceStop(), records[i]->getQueryStart() - qrstart, qrstop - records[i]->getQueryStop());
                                stopwatch_seqret.stop();
                                
                                score = -seqan::globalAlignmentScore(segments[i], segments[index_anchor], seqan::MyersBitVector());
                                logsink << std::setprecision(2) << "    +ALN " << i << " <=> " << index_anchor << tab << "qscore_loc = " << qlscore << "; qpid = " << qpid << "; score = " << score <<  std::endl;
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
                                
                                int score = -seqan::globalAlignmentScore(segments[index_anchor], qrseq, seqan::MyersBitVector());
                                large_unsigned_int matches = std::max(static_cast<large_unsigned_int>(std::max(seqan::length(segments[index_anchor]), seqan::length(qrseq)) - score), querymatches[index_anchor]);
                                double qpid = static_cast<double>(matches)/qrlength;
                                logsink << std::setprecision(2) << "    +ALN query <=> " << index_anchor << tab << "qscore_loc = " << records[index_anchor]->getScore() << "; qpid = " << qpid << "; score = " << score << "; matches = " << matches << std::endl;
                                queryscores[index_anchor] = score;
                                querymatches[index_anchor] = matches;
                                qscore_ex = score*bandfactor_max;
                                logsink << "      query: (" << qscore_ex << ") unknown" << std::endl;
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

//                 {
//                     int qscore_ex = qscore * bandfactor2.getFactor();
//                     logsink << std::endl << "    EXT\tqueryscore = " << qscore << "; threshold = " << qscore_ex << std::endl;
//                     std::list<uint>::iterator it = reevaluate.begin();
//                     while(it != reevaluate.end() && unode != lca_allnodes ) {
//                         const int score = queryscores[*it];
//                         const TaxonNode* cnode = records[*it]->getReferenceNode();
//                         if(score <= qscore_ex) {
//                             unode = this->taxinter_.getLCA(cnode, unode);
//                             logsink << "      current upper node: " << "("<< score <<") "<< unode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
//                         }
//                         ++it;
//                     }
//                 }

                logsink << std::endl;
            }
            logsink << "    NUMALN\t" << pass_2_counter << tab << pass_2_counter_naive - pass_2_counter << std::endl;
        }

//         lnode = this->taxinter_.getLCA(anchors_lnode);
//         if(! anchors_unode.empty()) unode = this->taxinter_.getRoot();  // root if no outgroup element
//         else unode = this->taxinter_.getLCA(anchors_unode);
        
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
    
    const seqan::Dna5String getSequence(const std::string& id, const large_unsigned_int start, const large_unsigned_int stop, const large_unsigned_int left_ext = 0, const large_unsigned_int right_ext = 0 ) {
        if(start <= stop) {
            large_unsigned_int newstart = left_ext < start ? start - left_ext : 1;
            large_unsigned_int newstop = stop + right_ext;
            return db_sequences_.getSequence(id, newstart, newstop); //TODO: can we avoid copying
        }
        large_unsigned_int newstart = right_ext < stop ? stop - right_ext : 1;
        large_unsigned_int newstop = start + left_ext;
        return db_sequences_.getSequenceReverseComplement(id, newstart, newstop); //TODO: can we avoid copying
    }

protected:
    typedef std::list<typename ContainerT::value_type> active_list_type_;
    QStorType& query_sequences_;
    const DBStorType& db_sequences_;
    SortByScoreFilter< active_list_type_ > sort_;
    compareTupleFirstLT< boost::tuple< int, uint >, 0 > tuple_1_cmp_le_;

private:
    const float exclude_alignments_factor_;
    const float reeval_bandwidth_factor_;
    StopWatchCPUTime measure_sequence_retrieval_;
    StopWatchCPUTime measure_pass_0_alignment_;
    StopWatchCPUTime measure_pass_1_alignment_;
    StopWatchCPUTime measure_pass_2_alignment_;
};

#endif // taxonpredictionmodelsequence_hh_
