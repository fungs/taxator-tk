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
        return bandfactor_;
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
class DoubleAnchorRPAPredictionModel : public TaxonPredictionModel< ContainerT > {
public:
    DoubleAnchorRPAPredictionModel(const Taxonomy* tax, QStorType& q_storage, const DBStorType& db_storage, float exclude_factor ,float reeval_bandwidth = .1) :
        TaxonPredictionModel< ContainerT >(tax),
        query_sequences_(q_storage),
        db_sequences_(db_storage),
        exclude_alignments_factor_(exclude_factor),
        reeval_bandwidth_factor_(1. - reeval_bandwidth), //TODO: check range
        measure_sequence_retrieval_("sequence retrieval using index"),
        measure_pass_0_alignment_("best reference re-evaluation alignments (pass 0)"),
        measure_pass_1_alignment_("best reference anchor alignments (pass 1)"),
        measure_pass_2_alignment_("distant anchor alignments (pass 2)"),
        measure_placement_algorithm_("placement algorithm")
    {};

    void predict(ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink) {
        this->initPredictionRecord(recordset, prec);  // just set querz name and length

        // cannot be part of class because it is not thread-safe
        boost::format db_seqname_fmt("%d:%d@%s tax=%s");
        boost::format query_seqname_fmt("%d:%d@%s len=%d");

        typename ContainerT::iterator rec_it = firstUnmaskedIter(recordset);

        if (rec_it == recordset.end()) {  // no alignments -> set unclassified
            TaxonPredictionModel< ContainerT >::setUnclassified(prec);
            return;
        }
        
        // set shortcuts
        const std::string& qid = (*rec_it)->getQueryIdentifier();
        large_unsigned_int qrstart = (*rec_it)->getQueryStart();
        large_unsigned_int qrstop = (*rec_it)->getQueryStop();
        const large_unsigned_int qlength = (*rec_it)->getQueryLength();
        float qmaxscore = .0;


        // determine position range of query to consider
        uint n = 0;  // number of alignments/reference segments
        
        // determine maximum alignment score
        do if(!(*rec_it)->isFiltered()) qmaxscore = std::max(qmaxscore, (*rec_it)->getScore());
        while(++rec_it != recordset.end());

        rec_it = firstUnmaskedIter(recordset);
        assert(qrstart <= qrstop);

        do {  // screen alignments and determine query range
            assert((*rec_it)->getQueryStart() <= (*rec_it)->getQueryStop());  // TODO: output error and exit
            if (! (*rec_it)->isFiltered()) {
                if((*rec_it)->getScore() < qmaxscore*exclude_alignments_factor_) (*rec_it)->filterOut();
                else {
                    qrstart = std::min((*rec_it)->getQueryStart(), qrstart);
                    qrstop = std::max((*rec_it)->getQueryStop(), qrstop);
                    ++n;
                }
            }
        } while (++rec_it != recordset.end());

        // measure number of alignment calculations in each of the three passes
        uint gcounter = 0;
        uint pass_0_counter = 0;
        uint pass_1_counter = 0;
        uint pass_2_counter = 0;

        // data storage  TODO: use of C-style arrays
        const seqan::Dna5String qrseq = query_sequences_.getSequence(qid, qrstart, qrstop);
        const std::string qrseqname = boost::str(query_seqname_fmt % qrstart % qrstop % qid % qlength);
        std::vector< AlignmentRecordTaxonomy* > records_ordered;
        std::vector< seqan::Dna5String > rrseqs_ordered; //TODO: boost ptr_container/seqan::StringSet/set to detect equal sequences
        std::vector< std::string > rrseqs_ordered_names;
        std::vector< int > rrseqs_qscores;
        std::vector< large_unsigned_int > rrseqs_matches;
        records_ordered.reserve(n);
        rrseqs_ordered.reserve(n);
        rrseqs_matches.reserve(n);
        rrseqs_qscores.reserve(n);


        {   // retrieve segment sequences
            measure_sequence_retrieval_.start();
            rec_it = recordset.begin();
            while (rec_it != recordset.end()) {
                if (! (*rec_it)->isFiltered()) {
                    records_ordered.push_back(*rec_it);
                    AlignmentRecordTaxonomy& rec = **rec_it;
                    const std::string& rid = rec.getReferenceIdentifier();
                    large_unsigned_int rstart = rec.getReferenceStart();
                    large_unsigned_int rstop = rec.getReferenceStop();
                    large_unsigned_int left_ext = rec.getQueryStart() - qrstart;
                    large_unsigned_int right_ext = qrstop - rec.getQueryStop();

                    // cut out reference region
                    if(rstart <= rstop) {
                        large_unsigned_int start = left_ext < rstart ? rstart - left_ext : 1;
                        large_unsigned_int stop = rstop + right_ext;
                        rrseqs_ordered.push_back(db_sequences_.getSequence(rid, start, stop)); //TODO: avoid copying
                        rrseqs_ordered_names.push_back(boost::str(db_seqname_fmt % start % stop % rid % rec.getReferenceNode()->data->taxid));
                    } else {
                        large_unsigned_int start = right_ext < rstop ? rstop - right_ext : 1;
                        large_unsigned_int stop = rstart + left_ext;
                        rrseqs_ordered.push_back(db_sequences_.getSequenceReverseComplement(rid, start, stop)); //TODO: avoid copying
                        rrseqs_ordered_names.push_back(boost::str(db_seqname_fmt % stop % start % rid % rec.getReferenceNode()->data->taxid));
                    }
                }
                ++rec_it;
            }
            measure_sequence_retrieval_.stop();
        }

        // logging
        logsink << "ID\t" << qrseqname << std::endl;
        logsink << "  NUMREF\t" << n << std::endl << std::endl;

        std::set< uint > qgroup;
        large_unsigned_int anchors_support = 0;
        const TaxonNode* rtax = NULL;  // taxon of closest evolutionary neighbor(s)

        {   // pass 0 (re-alignment to most similar reference segments)
            logsink << "  PASS\t0" << std::endl;
            measure_pass_0_alignment_.start();
            float dbalignment_score_threshold = reeval_bandwidth_factor_*qmaxscore;
            uint index_best = 0;

            for (uint i = 0; i < n; ++i) { //calculate scores for best-scoring references
                if (records_ordered[i]->getScore() >= dbalignment_score_threshold) { //TODO: avoid equal sequences using local alignment score/pid
                    qgroup.insert(i); //TODO: move to if-statement
                    const int score = -seqan::globalAlignmentScore(rrseqs_ordered[i], qrseq, seqan::MyersBitVector());
                    ++pass_0_counter;
                    const large_unsigned_int matches = std::max(static_cast<large_unsigned_int>(std::max(seqan::length(rrseqs_ordered[i]), seqan::length(qrseq)) - score), records_ordered[i]->getIdentities());
                    logsink << "    ALN " << i << " <=> query\tscore = " << score << "; matches = " << matches << std::endl;
                    rrseqs_qscores.push_back(score);
                    rrseqs_matches.push_back(matches);
                    anchors_support = std::max(anchors_support, matches);

                    if (score < rrseqs_qscores[index_best] || (score == rrseqs_qscores[index_best] && matches > rrseqs_matches[index_best])) {
                        index_best = i;
                    }
                } else {  // not similar -> fill in some dummy values
                    rrseqs_qscores.push_back(std::numeric_limits< int >::max());
                    rrseqs_matches.push_back(0);
                }
            }

            // only keep and use the best-scoring reference sequences
            std::list<const TaxonNode*> rtaxa_list;
            for (std::set< uint >::iterator it = qgroup.begin(); it != qgroup.end();) {
                if (rrseqs_qscores[*it] != rrseqs_qscores[index_best] || rrseqs_matches[*it] != rrseqs_matches[index_best]) qgroup.erase(it++);
                else {
                    rtaxa_list.push_back(records_ordered[*it]->getReferenceNode());
                    ++it;
                }
            }
            rtax = this->taxinter_.getLCA(rtaxa_list);
            measure_pass_0_alignment_.stop();
            logsink << "    NUMALN\t" << pass_0_counter << std::endl << std::endl;
        }

        std::vector< const TaxonNode* > anchors_lnode;      // taxa closer to best ref than query
        std::vector< const TaxonNode* > anchors_unode;      // outgroup-related taxa
        float anchors_taxsig = 1.;                          // a measure of tree-like scores  
        float anchors_ival = 0.;
        const TaxonNode* lnode;
        const TaxonNode* unode;
        std::map<uint, small_unsigned_int> outgroup;                 // outgroup sequences

        {   // pass 1 (best reference alignment)
            logsink << "  PASS\t1" << std::endl;
            const std::size_t num_qnodes = qgroup.size();
            anchors_lnode.reserve(num_qnodes);  //TODO: performance check
            anchors_unode.reserve(num_qnodes + 1);  // minimum

            uint alignments_counter = 0;
            uint alignments_counter_naive = 0;
            small_unsigned_int lca_root_dist_min = std::numeric_limits<small_unsigned_int>::max();
            do {  // determine query taxon range
                const uint index_anchor = *qgroup.begin();
                qgroup.erase(qgroup.begin());
                const int qscore = rrseqs_qscores[index_anchor];
                const TaxonNode* rnode = records_ordered[index_anchor]->getReferenceNode();
                lnode = rnode;
                unode = NULL;
                int lscore(-1.);
                int uscore(-1.);

                std::list< boost::tuple< uint, int > > outgroup_tmp;

                // align all others <=> anchor TODO: adaptive cut-off
                measure_pass_1_alignment_.start();
                int min_upper_score = std::numeric_limits< int >::max();
                logsink << "      query: (" << qscore << ") unknown" << std::endl;
                for(uint i = 0; i < n; ++i) {
                    const TaxonNode* cnode = records_ordered[i]->getReferenceNode();
                    int score, matches;
                    if (i == index_anchor) score = 0;
                    else {
                        // use triangle relation to avoid alignment
                        if (rrseqs_qscores[i] == 0 && rrseqs_qscores[index_anchor] == 0 ) { //&& rrseqs_matches[i]) {
                            score = rrseqs_qscores[i];
                            matches = rrseqs_matches[i];
                        }
                        else {
                            score = -seqan::globalAlignmentScore(rrseqs_ordered[ i ], rrseqs_ordered[ index_anchor ], seqan::MyersBitVector());
                            ++pass_1_counter;
                            ++alignments_counter;
                            matches = std::max(seqan::length(rrseqs_ordered[ i ]), seqan::length(rrseqs_ordered[ index_anchor ])) - score;
                            logsink << "    ALN " << i << " <=> " << index_anchor << "\tscore = " << score << "; matches = " << matches << std::endl;
                            
                            // update query alignment scores using triangle relation
                            if (rrseqs_qscores[index_anchor] == 0 && rrseqs_matches[i]) {  // TODO: why rrseqs_matches[i]? 
                                rrseqs_qscores[i] = score;
                                rrseqs_matches[i] = matches;
                            }
                        }
                        ++alignments_counter_naive;
                    }

                    // place sequence
                    if(score <= qscore) {
                        lnode = this->taxinter_.getLCA(lnode, cnode);
                        if(score > lscore) lscore = score;
                        logsink << "      current lower node: " << "("<< score <<") "<<lnode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
                        if(lnode == this->taxinter_.getRoot()) {
                            unode = this->taxinter_.getRoot();
                            break;
                        }
                    }

                    else {
                        if(score == min_upper_score) {
                            unode = this->taxinter_.getLCA(cnode, this->taxinter_.getLCA(lnode, unode));
                            logsink << "      current upper node: " << "("<< score <<") "<< unode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
                        }
                        else if (score < min_upper_score) {
                            uscore = score;
                            min_upper_score = score;
                            unode = this->taxinter_.getLCA(cnode, lnode);
                            outgroup_tmp.clear();
                            logsink << "      current upper node: " << "("<< score <<") "<< unode->data->annotation->name << " (* " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
                        }
                        outgroup_tmp.push_back(boost::make_tuple(i,score));
                    }

                    if (score == 0) qgroup.erase(i); //remove this from list of qnodes because it is sequence-identical
                        // TODO: consider tax. signal calculation
                }

                float ival = .0;                                 //  TODO: initialize
                if(!unode) {
                    unode = this->taxinter_.getRoot();
                    uscore = -1;
                    ival = 1.;
                } else if(unode != lnode) {
                    unode = this->taxinter_.getLCA(lnode, unode);
                    if (lscore < qscore) ival = (qscore - lscore)/static_cast<float>(uscore - lscore);
                }

                measure_pass_1_alignment_.stop();

                logsink << "    SCORE\tlscore = " << lscore << "; uscore = " << uscore << "; qscore = " << qscore << "; ival = " << ival  << std::endl;
                const float taxsig = .0;  // TODO: placer.getTaxSignal(qscore);

                anchors_ival = std::max(ival, anchors_ival);  // combine interpolation values conservatively
                anchors_taxsig = std::min(taxsig, anchors_taxsig);  // combine taxonomic signal values conservatively
                anchors_lnode.push_back(lnode);
                anchors_unode.push_back(unode);
                
                // push elements from temporary to outgroup set
                for(std::list< boost::tuple<uint,int> >::iterator it = outgroup_tmp.begin(); it != outgroup_tmp.end(); ++it) {
                    if(it->get<1>() == min_upper_score) {
                        const small_unsigned_int lca_root_dist = this->taxinter_.getLCA(records_ordered[it->get<0>()]->getReferenceNode(), rnode)->data->root_pathlength;
                        bool smaller = lca_root_dist < lca_root_dist_min;
                        if(smaller) lca_root_dist_min = lca_root_dist;
                        std::map<uint, small_unsigned_int>::iterator find_it = outgroup.find(it->get<0>());
                        if(find_it != outgroup.end()) { if(smaller) find_it->second = lca_root_dist; }
                        else outgroup.insert(std::make_pair(it->get<0>(), lca_root_dist));
                    }
                }
                logsink << std::endl;
            } while (! qgroup.empty() && lnode != this->taxinter_.getRoot());  // TODO: use lca_allnodes, not root

            logsink << "    NUMALN\t" << alignments_counter << tab << alignments_counter_naive - alignments_counter << std::endl;
        
            { // retain only taxa in outgroup which are closest to root
                std::size_t tmp_size = outgroup.size();
                for (std::map<uint, small_unsigned_int>::iterator it = outgroup.begin(); it !=  outgroup.end();) {
                    if(it->second > lca_root_dist_min) outgroup.erase(it++);  // TODO: check if this works
                    else ++it;
                }
                logsink << "    NUMOUTGRP\t" << outgroup.size() << tab << tmp_size - outgroup.size() << std::endl;
            }
            
        anchors_unode.push_back(unode);
        }

        // avoid unneccessary computation if root node is reached TODO: where?
        unode = this->taxinter_.getLCA(anchors_unode);
        anchors_unode.clear();  // TODO: remove; code simplification
        logsink << "    RANGE\t" << lnode->data->annotation->name << tab << unode->data->annotation->name << std::endl << std::endl;
        const TaxonNode* lca_allnodes = this->taxinter_.getRoot();  // TODO: set this to lower for optimization

        {   // pass 2 (stable upper node estimation alignment)
            logsink << "  PASS\t2" << std::endl;
            uint alignments_counter = 0;
            uint alignments_counter_naive = 0;
            while ((! outgroup.empty()) && unode != this->taxinter_.getRoot()) {
                std::list<int> reevaluate;
                BandFactor bandfactor(this->taxinter_, n);
                const uint index_anchor = outgroup.begin()->first;
                outgroup.erase(outgroup.begin());
//                 std::cerr << "after erase: " << outgroup.size() << std::endl;
//                 unode = records_ordered[index_anchor]->getReferenceNode();
//                 lnode = NULL;
                bandfactor.addSequence(0, records_ordered[index_anchor]->getReferenceNode());

                if (!rrseqs_matches[index_anchor]) { //need to align query <=> anchor
                    int score = -seqan::globalAlignmentScore(rrseqs_ordered[index_anchor], qrseq, seqan::MyersBitVector());
                    int matches = std::max(seqan::length(rrseqs_ordered[index_anchor]), seqan::length(qrseq)) - score;
                    logsink << "    ALN query <=> " << index_anchor << "\tscore = " << score << "; matches = " << matches << std::endl;
                    rrseqs_qscores[index_anchor] = score;
                    ++pass_2_counter;
                    ++alignments_counter;
                    ++alignments_counter_naive;
                    rrseqs_matches[index_anchor] = matches;
                }

                const int qscore = rrseqs_qscores[index_anchor];
                const TaxonNode* rnode = records_ordered[index_anchor]->getReferenceNode();
                logsink << "      query: (" << qscore << ") unknown" << std::endl;

                // align all others <=> anchor TODO: heuristic
                measure_pass_2_alignment_.start();

                for (uint i = 0; i < n; ++i) {
                    const TaxonNode* cnode = records_ordered[i]->getReferenceNode();
                    int score;
                    if (i == index_anchor) score = 0;
                    else {
                        score = -seqan::globalAlignmentScore(rrseqs_ordered[ i ], rrseqs_ordered[ index_anchor ], seqan::MyersBitVector());
                        logsink << "    ALN " << i << " <=> " << index_anchor << "\tscore = " << score <<  std::endl;
                        ++pass_2_counter;
                        rrseqs_qscores[i] = score;
                        ++alignments_counter;
                        ++alignments_counter_naive;
                    }

                    if (score == 0) outgroup.erase(i);
                    
                    bandfactor.addSequence(score, records_ordered[i]->getReferenceNode());

                    if(score <= qscore) {
                        unode = this->taxinter_.getLCA(unode, records_ordered[i]->getReferenceNode());
                        logsink << "      current upper node: " << "("<< score <<") "<< unode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
                    }
                    else reevaluate.push_back(i);
                }

                {
                    int qscore_ex = qscore * bandfactor.getFactor();
                    logsink << std::endl << "    EXT\tqscore = " << qscore << "; threshold = " << qscore_ex << std::endl;
                    std::list< int >::iterator it = reevaluate.begin();
                    while(it != reevaluate.end() && unode != lca_allnodes ) {
                        const int score = rrseqs_qscores[*it];
                        const TaxonNode* cnode = records_ordered[*it]->getReferenceNode();
                        if(score <= qscore_ex) {
                            unode = this->taxinter_.getLCA(cnode, unode);
                            logsink << "      current upper node: " << "("<< score <<") "<< unode->data->annotation->name << " (+ " << cnode->data->annotation->name << " at " << static_cast<int>(this->taxinter_.getLCA(cnode, rnode)->data->root_pathlength) << " )" << std::endl;
                        }
                        ++it;
                    }
                }

                logsink << std::endl;
                
                measure_pass_2_alignment_.stop();
                measure_placement_algorithm_.start();
                anchors_unode.push_back(unode);
                measure_placement_algorithm_.stop();
            }
            
            logsink << "    NUMALN\t" << alignments_counter << tab << alignments_counter_naive - alignments_counter << std::endl;
        }

        lnode = this->taxinter_.getLCA(anchors_lnode);
        if(anchors_unode.empty()) unode = this->taxinter_.getRoot();  // root if no outgroup element
        else unode = this->taxinter_.getLCA(anchors_unode);
        
        if(unode == lnode) anchors_ival = 1.;
        
        logsink << "    RANGE\t" << lnode->data->annotation->name << tab << unode->data->annotation->name << std::endl << std::endl;

        prec.setSignalStrength(anchors_taxsig);
        prec.setQueryFeatureBegin(qrstart);
        prec.setQueryFeatureEnd(qrstop);
        prec.setInterpolationValue(anchors_ival);
        prec.setNodeRange(lnode, unode, anchors_support);
        prec.setRtax(rtax);
        gcounter = pass_0_counter + pass_1_counter + pass_2_counter;
        float normalised_rt = (float)gcounter/(float)n;
        logsink << "STATS \"" << qrseqname << "\"\t" << n << "\t" << pass_0_counter << "\t"<< pass_1_counter << "\t" << pass_2_counter << "\t" << gcounter << "\t" << std::setprecision(2) << std::fixed <<normalised_rt << std::endl << std::endl;
    }

protected:
    QStorType& query_sequences_;
    const DBStorType& db_sequences_;
    SortByBitscoreFilter< ContainerT > sort_;
    compareTupleFirstLT< boost::tuple< int, uint >, 0 > tuple_1_cmp_le_;

private:
    const float exclude_alignments_factor_;
    const float reeval_bandwidth_factor_;
    StopWatchCPUTime measure_sequence_retrieval_;
    StopWatchCPUTime measure_pass_0_alignment_;
    StopWatchCPUTime measure_pass_1_alignment_;
    StopWatchCPUTime measure_pass_2_alignment_;
    StopWatchCPUTime measure_placement_algorithm_;
};

#endif // taxonpredictionmodelsequence_hh_
