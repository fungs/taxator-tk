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
class BandFactor{
    public:
        BandFactor( TaxonomyInterface& taxinter, uint reserve = 0 ) :
                bandfactor_(-1),
                taxinter_( taxinter ){
            if( reserve ) data_.reserve(reserve);
        }

        void addSequence( const int score, const TaxonNode* node ) {
			data_.push_back( boost::make_tuple( score, node ) );
		}

		float getFactor(){
            if(bandfactor_ < 0){
                 sort();
                 setBandFactor();
            }
            return bandfactor_;
		}

    private:
        void setBandFactor( const float min_bandfactor = 1., const float max_bandfactor = std::numeric_limits<  int >::max() ) { //data_ must be sorted TODO: optimize
			float bandfactor = min_bandfactor;
			int score;
			const TaxonNode* anchor;
			const TaxonNode* node;
			std::map< small_unsigned_int, int > worstscore_per_rank;
			boost::tie( score, anchor ) = data_[0];
			small_unsigned_int last_rank = anchor->data->root_pathlength;
			worstscore_per_rank[ last_rank ] = score; //TODO: change name

			for ( uint i = 1; i < data_.size(); ++i ) {
				boost::tie( score, node ) = data_[i];
				node = taxinter_.getLCA( node, anchor );
				small_unsigned_int rank = node->data->root_pathlength;
				if ( rank == last_rank ) {
				} else if ( rank < last_rank ) {
					worstscore_per_rank[rank] = score;
					last_rank = rank;
				} else { //disorder
					int refscore;
					small_unsigned_int r = rank - 1;
					do {
						std::map< small_unsigned_int, int >::iterator it = worstscore_per_rank.find( r );
						if ( it != worstscore_per_rank.end() ) {
							refscore = it->second;
							if ( refscore ) bandfactor = std::max( bandfactor, float( score )/float( refscore ) );
						}
					} while ( r-- );
				}
			}
			bandfactor_ = std::min( bandfactor, max_bandfactor );
		}

		void sort() { //sort by increasing order
			data_type_::iterator start_it = data_.begin();
			std::sort( ++start_it, data_.end(), comparator_ );
		}

        compareTupleFirstLT< boost::tuple< int, const TaxonNode* >, 0 > comparator_;
        typedef std::vector< boost::tuple< int, const TaxonNode* > > data_type_;
        float bandfactor_;
        data_type_ data_;
        TaxonomyInterface taxinter_;

};

// TODO: timers should be thread safe!
template< typename ContainerT, typename QStorType, typename DBStorType >
class DoubleAnchorRPAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		DoubleAnchorRPAPredictionModel( const Taxonomy* tax, QStorType& q_storage, const DBStorType& db_storage, float filterout ,float reeval_bandwidth = .1 ) :
					TaxonPredictionModel< ContainerT >( tax ),
					query_sequences_( q_storage ),
					db_sequences_( db_storage ),
					filterout_(filterout),
					reeval_bandwidth_factor_( 1. - reeval_bandwidth ), //TODO: check range
					measure_placement_algorithm_( "placement algorithm" ),
					measure_phase2_alignment_( "best reference anchor alignments (step 2)" ),
					measure_phase3_alignment_( "distant anchor alignments (step 3)" ),
					measure_reeval_alignment_( "best reference re-evaluation alignments (step 1)" ),
					measure_sequence_retrieval_( "sequence retrieval using index" )
					{};

		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec ); //set name and length
			// cannot be part of class because it is not thread-safe
			boost::format db_seqname_fmt( "%d:%d@%s tax=%s" );
			boost::format query_seqname_fmt( "%d:%d@%s len=%d" );

			typename ContainerT::iterator rec_it = firstUnmaskedIter( recordset );
            //typename ContainerT::iterator rec_it2 = rec_it;

			if ( rec_it == recordset.end() ) {
				TaxonPredictionModel< ContainerT >::setUnclassified( prec );
				return;
			}
			// set some shortcuts
			const std::string& qid = (*rec_it)->getQueryIdentifier();
			large_unsigned_int qrstart = (*rec_it)->getQueryStart();
			large_unsigned_int qrstop = (*rec_it)->getQueryStop();
			large_unsigned_int qlength = (*rec_it)->getQueryLength();
			float qprevscore = (*rec_it)->getScore();
            float qmaxscore = .0;
            //float factor = filterout_; // good 0.25


			// determine position range of query to consider
			uint n = 0;
			do {
                  if(!(*rec_it)->isFiltered()) qmaxscore = std::max( qmaxscore, (*rec_it)->getScore());
            } while(++rec_it != recordset.end());

            rec_it = firstUnmaskedIter( recordset );

			assert( qrstart <= qrstop );

			do {
                assert( (*rec_it)->getQueryStart() <= (*rec_it)->getQueryStop() );

                if((*rec_it)->getScore()<qmaxscore*filterout_){(*rec_it)->filterOut();}

				if ( ! (*rec_it)->isFiltered() ) {
					qrstart = std::min( (*rec_it)->getQueryStart(), qrstart );
					qrstop = std::max( (*rec_it)->getQueryStop(), qrstop );
	 				qprevscore = std::max( qprevscore, (*rec_it)->getScore() );
					++n;
				}
			} while ( ++rec_it != recordset.end() );

			//measure
			uint gcounter = 0;
			uint phase1_counter = 0;
			uint phase2_counter = 0;
			uint phase3_counter = 0;

			// lots of objects TODO: make use of C-style arrays
// 			int extra_mismatches = qlength - qrstop + qrstart - 1;
			const seqan::Dna5String qrseq = query_sequences_.getSequence( qid, qrstart, qrstop );
			const std::string qrseqname = boost::str( query_seqname_fmt % qrstart % qrstop % qid % qlength );
			std::vector< AlignmentRecordTaxonomy* > records_ordered;
			std::vector< seqan::Dna5String > rrseqs_ordered; //TODO: boost ptr_container/seqan::StringSet/set to detect equal sequences
			std::vector< std::string > rrseqs_ordered_names;
			std::vector< int > rrseqs_qscores;
			std::vector< large_unsigned_int > rrseqs_matches;
			records_ordered.reserve( n );
			rrseqs_ordered.reserve( n );
			rrseqs_matches.reserve( n );
			rrseqs_qscores.reserve( n );


			{ // retrieving sequence blocks
				measure_sequence_retrieval_.start();
				rec_it = recordset.begin();
				while ( rec_it != recordset.end() ) {
					if ( ! (*rec_it)->isFiltered() ) {
						records_ordered.push_back( *rec_it );
						AlignmentRecordTaxonomy& rec = **rec_it;
						const std::string& rid = rec.getReferenceIdentifier();
						large_unsigned_int rstart = rec.getReferenceStart();
						large_unsigned_int rstop = rec.getReferenceStop();
						large_unsigned_int left_ext = rec.getQueryStart() - qrstart;
						large_unsigned_int right_ext = qrstop - rec.getQueryStop();

						// cut out reference region
						if( rstart <= rstop ) {
							large_unsigned_int start = left_ext < rstart ? rstart - left_ext : 1;
							large_unsigned_int stop = rstop + right_ext;
							rrseqs_ordered.push_back( db_sequences_.getSequence( rid, start, stop ) ); //TODO: avoid copying
							rrseqs_ordered_names.push_back( boost::str( db_seqname_fmt % start % stop % rid % rec.getReferenceNode()->data->taxid ) );
						} else {
							large_unsigned_int start = right_ext < rstop ? rstop - right_ext : 1;
							large_unsigned_int stop = rstart + left_ext;
							rrseqs_ordered.push_back( db_sequences_.getSequenceReverseComplement( rid, start, stop ) ); //TODO: avoid copying
							rrseqs_ordered_names.push_back( boost::str( db_seqname_fmt % stop % start % rid % rec.getReferenceNode()->data->taxid ) );
						}
					}
					++rec_it;
				}
				measure_sequence_retrieval_.stop();
			}

			{ // logging
				logsink << "id\t" << qrseqname << std::endl;
				logsink << "# references:\t" << n << std::endl;
			}

			std::set< uint > anchor_indices_p1;
			large_unsigned_int anchors_support = 0;
			const TaxonNode* rtax = records_ordered[0]->getReferenceNode();

			{ // phase 1 (DB re-alignment within band)
				measure_reeval_alignment_.start();
				float dbalignment_score_threshold = reeval_bandwidth_factor_*qprevscore;
				uint index_best = 0;
				uint counter = 0;

				for ( uint i = 0; i < n; ++i ) { //calculate scores for best-scoring references
					if ( records_ordered[i]->getScore() >= dbalignment_score_threshold ) { //TODO: avoid equal sequences using local alignment score/pid
						anchor_indices_p1.insert( i ); //TODO: move to if-statement
						const int score = -seqan::globalAlignmentScore( rrseqs_ordered[i], qrseq, seqan::MyersBitVector() );
						++phase1_counter;
						const large_unsigned_int matches = std::max( static_cast<large_unsigned_int>( std::max( seqan::length( rrseqs_ordered[i] ), seqan::length( qrseq ) ) - score ), records_ordered[i]->getIdentities() );
						rrseqs_qscores.push_back( score );
						rrseqs_matches.push_back( matches );
						anchors_support = std::max( anchors_support, matches );

						if ( score < rrseqs_qscores[index_best] || ( score == rrseqs_qscores[index_best] && matches > rrseqs_matches[index_best] ) ) {
							index_best = i;
						}
						++counter;
					} else { //fill in some dummy values
						rrseqs_qscores.push_back( std::numeric_limits< int >::max() );
						rrseqs_matches.push_back( 0 );
					}
				}

				for ( std::set< uint >::iterator it = anchor_indices_p1.begin(); it != anchor_indices_p1.end(); ) { //reduce number of anchors
					if ( rrseqs_qscores[*it] != rrseqs_qscores[index_best] || rrseqs_matches[*it] != rrseqs_matches[index_best] ) anchor_indices_p1.erase( it++ );
					else{
					    rtax = this->taxinter_.getLCA(rtax,records_ordered[*it]->getReferenceNode());
                        ++it;
					}
				}
				measure_reeval_alignment_.stop();
				logsink << "# alignments step 1\t" << counter << std::endl;
			}

			std::vector< const TaxonNode* > anchors_lnode;
			std::vector< const TaxonNode* > anchors_unode;
			float anchors_taxsig = 1;
			float anchors_ival = 0.;
			const TaxonNode* lnode;
			const TaxonNode* unode;
			std::set< uint > anchor_indices_p2;

			{ //phase 2 (best reference alignment)
				const std::size_t anchornum_p1 = anchor_indices_p1.size();
				anchors_lnode.reserve( anchornum_p1 ); //TODO: performance check
				anchors_unode.reserve( anchornum_p1 + 1 ); //minimum

				uint cycle_count = 0;
				uint naive_cycle_count = anchornum_p1;
				do { // place using anchor
					++cycle_count;

					const uint index_anchor = *anchor_indices_p1.begin();
					const int qscore = rrseqs_qscores[index_anchor];
                    lnode = records_ordered[index_anchor]->getReferenceNode();
                    unode = NULL;
                    int lscore(-1.);
                    int uscore(-1.);

					std::list< boost::tuple< uint, int > > anchor_indices_p2_tmp;

					// align all others against anchor TODO: use cutoff for speedup
					measure_phase2_alignment_.start();
					int min_upper_score = std::numeric_limits< int >::max();
					for( uint i = 0; i < n; ++i ) {
						int score, matches;
						if ( i == index_anchor ) score = 0;
						else {
							if ( rrseqs_qscores[index_anchor] == 0 && rrseqs_matches[i] ) { //use triangle relation to avoid alignment
								score = rrseqs_qscores[i];
								matches = rrseqs_matches[i];
							}
							else {
								score = -seqan::globalAlignmentScore( rrseqs_ordered[ i ], rrseqs_ordered[ index_anchor ], seqan::MyersBitVector() ); //- extra_mismatches;
								++phase2_counter;
								matches = std::max( seqan::length( rrseqs_ordered[ i ] ), seqan::length( rrseqs_ordered[ index_anchor ] ) ) - score;
								if ( rrseqs_qscores[index_anchor] == 0 && rrseqs_matches[i] ) { //update using triangle relation
									rrseqs_qscores[i] = score;
									rrseqs_matches[i] = matches;
								}
							}
						}

						//place sequence
						if(score<=qscore){
						    lnode = this->taxinter_.getLCA(lnode,records_ordered[i]->getReferenceNode());
                            if(score > lscore) lscore = score;
                            logsink << "current lowernode : " << " ("<< score <<") "<<lnode->data->annotation->name << " " << records_ordered[i]->getReferenceNode()->data->annotation->name << "\n";
						    if(lnode == this->taxinter_.getRoot()){
                                unode = this->taxinter_.getRoot();
                                break;
                                }
                        }

                        else if(score > qscore && score<=min_upper_score){
                            if(score == min_upper_score) unode = this->taxinter_.getLCA(records_ordered[i]->getReferenceNode(),unode);
                            else{
                                uscore = score;
                                min_upper_score = score;
                                unode = this->taxinter_.getLCA(records_ordered[i]->getReferenceNode(),lnode);
                                anchor_indices_p2_tmp.clear();
                                }
                            anchor_indices_p2_tmp.push_back(boost::make_tuple(i,score));
                            logsink << "current uppernode : " << " ("<< score <<") "<< unode->data->annotation->name << " " << records_ordered[i]->getReferenceNode()->data->annotation->name << "\n";
                        }

                        if ( score == 0 ) anchor_indices_p1.erase( i ); //remove all identical sequences TODO: consider tax. signal calculation
					}

                    float ival ;

                    if(!unode){
                        unode = this->taxinter_.getRoot();
                        uscore=.0;
                        ival = 1.;
                    }else{unode = this->taxinter_.getLCA(lnode,unode);
                         if ( lscore >= qscore ) ival = .0;
                         else ival = (qscore - lscore)/static_cast<float>( uscore - lscore );
                    }

                    uint min_root_distance = std::numeric_limits< int >::max();

                    for ( std::list< boost::tuple< uint, int > >::iterator it =  anchor_indices_p2_tmp.begin(); it !=  anchor_indices_p2_tmp.end(); ) {
                        uint LCAroot = this->taxinter_.getLCA(records_ordered[it->get<0>()]->getReferenceNode(),lnode)->data->root_pathlength;
                        if( LCAroot < min_root_distance){ min_root_distance = LCAroot;}
                        ++it;
                    }

                    for ( std::list< boost::tuple< uint, int > >::iterator it =  anchor_indices_p2_tmp.begin(); it !=  anchor_indices_p2_tmp.end(); ) {
                        if(this->taxinter_.getLCA(records_ordered[it->get<0>()]->getReferenceNode(),lnode)->data->root_pathlength == min_root_distance){
                        anchor_indices_p2.insert( it->get<0>() );
                        }
                        ++it;
                    }

                    measure_phase2_alignment_.stop();

					logsink << "#ival: " << ival << "low: " << lscore << "| up: " << uscore << "\n";
					const float taxsig = .0;//placer.getTaxSignal( qscore );

					anchors_ival = std::max( ival, anchors_ival );
					anchors_taxsig = std::min( taxsig, anchors_taxsig );
					anchors_lnode.push_back( lnode );
					anchors_unode.push_back( unode );


				} while ( ! anchor_indices_p1.empty() && lnode != this->taxinter_.getRoot() );

				logsink << "# step 2 alignment runs\t" << cycle_count << std::endl;
				if ( naive_cycle_count ) logsink << "# step 2 saved alignment runs\t" << naive_cycle_count - cycle_count << std::endl;
			}

			// avoid unneccessary computation if root node is reached
			lnode = this->taxinter_.getLCA( anchors_unode );
			anchors_unode.clear();
			anchors_unode.push_back( lnode );

			{ //phase 3 (stable upper node estimation alignment)
			    logsink << "phase3\n";
				uint cycle_count = 0;
				uint naive_cycle_count = anchor_indices_p2.size();
				while ( (! anchor_indices_p2.empty()) && lnode != this->taxinter_.getRoot() ) {
					++cycle_count;
                    std::list<int> reevaluate;
					BandFactor bandfactor( this->taxinter_, n );
					const uint index_anchor = *anchor_indices_p2.begin();
					unode = records_ordered[index_anchor]->getReferenceNode();
					lnode = NULL;
					bandfactor.addSequence( 0, unode );

					if ( !rrseqs_matches[index_anchor] ) { //need to align query against anchor
						rrseqs_qscores[index_anchor] = -seqan::globalAlignmentScore( rrseqs_ordered[index_anchor], qrseq, seqan::MyersBitVector() ); //- extra_mismatches;;
						++phase3_counter;
						rrseqs_matches[index_anchor] = std::max( seqan::length( rrseqs_ordered[index_anchor] ), seqan::length( qrseq ) ) - rrseqs_qscores[index_anchor];
					}

					const int qscore = rrseqs_qscores[index_anchor];

					// align all others against anchor TODO: heuristic
					measure_phase3_alignment_.start();

					for ( uint i = 0; i < n; ++i ) {
						int score;
						if ( i == index_anchor ) score = 0;
						else {
							score = -seqan::globalAlignmentScore( rrseqs_ordered[ i ], rrseqs_ordered[ index_anchor ], seqan::MyersBitVector() ); //- extra_mismatches;
							++phase3_counter;
							rrseqs_qscores[i] = score;

						}

						if ( score == 0 ) anchor_indices_p2.erase( i );
						bandfactor.addSequence( score, records_ordered[i]->getReferenceNode() );

						if(score<=qscore)unode = this->taxinter_.getLCA(unode,records_ordered[i]->getReferenceNode());
                        else reevaluate.push_back(i);
					}

					int qscore_ex = qscore * bandfactor.getFactor();

                    for(std::list< int >::iterator it = reevaluate.begin(); it != reevaluate.end(); ){
                        if(rrseqs_qscores[*it] <= qscore_ex) unode = this->taxinter_.getLCA(records_ordered[*it]->getReferenceNode(),unode);
                        it++;
                    }

                    logsink<< "current uppernode " << unode->data->annotation->name << "\n";
					measure_phase3_alignment_.stop();
					measure_placement_algorithm_.start();
					anchors_unode.push_back( unode );
					measure_placement_algorithm_.stop();
				}
				logsink << "# step 3 alignment runs\t" << cycle_count << std::endl;
				if ( naive_cycle_count ) logsink << "# step 3 saved alignment runs\t" << naive_cycle_count - cycle_count << std::endl;
			}

			lnode = this->taxinter_.getLCA( anchors_lnode );
			unode = this->taxinter_.getLCA( anchors_unode );

			prec.setSignalStrength( anchors_taxsig );
			prec.setQueryFeatureBegin( qrstart );
			prec.setQueryFeatureEnd( qrstop );
			prec.setInterpolationValue( anchors_ival );
			prec.setNodeRange( lnode, unode, anchors_support );
			prec.setRtax( rtax );
			gcounter = phase1_counter + phase2_counter + phase3_counter;
			float normalised_rt = (float)gcounter/(float)n;
			logsink << "# ratio;" << "id: " << qrseqname << "\t" << n << "\t" << phase1_counter << "\t"<< phase2_counter << "\t" << phase3_counter << "\t" << gcounter << "\t" << std::setprecision(2) << std::fixed <<normalised_rt << std::endl;
 			logsink << "prediction:\t" << unode->data->annotation->name << " (" << static_cast<uint>( unode->data->root_pathlength ) << ")" << std::endl << std::endl;
		}

	protected:
		QStorType& query_sequences_;
		const DBStorType& db_sequences_;
		SortByBitscoreFilter< ContainerT > sort_;
		compareTupleFirstLT< boost::tuple< int, uint >, 0 > tuple_1_cmp_le_;

		private:
            const float filterout_;
			const float reeval_bandwidth_factor_;
			StopWatchCPUTime measure_placement_algorithm_;
			StopWatchCPUTime measure_phase2_alignment_;
			StopWatchCPUTime measure_phase3_alignment_;
			StopWatchCPUTime measure_reeval_alignment_;
			StopWatchCPUTime measure_sequence_retrieval_;
};

#endif // taxonpredictionmodelsequence_hh_
