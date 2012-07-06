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

#include <locale>

static std::locale loc;
static const std::collate< char >& coll = std::use_facet< std::collate< char > >( loc );

// helper class
template< bool smallerscore = true >
class PathScoreCollection {
	public:
		void add( const TaxonNode* node, int score ) {
			sorted_.insert( boost::make_tuple( node->data->root_pathlength, score, node ) );
		}
		
		void getMajority( const TaxonNode*& node_ret, int& score_ret, small_unsigned_int support = 100 ) { //optimized for percentage > 50
			std::multiset< boost::tuple< small_unsigned_int, int, const TaxonNode* > >::iterator it = sorted_.begin();
			
			if ( support == 100 || sorted_.size() == 1 ) {
				score_ret = it->get<1>();
				node_ret = it->get<2>();
				return;
			}
			
			std::size_t index_stop = sorted_.size()*(100 - support)/100;
			
			// go to index position
			for ( std::size_t i = 0; i < index_stop; ++i ) ++it;
			
			small_unsigned_int l = it->get<0>();
			const boost::tuple< small_unsigned_int, int, const TaxonNode* >* tmp = &*it;
			if ( smallerscore ) for ( std::multiset< boost::tuple< small_unsigned_int, int, const TaxonNode* > >::reverse_iterator it2( it ); it2 != sorted_.rend() && it2->get<0>() == l; ++it2 ) tmp = &*it2;
			else for ( std::multiset< boost::tuple< small_unsigned_int, int, const TaxonNode* > >::iterator it2( it ); it2 != sorted_.end() && it2->get<0>() == l; ++it2 ) tmp = &*it2;
			
			score_ret = tmp->get<1>();
			node_ret = tmp->get<2>();
		}
		
		small_unsigned_int getSupport( const TaxonNode* node ) {
			small_unsigned_int target_dist = node->data->root_pathlength;
			std::size_t count = sorted_.size();
			for ( std::multiset< boost::tuple< small_unsigned_int, int, const TaxonNode* > >::iterator it = sorted_.begin(); it != sorted_.end(); ++it) {
				small_unsigned_int tmp = it->get<0>();
				if ( tmp >= target_dist ) return count*100/sorted_.size();
				--count;
			};
			return 0;
		}
		
		void clear() {
			sorted_.clear();
		}
		
		std::size_t size() { return sorted_.size(); }
		
	private:
		std::multiset< boost::tuple< small_unsigned_int, int, const TaxonNode* > > sorted_;
};



class SandwichPlacerAdaptive {
	public:
		SandwichPlacerAdaptive( TaxonomyInterface& taxinter, uint reserve = 0 ) :
				prepared_( false ),
				taxinter_( taxinter ),
				lsupport_( 100 ),
				usupport_( 100 ), //TODO: see if we can lower value
				lsupport_stage2_( 80 ) { //TODO: study effects
			if ( reserve ) data_.reserve( reserve );
		}
		
		void addSequence( const int score, const TaxonNode* node ) {
			data_.push_back( boost::make_tuple( score, node ) );
			prepared_ = false;
		}
		
		float placeSequence( const int qscore, const TaxonNode*& lower_node, const TaxonNode*& upper_node ) { //TODO: only unique nodes
			assert( data_.size() ); //TODO: exception
			if ( ! prepared_ ) prepare();
			
			PathScoreCollection< false > lgroup;
			int score;
			const TaxonNode* node = NULL;
			
			// put asymmetric band around pairwise score to compensate the optimistic value
			const int qscore_ext = qscore*bandfactor_ + 1;
			
			// determine  lower node
			data_type_::const_iterator it = data_.begin();
			boost::tie( score, node ) = *it;
			assert( ! score ); //must be zero
			int lscore;
			const TaxonNode* lnode = NULL;
			const TaxonNode* anchor = node;
			
			do {
				lgroup.add( node, score );
				if ( ++it == data_.end() ) { //missing upper
					lgroup.getMajority( lnode, lscore, lsupport_ ); //TODO: changed from stage2 to normal
					lower_node = lnode;
					upper_node = taxinter_.getRoot();
					return 1.;
				};
				boost::tie( score, node ) = *it;
				node = taxinter_.getLCA( node, anchor );
			} while ( score <= qscore_ext );
			
			lgroup.getMajority( lnode, lscore, lsupport_ );
			
			// determine upper node
			const int critical_ext = score*bandfactor_ + 1;
			
			PathScoreCollection< false > ugroup;
			do {
				ugroup.add( node, score );
				if ( ++it == data_.end() ) break;
				boost::tie( score, node ) = *it;
				node = taxinter_.getLCA( node, anchor );
			} while ( score <= critical_ext );
			
			int uscore;
			const TaxonNode* unode = NULL;
			ugroup.getMajority( unode, uscore, usupport_ );
			
// 			std::cerr << "upper node is: " << unode->data->annotation->name << std::endl;
			
			// noise detection
			float ival;
			if ( taxinter_.isParentOf( lnode, unode ) ) {
				if ( lgroup.getSupport( unode ) >= lsupport_stage2_ ) {
					lnode = unode;
					lscore = uscore;
					ival = 1.;
				} else {
					unode = lnode;
					uscore = lscore;
					ival = 0.;
				}
			} else {
				if ( lscore >= qscore ) ival = .0;
				else ival = (qscore - lscore)/static_cast<float>( uscore - lscore );
			}
			
			lower_node = lnode;
			upper_node = unode;
			return ival;
		}
		
		float getTaxSignal( int placescore = std::numeric_limits<  int >::max() ) {
			assert( data_.size() ); //TODO: exception
			if ( ! prepared_ ) prepare();
			
			// create ordering of scores according to taxonomy
			typedef std::vector< boost::tuple< int, int > > ordertype;
			std::vector< int > scores, dists;
			scores.reserve( data_.size() );
			dists.reserve( data_.size() );
			
			{ //construction of arrays
				int score;
				const TaxonNode* node;
				const TaxonNode* anchor;
				data_type_::const_iterator it = data_.begin();
				boost::tie( score, anchor ) = *it;
				node = anchor;
				bool uppermode = false;
				while ( true ) {
					if ( score > placescore ) {
						if ( uppermode ) break;
						placescore = score*bandfactor_;
						uppermode = true;
					}
					scores.push_back( score );
					dists.push_back( -(node->data->root_pathlength) ); //negative pathlength for sorting
					if ( ++it == data_.end() ) break; //single exit condition
					boost::tie( score, node ) = *it;
					node = taxinter_.getLCA( node, anchor );
				}
			}
			
			double tau = kendall( scores, dists );
			return static_cast< float >( tau );
		};
		
	private:
		void setBandFactor( const float min_bandfactor = 1., const float max_bandfactor = std::numeric_limits<  int >::max() ) { //data_ must be sorted TODO: optimize
			float bandfactor = min_bandfactor;
			int score;
			const TaxonNode* anchor;
			const TaxonNode* node;
// 			small_unsigned_int deepest_rank = last_rank;
			std::map< small_unsigned_int, int > worstscore_per_rank;
			boost::tie( score, anchor ) = data_[0];
			small_unsigned_int last_rank = anchor->data->root_pathlength;
			worstscore_per_rank[ last_rank ] = score; //TODO: change name
			
			for ( uint i = 1; i < data_.size(); ++i ) {
				boost::tie( score, node ) = data_[i];
				node = taxinter_.getLCA( node, anchor );
				small_unsigned_int rank = node->data->root_pathlength;
				if ( rank == last_rank ) {
// 					std::cerr << "       rank " << static_cast< int >( rank ) << " with score " << score << std::endl;
// 					worstscore_per_rank[rank] = score;
// 					last_rank = rank;
				} else if ( rank < last_rank ) {
// 					std::cerr << "ok for rank " << static_cast< int >( rank ) << " with score " << score << std::endl;
					worstscore_per_rank[rank] = score;
					last_rank = rank;
				} else { //disorder
// 					std::cerr << " false rank " << static_cast< int >( rank ) << " with score " << score << std::endl;
// 					std::cerr << "old rank: " << static_cast< int >( last_rank ) << std::endl;
					int refscore;
					small_unsigned_int r = rank - 1;
					do {
// 						std::cerr << "searching for score with rank " << static_cast< int >( r ) << " ...";
						std::map< small_unsigned_int, int >::iterator it = worstscore_per_rank.find( r );
						if ( it != worstscore_per_rank.end() ) {
							refscore = it->second;
// 							std::cerr << " found with score " << refscore << std::endl;
							if ( refscore ) bandfactor = std::max( bandfactor, float( score )/float( refscore ) );
// 							std::cerr << "growing band to " << bandfactor << std::endl;
						} else std::cerr << std::endl;
					} while ( r-- );
				}
			}
			bandfactor_ = std::min( bandfactor, max_bandfactor );
// 			std::cerr << "final bandfactor is " << bandfactor << std::endl << std::endl;
		}
		
		void sort() { //sort by increasing order
			data_type_::iterator start_it = data_.begin();
			std::sort( ++start_it, data_.end(), comparator_ );
		}
		
		void prepare() {
			sort();
			setBandFactor();
			prepared_ = true;
		}
		
		double kendall( const std::vector<int>& X1, const std::vector<int>& X2 ) { //TODO: acknowledge contribution!
			assert( X1.size() == X2.size() );
			std::size_t len = X1.size();
			int m1 = 0, m2 = 0, s = 0, nPair;
			std::size_t i, j;
			double denominator1, denominator2;

			for(i = 0; i < len; i++) {
				for(j = i + 1; j < len; j++) {
					if(X2[i] > X2[j]) {
						if (X1[i] > X1[j]) {
								s++;
						} else if(X1[i] < X1[j]) {
								s--;
						} else {
								m1++;
						}
					} else if(X2[i] < X2[j]) {
						if (X1[i] > X1[j]) {
								s--;
						} else if(X1[i] < X1[j]) {
								s++;
						} else {
								m1++;
						}
					} else {
						m2++;

						if(X1[i] == X1[j]) {
								m1++;
						}
					}
				}
			}

			nPair = len * (len - 1) / 2;
			denominator1 = nPair - m1;
			denominator2 = nPair - m2;
			return s / sqrt(denominator1) / sqrt(denominator2);
		}
		
		compareTupleFirstLT< boost::tuple< int, const TaxonNode* >, 0 > comparator_;
		typedef std::vector< boost::tuple< int, const TaxonNode* > > data_type_;
		bool prepared_;
		float bandfactor_;
		data_type_ data_;
		TaxonomyInterface taxinter_;
		const small_unsigned_int lsupport_;
		const small_unsigned_int usupport_;
		const small_unsigned_int lsupport_stage2_;
};



class SandwichPlacer {
	public:
		SandwichPlacer( TaxonomyInterface& taxinter, const float bandwidth_factor, uint reserve = 0 ) :
				unsorted_( false ),
				taxinter_( taxinter ),
				bandwidth_factor_( bandwidth_factor ),
				lsupport_( 100 ),
				usupport_( 100 ), //TODO: see if we can lower value
				lsupport_stage2_( 80 ) { //TODO: study effects
			if ( reserve ) data_.reserve( reserve );
		}
		
		void addSequence( const int score, const TaxonNode* node ) {
			data_.push_back( boost::make_tuple( score, node ) );
			unsorted_ = true;
		}
		
		float placeSequence( const int qscore, const TaxonNode*& lower_node, const TaxonNode*& upper_node ) { //TODO: only unique nodes
			assert( data_.size() ); //TODO: exception
			if ( unsorted_ ) sort();
			
			PathScoreCollection< false > lgroup;
			int score;
			const TaxonNode* node = NULL;
			
			// put asymmetric band around pairwise score to compensate the optimistic value
			const int qscore_ext = qscore*bandwidth_factor_ + 1;
			
			// determine  lower node
			data_type_::const_iterator it = data_.begin();
			boost::tie( score, node ) = *it;
			assert( ! score ); //must be zero
			int lscore;
			const TaxonNode* lnode = NULL;
			const TaxonNode* anchor = node;
			
			do {
				lgroup.add( node, score );
				if ( ++it == data_.end() ) { //missing upper
					lgroup.getMajority( lnode, lscore, lsupport_ ); //TODO: changed from stage2 to normal
					//prec.setNodeRange( lnode, taxinter_.getRoot(), support );
					//prec.setInterpolationValue( 1. ); //needed because ival could be close to zero
					//return true;
					lower_node = lnode;
					upper_node = taxinter_.getRoot();
					return 1.;
				};
				boost::tie( score, node ) = *it;
				node = taxinter_.getLCA( node, anchor );
			} while ( score <= qscore_ext );
			
			lgroup.getMajority( lnode, lscore, lsupport_ );
			
// 			std::cerr << "lower node is: " << lnode->data->annotation->name << std::endl;
			
			// determine upper node
			const int critical_ext = score*bandwidth_factor_ + 1;
			
			PathScoreCollection< false > ugroup;
			do {
				ugroup.add( node, score );
				if ( ++it == data_.end() ) break;
				boost::tie( score, node ) = *it;
				node = taxinter_.getLCA( node, anchor );
			} while ( score <= critical_ext );
			
			int uscore;
			const TaxonNode* unode = NULL;
			ugroup.getMajority( unode, uscore, usupport_ );
			
// 			std::cerr << "upper node is: " << unode->data->annotation->name << std::endl;
			
			// noise detection
			float ival;
			if ( taxinter_.isParentOf( lnode, unode ) ) {
				if ( lgroup.getSupport( unode ) >= lsupport_stage2_ ) {
					lnode = unode;
					lscore = uscore;
					ival = 1.;
				} else {
					unode = lnode;
					uscore = lscore;
					ival = 0.;
				}
			} else {
				if ( lscore >= qscore ) ival = .0;
				else ival = (qscore - lscore)/static_cast<float>( uscore - lscore );
			}
			
			lower_node = lnode;
			upper_node = unode;
			return ival;
			//prec.setNodeRange( lnode, unode, support );
			//prec.setInterpolationValue( ival );
			//return true;
		}
		
		float getTaxSignal( int placescore = std::numeric_limits<  int >::max() ) {
			assert( data_.size() ); //TODO: exception
			if ( unsorted_ ) sort();
			
			// create ordering of scores according to taxonomy
			typedef std::vector< boost::tuple< int, int > > ordertype;
// 			ordertype ordering; //TODO: remove
// 			ordering.reserve( data_.size() ); //TODO: remove
			std::vector< int > scores, dists;
			scores.reserve( data_.size() );
			dists.reserve( data_.size() );
			
			{ //construction of arrays
				int score;
				const TaxonNode* node;
				const TaxonNode* anchor;
				data_type_::const_iterator it = data_.begin();
				boost::tie( score, anchor ) = *it;
				node = anchor;
				bool uppermode = false;
				while ( true ) {
// 					ordering.push_back( boost::make_tuple( -(node->data->root_pathlength), score ) ); //TODO: remove
					if ( score > placescore ) {
						if ( uppermode ) break;
						placescore = score*bandwidth_factor_;
						uppermode = true;
					}
					scores.push_back( score );
					dists.push_back( -(node->data->root_pathlength) ); //negative pathlength for sorting
					if ( ++it == data_.end() ) break; //single exit condition
					boost::tie( score, node ) = *it;
					node = taxinter_.getLCA( node, anchor );
				}
			}
			
			double tau = kendall( scores, dists );
			
// 			{ //sorting TODO: remove
// 				ordertype::iterator start_it = ordering.begin();
// 				std::sort( start_it, ordering.end() );
// 			}
			
// 			{ //DEBUG OUTPUT TODO: remove
// 				std::cerr << std::endl << "Taxonomic score array: " << std::endl;
// 				int score, depth;
// 				for ( ordertype::const_iterator it = ordering.begin(); it != ordering.end(); ++it ) {
// 					boost::tie( depth, score ) = *it;
// 					std::cerr << score << " at depth " << depth << std::endl;
// 				}
// 				std::cerr << "Kendall's tau: " << tau << std::endl;
// 			}
			
			return static_cast< float >( tau );
		};
		
	private:
		void sort() { //sort by increasing order
			data_type_::iterator start_it = data_.begin();
			std::sort( ++start_it, data_.end(), comparator_ );
			unsorted_ = false;
		}
		
		double kendall( const std::vector<int>& X1, const std::vector<int>& X2 ) { //TODO: acknowledge contribution!
			assert( X1.size() == X2.size() );
			std::size_t len = X1.size();
			int m1 = 0, m2 = 0, s = 0, nPair;
			std::size_t i, j;
			double denominator1, denominator2;

			for(i = 0; i < len; i++) {
				for(j = i + 1; j < len; j++) {
					if(X2[i] > X2[j]) {
						if (X1[i] > X1[j]) {
								s++;
						} else if(X1[i] < X1[j]) {
								s--;
						} else {
								m1++;
						}
					} else if(X2[i] < X2[j]) {
						if (X1[i] > X1[j]) {
								s--;
						} else if(X1[i] < X1[j]) {
								s++;
						} else {
								m1++;
						}
					} else {
						m2++;

						if(X1[i] == X1[j]) {
								m1++;
						}
					}
				}
			}

			nPair = len * (len - 1) / 2;
			denominator1 = nPair - m1;
			denominator2 = nPair - m2;
			return s / sqrt(denominator1) / sqrt(denominator2);
		}
		
		compareTupleFirstLT< boost::tuple< int, const TaxonNode* >, 0 > comparator_;
		typedef std::vector< boost::tuple< int, const TaxonNode* > > data_type_;
		bool unsorted_;
		data_type_ data_;
		TaxonomyInterface taxinter_;
		const float bandwidth_factor_;
		const small_unsigned_int lsupport_;
		const small_unsigned_int usupport_;
		const small_unsigned_int lsupport_stage2_;
};


// TODO: timers should be thread safe!
template< typename ContainerT, typename QStorType, typename DBStorType >
class RPAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		RPAPredictionModel( const Taxonomy* tax, QStorType& q_storage, const DBStorType& db_storage, float reeval_bandwidth = .1, std::ostream& debug_output = std::cerr ) :
					TaxonPredictionModel< ContainerT >( tax ),
					query_sequences_( q_storage ),
					db_sequences_( db_storage ),
					debug_output_( debug_output ),
					reeval_bandwidth_factor_( 1. - reeval_bandwidth ), //TODO: check range
					measure_placement_algorithm( "placement algorithm" ),
					measure_phase2_alignment( "phase 2 alignments" ),
					measure_reeval_alignment( "re-evaluation alignments" ),
					measure_sequence_retrieval( "sequence retrieval" )
					{};

		void predict( ContainerT& recordset, PredictionRecord& prec ) {
			this->initPredictionRecord( recordset, prec ); //set name and length
			
			// cannot be part of class because it is not thread-safe
			boost::format db_seqname_fmt( "%d:%d@%s tax=%s" );
			boost::format query_seqname_fmt( "%d:%d@%s len=%d" );

			typename ContainerT::iterator rec_it = firstUnmaskedIter( recordset );
			
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
			
			// determine position range of query to consider
			uint n = 0;
			assert( qrstart <= qrstop );			++n;
			while ( ++rec_it != recordset.end() ) {
				
				assert( (*rec_it)->getQueryStart() <= (*rec_it)->getQueryStop() );
				
				if ( ! (*rec_it)->isFiltered() ) {
					qrstart = std::min( (*rec_it)->getQueryStart(), qrstart );
					qrstop = std::max( (*rec_it)->getQueryStop(), qrstop );
					qprevscore = std::max( qprevscore, (*rec_it)->getScore() );
					++n;
				}
			}
			
			// lots of objects TODO: make use of C-style arrays
			int extra_mismatches = qlength - qrstop + qrstart - 1;
			const seqan::Dna5String qrseq = query_sequences_.getSequence( qid, qrstart, qrstop );
			const std::string qrseqname = boost::str( query_seqname_fmt % qrstart % qrstop % qid % qlength );
			const std::string debug_log_filename = qrseqname + ".log";
			const std::string debug_seqs_filename = qrseqname + ".fna";
			const std::string debug_dist_filename = qrseqname + ".dist";
			std::ofstream debug_log( debug_log_filename.c_str() );
			std::ofstream debug_seqs( debug_seqs_filename.c_str() );
			std::ofstream debug_dist( debug_dist_filename.c_str() );
			std::vector< AlignmentRecordTaxonomy* > records_ordered;
			std::vector< seqan::Dna5String > rrseqs_ordered; //TODO: boost ptr_container/seqan::StringSet/set to detect equal sequences
			std::vector< std::string > rrseqs_ordered_names;
			std::vector< int > rrseqs_qscores;
			std::vector< int > rrseqs_matches;
			records_ordered.reserve( n );
			rrseqs_ordered.reserve( n );
			rrseqs_matches.reserve( n );
			rrseqs_qscores.reserve( n );
			seqan::Align< seqan::Dna5String > aln1, aln2;
			seqan::resize( seqan::rows( aln1 ), 2);
			seqan::resize( seqan::rows( aln2 ), 2);
			seqan::assignSource( seqan::row( aln1, 0 ), qrseq );
			seqan::assignSource( seqan::row( aln2, 0 ), qrseq );
			seqan::Align< seqan::Dna5String >* bestaln = &aln1;
			seqan::Align< seqan::Dna5String >* aln = &aln2;
			
			{ // retrieving sequence blocks
				measure_sequence_retrieval.start();
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
							rrseqs_ordered.push_back( db_sequences_.getSequenceComplement( rid, start, stop ) ); //TODO: avoid copying
							rrseqs_ordered_names.push_back( boost::str( db_seqname_fmt % stop % start % rid % rec.getReferenceNode()->data->taxid ) );
						}
					}
					++rec_it;
				}
				measure_sequence_retrieval.stop();
			}

			std::list< uint > anchor_indices;
			{ // determine best hit(s) by re-aligning the extended query against each of the (extended) DB sequence ranges (within band)
				measure_reeval_alignment.start();
				float dbalignment_score_threshold = reeval_bandwidth_factor_*qprevscore;
				const int maxint = std::numeric_limits< int >::max();
				uint index_best = 0;
				
				for ( uint i = 0; i < n; ++i ) { //calculate scores for best-scoring references
					if ( records_ordered[i]->getScore() >= dbalignment_score_threshold ) { //TODO: avoid equal sequences
						anchor_indices.push_back( i );
						seqan::assignSource( seqan::row( *aln, 1 ), rrseqs_ordered[i] );
						const int score = -seqan::globalAlignment( *aln, seqan::SimpleScore(), seqan::MyersHirschberg() );
						rrseqs_qscores.push_back( score );
						rrseqs_matches.push_back( seqan::length( seqan::row( *aln, 0 ) ) - score );
						
						if ( score <= rrseqs_qscores[index_best] ) {
							index_best = i;
							std::swap( aln, bestaln ); //swap trick
						}
					} else { //fill in some dummy values
						rrseqs_qscores.push_back( maxint );
						rrseqs_matches.push_back( 0 );
					}
				}
				
				int threshold = 1.*rrseqs_qscores[index_best]; //TODO: make parameter
				std::list< uint >::iterator it = anchor_indices.begin();
				while ( it != anchor_indices.end() ) { //restrict anchors
					if ( rrseqs_qscores[*it] > threshold ) it = anchor_indices.erase( it );
					else ++it;
				}
				measure_reeval_alignment.stop();
			}

			SandwichPlacerAdaptive placer( this->taxinter_, n );
			{ // place for each anchor
				std::size_t anchornum = anchor_indices.size();
				std::vector< const TaxonNode* > anchors_lnode;
				anchors_lnode.reserve( anchornum );
				std::vector< const TaxonNode* > anchors_unode;
				anchors_unode.reserve( anchornum );
				float anchors_taxsig = 1;
				medium_unsigned_int anchors_support = std::numeric_limits< medium_unsigned_int >::max();
				float anchors_ival = 0.;
				
				while ( ! anchor_indices.empty() ) { // place using anchor
					uint& index_anchor = anchor_indices.front();
				
					// set anchor as reference in alignment
					placer.addSequence( 0, records_ordered[index_anchor]->getReferenceNode() );
					seqan::assignSource( seqan::row( *aln, 1 ), rrseqs_ordered[ index_anchor ] ); // reuse *aln
					
					// align all others against anchor
					measure_phase2_alignment.start();
					debug_dist << rrseqs_qscores[index_anchor] << '\t' << qrseqname << std::endl;
					debug_seqs << '>' << qrseqname << std::endl << qrseq << std::endl;
					for( uint i = 0; i < n; ++i ) { //TODO: avoid equal sequences
						seqan::assignSource( seqan::row( *aln, 0 ), rrseqs_ordered[ i ] );
						int score;
						if ( i == index_anchor ) score = 0;
						else score = -seqan::globalAlignment( *aln, seqan::SimpleScore(), seqan::MyersHirschberg() ); //- extra_mismatches;

						debug_dist << score << '\t' << rrseqs_ordered_names[i] << std::endl;
						debug_seqs << '>' << rrseqs_ordered_names[i] << std::endl << rrseqs_ordered[i] << std::endl;

						placer.addSequence( score, records_ordered[i]->getReferenceNode() );
					}
					measure_phase2_alignment.stop();
				
					// taxonomic placement of query
					measure_placement_algorithm.start();
					const TaxonNode* lnode;
					const TaxonNode* unode;
					const int& qscore = rrseqs_qscores[index_anchor];
					const medium_unsigned_int& support = rrseqs_matches[index_anchor];
					const float ival = placer.placeSequence( qscore, lnode, unode );
					const float taxsig = placer.getTaxSignal( qscore );
					anchors_ival = std::max( ival, anchors_ival );
					anchors_support = std::min( support, anchors_support );
					anchors_taxsig = std::min( taxsig, anchors_taxsig );
					anchors_lnode.push_back( lnode );
					anchors_unode.push_back( unode );
					
					measure_placement_algorithm.stop();

					debug_log << "query = " << qrseqname << std::endl;
					debug_log << "flanking_mismatches = " << extra_mismatches << std::endl;
					debug_log << "original score: " << qprevscore << std::endl;
					debug_log << "realignment factor: " << reeval_bandwidth_factor_ << std::endl << std::endl;

					debug_log << " PREDICTION:" << std::endl;
					debug_log << " lower node: " << lnode->data->annotation->name << " (" << lnode->data->taxid << ")" << std::endl;
					debug_log << "      query: " << qid << " with score " << qscore << std::endl;
					debug_log << " upper node: " << unode->data->annotation->name << " (" << unode->data->taxid << ")" << std::endl;
					debug_log << "    support: " << support << std::endl;
					debug_log << "       ival: " << ival << std::endl << std::endl;
					debug_log << "     taxsig: " << taxsig  << std::endl << std::endl;
				
					anchor_indices.pop_front();
				}
				prec.setSignalStrength( anchors_taxsig );
				prec.setQueryFeatureBegin( qrstart );
				prec.setQueryFeatureEnd( qrstop );
				prec.setInterpolationValue( anchors_ival );
				prec.setNodeRange( this->taxinter_.getLCA( anchors_lnode ), this->taxinter_.getLCA( anchors_unode ), anchors_support );
			}
		}

	protected:
		QStorType& query_sequences_;
		const DBStorType& db_sequences_;
		SortByBitscoreFilter< ContainerT > sort_;
		compareTupleFirstLT< boost::tuple< int, uint >, 0 > tuple_1_cmp_le_;

		private:
			std::ostream& debug_output_;
			const float reeval_bandwidth_factor_;
			StopWatchCPUTime measure_placement_algorithm;
			StopWatchCPUTime measure_phase2_alignment;
			StopWatchCPUTime measure_reeval_alignment;
			StopWatchCPUTime measure_sequence_retrieval;
};


#endif // taxonpredictionmodelsequence_hh_
