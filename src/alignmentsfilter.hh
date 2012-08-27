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

#ifndef alignmentsfilter_hh_
#define alignmentsfilter_hh_

#include <boost/regex.hpp>
#include <cmath>
#include <string>
#include "types.hh"
#include "alignmentrecord.hh"
#include "taxonomyinterface.hh"
#include "constants.hh"


// abstract base class
template< typename ContainerT >
class AlignmentsFilter {
	public:
		typedef ContainerT AlignmentRecordSetType;
		
		virtual ~AlignmentsFilter() {};
		virtual void filter( ContainerT& recordset ) = 0;
		virtual const std::string getInfo() { return description; };
	private:
		static const std::string description;
};

template< typename ContainerT >
const std::string AlignmentsFilter< ContainerT >::description = "abstract AlignmentsFilter";


// pseudo filters to get some information (needs redesign)

template< typename ContainerT >
class MaxBitscoreAlignmentFilter : public AlignmentsFilter< ContainerT > {
	public:
		typedef typename ContainerT::value_type AlignmentRecordPtrType;
		
		void filter( ContainerT& recordset ) {
			best_records.clear();
			if( ! recordset.empty() ) {

				// go to first valid alignment
				typename ContainerT::iterator record_it = recordset.begin();
				while( true ) {
					if( record_it == recordset.end() ) {
						return;
					}
					if( ! (*record_it)->isFiltered() ) {
						break;
					}
					++record_it;
				}

				float max_bs = (*record_it++)->getScore();

				for( ; record_it != recordset.end(); ++record_it ) {
					if( ! (*record_it)->isFiltered() && (*record_it)->getScore() > max_bs ) {
						max_bs = (*record_it)->getScore();
					}
				}

				// scan for maximum set
				for( ; record_it != recordset.end(); ++record_it ) {
					if( ! (*record_it)->isFiltered() ) {
						if( (*record_it)->getScore() == max_bs ) {
						best_records.push_back( *record_it );
						} else {
						  if( (*record_it)->getScore() > max_bs ) {
						    best_records.clear();
						    best_records.push_back( *record_it );
						    max_bs = (*record_it)->getScore();
              }
						}
					}
				}
			}
		}

		AlignmentRecordPtrType getBest() {
			if( best_records.empty() ) {
				return NULL;
			}
			return best_records.front();
		}

		const std::list< AlignmentRecordPtrType >& getBests() {
			return best_records;
		}


	private:
		std::list< AlignmentRecordPtrType > best_records;
		static const std::string description;
};

template< typename ContainerT >
const std::string MaxBitscoreAlignmentFilter< ContainerT >::description = "MaxBitscoreAlignmentFilter";



template< typename ContainerT >
class MinMaxBitscoreFilter : public AlignmentsFilter< ContainerT > {
	public:
		void filter( ContainerT& recordset ) {
			if( ! recordset.empty() ) {
				typename ContainerT::iterator record_it = recordset.begin();
				float tmp_max, tmp_min;

				while( record_it != recordset.end() ) {
					if( ! (*record_it)->isFiltered() ) {
						tmp_max = tmp_min = (*record_it)->getScore();
						++record_it;
						break;
					}
					++record_it;
				}

				while( record_it != recordset.end() ) {
					if( ! (*record_it)->isFiltered() ) {
						float& bitscore = (*record_it)->getScore();
						tmp_min = std::min( tmp_min, bitscore );
						tmp_max = std::max( tmp_max, bitscore );
					}
					++record_it;
				}
				min_bitscore = tmp_min;
				max_bitscore = tmp_max;
			} else {
				min_bitscore = 0.0;
				max_bitscore = 0.0;
			}
		}

		float getMin() {
			return min_bitscore;
		}

		float getMax() {
			return max_bitscore;
		}

		std::pair< float, float > getMinMax() {
			return std::make_pair( min_bitscore, max_bitscore );
		}

	private:
		float min_bitscore;
		float max_bitscore;
		static const std::string description;
};

template< typename ContainerT >
const std::string MinMaxBitscoreFilter< ContainerT >::description = "MinMaxBitscoreFilter";


template< typename ContainerT, typename Comparator = std::greater< float > >
class SortByBitscoreFilter : public AlignmentsFilter< ContainerT > { //includes masked records
	public:
		typedef typename ContainerT::value_type AlignmentRecordPtrType;
		
		void filter( ContainerT& recordset ) {
			std::multimap< float, AlignmentRecordPtrType, Comparator > sorttable;

			typename ContainerT::iterator record_it = recordset.begin();
			while( ! recordset.empty() ) {
				AlignmentRecordPtrType record = recordset.back();
				sorttable.insert( std::make_pair( record->getScore(), record ) );
				recordset.pop_back();
			}

			for( typename std::multimap< float, AlignmentRecordPtrType, std::greater< float > >::iterator it = sorttable.begin(); it != sorttable.end(); ++it ) {
				recordset.push_back( it->second );
			}
		}
	private:
		static const std::string description;
};

template< typename ContainerT, typename Comparator >
const std::string SortByBitscoreFilter< ContainerT, Comparator >::description = "SortByBitscoreFilter";





// real filter that mask a subset or all of the alignments based on some criteria

// experimental filter that takes a core set of good alignments and a taxonomy-distance [0,1] cutoff for all remaining
template< typename ContainerT >
class CleanseFDistAlignmentFilter : public SortByBitscoreFilter< ContainerT > {
	public:
		CleanseFDistAlignmentFilter( Taxonomy* tax, const float t1, const float t2 ) : taxinter( tax ), coreset_threshold( 1.0 - t1 ), cutoff( t2 ) {
			//min_rel_bs( t1 )
// 			std::cerr << "distance cutoff is: " << cutoff << std::endl;
		};

		void filter( ContainerT& recordset ) {
			if( ! recordset.empty() ) {
				SortByBitscoreFilter< ContainerT >::filter( recordset );

				typename ContainerT::iterator record_it = recordset.begin();
				std::list< const TaxonNode* > bestnodes;

				typename ContainerT::iterator it = recordset.begin();
				while( it != recordset.end() && (*it)->isFiltered() ) { ++it; } //get best valid alignment
				float best_bs = (*it)->getScore();
				const TaxonNode* tmpnode = (*it++)->getReferenceNode();
				bestnodes.push_back( tmpnode );

				for( ; it != recordset.end() && (*it)->getScore() >= coreset_threshold*best_bs; ++it ) { //collate all best hits until cutoff
					if( ! (*it)->isFiltered() ) {
						const TaxonNode* tmpnode = (*it++)->getReferenceNode();
						bestnodes.push_back( tmpnode );
					}
				}

// 				std::cerr << "number of alignments: " << recordset.size() << std::endl;
// 				std::cerr << "size of core set: " << bestnodes.size() << std::endl;

				// wheight remaining alignments by combined distance
				for( ; it != recordset.end(); ++it ) {
					if( ! (*it)->isFiltered() ) {
						const TaxonNode* tmpnode = (*it)->getReferenceNode();
						float bs_dist = 1.0 - (*it)->getScore() / best_bs;
						float tree_dist = getNormDist( tmpnode, bestnodes );
						float comb_dist = ( bs_dist + tree_dist ) / 2.0;
// 						std::cerr << "bitscore_dist: " << bs_dist << " tree_dist: " << tree_dist << " combined average distance: " << comb_dist << std::endl;
						if( comb_dist > cutoff ) {
// 							std::cerr << comb_dist << " > " <<  cutoff << ", masking this alignment" << std::endl;
							(*it)->filterOut();
						} else {
// 							std::cerr << comb_dist << " <= " <<  cutoff << ", not masking this alignment" << std::endl;
						}
					}
				}
			}
		}

	private:
		float getNormDist( const TaxonNode* n, std::list< const TaxonNode* >& bestnds ) {
			int a, b, c;
			int c_sum = 0;
			for( std::list< const TaxonNode* >::iterator it = bestnds.begin(); it != bestnds.end(); ++it ) {
				boost::tie( a, b, c ) = this->taxinter.getInterDistances( n, *it );
				c_sum += c;
			}
			return c_sum/(float)( bestnds.size()*default_rank_number );
		};

    TaxonomyInterface taxinter;
		const float coreset_threshold;
		const float cutoff;
// 		const float min_rel_bs;
// 		float threshold_square;
		static const std::string description;
};

template< typename ContainerT >
const std::string CleanseFDistAlignmentFilter< ContainerT >::description = "CleanseFDistAlignmentFilter";



template< typename ContainerT >
class OutlierDetectionAlignmentFilter : public MaxBitscoreAlignmentFilter< ContainerT > {
	public:
		OutlierDetectionAlignmentFilter( Taxonomy* tax, const float t, const int m = 3 ) : taxinter( tax ), threshold( t ), minsize( m ) {};

		void filter( ContainerT& recordset ) {
			if( ! recordset.empty() ) {
				MaxBitscoreAlignmentFilter< ContainerT >::filter( recordset );

				// 				std::cerr << "creating data point from recordset of size " << recordset.size() << std::endl;
				std::list< DataPoint > datapoints;
				std::list< AlignmentRecord* > bests = MaxBitscoreAlignmentFilter< ContainerT >::getBests();
				float best_bs = bests.front()->getScore();
				int best_num = bests.size() - 1;
				for( typename ContainerT::iterator it = recordset.begin(); it != recordset.end(); ++it ) {
					if( ! (*it)->isFiltered() ) {
						const TaxonNode* tmp_node = (*it)->getReferenceNode();
						if( tmp_node ) {
							datapoints.push_back( boost::make_tuple( tmp_node, *it ) );
						}
					}
				}
// 				std::cerr << "Number of data points: " << datapoints.size() << std::endl;

				std::multimap< float, AlignmentRecord*, std::greater< float > > dist_to_others;

				// simple n^2 distance calculations without optimization
				float avg_sum = 0.0;
				float n = (float) datapoints.size();
				float nminusone = n - 1.0;
				std::vector< float > sums;
				std::list< DataPoint >::iterator outer_it = datapoints.begin();
				while( outer_it != datapoints.end() ) {
					float sum = 0.0;
					std::list< DataPoint >::iterator inner_it = datapoints.begin();
					while( inner_it != outer_it ) { // starting at first data point
						sum += getDistance( *outer_it, *inner_it++ );
					}
					++inner_it; // jump over self distance here
					while( inner_it != datapoints.end() ) { // starting one after self data point
						sum += getDistance( *outer_it, *inner_it++ );
					}
					if( sum ) {
						sum /= nminusone;
					}
// 					std::cerr << "avg sum of distances: " << sum << std::endl;
					avg_sum += sum;
					sums.push_back( sum );
					dist_to_others.insert( std::make_pair( sum, outer_it->get<1>() ) );
					// 					std::cerr << "point w. sum of  distances to others: " << sum << std::endl;
					++outer_it;
				}

				//save calculation time
				if( datapoints.size() < minsize ) {
					return;
				}

				//use median as robust estimate of mean
				avg_sum /= n;
				float median;
				std::sort( sums.begin(), sums.end() );
				if( sums.size() % 2 == 0 ) {
					int index = sums.size()/2;
					median = ( sums[ index - 1 ] + sums[ index ] )/2.0;
				} else {
					median = sums[ (sums.size()-1)/2 ];
				}

				// calculate std_dev using median instead of mean
				float variance = 0.0;
				std::vector< float > abs_median_dists;
				for( std::vector< float >::iterator it = sums.begin(); it != sums.end(); ++it ) {
					float tmp = *it - median;
// 					std::cerr << "point minus median = " << tmp << std::endl;
					abs_median_dists.push_back( std::abs( tmp ) );
					variance += tmp*tmp;
				}
				sums.clear();
				if( variance ) {
					variance /= n;
				}
				float std_dev = std::sqrt( variance );

				// calulate MAD as replacement for std_dev
// 				float MAD;
// 				std::sort( abs_median_dists.begin(), abs_median_dists.end() );
// 				if( abs_median_dists.size() % 2 == 0 ) {
// 					int index = abs_median_dists.size()/2;
// 					MAD = ( abs_median_dists[ index - 1 ] + abs_median_dists[ index ] )/2.0;
// 				} else {
// 					MAD = abs_median_dists[ (abs_median_dists.size()-1)/2 ];
// 				}

				if( ! std_dev ) {
					return;
				}

				float cutoff = threshold*std_dev + median;
				/*if( MAD ) {
					cutoff = threshold*MAD/0.6745 + median; //use formula by Iglewicz and Hoaglin
				} else {
					if( std_dev ) {
						cutoff = threshold*std_dev + median; //fallback on std_dev
					} else {
						//will not mask anything
						return;
					}
				}*/

// 				std::cerr << "average: " << avg_sum << " median: " << median << " variance: " << variance << " std dev: " << std_dev << " MAD_stddev_replacement: " << 1.4826*MAD << " cutoff: " << cutoff << std::endl;

				// mark outliers
				for( std::multimap< float, AlignmentRecord* >::iterator it = dist_to_others.begin(); it != dist_to_others.end() && it->first > cutoff; ++it ) {
// 					std::cerr << "value for exclusion is: " << it->first << std::endl;
					if( it->second->getScore() == best_bs  ) {
						if( best_num ) {
// 							std::cerr << "skipping data point because it is the only best alignment left" << std::endl;
							continue;
						}
// 						std::cerr << "removing one of the best alignments!" << std::endl;
						--best_num;
					}
					it->second->filterOut();
// 					std::cerr << "* ";
// 					std::cerr << it->first << " (bitscore: " << it->second->getScore()/best_bs << ")\t";
				}
// 				std::cerr << std::endl;
// 				std::cerr << " (of " << datapoints.size() << ")" << std::endl;
			}
		}

	private:
		typedef boost::tuple< const TaxonNode*, AlignmentRecord* > DataPoint;

		float getDistance( const DataPoint& p1, const DataPoint& p2 ) {
			int a, b, c;
			boost::tie( a, b, c ) = this->taxinter.getInterDistances( p2.get<0>(), p1.get<0>() );
			return c/(float) default_rank_number;
		}

		TaxonomyInterface taxinter;
		const float threshold;
		const unsigned int minsize;
		static const std::string description;
};

template< typename ContainerT >
const std::string OutlierDetectionAlignmentFilter< ContainerT >::description = "OutlierDetectionAlignmentFilter";



template< typename ContainerT >
class RemoveRedundantFilter : public AlignmentsFilter< ContainerT > { //expects list to be sorted decreasingly
	public:
		RemoveRedundantFilter( Taxonomy* tax ) : taxinter( tax ) {};

		void filter( ContainerT& recordset ) {
			if( ! recordset.empty() ) {
				typename ContainerT::iterator record_it = recordset.begin();

				// set lca to first valid alignment
				const TaxonNode* lca;
				while( record_it != recordset.end() ) {
					if( ! (*record_it)->isFiltered() ) {
						lca = (*record_it)->getReferenceNode();
						++record_it;
						break;
					}
					++record_it;
				}

				// see whether the other alignments contribute or not
				while( record_it != recordset.end() ) {
					if( ! (*record_it)->isFiltered() ) {
						const TaxonNode* tmp_node = (*record_it)->getReferenceNode();
						if( lca == tmp_node || taxinter.isParentOf( lca, tmp_node ) ) {
							(*record_it)->filterOut();
// 							std::cerr << "RemoveRedundantFilter: masking alignment..." << std::endl;
						} else {
							lca = taxinter.getLCA( lca, tmp_node );
						}
					}
					++record_it;
				}
// 				std::cerr << "RemoveRedundantFilter: LCA is " << lca->data->annotation->name << std::endl;
			}
		}

	private:
		TaxonomyInterface taxinter;
		static const std::string description;
};

template< typename ContainerT >
const std::string RemoveRedundantFilter< ContainerT >::description = "RemoveRedundantFilter";



template< typename ContainerT >
class MinScoreTopPercentFilter : public AlignmentsFilter< ContainerT > {
	public:
		MinScoreTopPercentFilter( const float ms, const float tp ) : minscore( ms ), toppercent( tp ) {};

		void filter( ContainerT& recordset ) {
// 			if( minscore <= 0.0 && toppercent >= 1.0 ) { return; }

			float max_bitscore = .0;
			for( typename ContainerT::iterator record_it = recordset.begin(); record_it != recordset.end(); ++record_it ) {
				if( ! (*record_it)->isFiltered() ) {
					if( (*record_it)->getScore() < minscore ) {
						(*record_it)->filterOut();
					} else {
						if( (*record_it)->getScore() > max_bitscore ) {
							max_bitscore = (*record_it)->getScore();
						}
					}
				}
			}

// 			if( toppercent >= 1.0 ) { return; }
			max_bitscore = ( 1.0 - toppercent ) * max_bitscore;

			for( typename ContainerT::iterator record_it = recordset.begin(); record_it != recordset.end(); ++record_it ) {
				if ( ! (*record_it)->isFiltered() && (*record_it)->getScore() < max_bitscore ) {
					(*record_it)->filterOut();
				}
			}
		}

	private:
		const float minscore;
		const float toppercent;
		static const std::string description;
};

template< typename ContainerT >
const std::string MinScoreTopPercentFilter< ContainerT >::description = "MinScoreTopPercentFilter";



template< typename ContainerT >
class MinScoreMaxEvalueTopPercentFilter : public AlignmentsFilter< ContainerT > {
	public:
		MinScoreMaxEvalueTopPercentFilter( const float ms, const float ev, const float tp ) : minscore( ms ), maxevalue( ev ), toppercent( tp ) {};

		void filter( ContainerT& recordset ) {
			float max_bitscore = .0;
			for( typename ContainerT::iterator record_it = recordset.begin(); record_it != recordset.end(); ++record_it ) {
				if( ! (*record_it)->isFiltered() ) {
					if( (*record_it)->getScore() < minscore || (*record_it)->getEValue() > maxevalue ) {
						(*record_it)->filterOut();
					} else {
						if( (*record_it)->getScore() > max_bitscore ) {
							max_bitscore = (*record_it)->getScore();
						}
					}
				}
			}

			max_bitscore = ( 1.0 - toppercent ) * max_bitscore;

			for( typename ContainerT::iterator record_it = recordset.begin(); record_it != recordset.end(); ++record_it ) {
				if( (*record_it)->getScore() < max_bitscore ) {
					(*record_it)->filterOut();
				}
			}
		}

	private:
		const float minscore;
		const float maxevalue;
		const float toppercent;
		static const std::string description;
};

template< typename ContainerT >
const std::string MinScoreMaxEvalueTopPercentFilter< ContainerT >::description = "MinScoreMaxEvalueTopPercentFilter";



template< typename ContainerT >
class MinPIDFilter : public AlignmentsFilter< ContainerT > {
	public:
		MinPIDFilter( const float pid ) : minpid( pid ) {};

		void filter( ContainerT& recordset ) {
			typename ContainerT::iterator record_it = recordset.begin();
			while( record_it != recordset.end() ) {
				if( (*record_it)->getPID() < minpid ) {
					(*record_it)->filterOut();
				}
				++record_it;
			}
		}

	private:
		const float minpid;
		static const std::string description;
};

template< typename ContainerT >
const std::string MinPIDFilter< ContainerT >::description = "MinPIDFilter";



template< typename ContainerT >
class MaxEvalueMinScoreTopPercentFilter : public AlignmentsFilter< ContainerT > {
	public:
		MaxEvalueMinScoreTopPercentFilter( const float ms, const float tp, const double me ) : minscore( ms ), toppercent( tp ), maxevalue( me ) {};

		void filter( ContainerT& recordset ) {
			typename ContainerT::iterator record_it = recordset.begin();
			float max_bitscore = .0;
			while( record_it != recordset.end() ) {
				max_bitscore = std::max( (*record_it)->getScore(), max_bitscore ); //use it as toppercent value, no matter if it will be mask in the second step

				if( (*record_it)->getEValue() > maxevalue || (*record_it)->getScore() < minscore ) {
					(*record_it)->filterOut();
				}
				++record_it;
			}

			if( toppercent >= 1.0 ) { return; }
			record_it = recordset.begin();
			max_bitscore = ( 1.0 - toppercent ) * max_bitscore;
			while( record_it != recordset.end() ) {
				if( (*record_it)->getScore() < max_bitscore ) {
					(*record_it)->filterOut();
				}
				++record_it;
			}
		}

// 		const std::string getInfo() { return "MaxEvalueMinScoreTopPercentFilter"; };

	private:
		const float minscore;
		const float toppercent;
		const double maxevalue;
		static const std::string description;
};

template< typename ContainerT >
const std::string MaxEvalueMinScoreTopPercentFilter< ContainerT >::description = "MaxEvalueMinScoreTopPercentFilter";



template< typename ContainerT >
class MinSupportFilter : public AlignmentsFilter< ContainerT > {
	public:
		MinSupportFilter( const int ms ) : minsupport( ms ) {};

		void filter( ContainerT& recordset ) {
			typename ContainerT::iterator record_it = recordset.begin();
			int count = 0;
			while( record_it != recordset.end() ) {
        count += ! (*record_it++)->isFiltered();
			}

			if ( count < minsupport ) {
			  record_it = recordset.begin();
        while( record_it != recordset.end() ) {
          (*record_it++)->filterOut();
        }
		  }
		}

	private:
		const int minsupport;
		static const std::string description;
};

template< typename ContainerT >
const std::string MinSupportFilter< ContainerT >::description = "MinSupportFilter";



template< typename ContainerT >
class NumBestBitscoreFilter : public AlignmentsFilter< ContainerT > {
	public:
		NumBestBitscoreFilter( const int nbb ) : numbestbitscore( nbb ) {};

		void filter( ContainerT& recordset ) {
			std::multimap< float, AlignmentRecord*, std::greater<float> > sorted_bitscores;

			for( typename ContainerT::iterator record_it = recordset.begin(); record_it != recordset.end(); ++record_it ) {
				if( ! (*record_it)->isFiltered() ) {
					sorted_bitscores.insert( std::make_pair( (*record_it)->getScore(), (*record_it) ) );
				}
			}
			
			if ( sorted_bitscores.empty() ) return;

			// mark records to be masked
			unsigned int count = numbestbitscore;

			// spare the first nbb records not masked yet
			std::multimap< float, AlignmentRecord*, std::greater<float> >::iterator sb_it = sorted_bitscores.begin();
			float lastvalue = sb_it++->first;
			for( ; sb_it != sorted_bitscores.end(); ++sb_it ) {
				if( sb_it->first != lastvalue ) {
					if( --count <= 0 ) {
						break;
					}
					lastvalue = sb_it->first;
				}
			}

			// mask the rest
			for( ; sb_it != sorted_bitscores.end(); ++sb_it ) {
				sb_it->second->filterOut();
			}
		}

// 		const std::string getInfo() { return "NumBestBitscoreFilter"; };

	private:
		const int numbestbitscore;
		static const std::string description;
};

template< typename ContainerT >
const std::string NumBestBitscoreFilter< ContainerT >::description = "NumBestBitscoreFilter";



template< typename ContainerT >
class BestScorePerReferenceSeqIDFilter : public AlignmentsFilter< ContainerT > {
	public:
// 		BestScorePerReferenceIDFilter(){};

		void filter( ContainerT& recordset ) {
			std::map< std::string, AlignmentRecord* > keep;
			std::map< std::string, AlignmentRecord* >::iterator keep_it;
			//mask all records having the same gi but a worse bitscore
			for( typename ContainerT::iterator record_it = recordset.begin(); record_it != recordset.end(); ++record_it ) {
				if( ! (*record_it)->isFiltered() ) {
					const std::string& seqid = (*record_it)->getReferenceIdentifier();
					keep_it = keep.find( seqid );
					if( keep_it != keep.end() ) {
						if( keep_it->second->getScore() < (*record_it)->getScore() ) {
							keep_it->second->filterOut();
							keep_it->second = (*record_it);
						} else {
							(*record_it)->filterOut();
						}
					} else {
						keep[ seqid ] = (*record_it);
					}
				}
			}
		}
	private:
		static const std::string description;
};

template< typename ContainerT >
const std::string BestScorePerReferenceSeqIDFilter< ContainerT >::description = "BestScorePerReferenceSeqIDFilter";



template< typename ContainerT >
class BestScorePerReferenceTaxIDFilter : public AlignmentsFilter< ContainerT > {
	public:
 		BestScorePerReferenceTaxIDFilter(){};

		void filter( ContainerT& recordset ) {
			std::map< unsigned int, AlignmentRecord* > keep;
			std::map< unsigned int, AlignmentRecord* >::iterator keep_it;
			//mask all records having the same gi but a worse bitscore
			for( typename ContainerT::iterator record_it = recordset.begin(); record_it != recordset.end(); ++record_it ) {
				if( ! (*record_it)->isFiltered() ) { //TODO: change taxid to node pointer content
					TaxonID taxid = (*record_it)->getReferenceNode()->data->reference_taxid;
					keep_it = keep.find( taxid );
					if( keep_it != keep.end() ) {
						if( keep_it->second->getScore() < (*record_it)->getScore() ) {
							keep_it->second->filterOut();
							keep_it->second = (*record_it);
						} else {
							(*record_it)->filterOut();
						}
					} else {
						keep[ taxid ] = (*record_it);
					}
				}
			}
		}
	private:
		static const std::string description;
};

template< typename ContainerT >
const std::string BestScorePerReferenceTaxIDFilter< ContainerT >::description = "BestScorePerReferenceTaxIDFilter";





// the following filters are supervised: they expect the taxid of the query sequence to be known and contained in the name (needs redesign)

template< typename ContainerT >
class RemoveUnclassifiedQueriesFilter : public AlignmentsFilter< ContainerT > {
	public:
		RemoveUnclassifiedQueriesFilter( StrIDConverter& accessconv, Taxonomy* tax, const std::string extract_re ) : identifier_regex_( extract_re ), seqid2taxid( accessconv ), taxinter( tax ) {
		}

		void filter( ContainerT& recordset ) {
			if( ! recordset.empty() ) {
				typename ContainerT::iterator record_it = recordset.begin();
				const std::string seqid = extractIdentifier( (*record_it)->getQueryIdentifier() );
				if( taxinter.getNode( seqid2taxid[ seqid ] )->data->is_unclassified ) { //TODO: case taxid not found
					for( ; record_it != recordset.end(); ++record_it ) {
						(*record_it)->filterOut();
					}
				}
			}
		}
	private:
		const std::string extractIdentifier( const std::string& whole_name ) {
			boost::cmatch re_results;
			assert( boost::regex_match( whole_name.c_str(), re_results, identifier_regex_ ) );
			assert(  re_results.size() > 1 );
			return std::string( re_results[1].first, re_results[1].second );
		}
		
		const boost::regex identifier_regex_;
		static const std::string description;
		StrIDConverter& seqid2taxid;
		TaxonomyInterface taxinter;
};

template< typename ContainerT >
const std::string RemoveUnclassifiedQueriesFilter< ContainerT >::description = "RemoveUnclassifiedQueriesFilter";



template< typename ContainerT >
class TaxonMaskingFilter : public AlignmentsFilter< ContainerT > {
	public:
		TaxonMaskingFilter( StrIDConverter& staxon, StrIDConverter& rtaxon ) : staxon_( staxon ), rtaxon_( rtaxon ) {};
		
		void filter( ContainerT& recordset ) {
			if( ! recordset.empty() ) {
				typename ContainerT::iterator record_it = recordset.begin();
        const TaxonID qtax = staxon_[ (*record_it)->getQueryIdentifier() ];
				while( record_it != recordset.end() ) {
					if( qtax == rtaxon_[ (*record_it)->getReferenceIdentifier() ] ) {
						(*record_it)->filterOut();
					}
					++record_it;
				}
			}
		}
		
	private:
		StrIDConverter& staxon_;
		StrIDConverter& rtaxon_;
		static const std::string description;
};

template< typename ContainerT >
const std::string TaxonMaskingFilter< ContainerT >::description = "TaxonMaskingFilter";



template< typename ContainerT >
class RemoveIdentSeqIDFilter : public AlignmentsFilter< ContainerT > {
	public:
		RemoveIdentSeqIDFilter( const std::string extract_re ) : identifier_regex_( extract_re ) {}; //TODO: catch boost::regex_error
		
		void filter( ContainerT& recordset ) {
			if( ! recordset.empty() ) {
				typename ContainerT::iterator record_it = recordset.begin();
        const std::string seqid = extractIdentifier( (*record_it)->getQueryIdentifier() ); //extractFastaCommentField( (*record_it)->getQueryIdentifier(), "gi" );
				while( record_it != recordset.end() ) {
					if( seqid == (*record_it)->getReferenceIdentifier() ) {
						(*record_it)->filterOut();
					}
					++record_it;
				}
			}
		}
	private:
		const std::string extractIdentifier( const std::string& whole_name ) {
			boost::cmatch re_results;
			assert( boost::regex_match( whole_name.c_str(), re_results, identifier_regex_ ) );
			assert(  re_results.size() > 1 );
			return std::string( re_results[1].first, re_results[1].second );
		}
		
		const boost::regex identifier_regex_;
		static const std::string description;
};

template< typename ContainerT >
const std::string RemoveIdentSeqIDFilter< ContainerT >::description = "RemoveIdentSeqIDFilter";



template< typename ContainerT >
class RemoveIdentTaxIDFilter : public AlignmentsFilter< ContainerT > {
	public:
		RemoveIdentTaxIDFilter( StrIDConverter& accessconv, const std::string extract_re ) : identifier_regex_( extract_re ), seqid2taxid( accessconv ) {}; //TODO: catch boost::regex_error
		void filter( ContainerT& recordset ) {
			if( ! recordset.empty() ) {
				typename ContainerT::iterator record_it = recordset.begin();
				const std::string seqid = extractIdentifier( (*record_it)->getQueryIdentifier() ); //extractFastaCommentField( (*record_it)->getQueryIdentifier(), "gi" );
				try {
					TaxonID taxid = seqid2taxid[ seqid ];
					while( record_it != recordset.end() ) {
						TaxonID reftaxid = seqid2taxid[ (*record_it)->getReferenceIdentifier() ]; //TODO: use (*record_it)->getReferenceNode and traverse...
						
						if( taxid == reftaxid ) {
							(*record_it)->filterOut();
						}
						++record_it;
					}
				} catch ( std::out_of_range e ) { //mask all records
					std::cerr << "RemoveIdentTaxIDFilter: Could not map sequence id " << seqid << " to TaxID";
					std::cerr << ", skipping all records for record set..." << std::endl;
					while( record_it != recordset.end() ) {
						(*record_it)->filterOut();
						++record_it;
					}
				}
			}
		}
	private:
		const std::string extractIdentifier( const std::string& whole_name ) {
			boost::cmatch re_results;
			assert( boost::regex_match( whole_name.c_str(), re_results, identifier_regex_ ) );
			assert(  re_results.size() > 1 );
			return std::string( re_results[1].first, re_results[1].second );
		}
		
		const boost::regex identifier_regex_;
		static const std::string description;
		StrIDConverter& seqid2taxid;
};

template< typename ContainerT >
const std::string RemoveIdentTaxIDFilter< ContainerT >::description = "RemoveIdentTaxIDFilter";


// template< typename ContainerT >
// class TagEssentialFilter : public AlignmentsFilter< ContainerT > { //TODO: remove this filter
// 	public:
// 		TagEssentialFilter( StrIDConverter& accessconv, Taxonomy* tax ) : seqid2taxid( accessconv ), taxinter( tax ) {};
// 
// 		void filter( ContainerT& recordset ) {
// 			if( ! recordset.empty() ) {
// 				typename ContainerT::iterator record_it = recordset.begin();
//         std::string seqid = extractFastaCommentField( (*record_it)->query_identifier, "gi" );
// 				try {
// 					std::list< int > lcadepths;
// 					TaxonID taxid = seqid2taxid[ seqid ];
// 					const TaxonNode* label_node = taxinter.getNode( taxid );
// 					int maxdepth = 0;
// 
// 					while( record_it != recordset.end() ) {
// 						const TaxonNode* node = taxinter.getNode( (*record_it)->getNode()->data->reference_taxid );
// 						const TaxonNode* pair_lca = taxinter.getLCA( node, label_node );
// 
// 						int depth = pair_lca->data->root_pathlength;
// 
// 						lcadepths.push_back( depth );
// 
// 						if( ! (*record_it)->isFiltered() ) { //obeys former filtering for cut calculation
// 							maxdepth = std::max( depth, maxdepth );
// 						}
// 
// 						++record_it;
// 					}
// 
// 					std::list< int >::iterator depth_it = lcadepths.begin();
// 					record_it = recordset.begin();
// 					while( record_it != recordset.end() ) {
// 						if( *depth_it < maxdepth ) {
// 							(*record_it)->raw_line->append( default_field_separator + "excess" );
// 						} else {
// 							(*record_it)->raw_line->append( default_field_separator + "essential" );
// 						}
// 						record_it++;
// 						depth_it++;
// 					}
// 				} catch ( std::out_of_range e ) { //mask all records
// 					std::cerr << "TagEssentialFilter: Could not map sequence identifier " << seqid << " to TaxID";
// 					std::cerr << ", skipping all records for record set..." << std::endl;
// 					while( record_it != recordset.end() ) {
// 						(*record_it)->filterOut();
// 						++record_it;
// 					}
// 				}
// 			}
// 		}
// 	private:
// 		static const std::string description;
// 		StrIDConverter& seqid2taxid;
// 		TaxonomyInterface taxinter;
// };

// template< typename ContainerT >
// const std::string TagEssentialFilter< ContainerT >::description = "TagEssentialFilter";



#endif // alignmentsfilter_hh_
