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

#ifndef taxonpredictionmodel_hh_
#define taxonpredictionmodel_hh_

#include <cmath>
#include <boost/tuple/tuple.hpp>
#include <ostream>
#include "types.hh"
#include "alignmentrecord.hh"
#include "taxonomyinterface.hh"
#include "alignmentsfilter.hh"
#include "predictionrecord.hh"
#include "constants.hh"


inline float distWeight( const int dist, const float gamma = 1.2, const int max_dist = default_ranks.size() ) {
	return 1.0 - std::pow( float (1.0/max_dist)*dist, gamma );
}



template< typename ContainerT >
class TaxonPredictionModel {
	public:
		TaxonPredictionModel( const Taxonomy* tax ) : taxinter_( tax ), root_( taxinter_.getRoot() ) {};
		virtual void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) = 0;
	protected:
		void initPredictionRecord( ContainerT& recordset, PredictionRecord& prec ) {
			prec.initialize( recordset.front()->getQueryIdentifier(), recordset.front()->getQueryLength() );
		};
		const TaxonNode* unclassified() { return root_; };
		void setUnclassified( PredictionRecord& prec ) {
			prec.setNodePoint( root_, 0 );
		}
		TaxonomyInterface taxinter_;
		const TaxonNode* root_;
};



//only for debugging!
template< typename ContainerT >
class DummyPredictionModel : public TaxonPredictionModel< ContainerT > { //TODO: include winscore
	public:
		DummyPredictionModel( const Taxonomy* tax ) : TaxonPredictionModel< ContainerT >( tax ) {}

		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
};



template< typename ContainerT >
class LCASimplePredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		LCASimplePredictionModel( const Taxonomy* tax ) : TaxonPredictionModel< ContainerT >( tax ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			std::list< const TaxonNode* > refnodes;
			records2Nodes( recordset, refnodes );
			if( ! refnodes.empty() ) {
				prec.setNodePoint( this->taxinter_.getLCA( refnodes ) );
				return;
			}
			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
};



template< typename ContainerT >
class MeganLCAPredictionModel : public TaxonPredictionModel< ContainerT > { //TODO: include winscore
	public:
		MeganLCAPredictionModel( const Taxonomy* tax, bool iuc, const float toppercent, const float minscore, const int minsupport, const double maxevalue, const float winscore = 0.0) : TaxonPredictionModel< ContainerT >( tax ),
			msp( minsupport ), wsc( winscore ), ms_me_toppercent( minscore, maxevalue, toppercent ), ignore_unclassified( iuc ) {};
		
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			ms_me_toppercent.filter( recordset );
			std::list< const TaxonNode* > refnodes;
			
			typename ContainerT::iterator rec_it = firstUnmaskedIter( recordset );
			if ( rec_it != recordset.end() ) {
				large_unsigned_int qrstart = (*rec_it)->getQueryStart();
				large_unsigned_int qrstop = (*rec_it)->getQueryStop();
// 				large_unsigned_int qlength = (*rec_it)->getQueryLength();
				
				// determine range of query used
				large_unsigned_int n = 0;
				if( qrstart > qrstop ) std::swap( qrstart, qrstop ); //TODO: is this allowed? should we through an error/exception

				refnodes.push_back( (*rec_it)->getReferenceNode() );
				++n;
				++rec_it;
				
				for ( ;rec_it != recordset.end(); ++rec_it ) { 
					if ( ! (*rec_it)->isFiltered() ) {
						if( (*rec_it)->getQueryStart() <= (*rec_it)->getQueryStop() ) { //TODO: is this allowed? should we through an error/exception
							qrstart = std::min( (*rec_it)->getQueryStart(), qrstart );
							qrstop = std::max( (*rec_it)->getQueryStop(), qrstop );
						} else {
							qrstart = std::min( (*rec_it)->getQueryStop(), qrstart );
							qrstop = std::max( (*rec_it)->getQueryStart(), qrstop );
						}
						refnodes.push_back( (*rec_it)->getReferenceNode() );
						++n;
					}
				}
				
				// set region coverage on query
				prec.setQueryFeatureBegin( qrstart );
				prec.setQueryFeatureEnd( qrstop );
			}
			
			if( ignore_unclassified ) { //TODO: avoid extra pass
				std::list< const TaxonNode* >::iterator it = refnodes.begin();
				while( it != refnodes.end() ) {
					if( (*it)->data->is_unclassified ) {
						it = refnodes.erase( it );
					} else {
						++it;
					}
				}
			}

			if( refnodes.size() >= msp ) {
				prec.setNodePoint( this->taxinter_.getLCA( refnodes ) );
				return;
			}

			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
	private:
		unsigned int msp;
		float wsc;
		MinScoreMaxEvalueTopPercentFilter< ContainerT > ms_me_toppercent;
		bool ignore_unclassified;
};



template< typename ContainerT > //makes use of incomplete knowledge
class ICMeganLCAPredictionModel : public TaxonPredictionModel< ContainerT > { //TODO: include winscore
	public:
		ICMeganLCAPredictionModel( const Taxonomy* tax, const float toppercent, const float minscore, const int minsupport, const double maxevalue, const float winscore = 0.0) : TaxonPredictionModel< ContainerT >( tax ),
		msp( minsupport ), wsc( winscore ), ms_me_toppercent( minscore, maxevalue, toppercent ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			ms_me_toppercent.filter( recordset );
			std::list< const TaxonNode* > refnodes;
			records2Nodes( recordset, refnodes );
			if( refnodes.size() >= msp ) {
				prec.setNodePoint( this->taxinter_.getICLCA( refnodes ) );
				return;
			}

			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
	private:
		unsigned int msp;
		float wsc;
		MinScoreMaxEvalueTopPercentFilter< ContainerT > ms_me_toppercent;
};



template< typename ContainerT >
class NBestLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		NBestLCAPredictionModel( const Taxonomy* tax, const int n = 1) : TaxonPredictionModel< ContainerT >( tax ), findnbest( n ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			findnbest.filter( recordset );
			std::list< const TaxonNode* > refnodes;
			records2Nodes( recordset, refnodes );
			if( ! refnodes.empty() ) {
				prec.setNodePoint( this->taxinter_.getLCA( refnodes ) );
				return;
			}
			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
	private:
		NumBestBitscoreFilter< ContainerT > findnbest;
};



template< typename ContainerT > //TODO: requires Container of AlignmentRecordTaxonomy
class ClosestNodePredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		ClosestNodePredictionModel( const Taxonomy* tax, StrIDConverter* sidc ) : TaxonPredictionModel< ContainerT >( tax ), seqid2taxid( *sidc ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			
			// determine query label node
			const TaxonNode* qnode = this->taxinter_.getNode( seqid2taxid[ extractFastaCommentField( recordset.front()->getQueryIdentifier(), "gi" ) ] );

			// determine closest node in LCA tree
			const TaxonNode* pairlca = this->root_;
			for( typename ContainerT::iterator it = recordset.begin(); it != recordset.end(); ++it ) {
				if( ! (*it)->isFiltered() ) {
					const TaxonNode* tmpnode = (*it)->getReferenceNode();
					tmpnode = this->taxinter_.getLCA( tmpnode, qnode );
					if( tmpnode->data->root_pathlength > pairlca->data->root_pathlength ) {
						pairlca = tmpnode;
					}
				}
			}
			
			prec.setNodePoint( pairlca );
		};
	private:
		StrIDConverter& seqid2taxid;
};


// predicts the node that is closest to the query taxon, a parent of the best scoring and below-equals the overall lca!
template< typename ContainerT > //TODO: requires Container of AlignmentRecordTaxonomy
class CorrectionPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		CorrectionPredictionModel( const Taxonomy* tax, StrIDConverter* sidc ) : TaxonPredictionModel< ContainerT >( tax ), seqid2taxid( *sidc ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			
			// determine query label node
			const TaxonNode* qnode = this->taxinter_.getNode( seqid2taxid[ extractFastaCommentField( recordset.front()->getQueryIdentifier(), "gi" ) ] );

			// determine best scoring alignment and corresponding node in taxonomy
			findbestscoring.filter( recordset );
			std::list< const TaxonNode* > nodes;
			records2Nodes( findbestscoring.getBests(), nodes );

// 			std::cerr << "Total number of nodes: " << recordset.size() << std::endl;
// 			std::cerr << "Number of best nodes: " << nodes.size() << std::endl;
// 			std::cerr << "Query is: " << qnode->data->annotation->name << " (" << qnode->data->taxid << ")" << std::endl;

			const TaxonNode* lowest_lca;
			unsigned int longest_pathlength = 0;
			if( ! nodes.empty() ) {
				while( ! nodes.empty() ) {
					const TaxonNode* tmp = this->taxinter_.getLCA( nodes.back(), qnode );
					if( tmp->data->root_pathlength >= longest_pathlength ) {
						longest_pathlength = tmp->data->root_pathlength;
						lowest_lca = tmp;
					}
					nodes.pop_back();
				}

				// make upper range check
				if( lowest_lca ) {
					records2Nodes( recordset, nodes );
					const TaxonNode* lca_all = this->taxinter_.getLCA( nodes );
					
					if( this->taxinter_.isParentOf( lowest_lca, lca_all ) ) {
						prec.setNodeRange( this->root_, lca_all ); //TODO: wrong order?
						return;
					}
					prec.setNodeRange( this->root_, lowest_lca ); //TODO: wrong order?
					return;
				}
			}

			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
	private:
		MaxBitscoreAlignmentFilter< ContainerT > findbestscoring;
		StrIDConverter& seqid2taxid;
};



template< typename ContainerT > //TODO: requires Container of AlignmentRecordTaxonomy
class ICCorrectionPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		ICCorrectionPredictionModel( const Taxonomy* tax, StrIDConverter* sidc ) : TaxonPredictionModel< ContainerT >( tax ), seqid2taxid( *sidc ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);

			if( ! recordset.empty() ) {
				std::list< const TaxonNode* > nodes;
				std::list< std::pair< float, const TaxonNode* > > unclassified_nodes;
// 				std::list< float > unclassified_bitscores;
				std::list< const TaxonNode* > best_nodes;
// 				std::pair< const TaxonNode*, float > best_bs( NULL, 0.0 );
				float best_bs = 0.0;

				// determine query label node
				const TaxonNode* qnode = this->taxinter_.getNode( seqid2taxid[ extractFastaCommentField( recordset.front()->getQueryIdentifier(), "gi" ) ] );

				// sort out unclassified nodes with corresponding bitscore
				// get best scoring node in normal set
				for( typename ContainerT::iterator it = recordset.begin(); it != recordset.end(); ++it ) {
					if( ! (*it)->isFiltered() ) {
						const TaxonNode* tmpnode = (*it)->getReferenceNode();
						if( tmpnode->data->is_unclassified ) {
							tmpnode = this->taxinter_.mapUnclassified( tmpnode );
							unclassified_nodes.push_back( std::make_pair( (*it)->getScore(), tmpnode ) );
// 							unclassified_bitscores.push_back( (*it)->getScore() );
						} else {
							if( (*it)->getScore() > best_bs ) {
								best_nodes.clear();
								best_nodes.push_back( tmpnode );
								best_bs = (*it)->getScore();
							} else {
								if( (*it)->getScore() == best_bs ) {
									best_nodes.push_back( tmpnode );
								}
							}
						}
						nodes.push_back( tmpnode );
					}
				}
				const TaxonNode* lca_all = this->taxinter_.getLCA( nodes );

// 				std::cerr << "Total number of nodes: " << recordset.size() << std::endl;
// 				std::cerr << "Number of unclassified nodes: " << unclassified_nodes.size() << std::endl;
// 				std::cerr << "Number of best nodes: " << best_nodes.size() << std::endl;
// 				std::cerr << "Query is: " << qnode->data->annotation->name << " (" << qnode->data->taxid << ")" << std::endl;

				// remove unclassified nodes with bitscore less than maximum in normal set
				std::list< std::pair< float, const TaxonNode* > >::iterator n_it = unclassified_nodes.begin();
// 				std::list< float >::iterator b_it = unclassified_bitscores.begin();
				while( n_it != unclassified_nodes.end() ) {
					if( n_it->first < best_bs ) {
						n_it = unclassified_nodes.erase( n_it );
					} else {
						++n_it;
					}
				}

				if( ! unclassified_nodes.empty() ) {
					// sort nodes and find out lowest possible prediction
					unclassified_nodes.sort( std::greater< std::pair< float, const TaxonNode* > >() );
					const TaxonNode* lcc = this->root_;
					for( std::list< std::pair< float, const TaxonNode* > >::iterator n_it = unclassified_nodes.begin(); n_it != unclassified_nodes.end(); ++n_it ) {
// 						std::cerr << "Bitscore: " << n_it->first << " for node " << n_it->second->data->annotation->name << std::endl;
						lcc = this->taxinter_.getLCC( n_it->second, lcc );
						best_nodes.push_back( lcc );
					}
// 					best_nodes.push_back( this->taxinter_.getLCC( unclassified_nodes ) );
				}

// 				std::cerr << "Number of unclassified nodes with better score: " << unclassified_nodes.size() << std::endl;

				unsigned int max_depth = 0;
				const TaxonNode* correction = NULL;
				while( ! best_nodes.empty() ) {
					const TaxonNode* tmpnode = this->taxinter_.getLCA( best_nodes.back(), qnode );
					if( tmpnode->data->root_pathlength >= max_depth ) {
						correction = tmpnode;
						max_depth = correction->data->root_pathlength;
					}
					best_nodes.pop_back();
				}
				
				

				// make upper range check
				if( correction ) {
					if( this->taxinter_.isParentOf( correction, lca_all ) ) {
						prec.setNodePoint( lca_all );
						return;
					}
					prec.setNodePoint( correction );
					return;
				}
// 				const TaxonNode* lcc_unclassified = ;
// 				const TaxonNode* qlca_normal = this->taxinter_.getLCA( best_bs.first, qnode );
// 				const TaxonNode* qlca_unclassified = this->taxinter_.getLCA( lcc_unclassified, qnode );
//
// 				if( qlca_unclassified->data->root_pathlength > qlca_normal->data->root_pathlength ) {
// 					return qlca_unclassified;
// 				} else {
// 					return qlca_normal;
// 				}
			}

			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
	private:
		StrIDConverter& seqid2taxid;
};



template< typename ContainerT > //TODO: requires Container of AlignmentRecordTaxonomy
class QueryBestLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		QueryBestLCAPredictionModel( const Taxonomy* tax, StrIDConverter* sidc ) : TaxonPredictionModel< ContainerT >( tax ), seqid2taxid( *sidc ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			
			// determine query label node
			const TaxonNode* qnode = this->taxinter_.getNode( seqid2taxid[ extractFastaCommentField( recordset.front()->getQueryIdentifier(), "gi" ) ] );

			// determine best scoring alignment and corresponding node in taxonomy
			findbestscoring.filter( recordset );
			std::list< const TaxonNode* > nodes;
			records2Nodes( findbestscoring.getBests(), nodes );
			const TaxonNode* lowest_lca;
			unsigned int longest_pathlength = 0;
			if( ! nodes.empty() ) {
				while( ! nodes.empty() ) {
					const TaxonNode* tmp = this->taxinter_.getLCA( nodes.back(), qnode );
					if( tmp->data->root_pathlength >= longest_pathlength ) {
						longest_pathlength = tmp->data->root_pathlength;
						lowest_lca = tmp;
					}
					nodes.pop_back();
				}

				// return CA of query and best scoring alignment node
				if( lowest_lca ) {
					
					prec.setNodePoint( lowest_lca );
					return;
				}
			}

			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
	private:
		MaxBitscoreAlignmentFilter< ContainerT > findbestscoring;
		StrIDConverter& seqid2taxid;
};



template< typename ContainerT > //TODO: requires Container of AlignmentRecordTaxonomy
class ExtLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		ExtLCAPredictionModel( const Taxonomy* tax, const float thresh, const float gam ) : TaxonPredictionModel< ContainerT >( tax ), threshold( thresh ), gamma( gam ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);

			findbestscoring.filter( recordset );
			AlignmentRecordTaxonomy* bestaln = findbestscoring.getBest();
			const TaxonNode* best_node = bestaln->getReferenceNode();
			float &best_bs = bestaln->getScore();
			std::list< const TaxonNode* > refnodes;

			for( typename ContainerT::iterator it = recordset.begin(); it != recordset.end(); ++it ) {
				if( ! (*it)->isFiltered() ) {
					const TaxonNode* tmpnode = (*it)->getReferenceNode();
					int a, b, c;
					boost::tie( a, b, c ) = this->taxinter_.getInterDistances( tmpnode, best_node );
					float weight = distWeight( c, gamma );
					float bs_rel = (*it)->getScore() / best_bs;
					// 					std::cerr << weight << " * " << bs_rel << " = " << weight*bs_rel << std::endl;
					if( weight * bs_rel >= threshold ) {
						refnodes.push_back( tmpnode );
					}
				}
			}

			if( ! refnodes.empty() ) {
				prec.setNodePoint( this->taxinter_.getLCA( refnodes ) );
				return;
			}

			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};

	private:

		const float gamma;
		const float threshold;
		MaxBitscoreAlignmentFilter< ContainerT > findbestscoring;
};


template< typename ContainerT > //TODO: requires Container of AlignmentRecordTaxonomy
class TestExtLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		TestExtLCAPredictionModel( const Taxonomy* tax, const float thresh, const float dlb , const float b ) : TaxonPredictionModel< ContainerT >( tax ), threshold( thresh ), distweight_factor( ( 1.0 - dlb )/float( default_ranks.size() ) ), beta( b ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);

			sort.filter( recordset );
			std::list< const TaxonNode* > bestnodes;
			std::list< const TaxonNode* > refnodes;
			int a, b, c;

			float best_bs;
			typename ContainerT::iterator it = recordset.begin();
			for( ; it != recordset.end(); ++it ) { //get best valid alignment
				if( ! (*it)->isFiltered() ) {
					best_bs = (*it)->getScore();
					const TaxonNode* tmpnode = (*it)->getReferenceNode();
					bestnodes.push_back( tmpnode );
					refnodes.push_back( tmpnode );
					break;
				}
			}

			for( ; it != recordset.end() && (*it)->getScore() == best_bs; ++it ) { //ties will be included, too
				if( ! (*it)->isFiltered() ) {
					const TaxonNode* tmpnode = (*it)->getReferenceNode();
					bestnodes.push_back( tmpnode );
					refnodes.push_back( tmpnode );
				}
			}

			for( ; it != recordset.end(); ++it ) {
				if( ! (*it)->isFiltered() ) {
					const TaxonNode* tmpnode = (*it)->getReferenceNode();
					float bs_rel = (*it)->getScore() / best_bs;
					float norm_dist = getNormDist( tmpnode, bestnodes );
// 					std::cerr << "bs_rel: " << bs_rel << " norm_dist: " << norm_dist << " combined: " << combinedFeature( bs_rel, norm_dist ) << std::endl;
					if( combinedFeature( bs_rel, norm_dist )  >= threshold ) {
						refnodes.push_back( tmpnode );
					}
				}
			}

			if( ! refnodes.empty() ) {
				prec.setNodePoint( this->taxinter_.getLCA( refnodes ) );
				return;
			}

			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};

	private:
		float getNormDist( const TaxonNode* n, std::list< const TaxonNode* >& bestnds ) {
			int a, b, c;
			int c_sum = 0;
			for( std::list< const TaxonNode* >::iterator it = bestnds.begin(); it != bestnds.end(); ++it ) {
				boost::tie( a, b, c ) = this->taxinter_.getInterDistances( n, *it );
				c_sum += c;
			}
			return 1.0 - distweight_factor * c_sum/float( bestnds.size() );
		};

		float combinedFeature( float feat_a, float feat_b ) {
			return std::pow( feat_a, 1.0 - beta ) * std::pow( feat_b, beta );
		};

		const float threshold;
		const float beta;
		const float distweight_factor;
		SortByBitscoreFilter< ContainerT > sort;
};


template< typename ContainerT >
class ExtendedLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		ExtendedLCAPredictionModel( const Taxonomy* tax, const float t1, const float t2 ) : TaxonPredictionModel< ContainerT >( tax ), cleanse( tax, t1, t2 ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);

			cleanse.filter( recordset );

			std::list< const TaxonNode* > refnodes;
			records2Nodes( recordset, refnodes );

			if( ! refnodes.empty() ) {
				prec.setNodePoint( this->taxinter_.getLCA( refnodes ) );
				return;
			}

			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};

	private:

		CleanseFDistAlignmentFilter< ContainerT > cleanse;
};



template< typename ContainerT >
class TopPercentOutlierLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		TopPercentOutlierLCAPredictionModel( const Taxonomy* tax, const float t1 = 1.0, const float t2 = 0.5, const unsigned int m = 3 ) : TaxonPredictionModel< ContainerT >( tax ), outlier( tax, t2, m ), min_top( 0.0, t1 ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);

			min_top.filter( recordset );
			outlier.filter( recordset );

			std::list< const TaxonNode* > refnodes;
			records2Nodes( recordset, refnodes );

			if( ! refnodes.empty() ) {
				prec.setNodePoint( this->taxinter_.getLCA( refnodes ) );
				return;
			}

			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};

	private:
		OutlierDetectionAlignmentFilter< ContainerT > outlier;
		MinScoreTopPercentFilter< ContainerT > min_top;
};



// following classes are experimental

class PredictionsParser {
	public:
		PredictionsParser( std::istream& istr ) : source( istr ), not_empty( true ) {
			parseLine();
		};

		PredictionsParser( const std::string& filename ) : filehandle( filename.c_str() ), source( filehandle ), not_empty( true ) {
			parseLine();
		};

		boost::tuple< std::string, unsigned int > getNext() {
			boost::tuple< std::string, unsigned int> ret_type = boost::make_tuple( qid, taxid_pred );
			parseLine();
			return ret_type;
		};

		bool notEmpty() {
			return not_empty;
		};

	private:
		void parseLine() {
			std::list< std::string > fields;
			std::list< std::string >::iterator field_it;

			while( std::getline( source, line ) ) {
				if( ignoreLine( line ) ) {
					std::cerr << "PredictionsParser: Ignoring comment line" << std::endl;
					continue;
				}

				tokenizeSingleCharDelim( line, fields, default_field_separator, 2 );
				field_it = fields.begin();

				try {
					qid = *field_it++;
					taxid_pred = boost::lexical_cast< unsigned int >( *field_it );
					return;
				} catch ( boost::bad_lexical_cast e ) {
					std::cerr << "PredictionsParser: Could not parse line '" << line << '\'' << std::endl;
					fields.clear();
				}
			}
			not_empty = false;
		}

    std::ifstream filehandle;
		std::istream& source;
		std::string qid;
		std::string line;
		unsigned int taxid_pred;
		bool not_empty;
};



#endif // taxonpredictionmodel_hh_
