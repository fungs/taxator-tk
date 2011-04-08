/*
The taxatorTK predicts the taxon for DNA sequences based on sequence alignment.

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

#ifndef taxonprediction_hh_
#define taxonprediction_hh_

#include <cmath>
#include <boost/tuple/tuple.hpp>
#include "types.hh"
#include "alignmentrecord.hh"
#include "taxontree.hh"
#include "alignmentsfilter.hh"
#include "constants.hh"


inline float distWeight( const int dist, const float gamma = 1.2, const int max_dist = default_ranks.size() ) {
	return 1.0 - std::pow( float (1.0/max_dist)*dist, gamma );
};



template< typename ContainerT >
class TaxonPredictionModel {
	public:
		TaxonPredictionModel( Taxonomy* tax ) : taxinter( tax ), taxonomy( tax ), root( taxinter.getRoot() ) {};
		virtual TaxonNode* predict( ContainerT& recordset ) = 0;
		TaxonNode* unclassified() { return root; };
	protected:
		TaxonomyInterface taxinter;
		Taxonomy* taxonomy;
		TaxonNode* root;
};



template< typename ContainerT >
class DummyPredictionModel : public TaxonPredictionModel< ContainerT > { //TODO: include winscore
	public:
		DummyPredictionModel( Taxonomy* tax ) : TaxonPredictionModel< ContainerT >( tax ) {}

		TaxonNode* predict( ContainerT& recordset ) {
			return this->unclassified();
		};
};



template< typename ContainerT >
class LCASimplePredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		LCASimplePredictionModel( Taxonomy* tax ) : TaxonPredictionModel< ContainerT >( tax ) {};
		TaxonNode* predict( ContainerT& recordset ) {
			std::list< TaxonNode* > refnodes;
			records2Nodes( recordset, &this->taxinter, refnodes );
			if( ! refnodes.empty() ) {
				return this->taxinter.getLCA( refnodes );
			}
			return this->unclassified();
		};
};



template< typename ContainerT >
class MeganLCAPredictionModel : public TaxonPredictionModel< ContainerT > { //TODO: include winscore
	public:
		MeganLCAPredictionModel( Taxonomy* tax, bool iuc, const float toppercent, const float minscore = 0.0, const int minsupport = 1, const float winscore = 0.0) : TaxonPredictionModel< ContainerT >( tax ),
			msp( minsupport ), wsc( winscore ), minscore_toppercent( minscore, toppercent ), ignore_unclassified( iuc ) {};
		TaxonNode* predict( ContainerT& recordset ) {
			minscore_toppercent.filter( recordset );
			std::list< TaxonNode* > refnodes;
			records2Nodes( recordset, &this->taxinter, refnodes );
			if( ignore_unclassified ) {
				std::list< TaxonNode* >::iterator it = refnodes.begin();
				while( it != refnodes.end() ) {
					if( (*it)->data->is_unclassified ) {
						it = refnodes.erase( it );
					} else {
						++it;
					}
				}
			}

			if( refnodes.size() >= msp ) {
				return this->taxinter.getLCA( refnodes );
			}

			return this->unclassified();
		};
	private:
		unsigned int msp;
		float wsc;
		MinScoreTopPercentFilter< ContainerT > minscore_toppercent;
		bool ignore_unclassified;
};



template< typename ContainerT > //makes use of incomplete knowledge
class ICMeganLCAPredictionModel : public TaxonPredictionModel< ContainerT > { //TODO: include winscore
	public:
		ICMeganLCAPredictionModel( Taxonomy* tax, const float toppercent, const float minscore = 0.0, const int minsupport = 1, const float winscore = 0.0) : TaxonPredictionModel< ContainerT >( tax ),
		msp( minsupport ), wsc( winscore ), minscore_toppercent( minscore, toppercent ) {};
		TaxonNode* predict( ContainerT& recordset ) {
			minscore_toppercent.filter( recordset );
			std::list< TaxonNode* > refnodes;
			records2Nodes( recordset, &this->taxinter, refnodes );
			if( refnodes.size() >= msp ) {
				return this->taxinter.getICLCA( refnodes );
			}

			return this->unclassified();
		};
	private:
		unsigned int msp;
		float wsc;
		MinScoreTopPercentFilter< ContainerT > minscore_toppercent;
};



template< typename ContainerT >
class NBestLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		NBestLCAPredictionModel( Taxonomy* tax, const int n = 1) : TaxonPredictionModel< ContainerT >( tax ), findnbest( n ) {};
		TaxonNode* predict( ContainerT& recordset ) {
			findnbest.filter( recordset );
			std::list< TaxonNode* > refnodes;
			records2Nodes( recordset, &this->taxinter, refnodes );
			if( ! refnodes.empty() ) {
				return this->taxinter.getLCA( refnodes );
			}
			return this->unclassified();
		};
	private:
		NumBestBitscoreFilter< ContainerT > findnbest;
};



template< typename ContainerT >
class ClosestNodePredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		ClosestNodePredictionModel( Taxonomy* tax, StrIDConverter* sidc ) : TaxonPredictionModel< ContainerT >( tax ), seqid2taxid( *sidc ) {};
		TaxonNode* predict( ContainerT& recordset ) {
			AlignmentRecord* firstrec = recordset.front();

			// determine query label node
			TaxonNode* qnode = this->taxinter.getNode( seqid2taxid[ extractFastaCommentField( firstrec->query_identifier, "gi" ) ] );

			// determine closest node in LCA tree
			TaxonNode* pairlca = this->root;
			for( typename ContainerT::iterator it = recordset.begin(); it != recordset.end(); ++it ) {
				if( ! (*it)->mask ) {
					TaxonNode* tmpnode = this->taxinter.getNode( (*it)->reference_taxid );
					tmpnode = this->taxinter.getLCA( tmpnode, qnode );
					if( tmpnode->data->root_pathlength > pairlca->data->root_pathlength ) {
						pairlca = tmpnode;
					}
				}
			}
			return pairlca;
		};
	private:
		StrIDConverter& seqid2taxid;
};


// predicts the node that is closest to the query taxon, a parent of the best scoring and below-equals the overall lca!
template< typename ContainerT >
class CorrectionPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		CorrectionPredictionModel( Taxonomy* tax, StrIDConverter* sidc ) : TaxonPredictionModel< ContainerT >( tax ), seqid2taxid( *sidc ) {};
		TaxonNode* predict( ContainerT& recordset ) {
			AlignmentRecord* firstrec = recordset.front();

			// determine query label node
			TaxonNode* qnode = this->taxinter.getNode( seqid2taxid[ extractFastaCommentField( firstrec->query_identifier, "gi" ) ] );

			// determine best scoring alignment and corresponding node in taxonomy
			findbestscoring.filter( recordset );
			std::list< TaxonNode* > nodes;
			records2Nodes( findbestscoring.getBests(), &this->taxinter, nodes );

// 			std::cerr << "Total number of nodes: " << recordset.size() << std::endl;
// 			std::cerr << "Number of best nodes: " << nodes.size() << std::endl;
// 			std::cerr << "Query is: " << qnode->data->annotation->name << " (" << qnode->data->taxid << ")" << std::endl;

			TaxonNode* lowest_lca;
			unsigned int longest_pathlength = 0;
			if( ! nodes.empty() ) {
				while( ! nodes.empty() ) {
					TaxonNode* tmp = this->taxinter.getLCA( nodes.back(), qnode );
					if( tmp->data->root_pathlength >= longest_pathlength ) {
						longest_pathlength = tmp->data->root_pathlength;
						lowest_lca = tmp;
					}
					nodes.pop_back();
				}

				// make upper range check
				if( lowest_lca ) {
					records2Nodes( recordset, &this->taxinter, nodes );
					TaxonNode* lca_all = this->taxinter.getLCA( nodes );

					if( this->taxinter.isParentOf( lowest_lca, lca_all ) ) {
						return lca_all;
					}
					return lowest_lca;
				}
			}

			return this->unclassified();
		};
	private:
		MaxBitscoreAlignmentFilter< ContainerT > findbestscoring;
		StrIDConverter& seqid2taxid;
};



template< typename ContainerT >
class ICCorrectionPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		ICCorrectionPredictionModel( Taxonomy* tax, StrIDConverter* sidc ) : TaxonPredictionModel< ContainerT >( tax ), seqid2taxid( *sidc ) {};
		TaxonNode* predict( ContainerT& recordset ) {

			if( ! recordset.empty() ) {
				std::list< TaxonNode* > nodes;
				std::list< std::pair< float, TaxonNode* > > unclassified_nodes;
// 				std::list< float > unclassified_bitscores;
				std::list< TaxonNode* > best_nodes;
// 				std::pair< TaxonNode*, float > best_bs( NULL, 0.0 );
				float best_bs = 0.0;

				// determine query label node
				AlignmentRecord* firstrec = recordset.front();
				TaxonNode* qnode = this->taxinter.getNode( seqid2taxid[ extractFastaCommentField( firstrec->query_identifier, "gi" ) ] );

				// sort out unclassified nodes with corresponding bitscore
				// get best scoring node in normal set
				for( typename ContainerT::iterator it = recordset.begin(); it != recordset.end(); ++it ) {
					if( ! (*it)->mask ) {
						TaxonNode* tmpnode = this->taxinter.getNode( (*it)->reference_taxid );
						if( tmpnode->data->is_unclassified ) {
							tmpnode = this->taxinter.mapUnclassified( tmpnode );
							unclassified_nodes.push_back( std::make_pair( (*it)->bitscore, tmpnode ) );
// 							unclassified_bitscores.push_back( (*it)->bitscore );
						} else {
							if( (*it)->bitscore > best_bs ) {
								best_nodes.clear();
								best_nodes.push_back( tmpnode );
								best_bs = (*it)->bitscore;
							} else {
								if( (*it)->bitscore == best_bs ) {
									best_nodes.push_back( tmpnode );
								}
							}
						}
						nodes.push_back( tmpnode );
					}
				}
				TaxonNode* lca_all = this->taxinter.getLCA( nodes );

// 				std::cerr << "Total number of nodes: " << recordset.size() << std::endl;
// 				std::cerr << "Number of unclassified nodes: " << unclassified_nodes.size() << std::endl;
// 				std::cerr << "Number of best nodes: " << best_nodes.size() << std::endl;
// 				std::cerr << "Query is: " << qnode->data->annotation->name << " (" << qnode->data->taxid << ")" << std::endl;

				// remove unclassified nodes with bitscore less than maximum in normal set
				std::list< std::pair< float, TaxonNode* > >::iterator n_it = unclassified_nodes.begin();
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
					unclassified_nodes.sort( std::greater< std::pair< float, TaxonNode* > >() );
					TaxonNode* lcc = this->root;
					for( std::list< std::pair< float, TaxonNode* > >::iterator n_it = unclassified_nodes.begin(); n_it != unclassified_nodes.end(); ++n_it ) {
// 						std::cerr << "Bitscore: " << n_it->first << " for node " << n_it->second->data->annotation->name << std::endl;
						lcc = this->taxinter.getLCC( n_it->second, lcc );
						best_nodes.push_back( lcc );
					}
// 					best_nodes.push_back( this->taxinter.getLCC( unclassified_nodes ) );
				}

// 				std::cerr << "Number of unclassified nodes with better score: " << unclassified_nodes.size() << std::endl;

				unsigned int max_depth = 0;
				TaxonNode* correction = NULL;
				while( ! best_nodes.empty() ) {
					TaxonNode* tmpnode = this->taxinter.getLCA( best_nodes.back(), qnode );
					if( tmpnode->data->root_pathlength >= max_depth ) {
						correction = tmpnode;
						max_depth = correction->data->root_pathlength;
					}
					best_nodes.pop_back();
				}

				// make upper range check
				if( correction ) {
					if( this->taxinter.isParentOf( correction, lca_all ) ) {
						return lca_all;
					}
					return correction;
				}
// 				TaxonNode* lcc_unclassified = ;
// 				TaxonNode* qlca_normal = this->taxinter.getLCA( best_bs.first, qnode );
// 				TaxonNode* qlca_unclassified = this->taxinter.getLCA( lcc_unclassified, qnode );
//
// 				if( qlca_unclassified->data->root_pathlength > qlca_normal->data->root_pathlength ) {
// 					return qlca_unclassified;
// 				} else {
// 					return qlca_normal;
// 				}
			}

			return this->unclassified();
		};
	private:
		StrIDConverter& seqid2taxid;
};



template< typename ContainerT >
class QueryBestLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		QueryBestLCAPredictionModel( Taxonomy* tax, StrIDConverter* sidc ) : TaxonPredictionModel< ContainerT >( tax ), seqid2taxid( *sidc ) {};
		TaxonNode* predict( ContainerT& recordset ) {
			AlignmentRecord* firstrec = recordset.front();

			// determine query label node
			TaxonNode* qnode = this->taxinter.getNode( seqid2taxid[ extractFastaCommentField( firstrec->query_identifier, "gi" ) ] );

			// determine best scoring alignment and corresponding node in taxonomy
			findbestscoring.filter( recordset );
			std::list< TaxonNode* > nodes;
			records2Nodes( findbestscoring.getBests(), &this->taxinter, nodes );
			TaxonNode* lowest_lca;
			unsigned int longest_pathlength = 0;
			if( ! nodes.empty() ) {
				while( ! nodes.empty() ) {
					TaxonNode* tmp = this->taxinter.getLCA( nodes.back(), qnode );
					if( tmp->data->root_pathlength >= longest_pathlength ) {
						longest_pathlength = tmp->data->root_pathlength;
						lowest_lca = tmp;
					}
					nodes.pop_back();
				}

				// return CA of query and best scoring alignment node
				if( lowest_lca ) {
					return lowest_lca;
				}
			}

			return this->unclassified();
		};
	private:
		MaxBitscoreAlignmentFilter< ContainerT > findbestscoring;
		StrIDConverter& seqid2taxid;
};



template< typename ContainerT >
class ExtLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		ExtLCAPredictionModel( Taxonomy* tax, const float thresh, const float gam ) : TaxonPredictionModel< ContainerT >( tax ), threshold( thresh ), gamma( gam ) {};
		TaxonNode* predict( ContainerT& recordset ) {

			findbestscoring.filter( recordset );
			AlignmentRecord* bestaln = findbestscoring.getBest();
			TaxonNode* best_node = this->taxinter.getNode( bestaln->reference_taxid );
			float &best_bs = bestaln->bitscore;
			std::list< TaxonNode* > refnodes;

			for( typename ContainerT::iterator it = recordset.begin(); it != recordset.end(); ++it ) {
				if( ! (*it)->mask ) {
					TaxonNode* tmpnode = this->taxinter.getNode( (*it)->reference_taxid );
					int a, b, c;
					boost::tie( a, b, c ) = this->taxinter.getInterDistances( tmpnode, best_node );
					float weight = distWeight( c, gamma );
					float bs_rel = (*it)->bitscore / best_bs;
					// 					std::cerr << weight << " * " << bs_rel << " = " << weight*bs_rel << std::endl;
					if( weight * bs_rel >= threshold ) {
						refnodes.push_back( tmpnode );
					}
				}
			}

			if( ! refnodes.empty() ) {
				return this->taxinter.getLCA( refnodes );
			}

			return this->unclassified();
		};

	private:

		const float gamma;
		const float threshold;
		MaxBitscoreAlignmentFilter< ContainerT > findbestscoring;
};


template< typename ContainerT >
class TestExtLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		TestExtLCAPredictionModel( Taxonomy* tax, const float thresh, const float dlb , const float b ) : TaxonPredictionModel< ContainerT >( tax ), threshold( thresh ), distweight_factor( ( 1.0 - dlb )/float( default_ranks.size() ) ), beta( b ) {};
		TaxonNode* predict( ContainerT& recordset ) {

			sort.filter( recordset );
			std::list< TaxonNode* > bestnodes;
// 			TaxonNode* best_node = this->taxinter.getNode( bestaln->reference_taxid );
			std::list< TaxonNode* > refnodes;
			int a, b, c;

			float best_bs;
			typename ContainerT::iterator it = recordset.begin();
			for( ; it != recordset.end(); ++it ) { //get best valid alignment
				if( ! (*it)->mask ) {
					best_bs = (*it)->bitscore;
					TaxonNode* tmpnode = this->taxinter.getNode( (*it)->reference_taxid );
					bestnodes.push_back( tmpnode );
					refnodes.push_back( tmpnode );
					break;
				}
			}

			for( ; it != recordset.end() && (*it)->bitscore == best_bs; ++it ) { //ties will be included, too
				if( ! (*it)->mask ) {
					TaxonNode* tmpnode = this->taxinter.getNode( (*it)->reference_taxid );
					bestnodes.push_back( tmpnode );
					refnodes.push_back( tmpnode );
				}
			}

			for( ; it != recordset.end(); ++it ) {
				if( ! (*it)->mask ) {
					TaxonNode* tmpnode = this->taxinter.getNode( (*it)->reference_taxid );
					float bs_rel = (*it)->bitscore / best_bs;
					float norm_dist = getNormDist( tmpnode, bestnodes );
// 					std::cerr << "bs_rel: " << bs_rel << " norm_dist: " << norm_dist << " combined: " << combinedFeature( bs_rel, norm_dist ) << std::endl;
					if( combinedFeature( bs_rel, norm_dist )  >= threshold ) {
						refnodes.push_back( tmpnode );
					}
				}
			}

			if( ! refnodes.empty() ) {
				return this->taxinter.getLCA( refnodes );
			}

			return this->unclassified();
		};

	private:
		float getNormDist( TaxonNode* n, std::list< TaxonNode* >& bestnds ) {
			int a, b, c;
			int c_sum = 0;
			for( std::list< TaxonNode* >::iterator it = bestnds.begin(); it != bestnds.end(); ++it ) {
				boost::tie( a, b, c ) = this->taxinter.getInterDistances( n, *it );
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
		ExtendedLCAPredictionModel( Taxonomy* tax, const float t1, const float t2 ) : TaxonPredictionModel< ContainerT >( tax ), cleanse( tax, t1, t2 ) {};
		TaxonNode* predict( ContainerT& recordset ) {

			cleanse.filter( recordset );

			std::list< TaxonNode* > refnodes;
			records2Nodes( recordset, &this->taxinter, refnodes );

			if( ! refnodes.empty() ) {
				return this->taxinter.getLCA( refnodes );
			}

			return this->unclassified();
		};

	private:

		CleanseFDistAlignmentFilter< ContainerT > cleanse;
};



template< typename ContainerT >
class TopPercentOutlierLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		TopPercentOutlierLCAPredictionModel( Taxonomy* tax, const float t1 = 1.0, const float t2 = 0.5, const unsigned int m = 3 ) : TaxonPredictionModel< ContainerT >( tax ), outlier( tax, t2, m ), min_top( 0.0, t1 ) {};
		TaxonNode* predict( ContainerT& recordset ) {

			min_top.filter( recordset );
			outlier.filter( recordset );

			std::list< TaxonNode* > refnodes;
			records2Nodes( recordset, &this->taxinter, refnodes );

			if( ! refnodes.empty() ) {
				return this->taxinter.getLCA( refnodes );
			}

			return this->unclassified();
		};

	private:
		OutlierDetectionAlignmentFilter< ContainerT > outlier;
		MinScoreTopPercentFilter< ContainerT > min_top;
};



// following classes are experimental


class AtomicTaxonPredictionModel {
	public:
		AtomicTaxonPredictionModel( Taxonomy* tax ) : taxinter( tax ), taxonomy( tax ) {};
		virtual void loadModel( const std::string& filename ) = 0;
		virtual TaxonNode* predict( AlignmentRecord* record ) = 0;
	protected:
		TaxonomyInterface taxinter;
		Taxonomy* taxonomy;
};



class MyTaxonPredictionModel : AtomicTaxonPredictionModel {
	public:
		MyTaxonPredictionModel( Taxonomy* tax ) : AtomicTaxonPredictionModel( tax ) {};
		void loadModel( const std::string& filename );
		TaxonNode* predict( AlignmentRecord* record );

	private:
		std::map< unsigned int, float > bitscore_cutoff;
};



class MyTaxonPredictionModel2 : AtomicTaxonPredictionModel {
	public:
		MyTaxonPredictionModel2( Taxonomy* tax ) : AtomicTaxonPredictionModel( tax ) {};
		void loadModel( const std::string& filename );
		TaxonNode* predict( AlignmentRecord* record );

	private:
		std::map< const std::string*, float > ranks_cutoff;
};



class MyTaxonPredictionModel3 : AtomicTaxonPredictionModel {
	public:
		MyTaxonPredictionModel3( Taxonomy* tax ) : AtomicTaxonPredictionModel( tax ), max_bitscore( 0.0 ) {};
		void loadModel( const std::string& filename );
		void setMaxBitscore( const float bitscore );
		TaxonNode* predict( AlignmentRecord* record );

	private:
		std::map< const std::string*, float > ranks_cutoff;
		float max_bitscore;
};



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

				tokenizeSingleCharDelim( line, fields, FSEP, 2 );
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



#endif // taxonprediction_hh_
