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

#ifndef taxonomyinterface_hh_
#define taxonomyinterface_hh_

#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <iterator>
#include "types.hh"
#include "taxontree.hh"

class TaxonomyInterface {
	public:
		TaxonomyInterface( const Taxonomy* taxtree ) : tax( taxtree ) {} //TODO: make all inline or define in header
		
		TaxonomyInterface( const TaxonomyInterface& taxinter ) : tax( taxinter.tax ) {}

		const TaxonNode* getNode( const TaxonID taxid ) const;
		const TaxonNode* getRoot() const;
		small_unsigned_int getMaxDepth() { return tax->max_depth_; }

		const std::string& getRank( const TaxonNode* node ) const;
		const std::string& getRank( const TaxonID taxid ) const;
		
		const std::string& getRankInternal( const std::string& rank ) const;
		
		const std::string& getName( const TaxonNode* node ) const;
		const std::string& getName( const TaxonID taxid ) const;

		bool isParentOf( const TaxonNode* A, const TaxonNode* B ) const;
		bool isParentOf( const TaxonID A_taxid, const TaxonID B_taxid ) const;

		// if you know that one node is probably closer to the LCA than the other, set it to A (less traversal steps)
		const TaxonNode* getLCA( const TaxonNode* A, const TaxonNode* B ) const;
		const TaxonNode* getLCA( const TaxonID A_taxid, const TaxonID B_taxid ) const;

		const TaxonNode* getLCC( const TaxonNode* A, const TaxonNode* B ) const;
		const TaxonNode* getLCC( const TaxonID A_taxid, const TaxonID B_taxid ) const;



		template < typename ContainerT >
		const TaxonNode* getLCA( const ContainerT& nodescontainer ) const { //TODO: handle taxid not found
			if( nodescontainer.empty() ) {
				return NULL;
			}

			typename ContainerT::const_iterator node_it = nodescontainer.begin();
			const TaxonNode* tmplca = *node_it++;
			while( node_it != nodescontainer.end() ) {
				tmplca = getLCA( tmplca, *node_it );
				++node_it;
			}
			return tmplca;
		}



		template < typename ContainerT >
		const TaxonNode* getLCC( const ContainerT& nodescontainer ) const { //TODO: handle taxid not found
			if( nodescontainer.empty() ) {
				return NULL;
			}
			// we need to start with lowest overall concept
			std::pair< const TaxonNode*, small_unsigned_int > lowest( NULL, 0 );
			for( typename ContainerT::const_iterator node_it = nodescontainer.begin(); node_it != nodescontainer.end(); ++node_it ) {
				if( (*node_it)->data->root_pathlength >= lowest.second ) {
					lowest.first = *node_it;
					lowest.second = (*node_it)->data->root_pathlength;
				}
			}
			const TaxonNode* tmplcc = lowest.first;
			for( typename ContainerT::const_iterator node_it = nodescontainer.begin(); node_it != nodescontainer.end(); ++node_it ) {
				tmplcc = getLCC( tmplcc, *node_it );
			}
			return tmplcc;
		}



		template < typename ContainerT >
		const TaxonNode* getICLCA( ContainerT nodescontainer ) const { //TODO: handle taxid not found
			if( nodescontainer.empty() ) {
				return NULL;
			}

			//const TaxonNode* root = getRoot();
			for( typename ContainerT::iterator node_it = nodescontainer.begin(); node_it != nodescontainer.end(); ++node_it ) {
				if( (*node_it)->data->is_unclassified ) {
					*node_it = mapUnclassified( *node_it );
				}
			}

			return getLCC( nodescontainer );
		}

		const TaxonNode* mapUnclassified( const TaxonNode* node ) const;
		const TaxonNode* mapUnclassified( TaxonID taxid ) const;

		std::pair< small_unsigned_int, small_unsigned_int > getPathLength( const TaxonNode* A, const TaxonNode* B ) const;
		std::pair< small_unsigned_int, small_unsigned_int > getPathLength( const TaxonID A_taxid, const TaxonID B_taxid ) const;

		boost::tuple< small_unsigned_int, small_unsigned_int, small_unsigned_int > getInterDistances( const TaxonNode* A, const TaxonNode* B ) const;
		boost::tuple< small_unsigned_int, small_unsigned_int, small_unsigned_int > getInterDistances( const TaxonID A_taxid, const TaxonID B_taxid ) const;

		small_unsigned_int getPathLengthToParent( const TaxonNode* A, const TaxonNode* B ) const;
		small_unsigned_int getPathLengthToParent( const TaxonID A_taxid, const TaxonID B_taxid ) const;

		const std::string& getNameAtRank( const TaxonNode* node, const std::string& rank ) const;
		const std::string& getNameAtRank( const TaxonNode* node, const std::string* internal_rank ) const;
		
		Taxonomy::PathUpIterator traverseUp( const TaxonNode* node ) const;
		
		template< typename T >
		T traverseDown( const TaxonNode* node ) const {
			return T( tax->begin().node, node );
		}
		
		bool isLeaf( const TaxonNode* node ) const;

	private:
		const Taxonomy* const tax;
};

// void printPath( TaxonomyInterface& taxinter, TaxonNode* node, TaxonNode* ancestor );
// void printPath( Taxonomy* tax, TaxonNode* node, TaxonNode* ancestor );

#endif // taxonomyinterface_hh_
