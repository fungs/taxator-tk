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
#include "types.hh"
#include "taxontree.hh"

class TaxonomyInterface {
	public:
		TaxonomyInterface( Taxonomy* taxtree ) : tax( taxtree ) {};

		TaxonNode* const getNode( const unsigned int taxid );
		TaxonNode* const getRoot();

		const TTPString getRank( const TaxonNode* node );
		const TTPString getRank( const unsigned int taxid );
		const TTPString getName( const TaxonNode* node );
		const TTPString getName( const unsigned int taxid );

		bool isParentOf( const TaxonNode* A, const TaxonNode* B );
		bool isParentOf( const unsigned int A_taxid, const unsigned int B_taxid );

		TaxonNode* const getLCA( const TaxonNode* A, const TaxonNode* B );
		TaxonNode* const getLCA( const unsigned int A_taxid, const unsigned int B_taxid );

		TaxonNode* const getLCC( TaxonNode* A, TaxonNode* B );
		TaxonNode* const getLCC( const unsigned int A_taxid, const unsigned int B_taxid );



		template < typename ContainerT >
		TaxonNode* const getLCA( const ContainerT& nodescontainer ) { //TODO: handle taxid not found
			if( nodescontainer.empty() ) {
				return NULL;
			}

			typename ContainerT::const_iterator node_it = nodescontainer.begin();
			TaxonNode* tmplca = *node_it++;
			while( node_it != nodescontainer.end() ) {
				tmplca = getLCA( tmplca, *node_it );
				++node_it;
			}
			return tmplca;
		};



		template < typename ContainerT >
		TaxonNode* const getLCC( const ContainerT& nodescontainer ) { //TODO: handle taxid not found
			if( nodescontainer.empty() ) {
				return NULL;
			}
			// we need to start with lowest overall concept
			std::pair< TaxonNode*, unsigned int > lowest( NULL, 0 );
			for( typename ContainerT::const_iterator node_it = nodescontainer.begin(); node_it != nodescontainer.end(); ++node_it ) {
				if( (*node_it)->data->root_pathlength >= lowest.second ) {
					lowest.first = *node_it;
					lowest.second = (*node_it)->data->root_pathlength;
				}
			}
			TaxonNode* tmplcc = lowest.first;
			for( typename ContainerT::const_iterator node_it = nodescontainer.begin(); node_it != nodescontainer.end(); ++node_it ) {
				tmplcc = getLCC( tmplcc, *node_it );
			}
			return tmplcc;
		};



		template < typename ContainerT >
		TaxonNode* const getICLCA( ContainerT nodescontainer ) { //TODO: handle taxid not found
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
		};

		TaxonNode* mapUnclassified( TaxonNode* node );
		TaxonNode* mapUnclassified( unsigned int taxid );

		std::pair< int, int > getPathLength( const TaxonNode* A, const TaxonNode* B );
		std::pair< int, int > getPathLength( const unsigned int A_taxid, const unsigned int B_taxid );

		boost::tuple< int, int, int > getInterDistances( const TaxonNode* A, const TaxonNode* B );
		boost::tuple< int, int, int > getInterDistances( const unsigned int A_taxid, const unsigned int B_taxid );

		int getPathLengthToParent( const TaxonNode* A, const TaxonNode* B );
		int getPathLengthToParent( const unsigned int A_taxid, const unsigned int B_taxid );

		const std::string& getNameAtRank( TaxonNode* node, const std::string& rank );
		const std::string& getNameAtRank( TaxonNode* node, const std::string* internal_rank );

	private:
		Taxonomy* tax;
};

void printPath( TaxonomyInterface& taxinter, TaxonNode* node, TaxonNode* ancestor );
void printPath( Taxonomy* tax, TaxonNode* node, TaxonNode* ancestor );

#endif // taxonomyinterface_hh_
