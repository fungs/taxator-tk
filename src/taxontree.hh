/*
	<one line to give the program's name and a brief idea of what it does.>
	Copyright (C) <year>  <name of author>

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

#ifndef taxontree_hh_
#define taxontree_hh_

#include <map>
#include <vector>
#include <iostream>
#include <boost/tuple/tuple.hpp>
#include "types.hh"
#include "tree.hh"



class TaxonAnnotation {
		// contains all information like name, synonyms, rank etc.
	public:
		TaxonAnnotation( const TTPString& rankname ) : rank( rankname ) {};
		TaxonAnnotation( const TTPString& rankname, const TTPString& taxonname ) : rank( rankname ), name( taxonname ) {};
		const TTPString& rank;
		TTPString name;
	};



class Taxon {
	public:
		Taxon() : annotation( NULL ), mark_special( false ), is_unclassified( false ) {};
		Taxon( TaxonAnnotation* taxanno ) : annotation( taxanno ), mark_special( false ), is_unclassified( false ) {};
		~Taxon() {
			TaxonAnnotation* delanno = annotation;
				if( delanno ) { //delete object on heap
					delete delanno;
				}
		};
// 		Taxon( const Taxon& taxon);
		unsigned int taxid;
		unsigned int root_pathlength;
		unsigned int leftvalue; //nested set value
		unsigned int rightvalue; //nested set value
		TaxonAnnotation* annotation;
		bool mark_special;
		bool is_unclassified;
	};



typedef tree_node_<Taxon*> TaxonNode;



class TaxonTree : public tree< Taxon* > {
	friend class TaxonomyInterface;
	public:
		TaxonTree() {};
 		~TaxonTree();
		typedef tree_node Node;
		int indexSize();
		const std::string& getRankInternal( const std::string& rankname, const bool insert = false );
		void deleteUnmarkedNodes();
// 		void addDummyRankNodes( const std::vector< std::string >& ranks );
		void setRankDistances( const std::vector< std::string >& ranks );
		void recalcNestedSetInfo();
		void addToIndex( unsigned int taxid, Node* node );
		void recreateNodeIndex();

	private:
		std::set< std::string > ranks;
		std::map< unsigned int, Node* > taxid2node; //use boost::ptr_map<> -> no destructor needed
};



typedef TaxonTree Taxonomy;



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

// 		template < typename ContainerT >
// 		TaxonNode* const getICLCA2( const ContainerT& nodescontainer ) { //TODO: handle taxid not found
// 			if( nodescontainer.empty() ) {
// 				return NULL;
// 			}
// 			// 			std::cerr << "running iclca" << std::endl;
// 			std::list< TaxonNode* > incomplete_nodes, normal_nodes;
// 			const TaxonNode* root = getRoot();
// 			for( typename ContainerT::const_iterator node_it = nodescontainer.begin(); node_it != nodescontainer.end(); ++node_it ) {
// 				TaxonNode* node = *node_it;
// 				if( node->data->is_unclassified ) {
// 					incomplete_nodes.push_back( mapUnclassified( node ) );
// 				} else {
// 					normal_nodes.push_back( node );
// 				}
// 			}
// 			std::cerr << "nodes total: " << nodescontainer.size() << std::endl;
// 			std::cerr << "separated nodes, incomplete: " << incomplete_nodes.size() << ", normal: " << normal_nodes.size() << std::endl;
//
// 			TaxonNode* lca;
// 			if( normal_nodes.empty() ) {
// 				if( incomplete_nodes.empty() ) {
// 					return NULL;
// 				}
// 				lca = incomplete_nodes.back();
// 				incomplete_nodes.pop_back();
// 			} else {
// 				lca = getLCA( normal_nodes );
// 			}
//
// 			for( std::list< TaxonNode* >::iterator it = incomplete_nodes.begin(); it != incomplete_nodes.end(); ++it ) {
// 				if( ! isParentOf( (*it), lca ) ) {
// 					std::cerr << "incomplete knowledge incorporation: " << (*it)->data->annotation->name << " + " << lca->data->annotation->name;
// 					lca = getLCA( (*it), lca );
// 					std::cerr << " = " << lca->data->annotation->name << std::endl;
// 				}
// 			}
//
// 			return lca;
// 		};

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



void printPath( Taxonomy* tax, TaxonNode* node, TaxonNode* ancestor );
void printPath( TaxonomyInterface& taxinter, TaxonNode* node, TaxonNode* ancestor );


// template<class T>
// void print_subtree_pointer_bracketed(const tree<T*>& t, typename tree<T>::iterator iRoot, std::ostream& str);
//
// // Iterate over all roots (the head) and print each one on a new line
// // by calling printSingleRoot.
// template<class T>
// void print_tree_pointer_bracketed(const tree<T*>& t, std::ostream& str = std::cout )
// {
// 	int headCount = t.number_of_siblings(t.begin());
// 	int headNum = 0;
// 	for(typename tree<T*>::sibling_iterator iRoots = t.begin(); iRoots != t.end(); ++iRoots, ++headNum) {
// 		print_subtree_pointer_bracketed(t,iRoots,str);
// 		if (headNum != headCount) {
// 			str << std::endl;
// 		}
// 	}
// };
//
//
// // Print everything under this root in a flat, bracketed structure.
// template<class T>
// void print_subtree_pointer_bracketed(const tree<T*>& t, typename tree<T*>::iterator iRoot, std::ostream& str)
// {
// 	if(t.empty()) return;
// 	if (t.number_of_children(iRoot) == 0) {
// 		str << (*iRoot)->annotation->name;
// 	}
// 	else {
// 		// parent
// 		str << (*iRoot)->annotation->name;
// 		str << "(";
// 		// child1, ..., childn
// 		int siblingCount = t.number_of_siblings(t.begin(iRoot));
// 		int siblingNum;
// 		typename tree<T*>::sibling_iterator iChildren;
// 		for (iChildren = t.begin(iRoot), siblingNum = 0; iChildren != t.end(iRoot); ++iChildren, ++siblingNum) {
// 			// recursively print child
// 			print_subtree_pointer_bracketed(t,iChildren,str);
// 			// comma after every child except the last one
// 			if (siblingNum != siblingCount ) {
// 				str << ", ";
// 			}
// 		}
// 		str << ")";
// 	}
// };

#endif // taxontree_hh_
