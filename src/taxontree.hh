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

class TaxonomyInterface;

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

#endif // taxontree_hh_
