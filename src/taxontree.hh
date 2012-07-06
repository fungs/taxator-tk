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

#include "types.hh"
#include <tree.hh>
#include <boost/tuple/tuple.hpp>
#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <stack>



class TaxonAnnotation {
		// contains all information like name, synonyms, rank etc.
	public:
		TaxonAnnotation( const std::string& rankname ) : rank( rankname ) {};
		TaxonAnnotation( const std::string& rankname, const std::string& taxonname ) : rank( rankname ), name( taxonname ) {};
		const std::string& rank;
		std::string name;
	};



class Taxon {
	public:
		Taxon() : annotation( NULL ), mark_special( false ), is_unclassified( false ) {};
		Taxon( TaxonAnnotation* taxanno ) : annotation( taxanno ), mark_special( false ), is_unclassified( false ) {};
		~Taxon() {
				if( annotation ) { //delete object on heap
					delete annotation;
				}
		};
		
		//default order is pre-order
		bool operator<( const Taxon& t ) {
			return this->leftvalue < t.leftvalue;
		}
		
// 		Taxon( const Taxon& taxon);
		TaxonID taxid;
		small_unsigned_int root_pathlength;
		large_unsigned_int leftvalue; //nested set value
		large_unsigned_int rightvalue; //nested set value
		TaxonAnnotation* annotation;
		bool mark_special;
		bool is_unclassified;
};



typedef tree_node_<Taxon*> TaxonNode;



class TaxonomyInterface;



class TaxonTree : public tree< Taxon* > {
	friend class TaxonomyInterface;
	public:
		TaxonTree() : rank_not_found_( *ranks_.insert( "" ).first ) {};
 		~TaxonTree();
		typedef tree_node Node;
		int indexSize() const;
		const std::string& insertRankInternal( const std::string& rankname );
		const std::string& getRankInternal( const std::string& rankname ) const;
		void deleteUnmarkedNodes();
// 		void addDummyRankNodes( const std::vector< std::string >& ranks );
		void setRankDistances( const std::vector< std::string >& ranks );
		void setMaxDepth( small_unsigned_int depth ) { max_depth_ = depth; };
		void setMaxDepth();
		void recalcNestedSetInfo();
		void addToIndex( TaxonID taxid, Node* node );
		void recreateNodeIndex();
		
		// base class for path iterators (only forward)
		class PathIteratorBase {
		public:
			typedef std::forward_iterator_tag iterator_category;
			typedef PathIteratorBase self_type;
			typedef Node value_type;
			typedef size_t difference_type;
			typedef Node* pointer;
			typedef const Node* const_pointer;
			typedef Node& reference;
			typedef const Node& const_reference;
			
			const_reference operator*() const {
				return *current;
			};
			
			const_pointer operator->() const {
				return current;
			};
			
			bool operator==( const PathIteratorBase &it ) const {
				return current == it.current;
			};
			
			bool operator==( const_pointer node ) const {
				return current == node;
			};
			
			bool operator!=( const PathIteratorBase &it ) const {
				return current != it.current;
			};
			
			bool operator!=( const_pointer node ) const {
				return current != node;
			};
			
			virtual self_type& operator++() = 0;
			
		protected:
			const Node* current;
		};
		
		// attention: PathDownIterator's increment operator doesn't run in constant time but is linear in number of children
		class PathDownIterator : public PathIteratorBase {
		public:
			PathDownIterator( const_pointer startnode, const_pointer stopnode ) : stop_v( stopnode->data->leftvalue ) {
				current = startnode;
			};

			PathDownIterator& operator++() {
				for( sibling_iterator it = current->first_child; it != it.end(); ++it ) {
					Taxon* t = *it;
					if( t->leftvalue <= stop_v && t->rightvalue > stop_v ) {
						current = it.node;
						break;
					}
				}
				return *this;
			};
			
			PathDownIterator operator++( int ) {
				PathDownIterator tmp( *this );
				operator++();
				return tmp;
			}
			
		private:
			const large_unsigned_int stop_v; //leading value
		};
		
		class PathUpIterator : public PathIteratorBase { //iterator will go forever, if not checked for target
		public:
			PathUpIterator( const_pointer startnode ) {
				current = startnode;
			};

			PathUpIterator& operator++() {
				current = current->parent;
				return *this;
			};
			
			PathUpIterator operator++( int ) {
				PathUpIterator tmp( *this );
				operator++();
				return tmp;
			}
		};
			
		//cached version of PathDownIterator using a stack as cache to do everything in constant time 
		class CPathDownIterator : public PathIteratorBase {
		public:
			CPathDownIterator( const_pointer startnode, const_pointer stopnode ) {
				current = startnode;
				for( PathUpIterator it( stopnode ); it != startnode; ++it ) {
					path.push( &*it );
				}
			};

			CPathDownIterator& operator++() {
				current = path.top();
				path.pop();
				return *this;
			};
			
			CPathDownIterator operator++( int ) {
				CPathDownIterator tmp( *this );
				operator++();
				return tmp;
			}
			
		private:
			std::stack< const_pointer > path;
		};

	private:
		std::set< std::string > ranks_;
		const std::string& rank_not_found_;
		std::map< unsigned int, Node* > taxid2node_; //use boost::ptr_map<> -> no destructor needed, hash map is faster
		small_unsigned_int max_depth_;
};



typedef TaxonTree Taxonomy;

#endif // taxontree_hh_
