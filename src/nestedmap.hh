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

#ifndef nestedmap_hh_
#define nestedmap_hh_


#include<map>

// abstract interface class
template< typename KeyT, typename ValueT, class MapT = std::map< KeyT, ValueT > >
class NestedMapBase {
	public:
		virtual int size() = 0;
		virtual typename MapT::iterator insert( KeyT key, ValueT val ) = 0;
		virtual bool find( KeyT key, typename MapT::iterator& it ) = 0;
		virtual void appendToMap( MapT& m ) = 0;
	private:
};



template< typename KeyT, typename ValueT, class MapT = std::map< KeyT, ValueT > >
class NestedMapSTLWrapper : public NestedMapBase< KeyT, ValueT, MapT > {
	public:
		NestedMapSTLWrapper() {};
		NestedMapSTLWrapper( MapT m ) : stlmap( m ) {};

		int size() {
			return stlmap.size();
		};

		typename MapT::iterator insert( KeyT key, ValueT val ) {
			std::pair< typename MapT::iterator, bool > returnpair = stlmap.insert( make_pair( key, val ) );
			if( ! returnpair.second ) {
				returnpair.first->second = val;
			}
			return returnpair.first;
		};

		bool find( KeyT key, typename MapT::iterator& it ) {
			it = stlmap.find( key );
			return it != stlmap.end();
		};

		void appendToMap( MapT& m ) {
// 			std::cerr << "Calling appendToMap() in NestedMapSTLWrapper" << std::endl;
			if( stlmap.size() ) {
// 				std::cerr << ".insert()" << std::endl;
// 				std::cerr << "inserting map of size " << stlmap.size() << " into map of size " << m.size() << std::endl;
// 				typename MapT::iterator it = stlmap.begin();
// 				while( it != stlmap.end() ) {
// 					std::cerr << "Entry___key: " << it->first << " ___value: " << it->second->bitscore << std::endl;
// 					m[ it->first ] = it->second;
// 					++it;
// 				}
// 				std::cerr << "Manual insertion finished" << std::endl;
				m.insert( stlmap.begin(), stlmap.end() );
// 				std::cerr << "Range re-insertion finished" << std::endl;
			}
		};

	private:
		MapT stlmap;
};



template< typename KeyT, typename ValueT, class MapT = std::map< KeyT, ValueT > >
class NestedMap : public NestedMapBase< KeyT, ValueT, MapT > {
	public:
// 		typedef NestedMapBase< MapT::key_type, MapT::value_type >::iterator iterator;

		NestedMap( NestedMapBase< KeyT, ValueT, MapT >& m1 ) : map1( m1 ), maplength( m1.size() ), increase_counter( false ), optinsert( false ) {};

		int size() {
			return maplength;
		};

		typename MapT::iterator insert( KeyT key, ValueT val ) { //distorts maplength counter
			return map2.insert( key, val );
		};

		typename MapT::iterator insertAfterFind( KeyT key, ValueT val ) {
			if( optinsert ) {
				optinsert_it->second = val;
				return optinsert_it;
			}

			if( increase_counter ) {
				++maplength;
				increase_counter = false;
			}

			return insert( key, val );
		};

		bool find( KeyT key, typename MapT::iterator& it ) {
			increase_counter = false;
			if ( map2.find( key, it ) ) {
				optinsert_it = it; //store location for faster insert
				optinsert = true;
				return true;
			}

			optinsert = false; //no optimized insert possible

			if ( map1.find( key, it ) ) {
				return true;
			}

			increase_counter = true; //totally new element
			return false;
		};

		void appendToMap( MapT& m ) {
// 			std::cerr << "Calling appendToMap()" << std::endl;
			map2.appendToMap( m );
// 			std::cerr << "Returned correctly from appendToMap()" << std::endl;
// 			std::cerr << "Calling parent map of size " << map1.size() << std::endl;
			map1.appendToMap( m );
		};

	private:
		NestedMapBase< KeyT, ValueT, MapT >& map1; //just readable
		NestedMapSTLWrapper< KeyT, ValueT, MapT > map2;
		int maplength;
		typename MapT::iterator optinsert_it;
		bool optinsert;
		bool increase_counter;
};



#endif // nestedmap_hh_
