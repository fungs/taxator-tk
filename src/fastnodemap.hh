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

#ifndef fastnodemap_hh_
#define fastnodemap_hh_

#include "taxontree.hh"
#include "types.hh"
#include <map>
#include <vector>

// meta map without iterators

template< typename ValueType >
class FastNodeMap {
	public:
		FastNodeMap( small_unsigned_int max_depth ) : map_at_level_( max_depth ) {};
		typedef std::map< const TaxonNode*, ValueType > BasicMapType;
		
		typename std::vector< BasicMapType >::size_type size() const {
			typename std::vector< BasicMapType >::size_type sum = 0;
			for ( typename std::vector< BasicMapType >::const_iterator it = map_at_level_.begin(); it != map_at_level_.end(); ++it ) sum += it->size();
			return sum;
		}
		
		ValueType& operator[]( const TaxonNode* node ) {
			return map_at_level_[ node->data->root_pathlength ][ node ];
		};
		
		ValueType* find( const TaxonNode* node ) {
			BasicMapType& directmap = map_at_level_[ node->data->root_pathlength ];
			typename BasicMapType::iterator it = directmap.find( node );
			if ( it != directmap.end() ) return &it->second;
			return NULL;
		};
		
	private:
		std::vector< BasicMapType > map_at_level_;
};

#endif //fastnodemap_hh_
