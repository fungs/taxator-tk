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

#ifndef types_hh_
#define types_hh_

#include <boost/cstdint.hpp>
#include <vector>
#include <list>
#include <map>
#include <string>

// some integer type definitions using boost
typedef int_least8_t small_int; // –128 to 127
typedef int_least16_t medium_int; // –32,768 to 32,767
typedef int_least32_t large_int; // –2,147,483,648 to 2,147,483,647
typedef uint_least8_t small_unsigned_int; // 0 to 255
typedef uint_least16_t medium_unsigned_int; // 0 to 65,535
typedef uint_least32_t large_unsigned_int; // 0 to 4,294,967,295
typedef uint_least64_t very_large_unsigned_int; // 0 to 18,446,744,073,709,551,615

// basic mappings
typedef large_unsigned_int TaxonID; //maximum number at time of writing was 1,050,856

#endif //types_hh_
