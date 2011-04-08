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

#ifndef ncbidata_hh_
#define ncbidata_hh_

#include "taxontree.hh"
#include "sqlite3pp.hh"



Taxonomy* parseNCBIFlatFiles( const std::string& names, const std::string& nodes, const std::vector< std::string >* ranks_to_mark = NULL );



Taxonomy* loadTaxonomyFromEnvironment( const std::vector< std::string >* ranks_to_mark = NULL );



const std::string extractFastaCommentField( const std::string& comment, const std::string& key );


#endif // ncbidata_hh_
