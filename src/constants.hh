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

#ifndef constants_hh_
#define constants_hh_

#include <string>
#include <vector>

const char endline = '\n';
const char tab = '\t';
const std::string tab_as_str = "\t";
const std::string default_field_separator = "\t";
const std::string default_ranks_carray[7] =  { "superkingdom", "phylum", "class", "order", "family", "genus", "species" };
const std::vector< std::string > default_ranks( default_ranks_carray, default_ranks_carray+7 );
const int default_rank_number = 7;
const int standard_max_support_stop = -1; // negative values have no effect
const char default_comment_symbol = '#';
const char default_mask_symbol = '*';
const std::string empty_string;

const std::string ENVVAR_TAXONOMY_ROOT = "TAXATORTK_NCBI_ROOT";

#endif //constants_hh_
