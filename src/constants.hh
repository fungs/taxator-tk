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
#include <boost/iterator/iterator_concepts.hpp>

const char endline = '\n';
const char tab = '\t';
const std::string tab_as_str = {tab};
const std::string default_field_separator = tab_as_str;
const std::vector< std::string > default_ranks = { "superkingdom", "phylum", "class", "order", "family", "genus", "species" };
const int standard_max_support_stop = -1; // negative values have no effect
const char default_comment_symbol = '#';
const char default_mask_symbol = '*';
const std::string empty_string;
const std::string ENVVAR_TAXONOMY_NCBI = "TAXATORTK_TAXONOMY_NCBI";

namespace newick {
	const std::string nodestart = "(";
	const std::string nodestop = ")";
	const std::string nodesep = ",";
	const std::string treestop = ";\n";
}

const std::string program_version = "1.3.2";
const std::string citation_note = u8R"(
J. Dröge, I. Gregor, and A. C. McHardy
Taxator-tk: precise taxonomic assignment of metagenomes by fast approximation of evolutionary neighborhoods
Bioinformatics 2015 31: 817-824.
doi: 10.1093/bioinformatics/btu745
)";

#endif //constants_hh_

