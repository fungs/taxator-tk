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

#ifndef predictionrecord_hh_
#define predictionrecord_hh_

#include <string>
#include <fstream>
#include <iostream>
#include "types.hh"
#include "constants.hh"
#include "taxontree.hh"
#include "taxonomyinterface.hh"



class PredictionRecord {
	public:
		const TaxonNode* prediction_node;
		std::string query_identifier;
		void print( std::ostream& strm = std::cout ) const {
			strm << query_identifier << default_field_separator << prediction_node->data->taxid << std::endl;
		};
};



class PredictionFileParser {
	public:
		PredictionFileParser( const std::string& filename, const Taxonomy* tax ) : filehandle( filename.c_str() ), handle( filehandle ), taxinter( tax ) {};
		PredictionFileParser( std::istream& strm, const Taxonomy* tax ) : handle( strm ), taxinter( tax ) {};
		PredictionRecord* next();
		void destroyRecord( const PredictionRecord* rec ) const;
		bool eof() const { return handle.eof(); };
	protected:
		std::ifstream filehandle;
		std::istream& handle;
		TaxonomyInterface taxinter;
};



#endif // predictionrecord_hh_
