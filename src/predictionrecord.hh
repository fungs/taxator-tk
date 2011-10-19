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
#include "types.hh"
#include "utils.hh"



class PredictionRecord {
	public:
		PredictionRecord( const Taxonomy* tax ) : query_coverage( 1. ), taxinter_( tax ) {};
		std::string query_identifier;
		const TaxonNode* lower_node;
		const TaxonNode* upper_node;
		float alpha;
		float query_coverage;
		
		bool parse( const std::string& line ) {
			if( ignoreLine( line ) || maskedLine( line ) ) return false;
			
			std::vector< std::string > fields;
			tokenizeSingleCharDelim( line, fields, default_field_separator, 5, false );
			query_identifier = fields[0];
			
			try {
				TaxonID taxid = boost::lexical_cast< TaxonID >( fields[1] );
				const TaxonNode* node = taxinter_.getNode( taxid );
				if( node ) {
					lower_node = node;
				} else {
					std::cerr << "Could not find node with taxonomic id " << taxid << " in taxonomy";
					std::cerr << ", skipping record..." << std::endl;
					return false;
				}
				
				taxid = boost::lexical_cast< TaxonID >( fields[2] );
				node = taxinter_.getNode( taxid );
				if( node ) {
					upper_node = node;
				} else {
					std::cerr << "Could not find node with taxonomic id " << taxid << " in taxonomy";
					std::cerr << ", skipping record..." << std::endl;
					return false;
				}
				
			} catch( boost::bad_lexical_cast e ) {
				std::cerr << "Could not parse taxonomic ID in line " << line << " from input";
				std::cerr << ", skipping record..." << std::endl;
				return false;
			}
			
			try {
				alpha = boost::lexical_cast< float >( fields[3] );
			} catch ( boost::bad_lexical_cast e ) {
				std::cerr << "Could not parse alpha field in line " << line << " from input";
				std::cerr << ", skipping record..." << std::endl;
				return false;
			}
			
			try {
				query_coverage = boost::lexical_cast< float >( fields[4] );
			} catch ( boost::bad_lexical_cast e ) {
				std::cerr << "Could not parse query coverage field in line " << line << " from input";
				std::cerr << ", skipping record..." << std::endl;
				return false;
			}
			
			return true;
		}
		
		void print( std::ostream& strm = std::cout ) const {
			strm << query_identifier << default_field_separator << lower_node->data->taxid << default_field_separator << upper_node->data->taxid << default_field_separator << alpha << default_field_separator << query_coverage << std::endl;
		};
	private:
		TaxonomyInterface taxinter_;
};



class PredictionFileParser {
	public:
		PredictionFileParser( const std::string& filename, const Taxonomy* tax ) : filehandle( filename.c_str() ), handle( filehandle ), tax_( tax ) {};
		PredictionFileParser( std::istream& strm, const Taxonomy* tax ) : handle( strm ), tax_( tax ) {};
		PredictionRecord* next();
		void destroyRecord( const PredictionRecord* rec ) const;
		bool eof() const { return handle.eof(); };
	protected:
		std::ifstream filehandle;
		std::istream& handle;
		const Taxonomy* tax_;
};



std::ostream& operator<<( std::ostream& strm, const PredictionRecord& prec );


std::istream& operator>>( std::istream& strm, PredictionRecord& prec );


#endif // predictionrecord_hh_
