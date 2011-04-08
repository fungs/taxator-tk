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

#ifndef trainingrecord_hh_
#define trainingrecord_hh_

#include <fstream>
#include <string>
#include "types.hh"
#include "utils.hh"



class TrainingRecord {
	public:
		TrainingRecord() : mask( false ), raw_line( NULL ) {};
		TrainingRecord( std::string* line ) : mask( false ), raw_line( line ) {};
		~TrainingRecord() { if( raw_line ) delete raw_line; };
		unsigned int reference_taxid;
		unsigned int query_taxid;
		std::string query_identifier;
		float bitscore;
		double evalue;
		std::string* raw_line;
		bool mask;
};



class TrainingFileParser {
	public:
		TrainingFileParser( const std::string& filename, const bool keep_raw_lines = false ) : filehandle( filename.c_str() ), handle( filehandle ), keep_lines( keep_raw_lines ) {};
		TrainingFileParser( std::istream& strm, const bool keep_raw_lines = false ) : handle( strm ), keep_lines( keep_raw_lines ) {};
		TrainingRecord* next();
		bool eof() { return handle.eof(); };
	private:
		std::istream& handle;
		std::ifstream filehandle;
		const bool keep_lines;
};



class TrainingRecordSetGenerator {
	public:
		TrainingRecordSetGenerator( TrainingFileParser& p ) : parser( p ), last_record( p.next() ) {
			last_query_id = last_record ? &(last_record->query_identifier) : NULL;
		}

		template< typename ContainerT >
		void getNext( ContainerT& recordset ) {
			if( last_record ) {
				TrainingRecord* record = last_record;
				std::string& query_id = *last_query_id;

				do {
					if( query_id == record->query_identifier ) { //still the same query
						recordset.push_back( record );
					} else {
						last_query_id = &(record->query_identifier);
						last_record = record;
						break;
					}
					record = parser.next();
				} while( record );
				last_record = record;
			}
		};

		bool notEmpty() { return last_record; };

	private:
		TrainingFileParser& parser;
		std::string* last_query_id;
		TrainingRecord* last_record;
};



#endif // trainingrecord_hh_
