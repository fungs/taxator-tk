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

#ifndef alignmentrecord_hh_
#define alignmentrecord_hh_

#include <fstream>
#include <string>
#include "types.hh"
#include "utils.hh"
#include "ncbidata.hh"
#include "accessconv.hh"
#include "taxonomyinterface.hh"



class AlignmentRecord {
	public:
		AlignmentRecord( bool m = false ) : mask( m ), raw_line( NULL ) {};
		AlignmentRecord( std::string* line , bool m = false ) : mask( m ), raw_line( line ) {};
		~AlignmentRecord() { if( raw_line ) { delete raw_line; } };
		void print( std::ostream& strm = std::cout ) const {
			//TODO: strm << ... << default_field_separator << ..<< std::endl;
		};

		unsigned int reference_taxid;
		std::string reference_identifier;
		std::string query_identifier;
		float pid;
		float bitscore;
		double evalue;
		bool mask;
		std::string* raw_line;
};



class AlignmentFileParser {
	public:
		AlignmentFileParser( const std::string& filename, StrIDConverter* seqid2taxid, const bool keep_raw_lines = false ) : accessconv( seqid2taxid ), filehandle( filename.c_str() ), handle( filehandle ), keep_lines( keep_raw_lines ) {};
		AlignmentFileParser( std::istream& strm, StrIDConverter* seqid2taxid, const bool keep_raw_lines = false ) : accessconv( seqid2taxid ), handle( strm ), keep_lines( keep_raw_lines ) {};
		AlignmentRecord* next();
		void destroyRecord( const AlignmentRecord* rec ) const;
		bool eof() { return handle.eof(); };
	protected:
		StrIDConverter* accessconv;
		std::ifstream filehandle;
		std::istream& handle;
		const bool keep_lines;
};



class RecordSetGenerator {
	public:
		RecordSetGenerator( AlignmentFileParser& p ) : parser( p ), last_record( p.next() ) {
			last_query_id = last_record ? &(last_record->query_identifier) : NULL;
		}

		template< typename ContainerT >
		void getNext( ContainerT& recordset ) {
			if( last_record ) {
				AlignmentRecord* record = last_record;
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
		AlignmentFileParser& parser;
		AlignmentRecord* last_record;
		std::string* last_query_id;
};



template< typename ContainerT1, typename ContainerT2 >
void records2Nodes( const ContainerT1& recordset, TaxonomyInterface* inter, ContainerT2& refnodes ) {
	typename ContainerT1::const_iterator record_it = recordset.begin();
	const TaxonNode* node;
	while( record_it != recordset.end() ) {
		if( !(*record_it)->mask ) {
			node = inter->getNode( (*record_it)->reference_taxid );
			if( node ) {
				refnodes.push_back( node );
			}
		}
		++record_it;
	}
};



template< typename ContainerT >
void deleteRecords( ContainerT& recordset ) {
	while( recordset.size() ) {
		delete recordset.back();
		recordset.pop_back();
	}
};



// template< typename ContainerT1 >
// AlignmentRecord* getBestBitscoreRecord( ContainerT1& recordset ) {
// 	if( recordset.empty() ) {
// 		return NULL;
// 	}
//
// 	typename ContainerT1::iterator record_it = recordset.begin();
// 	typename ContainerT1::value_type record = *record_it;
// 	while( record_it != recordset.end() ) {
// 		if( ! (*record_it)->mask && (*record_it)->bitscore > record->bitscore ) {
// 			record = *record_it;
// 		}
// 		++record_it;
// 	}
// 	return record;
// };



#endif // alignmentrecord_hh_
