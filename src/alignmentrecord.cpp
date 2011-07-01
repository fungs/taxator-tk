#include <boost/lexical_cast.hpp>
// #include <boost/algorithm/string/case_conv.hpp>
#include <vector>
#include "alignmentrecord.hh"
#include "utils.hh"

inline bool lineIsMasked( const std::string& line ) {
	return ! line.empty() && line[0] == '*';
}

AlignmentRecord* AlignmentFileParser::next() { //TODO: increase performance by avoiding unneeded tests (line deletion)

	std::vector< std::string > fields;
	std::string* line = new std::string;
	AlignmentRecord* record;

	while( std::getline( handle, *line ) ) {
		if( ignoreLine( *line ) ) { continue; }
		bool mask_record;
		if( lineIsMasked( *line ) ) {
			mask_record = true;
			*line->erase( line->begin() );
		} else {
			mask_record = false;
		}

		// create new record on heap
// 		if( keep_lines ) {
// 			record = new AlignmentRecord( line );
// 		} else {
		record = new AlignmentRecord( mask_record );
// 		}

		tokenizeSingleCharDelim( *line, fields, "\t", 9, false );
		record->query_identifier = fields[0];

		try {
			record->reference_identifier = fields[3];
		} catch( boost::bad_lexical_cast e ) {
			std::cerr << "Could not parse sequence id field in line " << " from input";
			std::cerr << ", skipping record..." << std::endl;
			fields.clear();
			delete record;
			continue;
		}
		try {
			if( accessconv ) {
				record->reference_taxid = (*accessconv)[ record->reference_identifier ];
			}
		} catch( std::out_of_range e ) {
			std::cerr << "Alignment File Parser: Could not map reference sequence id to taxonomy " << record->reference_identifier << " to TaxID";
			std::cerr << ", skipping record..." << std::endl;
			fields.clear();
			delete record;
			continue;
		}
		try {
			record->pid = boost::lexical_cast< float >( fields[6] );
		} catch( boost::bad_lexical_cast e ) {
			std::cerr << "Could not parse PID field in line " << " from input";
			std::cerr << ", skipping record..." << std::endl;
			fields.clear();
			delete record;
			continue;
		}
		try {
			record->bitscore = boost::lexical_cast< float >( fields[7] );
		} catch( boost::bad_lexical_cast e ) {
			 std::cerr << "Could not parse BITSCORE field in line " << " from input";
			 std::cerr << ", skipping record..." << std::endl;
			fields.clear();
			delete record;
			continue;
		}
		try {
			record->evalue = boost::lexical_cast< double >( fields[8] );
		} catch( boost::bad_lexical_cast e ) {
			std::cerr << "Could not parse EVALUE field in line " << " from input";
			std::cerr << ", skipping record..." << std::endl;
			fields.clear();
			delete record;
			continue;
		}

		// clean up memory
		if( !keep_lines ) {
			delete line;
		} else {
			record->raw_line = line;
		}

		return record;
	}

	delete line;
	return NULL;
}



void AlignmentFileParser::destroyRecord(const AlignmentRecord* rec) const {
	delete rec;
}

