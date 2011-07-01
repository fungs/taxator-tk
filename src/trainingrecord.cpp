#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>
#include "trainingrecord.hh"
#include "utils.hh"
#include "constants.hh"

TrainingRecord* TrainingFileParser::next() { //TODO: increase performance by avoiding unneeded tests (line deletion)

	std::list< std::string > fields;
	std::list< std::string >::iterator field_it;
	std::string* line = new std::string;
	TrainingRecord* record;

	while( std::getline( handle, *line ) ) {
		if( ignoreLine( *line ) ) { continue; }

		// create new record on heap
		if( keep_lines ) {
			record = new TrainingRecord( line );
		} else {
			record = new TrainingRecord;
		}

		tokenizeSingleCharDelim( *line, fields, default_field_separator, 4 );
		field_it = fields.begin();
		
		try {
			record->reference_taxid = boost::lexical_cast< unsigned int >( *field_it++ );
		} catch( boost::bad_lexical_cast e ) {
			std::cerr << "Could not parse reference taxid field from input in line \"" << line << "\", skipping record..." << std::endl;
			fields.clear();
			delete record;
			continue;
		}
		
		try {
			record->query_taxid = boost::lexical_cast< unsigned int >( *field_it++ );
		} catch( boost::bad_lexical_cast e ) {
			std::cerr << "Could not parse query taxid field from input in line \"" << line << "\", skipping record..." << std::endl;
			fields.clear();
			delete record;
			continue;
		}
		
		record->query_identifier = *field_it++;
		
		try {
			record->bitscore = boost::lexical_cast< float >( *field_it++ );
		} catch( boost::bad_lexical_cast e ) {
			std::cerr << "Could not parse bitscore field from input in line \"" << line << "\", skipping record..." << std::endl;
			fields.clear();
			delete record;
			continue;
		}
		
		try {
			record->evalue = boost::lexical_cast< double >( *field_it++ );
		} catch( boost::bad_lexical_cast e ) {
			std::cerr << "Could not parse evalue field from input in line \"" << line << "\", skipping record..." << std::endl;
			fields.clear();
			delete record;
			continue;
		}

		// clean up memory
		if( !keep_lines ) {
			delete line;
		}

		return record;
	}

	delete line;
	return NULL;
}
