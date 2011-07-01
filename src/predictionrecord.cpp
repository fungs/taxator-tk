#include <boost/lexical_cast.hpp>
#include "predictionrecord.hh"
#include "types.hh"
#include "utils.hh"

PredictionRecord* PredictionFileParser::next() {

	std::vector< std::string > fields;
	std::string line;
	PredictionRecord* record = new PredictionRecord();

	while( std::getline( handle, line ) ) {
		if( ignoreLine( line ) || maskedLine( line ) ) { continue; }

		tokenizeSingleCharDelim( line, fields, default_field_separator, 3, false );
		record->query_identifier = fields[0];
		
		try {
			TaxonID taxid = boost::lexical_cast< TaxonID >( fields[1] );
			const TaxonNode* node = taxinter.getNode( taxid );
			if( node ) {
				record->prediction_node = node;
				return record;
			} else {
				std::cerr << "Could not find node with taxonomic id " << taxid << " in taxonomy";
				std::cerr << ", skipping record..." << std::endl;
				std::cerr << node << " | " << node->data->taxid << std::endl;
				fields.clear();
			}
		} catch( boost::bad_lexical_cast e ) {
			std::cerr << "Could not parse sequence id field in line " << line << " from input";
			std::cerr << ", skipping record..." << std::endl;
			fields.clear();
		}
	}
	
	delete record;
	return NULL;
};



void PredictionFileParser::destroyRecord( const PredictionRecord* rec ) const {
	delete rec;
};
