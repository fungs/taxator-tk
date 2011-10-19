#include <boost/lexical_cast.hpp>
#include "predictionrecord.hh"

PredictionRecord* PredictionFileParser::next() {

		std::string line;
	PredictionRecord* record = new PredictionRecord( tax_ );

	while( std::getline( handle, line ) ) {
		if( record->parse( line ) ) return record;
	}
	delete record;
	return NULL;
}



void PredictionFileParser::destroyRecord( const PredictionRecord* rec ) const {
	delete rec;
}



//overload ostream operator for class AlignmentRecord ->print()
std::ostream& operator<<( std::ostream& strm, const PredictionRecord& prec ) {
	prec.print( strm );
	return strm;
}



//overload istream operator for class AlignmentRecord ->parse()
std::istream& operator>>( std::istream& strm, PredictionRecord& prec ) { //TODO: what to do if return of parse() is false
	std::string line;
	std::getline( strm, line );
	prec.parse( line );
	return strm;
}