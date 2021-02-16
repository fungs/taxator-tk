#include <boost/lexical_cast.hpp>
#include "predictionrecord.hh"

//overload ostream operator for class PredictionRecordBase ->print()
std::ostream& operator<<( std::ostream& strm, const PredictionRecordBase& prec ) {
	prec.print( strm );
	return strm;
}



//overload istream operator for class PredictionRecordBase ->parse()
std::istream& operator>>( std::istream& strm, PredictionRecordBase& prec ) { //TODO: what to do if return of parse() is false
	std::string line;
	std::getline( strm, line );
	prec.parse( line );
	return strm;
}



std::ostream& operator<<( std::ostream& strm, const GFF3Header& ) {
	strm << "##gff-version 3" << std::endl;
	return strm;
}
