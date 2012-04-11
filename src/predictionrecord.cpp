#include <boost/lexical_cast.hpp>
#include "predictionrecord.hh"

//overload ostream operator for class AlignmentRecord ->print()
std::ostream& operator<<( std::ostream& strm, const PredictionRecordBase& prec ) {
	prec.print( strm );
	return strm;
}



//overload istream operator for class AlignmentRecord ->parse()
std::istream& operator>>( std::istream& strm, PredictionRecordBase& prec ) { //TODO: what to do if return of parse() is false
	std::string line;
	std::getline( strm, line );
	prec.parse( line );
	return strm;
}