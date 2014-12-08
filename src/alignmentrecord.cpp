#include "alignmentrecordRE.hh"

//overload ostream operator for class AlignmentRecord ->print()
std::ostream& operator<<( std::ostream& strm, const AlignmentRecord& rec ) {
	rec.print( strm );
	return strm;
}



//overload istream operator for class AlignmentRecord ->parse()
std::istream& operator>>( std::istream& strm, AlignmentRecord& rec ) {
	std::string line;
	std::getline( strm, line );
	rec.parse( line );
	return strm;
}
