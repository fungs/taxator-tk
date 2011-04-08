#include "accessconv.hh"

StrIDConverter* loadStrIDConverterFromFile( const std::string& filename, unsigned int cachesize ) {
  return loadAccessIDConverterFromFile< std::string >( filename, cachesize );
};
