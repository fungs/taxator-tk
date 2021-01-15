#include "accessconv.hh"

StrIDConverter* loadStrIDConverterFromFile( const std::string& filename, unsigned int cachesize ) { //TODO: remove depricated
  return loadAccessIDConverterFromFile< std::string >( filename, cachesize );
}
