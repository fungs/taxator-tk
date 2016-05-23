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

#ifndef accessconv_hh_
#define accessconv_hh_

#include <string>
#include <map>
#include <queue>
#include <iostream>
#include <list>
#include <fstream>
#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include "constants.hh"
#include "types.hh"
#include "utils.hh"
#include "exception.hh"



// converts from access identifier to taxonomic id
template< typename TypeT > //TODO: add operator[] const (avoid caching, history etc.)
class AccessIDConverter {
public:
    virtual ~AccessIDConverter() {};
    virtual TaxonID operator[]( const TypeT& acc ) /*throw( std::out_of_range )*/ = 0;
};



template< typename TypeT >
class AccessIDConverterFlatfileMemory : public AccessIDConverter< TypeT > {
public:
    AccessIDConverterFlatfileMemory( const std::string& flatfile_filename ) : filename_(flatfile_filename) {
        parse(flatfile_filename);
    }

    TaxonID operator[](const TypeT& acc) { /*throw( std::out_of_range )*/
        typename std::map< TypeT, TaxonID >::iterator it = accessidconv.find(acc);
        if(it == accessidconv.end()) BOOST_THROW_EXCEPTION(TaxonMappingNotFound{} << seqid_info{acc} << file_info{filename_});
        return it->second;
    }

private:
    void parse( const std::string& flatfile_filename ) {
        std::list< std::string > fields;
        std::list< std::string >::iterator field_it;
        std::string line;
        std::ifstream flatfile( flatfile_filename.c_str() );
        TypeT acc;
        while( std::getline( flatfile, line ) ) {
            if( ignoreLine( line ) ) {
                continue;
            }
            fields.clear();
            tokenizeSingleCharDelim( line, fields, default_field_separator, 2 );
            field_it = fields.begin();

            try {
                acc = boost::lexical_cast< TypeT >( *field_it );
                ++field_it;
                TaxonID taxid = boost::lexical_cast< TaxonID >( *field_it );
                accessidconv[ acc ] = taxid;
            } catch( boost::bad_lexical_cast &e ) {  // TODO: pass info via exception
                std::cerr << "Could not parse line: " << line << ", skipping alignment." << std::endl;
                std::cerr << "key:" << acc << std::endl;
                std::cerr << "error parsing taxonomic ID: " << *field_it << std::endl;
                throw e;
            }
        }
        flatfile.close();
    };

    typename std::map< TypeT, TaxonID > accessidconv; //TODO: hash_map aka unordered_map would be better
    const std::string filename_;
};



template< typename TypeT >
AccessIDConverter< TypeT >* loadAccessIDConverterFromFile( const std::string& filename, unsigned int cachesize = 0 ) {
    AccessIDConverter< TypeT >* accidconv = NULL;
    if( ! boost::filesystem::exists( filename ) ) BOOST_THROW_EXCEPTION(FileNotFound {} << file_info {filename});
    accidconv = new AccessIDConverterFlatfileMemory< TypeT >( filename );
    return accidconv;
}


// converts general string sequence identifier to taxonomic id
typedef AccessIDConverter< std::string > StrIDConverter;
typedef AccessIDConverterFlatfileMemory< std::string > StrIDConverterFlatfileMemory;


// alias function (TODO: a function pointer might be better)
StrIDConverter* loadStrIDConverterFromFile( const std::string& filename, unsigned int cachesize = 0 );



#endif // accessconv_hh_
