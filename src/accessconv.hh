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
#include "sqlite3pp.hh"
#include "constants.hh"
#include "types.hh"
#include "utils.hh"



// converts from access identifier to taxonomic id
template< typename TypeT >
class AccessIDConverter {
	public:
		virtual ~AccessIDConverter() {};
		virtual TaxonID operator[]( const TypeT& acc ) /*throw( std::out_of_range )*/ = 0;
};



template< typename TypeT >
class AccessIDConverterSQLite : public AccessIDConverter< TypeT > {
	public:
		AccessIDConverterSQLite( const std::string& database_filename ) {
      db = new sqlite3pp::database( database_filename.c_str() );
      initialization_error = db->error_code();
    };

		~AccessIDConverterSQLite() {
		  delete db;
    };

		TaxonID operator[]( const TypeT& acc ) /*throw( std::out_of_range)*/ {
      const std::string sql_statement = boost::str( boost::format("select ncbi_taxon_id from ncbi_gitaxid_mapping where ncbi_gi = %1%" ) % acc );
      sqlite3pp::query qry( *db, sql_statement.c_str() );
      sqlite3pp::query::iterator q_it = qry.begin();
      if( q_it == qry.end() ) {
        throw std::out_of_range( sql_statement );
      }

      return q_it->get< int >( 0 );
    };

		bool fail( bool print_message = false ) {
		  bool code = db->error_code();
      if( code && print_message ) {
        std::cerr << db->error_msg() << std::endl;
        return true;
      }
      return false;
    };

 	private:
		sqlite3pp::database* db;
		bool initialization_error;
};



template< typename TypeT >
class AccessIDConverterSQLiteCache : public AccessIDConverterSQLite< TypeT > {
	public:
		AccessIDConverterSQLiteCache( const std::string& database_filename, unsigned int cachesize = 0 ) : AccessIDConverterSQLite< TypeT >( database_filename ), last_lookup( std::make_pair(TypeT(),0) ) {}; //TODO: initialize differently

		TaxonID operator[]( const TypeT& acc ) /*throw( std::out_of_range )*/ {
      // check for repeated lookup
      if( acc == last_lookup.first ) {
        return last_lookup.second;
      }

      // search in cache
      typename std::map< TypeT, TaxonID >::iterator cache_it = cache.find( acc );
      if( cache_it != cache.end() ) { //if found in cache
        return cache_it->second;
      }

      // fall back on subclass operator[]
      TaxonID taxid = AccessIDConverterSQLite< TypeT >::operator[]( acc ); //can throw exception if not found

      // update cache
      if( max_cache_size && history.size() >= max_cache_size ) {
        cache.erase( history.back() );
        history.pop();
      }

      history.push( cache.insert( std::make_pair( acc, taxid ) ).first );
      return taxid;
		}

	private:
		unsigned int max_cache_size;
		typename std::pair< TypeT, TaxonID > last_lookup; //constant lookup time
		typename std::map< TypeT, TaxonID > cache; //logarithmic lookup time (in memory)
		typename std::queue< typename std::map< TypeT, TaxonID >::iterator > history;
};



template< typename TypeT >
class AccessIDConverterFlatfileMemory : public AccessIDConverter< TypeT > {
	public:
		AccessIDConverterFlatfileMemory( const std::string& flatfile_filename ) {
		  parse( flatfile_filename );
    }

		TaxonID operator[]( const TypeT& acc ) /*throw( std::out_of_range )*/ {
      typename std::map< TypeT, TaxonID >::iterator it = accessidconv.find( acc );
      if( it == accessidconv.end() ) {
        throw std::out_of_range( boost::lexical_cast<std::string>( acc ) );
      }
      return it->second;
		}

	private:
		void parse( const std::string& flatfile_filename ) {
      std::list< std::string > fields;
      std::list< std::string >::iterator field_it;
      std::string line;
      std::ifstream flatfile( flatfile_filename.c_str() );
      TypeT acc;
      TaxonID taxid;
      while( std::getline( flatfile, line ) ) {
        if( ignoreLine( line ) ) { continue; }
        fields.clear();
        tokenizeSingleCharDelim( line, fields, default_field_separator, 2 );
        field_it = fields.begin();

        try {
          acc = boost::lexical_cast< TypeT >( *field_it );
          ++field_it;
          taxid = boost::lexical_cast< TaxonID >( *field_it );
        } catch( boost::bad_lexical_cast e ) {
          std::cerr << "Could not parse line: " << line << ", skipping..." << std::endl;
          std::cerr << "key:" << acc << std::endl;
          std::cerr << "taxid:" << taxid << std::endl;
          std::cerr << "error parsing: " << *field_it << std::endl;
          throw e;
          continue;
        }

        accessidconv[ acc ] = taxid;
      }
      flatfile.close();
		};

		typename std::map< TypeT, TaxonID > accessidconv;
};



template< typename TypeT >
AccessIDConverter< TypeT >* loadAccessIDConverterFromFile( const std::string& filename, unsigned int cachesize = 0 ) {
  AccessIDConverter< TypeT >* accidconv;
  if( filename.substr( filename.size() - 9 ) == ".sqlitedb" ) {
    accidconv = new AccessIDConverterSQLiteCache< TypeT >( filename, cachesize );
  } else {
    accidconv = new AccessIDConverterFlatfileMemory< TypeT >( filename );
  }
  return accidconv;
}


// converts general string sequence identifier to taxonomic id
typedef AccessIDConverter< std::string > StrIDConverter;

typedef AccessIDConverterSQLite< std::string > StrIDConverterSQLite;

typedef AccessIDConverterSQLiteCache< std::string > StrIDConverterSQLiteCache;

typedef AccessIDConverterFlatfileMemory< std::string > StrIDConverterFlatfileMemory;


// alias function (TODO: a function pointer might be better)
StrIDConverter* loadStrIDConverterFromFile( const std::string& filename, unsigned int cachesize = 0 );



#endif // accessconv_hh_
