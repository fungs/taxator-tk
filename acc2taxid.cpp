#include <iostream>
#include <string>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include "src/ncbidata.hh"
#include "src/alignmentrecord.hh"
#include "src/utils.hh"
#include "src/constants.hh"



using namespace std;



int main( int argc, char** argv ) {
//  std::string access_id_type;
  std::string conversion_filename;
  unsigned int cache_size; //standard cache size (number of elements)
  unsigned int field_pos;
  stringstream buffer;

  namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
	( "help,h", "show help message")
//	( "type,t", po::value< string >( &access_id_type )->default_value( "ncbi" ), "specify identifier type (NCBI or EMBL/GenBank)" )
  ( "field,f", po::value< unsigned int >( &field_pos )->default_value( 1 ), "column number to use" )
  ( "mapping-file,m", po::value< string >( &conversion_filename ), "give name to mappings file (>=2-col-tab-separated or specific SQLite DB \"*.sqlite\" file)" )
  ( "cache,c", po::value< unsigned int >( &cache_size )->default_value( 100 ), "number of records to cache (effective for SQLite files only)" );

	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( desc ).run(), vm);
	po::notify(vm);

	//std::transform( access_id_type.begin(), access_id_type.end(), access_id_type.begin(), ::tolower );

	// command line arguments
  if( ! vm.count( "field" ) ) {
		cout << desc << endl;
		return EXIT_FAILURE;
	}

  if( ! vm.count( "mapping-file" ) ) {
		cout << desc << endl;
		return EXIT_FAILURE;
	}

	if( field_pos < 1 ) {
	  cerr << "Field number index is 1-based" << endl;
	  return EXIT_FAILURE;
  }

	// parse line by line
	string line;
	list< string > fields;
	list< string >::iterator field_it;

  StrIDConverter* accessconverter = loadStrIDConverterFromFile( conversion_filename, cache_size );
  StrIDConverter& seqid2taxid = *accessconverter;

  while( getline( cin, line ) ) {
    if( ignoreLine( line ) ) { continue; }

    tokenizeSingleCharDelim( line, fields, FSEP, field_pos );

    if( fields.size() < field_pos ) {
		  cout << line << endl; //leave line unchanged
		} else {
      field_it = fields.begin();
      for( unsigned int i = 1; i < field_pos; ++i ) {
          buffer << *field_it++ << FSEP;
      }

      std::string acc = *field_it++;
      try {
        cout << buffer.str() << seqid2taxid[ acc ];
        if( ! field_it->empty() ) { //don't print empty field at the end
          cout << FSEP << *field_it;
        }
        cout << endl;
      } catch ( out_of_range e ) {
        cerr << "acc2taxid: Could not map sequence identifier " << acc << " to taxonomic identifier";
        cerr << ", skipping record..." << endl;
      }
		}
    fields.clear();
    buffer.str("");
    buffer.clear();
  }

  delete accessconverter; //tidy up
 // }
	return EXIT_SUCCESS;
}
