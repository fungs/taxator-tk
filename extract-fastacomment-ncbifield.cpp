#include <iostream>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include "src/ncbidata.hh"
#include "src/utils.hh"
#include "src/constants.hh"



using namespace std;



int main( int argc, char** argv ) {

	string key;
	unsigned int field_pos;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		("help,h", "show help message")
		("key,k", po::value<string>(&key)->default_value( "gi" ), "set key for which to extract the value")
		( "field,f", po::value< unsigned int >( &field_pos )->default_value( 4 ), "column number to use" );

	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( desc ).run(), vm);
	po::notify(vm);

	if( vm.count( "help" ) ) {
		cout << desc << endl;
		return EXIT_SUCCESS;
	}

	if( field_pos < 1 ) {
	  cerr << "Field number index is 1-based" << endl;
	  return EXIT_FAILURE;
  }

	// parse line by line
	string line;
	list< string > fields;
	list< string >::iterator field_it;

	while( getline( cin, line ) ) {
		if( ignoreLine( line ) ) { continue; }

		tokenizeSingleCharDelim( line, fields, FSEP, field_pos );

		if( fields.size() < field_pos ) {
		  cout << line << endl; //leave line unchanged
		} else {
      field_it = fields.begin();
      for( unsigned int i = 1; i < field_pos; ++i ) {
          cout << *field_it++ << FSEP;
      }

      cout << extractFastaCommentField( *field_it++, key );

      if( ! field_it->empty() ) { //don't print empty field at the end
        cout << FSEP << *field_it;
      }
      cout << endl;
		}
    fields.clear();
	}

	return EXIT_SUCCESS;
}
