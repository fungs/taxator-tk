#include <iostream>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <assert.h>
#include "src/taxontree.hh"
#include "src/ncbidata.hh"
#include "src/utils.hh"
#include "src/constants.hh"



using namespace std;



int main( int argc, char** argv ) {
	std::string rank_name, invalid_replace_value;
	unsigned int field_pos;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
	( "help,h", "show help message")
	( "keep-not-rank,k", "unmappable taxids stay the same (otherwise mapped to root node = 1)" )
	( "keep-not-taxid,t", "taxonomic IDs not found will be kept (otherwise line is skipped)" )
	( "set-invalid-value,b", po::value< string >( &invalid_replace_value ),"replace all taxids that are invalid by this given value" )
	( "field,f", po::value< unsigned int >( &field_pos )->default_value( 1 ), "input column number to use" )
	( "rank,r", po::value< string >( &rank_name ), "traverse taxonomy up to this rank" );

	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( desc ).run(), vm);
	po::notify(vm);

	// command line arguments
	if( ! vm.count( "rank" ) ) {
		cout << desc << endl;
		return EXIT_SUCCESS;
	}

  if( ! vm.count( "rank" ) ) {
		cout << desc << endl;
		return EXIT_FAILURE;
	}

  if( field_pos < 1 ) {
	  cerr << "Field number index is 1-based" << endl;
	  return EXIT_FAILURE;
  }

	bool keep_not_rank = vm.count( "keep-not-rank" );
  bool keep_not_taxid = vm.count( "keep-not-taxid" );
  bool invalid_replace = vm.count( "set-invalid-value" );

	// create taxonomy
	Taxonomy* tax = loadTaxonomyFromEnvironment();
	TaxonomyInterface interface( tax );

	// string reference comparison is fastest
	const string& rank = tax->getRankInternal( rank_name );

	if( rank.empty() ) {
		cerr << "Rank '" << rank_name << "' not found in taxonomy, expect identity mapping..." << endl;
	}

	// parse line by line
	string line;
	stringstream buffer;
	list< string > fields;
	list< string >::iterator field_it;
	unsigned int taxid;
	TaxonNode* rootnode = interface.getRoot();
	TaxonNode* node;

	while( getline( cin, line ) ) {
		if( ignoreLine( line ) ) { continue; }

		tokenizeSingleCharDelim( line, fields, FSEP, field_pos );
		field_it = fields.begin();
		unsigned int i = 1;
		while( field_it != fields.end() ) {
      if( i < field_pos ) {
        buffer << *field_it++ << FSEP;
        ++i;
      } else {
        try {
          //cerr << "Converting " << *field_it << " to unsigned int" << endl;
          taxid = boost::lexical_cast< unsigned int >( *field_it );
          //cerr << "Converted without exception!" << endl;
          node = interface.getNode( taxid );
          if( node ) {
            while( ! node->data->annotation || ( &(node->data->annotation->rank) != &rank && node != rootnode ) ) {
              node = node->parent;
            }
            if( keep_not_rank && node == rootnode ) {
              cout << buffer.str();
              if( invalid_replace ) {
                cout << taxid;
              } else {
                cout << invalid_replace_value;
              }
            } else {
              cout << buffer.str() << node->data->taxid;
            }
            if( (++field_it)->empty() ) {
              cout << endl;
            } else {
              cout << FSEP << *field_it << endl;
            }
          } else {
            cerr << "rank-filter: Could not find node with taxid " << *field_it << " in the taxonomy";
            if( keep_not_taxid ) {
              cerr << endl;
              cout << buffer.str();
              if( invalid_replace ) {
                cout << invalid_replace_value;
              } else {
                cout << taxid;
              }
              if( ! (++field_it)->empty() ) {
                cout << FSEP << *field_it;
              }
              cout << endl;
            } else {
              cerr << ", skipping record..." << endl;
            }
          }
        } catch( boost::bad_lexical_cast e ) {
          cerr << "rank-filter: Could not parse taxid " << *field_it << " in line \"" << line << "\", skipping record..." << endl;
        }
        break;
      }
    }
		fields.clear();
		buffer.str("");
		buffer.clear();
	}

	// tidy up and quit
	delete tax;
	return EXIT_SUCCESS;
}
