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

	string show_what, invalid_replace_value;
	unsigned int field_pos;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
	( "help,h", "show help message")
	( "field,f", po::value< unsigned int >( &field_pos )->default_value( 1 ), "column number to use" )
	( "set-invalid-value,b", po::value< string >( &invalid_replace_value ),"replace all taxids that are invalid by this given value" )
	( "show,s", po::value< string >( &show_what )->default_value( "name" ), "either 'name' or 'rank'" );

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

	bool invalid_replace = vm.count( "set-invalid-value" );

	// create taxonomy
	Taxonomy* tax = loadTaxonomyFromEnvironment();
	if( !tax ) {
		return EXIT_FAILURE;
	}

	TaxonomyInterface interface( tax );

	// parse line by line
	string line;
	list< string > fields;
	list< string >::iterator field_it;
	unsigned int taxid;
	TaxonNode* node;
	stringstream buffer;

	if ( show_what == "name" ) {
		while( getline( cin, line ) ) {
			tokenizeSingleCharDelim( line, fields, FSEP, field_pos );
			field_it = fields.begin();
			unsigned int i = 1;
			while( field_it != fields.end() ) {
			  if( i < field_pos ) {
          buffer << *field_it++ << FSEP;
          ++i;
			  } else {
			    try {
			      taxid = boost::lexical_cast< unsigned int >( *field_it );
			      node = interface.getNode( taxid );
            if( node ) {
              if( node->data->annotation ) {
                cout << buffer.str() << node->data->annotation->name;
                if( ! (++field_it)->empty() ) {
                  cout << FSEP << *field_it;
                }
                cout << endl;
              } else {
                cout << "node_without_annotation";
              }
            } else {
              cerr << "no taxon with taxid " << taxid << " found in taxonomy" << endl;
              if ( invalid_replace ) {
                cout << buffer.str() << invalid_replace_value;
                if( ! (++field_it)->empty() ) {
                  cout << FSEP << *field_it;
                }
                cout << endl;
              }
            }
			    } catch( boost::bad_lexical_cast e ) {
            cerr << "Could not parse taxonomic id from field \"" << *field_it << '\"' << endl;
            if ( invalid_replace ) {
              cout << buffer.str() << invalid_replace_value;
              if( ! (++field_it)->empty() ) {
                cout << FSEP << *field_it;
              }
              cout << endl;
            }
			    }
			  break;
			  }
			}
			fields.clear();
			buffer.str("");
			buffer.clear();
		}
	} else {
		if( show_what == "rank" ) {
			while( getline( cin, line ) ) {
				tokenizeSingleCharDelim( line, fields, "\t", 1 );
				field_it = fields.begin();
        unsigned int i = 1;
        while( field_it != fields.end() ) {
          if( i < field_pos ) {
            buffer << *field_it++ << FSEP;
            ++i;
          } else {
            try {
              taxid = boost::lexical_cast< unsigned int >( *field_it++ );
              node = interface.getNode( taxid );
              if( node ) {
                if( node->data->annotation ) {
                  cout << buffer.str() << node->data->annotation->rank;
                  if( ! field_it->empty() ) {
                    cout << FSEP << *field_it;
                  }
                  cout << endl;
                } else {
                  cout << "node_without_annotation";
                }
              } else {
                cerr << "no taxon with taxid " << taxid << " found in taxonomy" << endl;
               if ( invalid_replace ) {
                  cout << buffer.str() << invalid_replace_value;
                  if( ! (++field_it)->empty() ) {
                    cout << FSEP << *field_it;
                  }
                  cout << endl;
                }
              }
            } catch( boost::bad_lexical_cast e ) {
              cerr << "Could not parse taxonomic id from field \"" << *field_it << '\"' << endl;
              if ( invalid_replace ) {
                cout << buffer.str() << invalid_replace_value;
                if( ! (++field_it)->empty() ) {
                  cout << FSEP << *field_it;
                }
                cout << endl;
              }
			      }
            break;
          }
        }
        fields.clear();
        buffer.str("");
        buffer.clear();
			}
		} else {
			cerr << "unknown parameter for --show / -s" << endl;
			cout << desc << endl;
			return EXIT_FAILURE;
		}
	}

	// tidy up and quit
	delete tax;
	return EXIT_SUCCESS;
}
