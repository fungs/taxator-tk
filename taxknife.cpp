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
#include <iostream>
#include <sstream>
#include <vector>
#include <set>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/scoped_ptr.hpp>
#include <assert.h>
#include "src/taxonomyinterface.hh"
#include "src/ncbidata.hh"
#include "src/utils.hh"
#include "src/constants.hh"

using namespace std;

int main( int argc, char** argv ) {
	string show_what, invalid_replace_value_traversal, invalid_replace_value_annotation;
	vector< string > rank_names, ranks, operations;
	unsigned int field_pos;
	bool allnodes = false;

	namespace po = boost::program_options;

    po::positional_options_description mode;
    mode.add("mode",2);

	po::options_description desc("Allowed options");
	desc.add_options()
	( "help,h", "show help message")
	( "mode", po::value< vector <string> >( &operations),"choose mode:\n"
                                          "\"traversal\": follow nodes upwards in taxonomy\n\n"
                                          "\"annotation\": looks up metainformation attached to nodes (e.g. names)\n\n")
	( "field,f", po::value< unsigned int >( &field_pos )->default_value( 1 ), "input column\n" );

	// TODO: put option parsing in separate objects which can be chained
	po::options_description rankf("traversal mode");
	rankf.add_options()
	( "keep-not-rank,k", "unmappable taxids remain (otherwise mapped to root)" )
	( "keep-not-taxid,t", "unknown taxonomic IDs are kept (otherwise skipped)" )
	( "set-invalid-traversal,b", po::value< string >( &invalid_replace_value_traversal ),"replace unknown taxids by this given value" )
	( "traverse-ranks,r", po::value< vector <string> >( &rank_names)-> multitoken(),"traverse taxonomy up to one of these ranks (space separated list)");

	// TODO: put option parsing in separate objects which can be chained
	po::options_description namef("annotation mode");
	namef.add_options()
	( "allnodes,a", "if set, all nodes will be used, not only at selected ranks" )
	( "set-invalid-annotation,c", po::value< string >( &invalid_replace_value_annotation ),"replace all taxids that are invalid by this given value" )
	( "show,s", po::value< string >( &show_what )->default_value( "name" ),"either 'name', 'rank', 'path' or 'taxid-path'" )
	( "name-ranks,n", po::value< vector <string> >( &ranks)-> multitoken(),"select ranks to be considered; if not set, default ranks will be used");

    desc.add(rankf).add(namef);

    po::variables_map vm;
    po::store(po::command_line_parser( argc, argv ).options( desc ).positional(mode).run(), vm);
    po::notify(vm);

		// global sanity checks
    if( vm.count("help")) { cout << desc << endl; return EXIT_SUCCESS; }
    if(operations.empty()) { cout << "\n Please choose a mode.\n" << endl; cout << desc <<endl; return EXIT_FAILURE; }
    if( field_pos < 1 ) { cerr << "Field number index is 1-based" << endl; return EXIT_FAILURE; }

    // build taxonomy
    boost::scoped_ptr< Taxonomy > tax( loadTaxonomyFromEnvironment( &default_ranks ) );
    if( ! tax ) return EXIT_FAILURE;
    TaxonomyInterface interface( tax.get() );

    if(find(operations.begin(), operations.end(),"traversal") != operations.end()){
        // command line arguments

//         if( ! vm.count( "traversal" ) ) { //TODO: always false?
//             cout << "Required argument \"--traversal\" is not given.\n" << desc << endl;
//             return EXIT_FAILURE;
//         }

        bool keep_not_rank = vm.count( "keep-not-rank" );
        bool keep_not_taxid = vm.count( "keep-not-taxid" );
        bool replace_invalid = vm.count( "set-invalid-traversal" );

        // create taxonomy
        /*
        boost::scoped_ptr< Taxonomy > tax( loadTaxonomyFromEnvironment( &default_ranks ) );
        if( ! tax ) return EXIT_FAILURE;
        TaxonomyInterface interface( tax.get() );
        */

        //TODO: change code to use set of ranks, not single rank to traverse
        set< const string* > ranks;
        for (vector< string >::iterator it = rank_names.begin(); it != rank_names.end(); ++it ) {
            const string& rank = tax->getRankInternal( *it );
            if( rank.empty() ) cerr << "Rank '" << *it << "' not found in taxonomy, not using for mapping..." << endl;
            else ranks.insert( &rank );
        }

        // parse line by line
        string line;
        stringstream buffer;
        list< string > fields;
        list< string >::iterator field_it;
        unsigned int taxid;
        const TaxonNode* rootnode = interface.getRoot();
        const TaxonNode* node;

        while( getline( cin, line ) ) {
            if( ignoreLine( line ) ) continue;

            tokenizeSingleCharDelim( line, fields, default_field_separator, field_pos );
            field_it = fields.begin();
            unsigned int i = 1;
            while( field_it != fields.end() ) {
          if( i < field_pos ) {
            buffer << *field_it++ << default_field_separator;
            ++i;
          } else {
            try {
              //cerr << "Converting " << *field_it << " to unsigned int" << endl;
              taxid = boost::lexical_cast< unsigned int >( *field_it );
              //cerr << "Converted without exception!" << endl;
              node = interface.getNode( taxid );
              if( node ) {
                while( ! node->data->annotation || ( ! ranks.count( &(node->data->annotation->rank) ) && node != rootnode ) ) {
                  node = node->parent;
                }
                if( keep_not_rank && node == rootnode ) {
                  cout << buffer.str();
                  if( replace_invalid ) {
                    cout << invalid_replace_value_traversal;
                  } else {
                                    cout << taxid;
                  }
                } else {
                  cout << buffer.str() << node->data->taxid;
                }
                if( (++field_it)->empty() ) {
                  cout << endl;
                } else {
                  cout << default_field_separator << *field_it << endl;
                }
              } else {
                cerr << "traversal: Could not find node with taxid " << *field_it << " in the taxonomy";
                if( keep_not_taxid ) {
                  cerr << endl;
                  cout << buffer.str();
                  if( replace_invalid ) { //TODO: only works in combination with keep-invalid
                    cout << invalid_replace_value_traversal;
                  } else {
                    cout << taxid;
                  }
                  if( ! (++field_it)->empty() ) {
                    cout << default_field_separator << *field_it;
                  }
                  cout << endl;
                } else {
                  cerr << ", skipping record..." << endl;
                }
              }
            } catch( boost::bad_lexical_cast e ) {
              cerr << "traversal: Could not parse taxid " << *field_it << " in line \"" << line << "\", skipping record..." << endl;
            }
            break;
          }
        }
            fields.clear();
            buffer.str("");
            buffer.clear();
        }

        return EXIT_SUCCESS;
	}
    else if(find(operations.begin(), operations.end(),"annotation") != operations.end()){

        if( ! vm.count( "name-ranks" ) ) ranks = default_ranks;

        if( ! vm.count( "allnodes" ) ) allnodes = false;
        else allnodes = true;

        bool replace_invalid = vm.count( "set-invalid-annotation" );

        // parse line by line
        string line;
        list< string > fields;
        list< string >::iterator field_it;
        unsigned int taxid;
        const TaxonNode* node;
        stringstream buffer;

        if ( show_what == "name" ) {
            while( getline( cin, line ) ) {
                if ( ignoreLine( line ) ) continue;
                tokenizeSingleCharDelim( line, fields, default_field_separator, field_pos );
                field_it = fields.begin();
                unsigned int i = 1;
                while( field_it != fields.end() ) {
                  if( i < field_pos ) {
              buffer << *field_it++ << default_field_separator;
              ++i;
                  } else {
                    try {
                      taxid = boost::lexical_cast< unsigned int >( *field_it );
                      node = interface.getNode( taxid );
                if( node ) {
                  if( node->data->annotation ) {
                    cout << buffer.str() << node->data->annotation->name;
                    if( ! (++field_it)->empty() ) {
                      cout << default_field_separator << *field_it;
                    }
                    cout << endl;
                  } else {
                    cout << "node_without_annotation";
                  }
                } else {
                  cerr << "Could not find node with taxonomic id " << taxid << " in taxonomy" << endl;
                  if ( replace_invalid ) {
                    cout << buffer.str() << invalid_replace_value_annotation;
                    if( ! (++field_it)->empty() ) {
                      cout << default_field_separator << *field_it;
                    }
                    cout << endl;
                  }
                }
                    } catch( boost::bad_lexical_cast e ) {
                cerr << "Could not parse taxonomic id from field \"" << *field_it << '\"' << endl;
                if ( replace_invalid ) {
                  cout << buffer.str() << invalid_replace_value_annotation;
                  if( ! (++field_it)->empty() ) {
                    cout << default_field_separator << *field_it;
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
                    if ( ignoreLine( line ) ) continue;
                    tokenizeSingleCharDelim( line, fields, default_field_separator, field_pos );
                    field_it = fields.begin();
            unsigned int i = 1;
            while( field_it != fields.end() ) {
              if( i < field_pos ) {
                buffer << *field_it++ << default_field_separator;
                ++i;
              } else {
                try {
                  taxid = boost::lexical_cast< unsigned int >( *field_it );
                  node = interface.getNode( taxid );
                  if( node ) {
                    if( node->data->annotation ) {
                      cout << buffer.str() << node->data->annotation->rank;
                      if( ! (++field_it)->empty() ) {
                        cout << default_field_separator << *field_it;
                      }
                      cout << endl;
                    } else {
                      cout << "node_without_annotation";
                    }
                  } else {
                    cerr << "no taxon with taxid " << taxid << " found in taxonomy" << endl;
                   if ( replace_invalid ) {
                      cout << buffer.str() << invalid_replace_value_annotation;
                      if( ! (++field_it)->empty() ) {
                        cout << default_field_separator << *field_it;
                      }
                      cout << endl;
                    }
                  }
                } catch( boost::bad_lexical_cast e ) {
                  cerr << "Could not parse taxonomic id from field \"" << *field_it << '\"' << endl;
                  if ( replace_invalid ) {
                    cout << buffer.str() << invalid_replace_value_annotation;
                    if( ! (++field_it)->empty() ) {
                      cout << default_field_separator << *field_it;
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
                if ( show_what == "path" ) {
                    while( getline( cin, line ) ) {
                        if ( ignoreLine( line ) ) continue;
                        tokenizeSingleCharDelim( line, fields, default_field_separator, field_pos );
                        field_it = fields.begin();
                        unsigned int i = 1;
                        while( field_it != fields.end() ) {
                            if( i < field_pos ) {
                                buffer << *field_it++ << default_field_separator;
                                ++i;
                            } else {
                                try {
                                    taxid = boost::lexical_cast< unsigned int >( *field_it );
                                    node = interface.getNode( taxid );
                                    if( node ) {
                                        cout << buffer.str();
                                        const TaxonNode* root = interface.getRoot();
                                        for ( Taxonomy::CPathDownIterator it( root, node ); it != node; ++it ) {
                                            if ( allnodes || it->data->mark_special ) {
                                                if( it->data->annotation ) {
                                                    cout << it->data->annotation->name << ';';
                                                } else {
                                                    cout << "node_without_annotation;";
                                                }
                                            }
                                        }
                                        if ( allnodes || node->data->mark_special ) {
                                            cout << node->data->annotation->name << ';';
                                        }
                                        if( ! (++field_it)->empty() ) {
                                            cout << default_field_separator << *field_it;
                                        }
                                            cout << endl;
                                    } else {
                                        cerr << "no taxon with taxid " << taxid << " found in taxonomy" << endl;
                                    if ( replace_invalid ) {
                                            cout << buffer.str() << invalid_replace_value_annotation;
                                            if( ! (++field_it)->empty() ) {
                                                cout << default_field_separator << *field_it;
                                            }
                                            cout << endl;
                                        }
                                    }
                                } catch( boost::bad_lexical_cast e ) {
                                    cerr << "Could not parse taxonomic id from field \"" << *field_it << '\"' << endl;
                                    if ( replace_invalid ) {
                                        cout << buffer.str() << invalid_replace_value_annotation;
                                        if( ! (++field_it)->empty() ) {
                                            cout << default_field_separator << *field_it;
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
                    if ( show_what == "taxid-path" ) {
                            while( getline( cin, line ) ) {
                            if ( ignoreLine( line ) ) continue;
                            tokenizeSingleCharDelim( line, fields, default_field_separator, field_pos );
                            field_it = fields.begin();
                            unsigned int i = 1;
                            while( field_it != fields.end() ) {
                                if( i < field_pos ) {
                                    buffer << *field_it++ << default_field_separator;
                                    ++i;
                                } else {
                                    try {
                                        taxid = boost::lexical_cast< unsigned int >( *field_it );
                                        node = interface.getNode( taxid );
                                        if( node ) {
                                            cout << buffer.str();
                                            const TaxonNode* root = interface.getRoot();
                                            for ( Taxonomy::CPathDownIterator it( root, node ); it != node; ++it ) {
                                                if ( allnodes || it->data->mark_special ) {
                                                    cout << it->data->taxid << ';';
                                                }
                                            }
                                            if ( allnodes || node->data->mark_special ) {
                                                cout << node->data->taxid << ';';
                                            }
                                            if( ! (++field_it)->empty() ) {
                                                cout << default_field_separator << *field_it;
                                            }
                                                cout << endl;
                                        } else {
                                            cerr << "no taxon with taxid " << taxid << " found in taxonomy" << endl;
                                        if ( replace_invalid ) {
                                                cout << buffer.str() << invalid_replace_value_annotation;
                                                if( ! (++field_it)->empty() ) {
                                                    cout << default_field_separator << *field_it;
                                                }
                                                cout << endl;
                                            }
                                        }
                                    } catch( boost::bad_lexical_cast e ) {
                                        cerr << "Could not parse taxonomic id from field \"" << *field_it << '\"' << endl;
                                        if ( replace_invalid ) {
                                            cout << buffer.str() << invalid_replace_value_annotation;
                                            if( ! (++field_it)->empty() ) {
                                                cout << default_field_separator << *field_it;
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
            }
        }

        return EXIT_SUCCESS;
    }
}
