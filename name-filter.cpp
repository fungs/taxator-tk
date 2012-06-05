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

	string show_what, invalid_replace_value;
	unsigned int field_pos;
	vector< string > ranks;
	bool allnodes = false;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
	( "help,h", "show help message")
	( "field,f", po::value< unsigned int >( &field_pos )->default_value( 1 ), "column number to use" )
	( "ranks,r", po::value< vector< string > >( &ranks )->multitoken(), "select ranks to be considered; if not set, default ranks will be used" )
	( "allnodes,a", "if set, all nodes will be used, not only at selected ranks" )
	( "set-invalid-value,b", po::value< string >( &invalid_replace_value ),"replace all taxids that are invalid by this given value" )
	( "show,s", po::value< string >( &show_what )->default_value( "name" ), "either 'name', 'rank', 'path' or 'taxid-path'" );

	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( desc ).run(), vm);
	po::notify(vm);

	if( vm.count( "help" ) ) {
		cout << desc << endl;
		return EXIT_SUCCESS;
	}
	
	if( ! vm.count( "ranks" ) ) ranks = default_ranks;
	
	if( ! vm.count( "allnodes" ) ) allnodes = false;
	else allnodes = true;

	if( field_pos < 1 ) {
	  cerr << "Field number index is 1-based" << endl;
	  return EXIT_FAILURE;
  }

	bool invalid_replace = vm.count( "set-invalid-value" );

	// create taxonomy
	boost::scoped_ptr< Taxonomy > tax( loadTaxonomyFromEnvironment( &ranks ) );
	if( ! tax ) return EXIT_FAILURE;

	TaxonomyInterface interface( tax.get() );

	// parse line by line
	string line;
	list< string > fields;
	list< string >::iterator field_it;
	unsigned int taxid;
	const TaxonNode* node;
	stringstream buffer;

	if ( show_what == "name" ) {
		while( getline( cin, line ) ) {
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
              if ( invalid_replace ) {
                cout << buffer.str() << invalid_replace_value;
                if( ! (++field_it)->empty() ) {
                  cout << default_field_separator << *field_it;
                }
                cout << endl;
              }
            }
			    } catch( boost::bad_lexical_cast e ) {
            cerr << "Could not parse taxonomic id from field \"" << *field_it << '\"' << endl;
            if ( invalid_replace ) {
              cout << buffer.str() << invalid_replace_value;
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
               if ( invalid_replace ) {
                  cout << buffer.str() << invalid_replace_value;
                  if( ! (++field_it)->empty() ) {
                    cout << default_field_separator << *field_it;
                  }
                  cout << endl;
                }
              }
            } catch( boost::bad_lexical_cast e ) {
              cerr << "Could not parse taxonomic id from field \"" << *field_it << '\"' << endl;
              if ( invalid_replace ) {
                cout << buffer.str() << invalid_replace_value;
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
								if ( invalid_replace ) {
										cout << buffer.str() << invalid_replace_value;
										if( ! (++field_it)->empty() ) {
											cout << default_field_separator << *field_it;
										}
										cout << endl;
									}
								}
							} catch( boost::bad_lexical_cast e ) {
								cerr << "Could not parse taxonomic id from field \"" << *field_it << '\"' << endl;
								if ( invalid_replace ) {
									cout << buffer.str() << invalid_replace_value;
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
									if ( invalid_replace ) {
											cout << buffer.str() << invalid_replace_value;
											if( ! (++field_it)->empty() ) {
												cout << default_field_separator << *field_it;
											}
											cout << endl;
										}
									}
								} catch( boost::bad_lexical_cast e ) {
									cerr << "Could not parse taxonomic id from field \"" << *field_it << '\"' << endl;
									if ( invalid_replace ) {
										cout << buffer.str() << invalid_replace_value;
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
