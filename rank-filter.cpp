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
	string invalid_replace_value;
	vector< string > rank_names;
	unsigned int field_pos;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
	( "help,h", "show help message")
	( "keep-not-rank,k", "unmappable taxids remain (otherwise mapped to root)" )
	( "keep-not-taxid,t", "unknown taxonomic IDs are kept (otherwise skipped)" )
	( "set-invalid-value,b", po::value< string >( &invalid_replace_value ),"replace unknown taxids by this given value" )
	( "field,f", po::value< unsigned int >( &field_pos )->default_value( 1 ), "input column number to use" )
	( "ranks,r", po::value< vector< string > >( &rank_names )->multitoken(), "traverse taxonomy up to one of these rank (space separated list)" );

	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( desc ).run(), vm);
	po::notify(vm);

	// command line arguments
	if( ! vm.count( "ranks" ) ) {
		cout << desc << endl;
		return EXIT_SUCCESS;
	}

  if( ! vm.count( "ranks" ) ) {
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
	boost::scoped_ptr< Taxonomy > tax( loadTaxonomyFromEnvironment( &default_ranks ) );
	if( ! tax ) return EXIT_FAILURE;
	TaxonomyInterface interface( tax.get() );

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
              if( invalid_replace ) {
                cout << invalid_replace_value;
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
                cout << default_field_separator << *field_it;
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

	return EXIT_SUCCESS;
}
