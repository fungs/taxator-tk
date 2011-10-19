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
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/scoped_ptr.hpp>
#include <cmath>
#include "src/predictionrecord.hh"
#include "src/taxonomyinterface.hh"
#include "src/ncbidata.hh"
#include "src/utils.hh"
#include "src/constants.hh"
#include "src/types.hh"



using namespace std;



int main( int argc, char** argv ) {
	vector< string > ranks;
	bool delete_unmarked, rank_distances;
	string predictionsfile, correctionsfile;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
	( "help,h", "show help message")
	( "ranks,r", po::value< vector< string > >( &ranks )->multitoken(), "set node ranks at which to do predictions" )
	( "delete-notranks,d", po::value< bool >( &delete_unmarked )->default_value( true ), "delete all nodes that don't have any of the given ranks" )
	( "normalized-tree,n", po::value< bool >( &rank_distances )->default_value( true  ), "set distances to correct for a balanced tree by ncbi major ranks given" )
	( "predictions-file,p", po::value< string >( &predictionsfile ), "set two-column predictions file to convert into distances using prediction query-label" )
	( "corrections-file,c", po::value< string >( &correctionsfile ), "set two-column correction predictions file as corrected query-labels" );

	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( desc ).run(), vm);
	po::notify(vm);

	if( vm.count( "help" ) ) {
		cout << desc << endl;
		return EXIT_FAILURE;
	}

	if( ! vm.count( "ranks" ) ) ranks = default_ranks;
	
	if( ! vm.count( "corrections-file" ) ) {
		cout << "You need to specify a file with correct taxa for distance calculation." << endl;
		return EXIT_FAILURE;
	}

	istream* predictions = NULL;
	if( vm.count( "predictions-file" ) ) {
		predictions = new ifstream( predictionsfile.c_str() );
	} else {
		predictions = &cin;
	}

	// create taxonomy
 	boost::scoped_ptr< Taxonomy > tax( loadTaxonomyFromEnvironment( &ranks ) );
	if( ! tax ) return EXIT_FAILURE;

	if( delete_unmarked ) {
		tax->deleteUnmarkedNodes();
		if( rank_distances ) { //does only make sense with the deletion before!
			tax->setRankDistances( ranks );
		}
	}

	TaxonomyInterface interface( tax.get() );
	small_unsigned_int a1, b1, c1, a2, b2, c2;

	ifstream corrections( correctionsfile.c_str() );
	string line1, line2;
	list< string > fields;
	list< string >::iterator field_it;
	unsigned int correct_taxid;
	PredictionRecord prec( tax.get() );

	if ( getline( *predictions, line1 ) && getline( corrections, line2 ) ) {
		do {
			if ( ignoreLine( line1 ) ) {
				if ( getline( *predictions, line1 ) ) continue;
				break;
			}
			
			if ( ignoreLine( line2 ) ) {
				if ( getline( corrections, line2 ) ) continue;
				break;
			}
			
			if ( ! prec.parse( line1 ) ) {
				if ( getline( *predictions, line1 ) && getline( corrections, line2 ) ) continue;
				cerr << "could not parse prediction record" << endl;
				break;
			}
			
			tokenizeSingleCharDelim( line2, fields, default_field_separator, 2 );
			field_it = fields.begin();
			const string& qid = *field_it++;
			try {
				correct_taxid = boost::lexical_cast< unsigned int >( *field_it );
			} catch( boost::bad_lexical_cast e ) {
				cerr << "could not parse taxid in input line: '" << line2 << "' in file '" << correctionsfile << "'" << endl;
				fields.clear();
				continue;
			}
			
			const TaxonNode* correct_node = interface.getNode( correct_taxid ); //TODO: test for NULL
			if ( ! correct_node ) {
				cerr << "could not find taxid '" << correct_taxid << "' in taxonomy" << endl;
				if ( getline( *predictions, line1 ) && getline( corrections, line2 ) ) continue;
				break;
			}
			
			boost::tie( a1, b1, c1 ) = interface.getInterDistances( correct_node, prec.lower_node );
			boost::tie( a2, b2, c2 ) = interface.getInterDistances( correct_node, prec.upper_node );
			
			small_unsigned_int rangesize = prec.lower_node == prec.upper_node ? 0 : interface.getPathLengthToParent( prec.lower_node, prec.upper_node );
			
			cout << qid << default_field_separator << static_cast<int>( a1 ) << default_field_separator << static_cast<int>( c2 ) << default_field_separator << static_cast<int>( rangesize ) << endl;
			
			fields.clear();
			if ( !( getline( *predictions, line1 ) && getline( corrections, line2 ) ) ) break; 
		} while ( true );
	}
	
	// tidy up and quit
	if( predictions != &cin ) { //TODO: use auto_ptr?
		delete predictions;
	}
	return EXIT_SUCCESS;
}
