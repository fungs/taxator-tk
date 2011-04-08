#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <cmath>
#include "src/taxontree.hh"
#include "src/ncbidata.hh"
#include "src/utils.hh"
#include "src/constants.hh"



using namespace std;



int main( int argc, char** argv ) {
	vector< string > ranks;
	bool delete_unmarked, rank_distances;
	string predictionsfile, correctionsfile;
	bool single_mode = true;

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

	if( ! vm.count( "ranks" ) ) { //set to fallback if not given
		ranks = default_ranks;
	}

	istream* predictions = NULL;
	if( vm.count( "predictions-file" ) ) {
		predictions = new ifstream( predictionsfile.c_str() );
	} else {
		predictions = &cin;
	}

	if( vm.count( "corrections-file" ) ) {
		single_mode = false;
	}

	// create taxonomy
	Taxonomy* tax = loadTaxonomyFromEnvironment( &ranks );
	if( !tax ) {
		return EXIT_FAILURE;
	}

	if( delete_unmarked ) {
		tax->deleteUnmarkedNodes();
		if( rank_distances ) { //does only make sense with the deletion before!
			tax->setRankDistances( ranks );
		}
	}

	TaxonomyInterface interface( tax );
	int a, b, c;

	if( single_mode ) {
		// parse line by line
		string line;
		list< string > fields;
		list< string >::iterator field_it;
		unsigned int query_taxid, prediction_taxid;

		while( getline( *predictions, line ) ) {
			if( ignoreLine( line ) ) { continue; }

			tokenizeSingleCharDelim( line, fields, FSEP, 2 );
			field_it = fields.begin();

			try {
				const string& qid = *field_it;
				query_taxid = boost::lexical_cast< unsigned int >( *field_it++ );
				prediction_taxid = boost::lexical_cast< unsigned int >(  *field_it );
				boost::tie( a, b, c ) = interface.getInterDistances( query_taxid, prediction_taxid );
				cout << qid << FSEP << a << FSEP << b << FSEP << c << endl;
			} catch( boost::bad_lexical_cast e ) {
				cerr << "Could not parse line in input: " << line << endl;
			}
			fields.clear();
		}
	} else {
		ifstream corrections( correctionsfile.c_str() );
		string line1, line2;
		list< string > fields1, fields2;
		list< string >::iterator field_it1, field_it2;
		unsigned int query_taxid, prediction_taxid;

		while( getline( *predictions, line1 ) && getline( corrections, line2 ) ) {
			if( ignoreLine( line1 ) || ignoreLine( line2 ) ) {
				continue;
			}
			tokenizeSingleCharDelim( line1, fields1, FSEP, 2 );
			field_it1 = fields1.begin();
			try {
				prediction_taxid = boost::lexical_cast< unsigned int >(  *++field_it1 );
			} catch( boost::bad_lexical_cast e ) {
				cerr << "Could not parse line in input: " << line1 << " in file '" << predictionsfile << "'" << endl;
				fields1.clear();
				continue;
			}
			tokenizeSingleCharDelim( line2, fields2, FSEP, 2 );
			field_it2 = fields2.begin();
			try {
				const string& qid = *field_it2++;
				query_taxid = boost::lexical_cast< unsigned int >( *field_it2 );
				boost::tie( a, b, c ) = interface.getInterDistances( query_taxid, prediction_taxid );
				cout << qid << FSEP << a << FSEP << b << FSEP << c << endl;
			} catch( boost::bad_lexical_cast e ) {
				cerr << "Could not parse line in input: " << line2 << " in file '" << correctionsfile << "'" << endl;
			}
			fields1.clear();
			fields2.clear();
		}
	}

	// tidy up and quit
	if( predictions != &cin ) {
		delete predictions;
	}
	delete tax;
	return EXIT_SUCCESS;
}
