#include <iostream>
#include <sstream>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/tuple/tuple.hpp>
#include <map>
#include <list>
#include "src/utils.hh"
#include "src/constants.hh"
#include "src/losses.hh"



using namespace std;



int main( int argc, char** argv ) {
	vector< string > ranks;
	string loss_name;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
	( "help,h", "show help message")
	( "show-categories,s", "whether to show a summary with errors for each of the categories" )
	( "loss,l", po::value< string >( &loss_name ), "select loss to use (l1 and l2 supported)" );

	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( desc ).run(), vm);
	po::notify(vm);

	if( vm.count( "help" ) ) {
		cout << desc << endl;
		return EXIT_FAILURE;
	}


	bool show_categories = vm.count( "show-categories" );

	BaseLossFactory* factory = NULL;
	// set error object according to parameters
	if( loss_name == "l1" ) {
		factory = new LossFactory< L1Loss >;
	} else {
		factory = new LossFactory< L2Loss >;
	}

	// parse line by line
	string line;
	list< string > fields;
	list< string >::iterator field_it;
	LossMap losses( factory );
	BaseLoss* loss_all = (*factory)();

	while( getline( cin, line ) ) {
		if( ignoreLine( line ) ) { continue; }

		tokenizeSingleCharDelim( line, fields, FSEP, 4 );
		field_it = fields.begin();

		try {
			const string& id = *field_it;
			int a = boost::lexical_cast< int >( *++field_it );
			int b = boost::lexical_cast< int >( *++field_it );
			int c = boost::lexical_cast< int >( *++field_it );
			losses[ id ].add( a, b, c );
			loss_all->add( a, b, c );
		} catch( boost::bad_lexical_cast e ) {
			cerr << "Could not parse line in input: " << line << endl;
		}
		fields.clear();
	}

	std::stringstream summary_strm;

	{
		// output all losses for all ids and also the average TODO: minimum lines of input: 1
		vector< float > average( losses.front().valNum(), 0.0 );
		int count = 0;
		for( LossMap::iterator it = losses.begin(); it != losses.end(); ++it ) {
			vector< float > loss = it->second->get();
			summary_strm << it->first;
			for( unsigned int i = 0; i < loss.size(); ++i ) {
				summary_strm << FSEP << loss[i];
				average[i] += loss[i];
			}
			summary_strm << FSEP << it->second->support() << endl;
			++count;
		}
// 		summary_strm << "#======";
// 		for( int i = 0; i < average.size(); ++i ) {
// 			summary_strm << FSEP << "=======";
// 		}
// 		summary_strm << endl;

		if( show_categories ) {
			cout << summary_strm.str();
		}

		// output average
		cout << "NORM";
		for( unsigned int i = 0; i < average.size(); ++i ) {
			cout << FSEP << average[i] / float( count );
		}
		cout << endl;
	}

	{
		vector< float > all = loss_all->get();
		cout << "ALL";
		for( unsigned int i = 0; i < all.size(); ++i ) {
			cout << FSEP << all[i];
		}
		cout << endl;
	}

	delete loss_all;
	delete factory;
	return EXIT_SUCCESS;
}
