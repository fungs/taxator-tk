#include <iostream>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include "src/alignmentrecord.hh"
#include "src/alignmentsfilter.hh"



using namespace std;



int main( int argc, char** argv ) {

	float minscore, toppercent, minpid;
	double maxevalue;
	unsigned int numbestbitscore, minsupport;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		( "help,h", "show help message")
		( "min-score,m", po::value< float >( &minscore )->default_value( 0.0 ), "set min-score filter value" )
		( "min-pid,p", po::value< float >( &minpid )->default_value( 0.0 ), "set minimal PID to consider" )
		( "top-percent,t", po::value< float >( &toppercent )->default_value( 1.0 ), "set top-percent filter value" )
		( "max-evalue,e", po::value< double >( &maxevalue )->default_value( -1.0 ), "set maximum evalue for filtering" )
		( "best-alignments,b", po::value< unsigned int >( &numbestbitscore )->default_value( 0 ), "set number of top bitscore alignments to consider (after toppercent filter)" )
		( "sort-bitscore,s", "sort alignments by decreasing bitscore" )
		( "keep-best-per-ref,k", "for each combination of query and reference sequence id all but the best scoring alignment are removed" )
		( "min-support,c", po::value< unsigned int >( &minsupport )->default_value( 1 ), "set minimum number of hits an alignment needs to have (after filtering)" )
		( "mask-by-star,z", "instead of suppressing filtered alignments mask them by prefixing a star at the line start" );

	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( desc ).run(), vm);
	po::notify(vm);

	if( vm.count( "help" ) ) {
		cout << desc << endl;
		return EXIT_SUCCESS;
	}

	bool sort_by_bitscore = vm.count( "sort-bitscore" );
	bool keep_best_per_gi = vm.count( "keep-best-per-ref" );
	bool mask_by_star = vm.count( "mask-by-star" );

	AlignmentFileParser parser( cin, NULL, true );
	RecordSetGenerator recgen( parser );

	typedef list< AlignmentRecord* > RecordSetType;
	RecordSetType recordset;
	RecordSetType::iterator rec_it;
	AlignmentRecord* record;

	boost::ptr_list< AlignmentRecordSetFilter< RecordSetType > > filters; //takes care of object destruction by itself

	// put filters in queue
	if( keep_best_per_gi ) {
		filters.push_back( new BestScorePerReferenceSeqIDFilter< RecordSetType >() );
	}

	if( sort_by_bitscore ) {
		filters.push_back( new SortByBitscoreFilter< RecordSetType >() );
	}

	if( minpid > 0.0 ) {
		filters.push_back( new MinPIDFilter< RecordSetType >( minpid ) );
	}

	if( maxevalue > 0 ) {
		filters.push_back( new MaxEvalueMinScoreTopPercentFilter< RecordSetType >( minscore, toppercent, maxevalue ) );
	} else {
		if( minscore || toppercent != 1.0 ) {
			filters.push_back( new MinScoreTopPercentFilter< RecordSetType >( minscore, toppercent ) );
		}
	}

	if( numbestbitscore ) {
		filters.push_back( new NumBestBitscoreFilter< RecordSetType >( numbestbitscore ) );
	}

	if( minsupport ) {
	  filters.push_back( new MinSupportFilter< RecordSetType >( minsupport ) );
	}

	boost::ptr_list< AlignmentRecordSetFilter< RecordSetType > >::iterator filter_it;

	//apply list of filters to each record set

	while( recgen.notEmpty() ) {
		recgen.getNext( recordset );

		// apply all filters
		filter_it = filters.begin();
		while( filter_it != filters.end() ) {
// 			cerr << "applying filter: " << filter_it->getInfo() << endl;
			filter_it++->filter( recordset );
		}

		rec_it = recordset.begin();
		while( rec_it != recordset.end() ) {
			record = *rec_it;
			if( ! record->mask ) {
				cout << *record->raw_line << endl;
			} else {
				if( mask_by_star ) {
					cout << '*' << *record->raw_line << endl;
				}
			}
			++rec_it;
			delete record; //clear memory again
		}
		recordset.clear();
	}

	// delete filters (boost pointer list magic)
	filters.clear();

	return EXIT_SUCCESS;
}
