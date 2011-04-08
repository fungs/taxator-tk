#include <iostream>
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include "src/alignmentrecord.hh"
#include "src/alignmentsfilter.hh"
#include "src/accessconv.hh"
#include "src/utils.hh"
#include "src/constants.hh"



using namespace std;



int main( int argc, char** argv ) {

	string ident_level, accessconverter_filename;
	vector< string > ranks;
	bool delete_unmarked;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		( "help,h", "show help message")
		( "filter-ident-level,i", po::value< string >( &ident_level ), "remove alignments that are equivalent with the query label on the given level (either 'seqid' or 'taxid')")
		( "seqid-conv-file,g", po::value< string >( &accessconverter_filename ), "filename of seqid->taxid mappings" )
		( "filter-unclassified,u", "remove queries that were sampled from unclassified organisms (determined by sequence identifier in header)" )
		( "tag-essential,t", "tag alignment to be essential or in excess for the LCA classification algorithm" )
		( "delete-notranks,d", po::value< bool >( &delete_unmarked )->default_value( true ), "delete all nodes that don't have any of the given ranks" )
		( "ranks,r", po::value< vector< string > >( &ranks )->multitoken(), "set node ranks at which to do predictions (if not set it defaults to major NCBI ranks)" )
		( "mask-by-star,z", "instead of suppressing filtered alignments mask them by prefixing a star at the line start" );

	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( desc ).run(), vm);
	po::notify(vm);

	if( vm.count( "help" ) ) {
		cout << desc << endl;
		return EXIT_SUCCESS;
	}

	//parameter consistancy

	if( ! vm.count( "ranks" ) ) { //set to fallback if not given
		ranks = default_ranks;
	}

	bool tag_essential = vm.count( "tag-essential" );
	bool filter_ident = vm.count( "filter-ident-level" );
	bool mask_by_star = vm.count( "mask-by-star" );
	bool filter_unclassified = vm.count( "filter-unclassified" );

	if( filter_unclassified && ! vm.count( "seqid-conv-file" ) ) {
		cerr << "Filtering for unclassified queries requires seqid->taxid mapping" << endl;
		return EXIT_FAILURE;
	}

	if( filter_ident && ident_level != "taxid" && ident_level != "seqid" ) {
		cerr << "Level for filtering can be either 'seqid' or 'taxid'" << endl;
		return EXIT_FAILURE;
	}

	if( ( tag_essential || ( filter_ident && ident_level == "taxid" ) ) && ! vm.count( "seqid-conv-file" ) ) {
		cerr << "For taxid filtering or essential tagging you must give a seqid->taxid mapping" << endl;
		return EXIT_FAILURE;
	}

	StrIDConverter* accessconverter = NULL;
	AlignmentFileParser* parser = NULL;
	boost::ptr_list< AlignmentRecordSetFilter< list< AlignmentRecord* > > > filters;

	if( ident_level == "taxid" || tag_essential || filter_unclassified ) {
		accessconverter = loadStrIDConverterFromFile( accessconverter_filename, 1000 );
		parser = new AlignmentFileParser( cin, accessconverter, true );
	}

	Taxonomy* tax = NULL;
	if( filter_unclassified ) {
		tax = loadTaxonomyFromEnvironment( &ranks );
		if( !tax ) {
			return EXIT_FAILURE;
		}
		filters.push_back( new RemoveUnclassifiedQueriesFilter< list< AlignmentRecord*> >(*accessconverter, tax ) );
	}

	if( filter_ident ) {
		if( ident_level == "taxid" ) {
// 			cerr << "adding RemoveIdentTaxIDFilter" << endl;
			filters.push_back( new RemoveIdentTaxIDFilter< list< AlignmentRecord* > >( *accessconverter ) );
		} else { //this means filter on seqid level only
			filters.push_back( new RemoveIdentSeqIDFilter< list< AlignmentRecord* > >() );
		}
	}

	if( tag_essential ) {
		Taxonomy* tax = loadTaxonomyFromEnvironment( &ranks );
		if( ! tax ) {
			filters.clear();
			if( accessconverter ) {
				delete accessconverter;
			}
			if( parser ) {
				delete parser;
			}
			return EXIT_FAILURE;
		}
		if( delete_unmarked ) {
			tax->deleteUnmarkedNodes(); //should speed up the whole process of LCA finding
		}
// 		cerr << "adding TagEssentialFilter" << endl;
		filters.push_back( new TagEssentialFilter< list< AlignmentRecord* > >( *accessconverter, tax ) );
	}

	if( ! parser ) { //in case not set, yet
		parser = new AlignmentFileParser( cin, NULL, true );
	}

	//start applying the filters

	RecordSetGenerator recgen( *parser );
	list< AlignmentRecord* > recordset;
	list< AlignmentRecord* >::iterator rec_it;
	boost::ptr_list< AlignmentRecordSetFilter< list< AlignmentRecord* > > >::iterator filter_it;
	AlignmentRecord* record;

	while( recgen.notEmpty() ) {
		recgen.getNext( recordset );

		// apply all filters
		filter_it = filters.begin();
		while( filter_it != filters.end() ) {
			filter_it++->filter( recordset );
		}

		// output lines not filtered
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


	//tidy up
	if( accessconverter ) {
		delete accessconverter;
	}
	delete parser;

	return EXIT_SUCCESS;
}
