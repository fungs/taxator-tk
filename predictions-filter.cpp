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
#include <stack>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>
#include "boost/filesystem.hpp"
#include "src/taxontree.hh"
#include "src/ncbidata.hh"
#include "src/constants.hh"
#include "src/predictionrecordbinning.hh"
#include "src/taxonomyinterface.hh"
#include "src/predictionranges.hh"
#include "src/fastnodemap.hh"



using namespace std;



int main ( int argc, char** argv ) {

	vector< string > ranks, files;
	bool delete_unmarked;
	large_unsigned_int min_support_in_sample;
	float signal_majority_per_sequence;
	medium_unsigned_int min_support_per_sequence;
	boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::size_type num_queries_preallocation;

	namespace po = boost::program_options;
	po::options_description desc ( "Allowed options" );
	desc.add_options()
	( "help,h", "show help message" )
	( "ranks,r", po::value< vector< string > >( &ranks )->multitoken(), "set node ranks at which to do predictions" )
	( "delete-notranks,d", po::value< bool >( &delete_unmarked )->default_value( false ), "delete all nodes that don't have any of the given ranks (make sure that input taxons are at those ranks)" )
	( "files,f", po::value< vector< string > >( &files )->multitoken(), "arbitrary number of prediction files" )
	( "sample-min-support,m", po::value< large_unsigned_int >( &min_support_in_sample )->default_value( 0 ), "minimum support (in positions) for taxon for noise pruning" )
	( "sequence-min-support,s", po::value< medium_unsigned_int >( &min_support_per_sequence )->default_value( 50 ), "minimum number of positions supporting a taxonomic signal on any sequence, if not met a fall-back on a more robust algorthm will be used" )
	( "signal-majority,j", po::value< float >( &signal_majority_per_sequence )->default_value( .6 ), "majority fraction of support of any single sequence required for combination of different signals (> 0.5 to be stable)" )
	( "preallocate-num-queries", po::value< boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::size_type >( & num_queries_preallocation )->default_value( 5000 ), "advanced parameter for better memory allocation, set to number of query sequences or similar (no need to be set)" );

	po::variables_map vm;
	po::store ( po::command_line_parser ( argc, argv ).options ( desc ).run(), vm );
	po::notify ( vm );

	if ( vm.count ( "help" ) ) {
		cout << desc << endl;
		return EXIT_FAILURE;
	}

	if ( ! vm.count ( "ranks" ) ) ranks = default_ranks;
	
	set< string > additional_files;
	
	// create taxonomy
	boost::scoped_ptr< Taxonomy > tax( loadTaxonomyFromEnvironment( &ranks ) );
	if( ! tax ) return EXIT_FAILURE;
	if( delete_unmarked ) tax->deleteUnmarkedNodes();
	TaxonomyInterface taxinter ( tax.get() );
	
	// setup parser for primary input file (that determines the output order)
	boost::scoped_ptr< PredictionFileParser< PredictionRecordBinning > > parse;
	if ( files.empty() ) {
		parse.reset( new PredictionFileParser< PredictionRecordBinning > ( std::cin, tax.get() ) );
	} else {
		vector< string >::iterator file_it = files.begin();
		while( file_it != files.end() ) {
			if( *file_it == "-" ) {
				parse.reset( new PredictionFileParser< PredictionRecordBinning > ( std::cin, tax.get() ) );
				++file_it;
				break;
			} else {
				if( boost::filesystem::exists( *file_it ) ) {
					parse.reset( new PredictionFileParser< PredictionRecordBinning > ( *file_it, tax.get() ) );
					break;
				} else {
					cerr << "Could not read file \"" << *file_it++ << "\"" << endl;
				}
			}
		}
		
		if ( ! parse ) {
			cerr << "There was no valid input file" << endl;
			return EXIT_FAILURE;
		}
		
		// define additional input files
		if ( file_it != files.end() ) {
			const std::string& primary_file = *file_it;
			do {
				if( boost::filesystem::exists( *file_it ) ) {
					additional_files.insert( *file_it );
				} else {
					cerr << "Could not read file \"" << *file_it << "\"" << endl;
				}
			} while ( ++file_it != files.end() );
			additional_files.erase( primary_file );
		}
	}

	// parse primary input
	// default output order corresponds to the first input file with additional records appended at the end
	boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > > predictions_per_query; //future owner of all dynamically allocated objects
 	predictions_per_query.reserve( num_queries_preallocation ); //avoid early re-allocation TODO: command line parameter
	
	{
		if ( additional_files.empty() ) { //parse only primary file (no lookup is done here! this means predictions for same sequences must be consecutive)
			const std::string* prev_name = &empty_string;
			boost::ptr_list< PredictionRecordBinning >* last_added_rec_list = NULL;
			for ( PredictionRecordBinning* rec = parse->next(); rec; rec = parse->next() ) {
				if ( rec->getQueryIdentifier() != *prev_name ) {
					prev_name = &rec->getQueryIdentifier();
// 					std::cerr << "new query: " << rec->getQueryIdentifier() << std::endl;
// 					std::cerr << "entry output is: " << *rec;
					last_added_rec_list = new boost::ptr_list< PredictionRecordBinning >();
					predictions_per_query.push_back( last_added_rec_list );
				}
				last_added_rec_list->push_back( rec ); //will take ownership of the record
			}
		} else { //parse additional
			
			std::map< string, boost::ptr_list< PredictionRecordBinning >& > records_by_queryid; //TODO: use save mem trick
			
			{ //parse primary in case of multiple files (with lookup!)
				const std::string* prev_name = &empty_string;
				boost::ptr_list< PredictionRecordBinning >* last_added_rec_list = NULL;
				for ( PredictionRecordBinning* rec = parse->next(); rec; rec = parse->next() ) {
					if ( rec->getQueryIdentifier() == *prev_name ) last_added_rec_list->push_back( rec ); //will take ownership of the record_container
					else {
						prev_name = &rec->getQueryIdentifier();
						std::map< string, boost::ptr_list< PredictionRecordBinning >& >::iterator find_it = records_by_queryid.find( rec->getQueryIdentifier() );
						if ( find_it != records_by_queryid.end() ) {
							find_it->second.push_back( rec ); //transfer ownership
						} else {
							last_added_rec_list = new boost::ptr_list< PredictionRecordBinning >();
							predictions_per_query.push_back( last_added_rec_list ); //transfer ownership
							records_by_queryid[ rec->getQueryIdentifier() ] = *last_added_rec_list;
							last_added_rec_list->push_back( rec ); //transfer ownership
						}
					}
				}
			}
			
			// parse additional files
			for (set< string >::const_iterator file_it = additional_files.begin(); file_it != additional_files.end(); ++file_it ) {
				PredictionFileParser< PredictionRecordBinning > parse( *file_it, tax.get() );
				boost::ptr_list< PredictionRecordBinning >* last_added_rec_list = NULL;
				for ( PredictionRecordBinning* rec = parse.next(); rec; rec = parse.next() ) {
					std::map< string, boost::ptr_list< PredictionRecordBinning >& >::iterator find_it = records_by_queryid.find( rec->getQueryIdentifier() );
					if ( find_it == records_by_queryid.end() ) {
						last_added_rec_list = new boost::ptr_list< PredictionRecordBinning >();
						predictions_per_query.push_back( last_added_rec_list ); //transfer ownership
						records_by_queryid[ rec->getQueryIdentifier() ] = *last_added_rec_list;
						last_added_rec_list->push_back( rec ); //transfer ownership
					} else {
						find_it->second.push_back( rec ); //transfer ownership
					}
				}
			}
		}
	}
	
	// range pruning
	// in this step the overall sample support for each node is recorded and each
	// range is shrunk such that the remaining nodes have a minimum support (unit is bp)

	//counting support of nodes
	std::cerr << "analyzing sample composition by signal counting...";
	medium_unsigned_int minimum_support_found = std::numeric_limits< medium_unsigned_int >::max();
	const TaxonNode* const root_node = taxinter.getRoot();
	FastNodeMap< large_unsigned_int > support( taxinter.getMaxDepth() );
	large_unsigned_int& root_support = support[ root_node ];
	for ( boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::iterator query_it = predictions_per_query.begin(); query_it != predictions_per_query.end(); ++query_it ) {
		for ( boost::ptr_list< PredictionRecordBinning >::iterator prec_it = query_it->begin(); prec_it != query_it->end(); ++prec_it ) {
			Taxonomy::PathUpIterator pit = taxinter.traverseUp( prec_it->getLowerNode() );
			
			// process lowest node
			medium_unsigned_int total_node_support = prec_it->getSupportAt( &*pit );
			minimum_support_found = std::min( minimum_support_found, total_node_support );
			large_unsigned_int* value_found = support.find( &*pit );
			if ( value_found ) *value_found += total_node_support;
			else support[ &*pit ] = total_node_support;
			
			// process rest
			if ( pit != root_node ) {
				while ( ++pit != root_node ) {
					total_node_support = std::max( total_node_support, prec_it->getSupportAt( &*pit ) );
					large_unsigned_int* value_found = support.find( &*pit );
					if ( value_found ) *value_found += total_node_support;
					else support[ &*pit ] = total_node_support;
	// 				std::cerr << "after adding, node " << pit->data->annotation->name << " it has support " << support[ &*pit ] << std::endl;
				}
				total_node_support = std::max( total_node_support, prec_it->getSupportAt( root_node ) );
				root_support += total_node_support;
			}
		}
	}
	std::cerr << " done: sample contains " << support.size() << " (nested) taxa with total support of " << support[ root_node ] << " bp" << std::endl;
	
	// shrink ranges from lower end
	std::cerr << "noise removal...";
	std::set< const TaxonNode* > pruned_nodes;
	if ( minimum_support_found < min_support_in_sample ) {
		for ( boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::iterator query_it = predictions_per_query.begin(); query_it != predictions_per_query.end(); ++query_it ) {
			for ( boost::ptr_list< PredictionRecordBinning >::iterator prec_it = query_it->begin(); prec_it != query_it->end(); ) {
				const TaxonNode* lower_node = prec_it->getLowerNode();
				const TaxonNode* upper_node = prec_it->getUpperNode();
				
				Taxonomy::PathUpIterator pit = taxinter.traverseUp( lower_node );
				while ( pit != upper_node && support[ &*pit ] < min_support_in_sample ) {
					pruned_nodes.insert( &*pit );	
					++pit;
				}
				
				if ( pit == upper_node && support[ &*pit ] < min_support_in_sample ) { //remove whole range
					pruned_nodes.insert( &*pit );
					prec_it = query_it->erase( prec_it );
					continue;
				}
				
				if ( pit != lower_node ) prec_it->pruneLowerNode( &*pit ); //prune
				++prec_it;
			}
		}
	}
	std::cerr << " done: " << pruned_nodes.size() << " taxa were removed" << std::endl;
	
	// binning
	// in this step multiple ranges are combined into a single range by combining
	// evidence for sub-ranges. This algorithm should consider support, signal
	// strength and possibly the interpolation value to come up with a good
	// prediction range
	std::cerr << "binning step... ";
	std::ofstream binning_debug_output( "binning.log" );
	for ( boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::iterator it = predictions_per_query.begin(); it != predictions_per_query.end(); ++it ) {
		if ( it->size() == 1 ) { // pass-through
			PredictionRecordBinning& prec = it->front();
			prec.setQueryFeatureBegin( 1 );
			prec.setQueryFeatureEnd( prec.getQueryLength() );
			prec.setBinningType( PredictionRecordBinning::single );
			cout << prec;
		}
		
		else if ( it->size() > 1 ) { //run combination algo
			boost::scoped_ptr< PredictionRecordBinning > prec( combinePredictionRanges( *it, tax.get(), signal_majority_per_sequence, min_support_per_sequence, binning_debug_output ) );
			std::cout << *prec;
		}
	}
	std::cerr << " done" << std::endl;

	return EXIT_SUCCESS;
}
