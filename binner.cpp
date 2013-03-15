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
	float signal_majority_per_sequence, min_support_in_sample_percentage( 0. );
	string min_support_in_sample_str, log_filename;
	medium_unsigned_int min_support_per_sequence;
	boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::size_type num_queries_preallocation;

	namespace po = boost::program_options;
	po::options_description visible_options ( "Allowed options" );
	visible_options.add_options()
	( "help,h", "show help message" )
	( "sample-min-support,m", po::value< std::string >( &min_support_in_sample_str )->default_value( "0.01" ), "minimum support in positions (>=1) or fraction of total support (<1) for any taxon" )
	( "sequence-min-support,s", po::value< medium_unsigned_int >( &min_support_per_sequence )->default_value( 50 ), "minimum number of positions supporting a taxonomic signal for any single sequence. If not reached, a fall-back on a more robust algorthm will be used" )
	( "signal-majority,j", po::value< float >( &signal_majority_per_sequence )->default_value( .7 ), "minimum combined fraction of support for any single sequence (> 0.5 to be stable)" )
	( "identity-constrain,i", po::value< vector< string > >(), "minimum required identity for this rank (e.g. -i species:0.8 -i genus:0.7)")
	( "ranks,r", po::value< vector< string > >( &ranks )->multitoken(), "set ranks at which to do predictions" )
	( "files,f", po::value< vector< string > >( &files )->multitoken(), "arbitrary number of prediction files (replaces standard input, use \"-\" to specify a combination of both)" )
	( "logfile,l", po::value< std::string >( &log_filename )->default_value( "binning.log" ), "specify name of file for logging (appending lines)" );

	po::options_description hidden_options("Hidden options");
	hidden_options.add_options()
	( "preallocate-num-queries", po::value< boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::size_type >( & num_queries_preallocation )->default_value( 5000 ), "advanced parameter for better memory allocation, set to number of query sequences or similar (no need to be set)" )
	( "delete-notranks,d", po::value< bool >( &delete_unmarked )->default_value( true ), "delete all nodes that don't have any of the given ranks (make sure that input taxons are at those ranks)" );
	
	po::options_description all_options;
	all_options.add( visible_options ).add( hidden_options );
	
	po::variables_map vm;
	po::store ( po::command_line_parser ( argc, argv ).options ( all_options ).run(), vm );
	po::notify ( vm );

	if ( vm.count ( "help" ) ) {
		cout << visible_options << endl;
		return EXIT_FAILURE;
	}

	if ( ! vm.count ( "ranks" ) ) ranks = default_ranks;
	
	// interpret given sample support
	if ( min_support_in_sample_str.find( '.' ) == std::string::npos ) min_support_in_sample = boost::lexical_cast< large_unsigned_int >( min_support_in_sample_str );
	else min_support_in_sample_percentage = boost::lexical_cast< float >( min_support_in_sample_str );
	
	set< string > additional_files;
	
	// create taxonomy
	boost::scoped_ptr< Taxonomy > tax( loadTaxonomyFromEnvironment( &ranks ) );
	if( ! tax ) return EXIT_FAILURE;
	if( ! ranks.empty() && delete_unmarked ) tax->deleteUnmarkedNodes(); //collapse taxonomy to contain only specified ranks
	TaxonomyInterface taxinter ( tax.get() );
	
	map< const string*, float > pid_per_rank;
	if ( vm.count ( "identity-constrain" ) ) {
		vector< string > fields;
		const vector< string >& values = vm["identity-constrain"].as< const vector< string > >();
		for( vector< string >::const_iterator it = values.begin(); it != values.end(); ++it ) {
			tokenizeSingleCharDelim<>( *it, fields, ":", 1 );
			try {
				if ( fields[0].empty() ) {
					cerr << "Could not read identity constrain: rank cannot be empty string, use e.g. \"-i species:0.8\"" << endl;
					return EXIT_FAILURE;
				}
				pid_per_rank[&(tax->getRankInternal( fields[0]))] = boost::lexical_cast< float >( fields[1] );
			} catch ( const boost::bad_lexical_cast& ) {
				cerr << "Could not read identity constrain: \"" << fields[1] << "\" for rank \"" << fields[0] << "\" as float, use e.g. \"-i species:0.8\"" << endl;
				return EXIT_FAILURE;
			}
			fields.clear();
		}
	}
	
	//STEP 0: PARSING INPUT
	
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
 	predictions_per_query.reserve( num_queries_preallocation ); //avoid early re-allocation
	
	{
		if ( additional_files.empty() ) { //parse only primary file (predictions for same sequences must be consecutive!)
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
			
			std::map< string, boost::ptr_list< PredictionRecordBinning >* > records_by_queryid; //TODO: use save mem trick
			
			{ //parse primary in case of multiple files (with lookup!)
				const std::string* prev_name = &empty_string;
				boost::ptr_list< PredictionRecordBinning >* last_added_rec_list = NULL;
				for ( PredictionRecordBinning* rec = parse->next(); rec; rec = parse->next() ) {
					if ( rec->getQueryIdentifier() == *prev_name ) last_added_rec_list->push_back( rec ); //transfer ownership of record_container
					else {
						prev_name = &rec->getQueryIdentifier();
						std::map< string, boost::ptr_list< PredictionRecordBinning >* >::iterator find_it = records_by_queryid.find( rec->getQueryIdentifier() );
						if ( find_it != records_by_queryid.end() ) {
							find_it->second->push_back( rec ); //transfer ownership
						} else {
							last_added_rec_list = new boost::ptr_list< PredictionRecordBinning >();
							predictions_per_query.push_back( last_added_rec_list ); //transfer ownership
							records_by_queryid[ rec->getQueryIdentifier() ] = last_added_rec_list;
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
					std::map< string, boost::ptr_list< PredictionRecordBinning >* >::iterator find_it = records_by_queryid.find( rec->getQueryIdentifier() );
					if ( find_it == records_by_queryid.end() ) {
						last_added_rec_list = new boost::ptr_list< PredictionRecordBinning >();
						predictions_per_query.push_back( last_added_rec_list ); //transfer ownership
						records_by_queryid[ rec->getQueryIdentifier() ] = last_added_rec_list;
						last_added_rec_list->push_back( rec ); //transfer ownership
					} else {
						find_it->second->push_back( rec ); //transfer ownership
					}
				}
			}
		}
	}
	
	
	
	// STEP 1: RANGE PRUNING
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
	std::cerr << " done: " << support.size() << " nested taxa with total support of " << support[ root_node ] << " bp" << std::endl;
	
	// if min_support_in_sample was given as fraction
	if ( min_support_in_sample_percentage ) min_support_in_sample = support[ root_node ]*min_support_in_sample_percentage;
	
	// shrink ranges from lower end if support is smaller than the minimum required or if it does not comply with user-defined PID per rank.
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
					prec_it = query_it->erase( prec_it ); //TODO: mask instead of delete
					continue;
				}
				
				if ( pit != lower_node ) prec_it->pruneLowerNode( &*pit ); //prune
				++prec_it;
			}
		}
	}
	std::cerr << " done: " << pruned_nodes.size() << " taxa were removed" << std::endl;
	
	// STEP 2: BINNING
	// in this step multiple ranges are combined into a single range by combining
	// evidence for sub-ranges. This algorithm considers only support. Signal
	// strength and interpolation values are ignored. This heuristic seems quite
	// robust
	
	std::cerr << "binning step... ";
	std::ofstream binning_debug_output( log_filename.c_str() );
	for ( boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::iterator it = predictions_per_query.begin(); it != predictions_per_query.end(); ++it ) {
		if( it->empty() ) continue;
		boost::scoped_ptr< PredictionRecordBinning > prec_sptr;
		const PredictionRecordBinning* prec;
		if ( it->size() > 1 ) { //run combination algo for sequence segments
			prec_sptr.reset( combinePredictionRanges( *it, tax.get(), signal_majority_per_sequence, min_support_per_sequence, binning_debug_output ) );
			prec = prec_sptr.get();
		} else { // pass-through segment prediction for whole sequence
			prec = &it->front();
// 			prec->setBinningType( PredictionRecordBinning::single );
		}
		// apply user-defined constrain
		if ( prec->getUpperNode() != root_node && ! pid_per_rank.empty() ) {
			const double seqlen = static_cast< double >( prec->getQueryLength() );
			float min_pid = 0.; //enforce consistency when walking down
			map< const string*, float >::const_iterator find_it;
			const TaxonNode* predict_node = root_node;
			const TaxonNode* target_node = prec->getUpperNode();
			Taxonomy::CPathDownIterator pit = taxinter.traverseDown<Taxonomy::CPathDownIterator>( target_node );
			do {
				pit++;
				find_it = pid_per_rank.find( &(pit->data->annotation->rank) );
				if ( find_it != pid_per_rank.end() ) min_pid = max( min_pid, find_it->second );
				if ( prec->getSupportAt( &*pit )/seqlen < min_pid ) break;
				predict_node = &*pit;
			} while ( pit != target_node );
			std::cout << prec->getQueryIdentifier() << tab << predict_node->data->taxid << endline;
		} else {
			std::cout << prec->getQueryIdentifier() << tab << prec->getUpperNode()->data->taxid << endline;
		}
	}
	std::cerr << " done" << std::endl;

	return EXIT_SUCCESS;
}
