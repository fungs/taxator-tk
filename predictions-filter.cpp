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
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>
#include "boost/filesystem.hpp"
#include "src/taxontree.hh"
#include "src/ncbidata.hh"
#include "src/constants.hh"
#include "src/predictionrecord.hh"
#include "src/taxonomyinterface.hh"



using namespace std;



int main ( int argc, char** argv ) {

	vector< string > ranks, files;
	string cmb_mode;
	float minval;
	bool delete_unmarked;

	namespace po = boost::program_options;
	po::options_description desc ( "Allowed options" );
	desc.add_options()
	( "help,h", "show help message" )
	( "ranks,r", po::value< vector< string > >( &ranks )->multitoken(), "set node ranks at which to do predictions" )
	( "delete-notranks,d", po::value< bool >( &delete_unmarked )->default_value( false ), "delete all nodes that don't have any of the given ranks (make sure that input taxons are at those ranks)" )
	( "files,f", po::value< vector< string > >( &files )->multitoken(), "arbitrary number of prediction files" )
	( "combination-mode,c", po::value< string >( &cmb_mode ), "define way to combine the files" )
	( "group-min,m", po::value< float >( &minval )->default_value( 1.0 ), "minimum fraction (<1) of sample or absolute count (>=1) for group" );

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
	PredictionFileParser* parse = NULL;
	if ( files.empty() ) {
		parse = new PredictionFileParser ( std::cin, tax.get() );
	} else {
		vector< string >::iterator file_it = files.begin();
		while( file_it != files.end() ) {
			if( *file_it == "-" ) {
				parse = new PredictionFileParser ( std::cin, tax.get() );
				++file_it;
				break;
			} else {
				if( boost::filesystem::exists( *file_it ) ) {
					parse = new PredictionFileParser ( *file_it, tax.get() );
					break;
				} else {
					cerr << "Could not read file \"" << *file_it++ << "\"" << endl;
				}
			}
		}
		
		if( ! parse ) {
			cerr << "There was no valid input file" << endl;
			return EXIT_FAILURE;
		}
		
		if ( file_it != files.end() ) { //process all others
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
	boost::ptr_vector< PredictionRecord > record_container ( 1000 ); //assume there are at least 1000 predictions in the file
	
	{
		if ( additional_files.empty() ) { //parse only primary file
			for ( PredictionRecord* rec = parse->next(); rec; rec = parse->next() ) {
				record_container.push_back ( rec ); //it will take ownership of the records
			}
		} else { //parse additional
			const TaxonNode* (TaxonomyInterface::*resolveNodePair)( const TaxonNode* A, const TaxonNode* B ) const = NULL;
			if ( cmb_mode == "lca" ) {
				resolveNodePair = &TaxonomyInterface::getLCA;
			} else {
				if ( cmb_mode == "lcc" ) {
					resolveNodePair = &TaxonomyInterface::getLCC;
				} else {
					cerr << "Parameter for --combination-mode, -c must be either \"lca\" or \"lcc\"" << endl;
					delete parse;
					return EXIT_FAILURE;
				}
			}
			
			map< string, PredictionRecord* > records_by_queryid;
			for ( PredictionRecord* rec = parse->next(); rec; rec = parse->next() ) {
				record_container.push_back ( rec ); //it will take ownership of the records
				records_by_queryid[ rec->query_identifier ] = rec; //overwritten if it exists twice in the file
			}
			
			// parse additional files
			for (set< string >::const_iterator file_it = additional_files.begin(); file_it != additional_files.end(); ++file_it ) {
				PredictionFileParser parse( *file_it, tax.get() );
				for ( PredictionRecord* rec = parse.next(); rec; rec = parse.next() ) {
					map< string, PredictionRecord* >::const_iterator recfind_it = records_by_queryid.find( rec->query_identifier );
					if ( recfind_it == records_by_queryid.end() ) {
						record_container.push_back( rec );
						records_by_queryid[ rec->query_identifier ] = rec;
					} else {
						PredictionRecord* primary_record = recfind_it->second;
						primary_record->prediction_node = (taxinter.*resolveNodePair)( primary_record->prediction_node, rec->prediction_node );
						parse.destroyRecord( rec );
					}
				}
			}
			delete parse;
		}
	}
	
	// here we are left with only one prediction per query sequence and perform the clean-up step
	// easy implementation using hash maps
	map< const TaxonNode*, unsigned int > leaves;
	map< const TaxonNode*, unsigned int > internals;

	const TaxonNode* const root = taxinter.getRoot();
	internals[ root ] = record_container.size(); //implicit

	// first pass: count
	for ( boost::ptr_vector< PredictionRecord >::const_iterator it = record_container.begin(); it != record_container.end(); ++it ) {

		Taxonomy::PathUpIterator path_it = taxinter.traverseUp ( it->prediction_node );

		if ( taxinter.isLeaf( &*path_it ) ) {
			map< const TaxonNode*, unsigned int >::iterator find_it = leaves.find ( &*path_it );

			if ( find_it != leaves.end() ) {
				find_it->second += 1;

			} else {
				leaves[ &*path_it ] = 1;
			}

			path_it++;
		}

		for ( ; path_it != root; ++path_it ) {
			map< const TaxonNode*, unsigned int >::iterator find_it = internals.find ( &*path_it );

			if ( find_it != internals.end() ) {
				find_it->second += 1;

			} else {
				internals[ &*path_it ] = 1;
			}
		}
	}

	unsigned int count_critical;
	if( minval ) { //set critical count value according to command line parameter
		if( minval < 1.0 ) {
			count_critical = minval * (float) record_container.size();
		} else {
			count_critical = (unsigned int) minval;
		}
	}
	
	// second pass: consolidation
	for ( boost::ptr_vector< PredictionRecord >::iterator it = record_container.begin(); it != record_container.end(); ++it ) {
		Taxonomy::PathUpIterator path_it = taxinter.traverseUp ( it->prediction_node );
		unsigned int count;

		if( path_it != root ) {
			if ( taxinter.isLeaf ( &*path_it ) ) {
				count = leaves[ &*path_it ];
			} else {
				count = internals[ &*path_it ];
			}
			
			if( count < count_critical ) {
				do {
					count = internals[ &*++path_it ];
				} while ( count < count_critical && path_it != root );
			}
		}

		it->prediction_node = &*path_it;

		// print result
		std::cout << *it;
	}

	return EXIT_SUCCESS;
}
