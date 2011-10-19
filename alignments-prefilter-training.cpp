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

#include "src/alignmentrecord.hh"
#include "src/alignmentsfilter.hh"
#include "src/accessconv.hh"
#include "src/utils.hh"
#include "src/constants.hh"
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/scoped_ptr.hpp>
#include <iostream>
#include <string>

using namespace std;



template< typename ALType >
void doFiltering( AlignmentRecordFactory< ALType >& fac, boost::ptr_list< AlignmentsFilter< list< ALType* > > >& filters, bool mask_lines = false ) {
	AlignmentFileParser< ALType > parser( cin, fac );
	RecordSetGenerator< ALType > recgen( parser );
	list< ALType* > recordset;
	typename list< ALType* >::iterator rec_it;
	typename boost::ptr_list< AlignmentsFilter< list< ALType* > > >::iterator filter_it;
	ALType* record;
	
	if( mask_lines ) {	
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
				cout << **rec_it;
				fac.destroy( *rec_it++ );
			}
			recordset.clear();
		}
	} else { //a little code duplication...
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
				if ( ! record->isFiltered() ) cout << *record;
				++rec_it;
				fac.destroy( record );
			}
			recordset.clear();
		}
	}
}



int main( int argc, char** argv ) {

	string ident_level, accessconverter_filename, regex_identifier;
	vector< string > ranks;
	bool delete_unmarked;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		( "help,h", "show help message")
		( "filter-ident-level,i", po::value< string >( &ident_level )->default_value( "taxid" ), "remove alignments that are equivalent with the query label on the given level (either 'seqid' or 'taxid')")
		( "seqid-conv-file,g", po::value< string >( &accessconverter_filename ), "filename of seqid->taxid mappings" )
		( "regex-identifier,x", po::value< string >( &regex_identifier ), "XXX-style regular expression for extracting the identifier from artificial sequence name" )
		( "filter-unclassified,u", "remove queries that were sampled from unclassified organisms (determined by sequence identifier in header)" )
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
	
	boost::scoped_ptr< Taxonomy > tax;
	boost::scoped_ptr< StrIDConverter > seqid2taxid;
	boost::ptr_list< AlignmentsFilter< list< AlignmentRecordTaxonomy* > > > filters_taxonomy;
	boost::ptr_list< AlignmentsFilter< list< AlignmentRecord* > > > filters_regular;
	
	if ( filter_unclassified ) {
		seqid2taxid.reset( loadStrIDConverterFromFile( accessconverter_filename, 1000 ) ); //TODO: add caching parameter
		tax.reset( loadTaxonomyFromEnvironment( &ranks ) );
		
		if( ! tax ) return EXIT_FAILURE;
		
		if( delete_unmarked ) tax->deleteUnmarkedNodes(); //should speed up the whole process of LCA finding
		
		filters_taxonomy.push_back( new RemoveUnclassifiedQueriesFilter< list< AlignmentRecordTaxonomy* > >( *seqid2taxid, tax.get(), regex_identifier ) );
	}
	
	if( filter_ident ) {
		if ( ident_level == "taxid" ) {
			if( ! seqid2taxid ) {
				seqid2taxid.reset( loadStrIDConverterFromFile( accessconverter_filename, 1000 ) ); //TODO: add caching parameter
			}
			if ( filters_taxonomy.empty() ) {
				filters_regular.push_back( new RemoveIdentTaxIDFilter< list< AlignmentRecord* > >( *seqid2taxid, regex_identifier ) );
			} else {
				filters_taxonomy.push_back( new RemoveIdentTaxIDFilter< list< AlignmentRecordTaxonomy* > >( *seqid2taxid, regex_identifier ) );
			}
		} else { //ident_level == "seqid"
			if ( filters_taxonomy.empty() ) {
				filters_regular.push_back( new RemoveIdentSeqIDFilter< list< AlignmentRecord* > >( regex_identifier ) );
			} else {
				filters_taxonomy.push_back( new RemoveIdentSeqIDFilter< list< AlignmentRecordTaxonomy* > >( regex_identifier ) );
			}
		}
	}

	if( filters_taxonomy.empty() ) {
		AlignmentRecordFactory< AlignmentRecord > fac;
		doFiltering( fac, filters_regular, mask_by_star );
	} else {
		assert( tax ); //check to be sure
		assert( seqid2taxid ); //check to be sure
		AlignmentRecordFactory< AlignmentRecordTaxonomy > fac( *seqid2taxid, tax.get() );
		doFiltering( fac, filters_taxonomy, mask_by_star );
	}
	
	return EXIT_SUCCESS;
}
