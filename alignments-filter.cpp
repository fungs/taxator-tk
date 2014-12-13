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
#include <string>
#include <boost/lexical_cast.hpp>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#include <boost/scoped_ptr.hpp>
#include "src/alignmentrecord.hh"
#include "src/alignmentsfilter.hh"



using namespace std;

template< typename AlignmentsFilterListType >
void parseAndFilter( AlignmentsFilterListType& filters, bool mask = true ) {
// 	template <template <typename> class Container, typename AlignmentRecordType>
// 	void parseAndFilter( boost::ptr_list< AlignmentsFilter< Container::< AlignmentRecordType* > > >& filters, bool mask = true ) {

	// some type tricks
	typedef typename boost::remove_pointer< typename AlignmentsFilterListType::value_type >::type AlignmentsFilterType; //expect stdcontainer
	typedef typename AlignmentsFilterType::AlignmentRecordSetType AlignmentRecordSetType;
	typedef typename boost::remove_pointer< typename AlignmentRecordSetType::value_type >::type AlignmentRecordType; //expect stdcontainer

	typename AlignmentsFilterListType::iterator filter_it;
	AlignmentRecordFactory< AlignmentRecordType > recfac;
	AlignmentFileParser< AlignmentRecordType > parser( cin, recfac );
	RecordSetGeneratorSort< AlignmentRecordType, AlignmentRecordSetType, false > recgen( parser ); // Eik geaendert

	AlignmentRecordSetType recordset;
	typename AlignmentRecordSetType::iterator rec_it;
	AlignmentRecordType* record;

	//apply list of filters to each record set
	while( recgen.notEmpty() ) {
	    //std::cerr << "iterate ov recgen\n";
		recgen.getNext( recordset );
		// apply all filters
		filter_it = filters.begin();
		while( filter_it != filters.end() ) {
 			//cerr << "applying filter: " << filter_it->getInfo() << endl;
			filter_it++->filter( recordset );
		}

		rec_it = recordset.begin();
		while( rec_it != recordset.end() ) {
			record = *rec_it;
			if( ( record->isFiltered() && mask ) || ! record->isFiltered() ) {
				cout << *record;
			}
			++rec_it;
			recfac.destroy( record ); //clear memory again
		}
		recordset.clear();
	}
}


int main( int argc, char** argv ) {

	float minscore, toppercent, minpid;
	double maxevalue;
	unsigned int numbestscore, minsupport;

	std::string tax_map1_filename, tax_map2_filename;

	namespace po = boost::program_options;
	po::options_description desc("Allowed options");
	desc.add_options()
		( "help,h", "show help message")
		( "min-score,m", po::value< float >( &minscore )->default_value( 0.0 ), "set min-score filter value" )
		( "min-pid,p", po::value< float >( &minpid )->default_value( 0.0 ), "set minimal PID to consider" )
		( "top-percent,t", po::value< float >( &toppercent )->default_value( 1.0 ), "set top-percent filter value" )
		( "max-evalue,e", po::value< double >( &maxevalue )->default_value( -1.0 ), "set maximum evalue for filtering" )
		( "best-alignments,b", po::value< unsigned int >( &numbestscore )->default_value( 0 ), "set number of top score alignments to consider (after toppercent filter)" )
		( "sort-score,s", "sort alignments by decreasing score" )
		( "keep-best-per-ref,k", "for each combination of query and reference sequence id all but the best scoring alignment are removed" )
		( "min-support,c", po::value< unsigned int >( &minsupport )->default_value( 1 ), "set minimum number of hits an alignment needs to have (after filtering)" )
		( "remove-ref-from-query-taxon,r", "remove alignments for labeled data to test different degrees of taxonomic distance" )
		( "taxon-mapping-sample,x", po::value< std::string >( &tax_map1_filename ), "map sample identifier to taxon" )
		( "taxon-mapping-reference,y", po::value< std::string >( &tax_map2_filename ), "map reference identifier to taxon" )
		( "mask-by-star,z", "instead of suppressing filtered alignments mask them by prefixing a star at the line start" );

	po::variables_map vm;
	po::store(po::command_line_parser( argc, argv ).options( desc ).run(), vm);
	po::notify(vm);

	if( vm.count( "help" ) ) {
		cout << desc << endl;
		return EXIT_SUCCESS;
	}

	bool sort_by_score = vm.count( "sort-score" );
	bool keep_best_per_gi = vm.count( "keep-best-per-ref" );
	bool mask_by_star = vm.count( "mask-by-star" );
	bool remove_same_taxon = vm.count( "remove-ref-from-query-taxon" );

	typedef list< AlignmentRecord* > RecordSetType;
	boost::ptr_list< AlignmentsFilter< RecordSetType > > filters; //takes care of object destruction by itself

	boost::scoped_ptr< StrIDConverter > seqid2taxid_sample;
	boost::scoped_ptr< StrIDConverter > seqid2taxid_reference;

	// put filters in queue
	if( remove_same_taxon ) {
		if( tax_map1_filename.empty() || tax_map2_filename.empty() ) {
			cout << "'--remove-ref-from-query-taxon' requires two mapping files: '--taxon-mapping-sample' and '--taxon-mapping-reference'" << endl;
			return EXIT_SUCCESS;
		}

		seqid2taxid_sample.reset( loadStrIDConverterFromFile( tax_map1_filename ) );
		seqid2taxid_reference.reset( loadStrIDConverterFromFile( tax_map2_filename, 1000 ) );

		filters.push_back( new TaxonMaskingFilter< RecordSetType >( *seqid2taxid_sample, *seqid2taxid_reference ) );
	}
	if( keep_best_per_gi ) {
		filters.push_back( new BestScorePerReferenceSeqIDFilter< RecordSetType >() );
	}
	if( sort_by_score ) {
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
	if( numbestscore ) {
		filters.push_back( new NumBestBitscoreFilter< RecordSetType >( numbestscore ) );
	}
	if( minsupport ) {
	  filters.push_back( new MinSupportFilter< RecordSetType >( minsupport ) );
	}

	parseAndFilter( filters, mask_by_star );

	// delete filters (boost pointer list magic)
	filters.clear();

	return EXIT_SUCCESS;
}
