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

#ifndef predictionranges_hh_
#define predictionranges_hh_

#include <boost/tuple/tuple.hpp>
#include "predictionrecordbinning.hh"
#include <iomanip>



namespace details {
	typedef boost::tuple< Taxonomy::CPathDownIterator, std::vector<medium_unsigned_int>, std::vector<medium_unsigned_int>, bool > TupleRangeCombine;
	typedef std::list< TupleRangeCombine > TupleRangeCombineList;

	template < int i >
	bool all( const TupleRangeCombineList& tlist ) {
		for ( TupleRangeCombineList::const_iterator it = tlist.begin(); it != tlist.end(); ++it ) if ( ! it->get<i>() ) return false;
		return true;
	}

	template < int i >
	bool any( const TupleRangeCombineList& tlist ) {
		for ( TupleRangeCombineList::const_iterator it = tlist.begin(); it != tlist.end(); ++it ) if ( it->get<i>() ) return true;
		return false;
	}

	template < int i >
	bool removeIf( TupleRangeCombineList& tlist ) {
		bool tmp = false;
		for ( TupleRangeCombineList::iterator it = tlist.begin(); it != tlist.end(); ) {
			if ( it->get<i>() ) {
				it = tlist.erase( it );
				tmp = true;
			}	else ++it;
		}
		return tmp;
	}

	void getSupport( const TupleRangeCombineList& tlist, medium_unsigned_int& direct_support, medium_unsigned_int& total_support ) {
		medium_unsigned_int direct = 0;
		medium_unsigned_int indirect = 0;
		for ( TupleRangeCombineList::const_iterator it = tlist.begin(); it != tlist.end(); ++it ) {
			direct += it->get<1>()[ it->get<0>()->data->root_pathlength ];
			indirect += it->get<2>()[ it->get<0>()->data->root_pathlength ];
		}
		direct_support = direct;
		total_support = indirect;
// 		std::cerr << "return support values: (" << static_cast<int>( direct_support ) << "," << static_cast<int>( total_support ) << ")" << std::endl;
	}

	void setPathEndState( TupleRangeCombineList& tlist ) {
		for ( TupleRangeCombineList::iterator it = tlist.begin(); it != tlist.end(); ++it ) it->get<3>() = it->get<0>()->data->root_pathlength == (it->get<1>().size() - 1);
	}

	void stepDown( TupleRangeCombineList& tlist ) {
		for ( TupleRangeCombineList::iterator it = tlist.begin(); it != tlist.end(); ++it ) it->get<0>()++;
	}

	bool reduceToMajority( TupleRangeCombineList& tlist ) { //2-pass function operating on small lists
// 		std::cerr << "reduceToMajority: list contains " << tlist.size() << " elements before" << std::endl;

		if ( tlist.size() < 2 ) return false;

		std::map< const TaxonNode*, float > supports;
		const TaxonNode* max_node = NULL;
		float max_support = .0;

		// pass 1: count and remember highest score
		for ( TupleRangeCombineList::iterator it = tlist.begin(); it != tlist.end(); ++it ) {
			const TaxonNode* node = &*it->get<0>();
			std::map< const TaxonNode*, float >::iterator find_it = supports.find( node );
			if ( find_it != supports.end() ) find_it->second += it->get<2>()[ it->get<0>()->data->root_pathlength ]; //count total support of node
			else find_it = supports.insert( std::make_pair( node, it->get<2>()[ it->get<0>()->data->root_pathlength ] ) ).first;

// 			std::cerr << "node partial total support: " << find_it->second << std::endl;

			if ( find_it->second > max_support ) { //TODO: if tie, stop?
				max_support = find_it->second;
				max_node = node;
			}
		}

		if ( supports.size() == 1 ) return false;

		// pass 2: keep majority, remove rest
		for ( TupleRangeCombineList::iterator it = tlist.begin(); it != tlist.end(); ) {
			if ( it->get<0>() != max_node ) it = tlist.erase( it );
			else ++it;
		}

// 		std::cerr << "reduceToMajority: list contains " << tlist.size() << " elements after" << std::endl;
		return true;
	}
}



template < class ContainerT >
PredictionRecordBinning* combinePredictionRanges( const ContainerT& predictions, const std::string& identifier, const Taxonomy* tax, float min_signal_percentage, medium_unsigned_int min_support, std::ostream& debug_output = std::cerr ) {

	assert( predictions.size() > 1 ); //TODO: debug mode

	using namespace details;
	TupleRangeCombineList tlist;
	TaxonomyInterface taxinter( tax );

	// initialize temporary data structure (list of tuples)
	medium_unsigned_int summed_support = 0;
	large_unsigned_int summed_length = 0;
	{
		std::set<std::string> processed_identifiers;
		for ( typename ContainerT::const_iterator it = predictions.begin(); it != predictions.end(); ++it ) {
			const TaxonNode* lower_node = it->getLowerNode();
			int i = lower_node->data->root_pathlength;
			medium_unsigned_int support = it->getSupportAt( i );
			summed_support += support;
			if(processed_identifiers.emplace(it->getQueryIdentifier()).second) {
				summed_length += it->getQueryLength();
			}

			const std::size_t depth = i + 1;
			tlist.push_back( boost::make_tuple( taxinter.traverseDown< Taxonomy::CPathDownIterator >( lower_node ), std::vector<medium_unsigned_int>( depth ), std::vector<medium_unsigned_int>( depth ), false ) );

			std::vector< medium_unsigned_int >& direct_support = tlist.back().get<1>(); //TODO: initialize direct_support with taxon_support
			std::vector< medium_unsigned_int >& total_support = tlist.back().get<2>();
			total_support[i] = direct_support[i] = support;

	// 		debug_output << *it;
	// 		debug_output << " summed support at level " << i << " is " << direct_support[i] << " (" << total_support[i] << ")" << std::endl;
			while ( --i >= 0 ) {
				direct_support[i] = it->getSupportAt( i );
				total_support[i] = std::max( total_support[i+1], direct_support[i] );
	// 			debug_output << "summed support at level " << i << " is " << direct_support[i] << " (" << total_support[i] << ")" << std::endl;
			}
		}
	}

	PredictionRecordBinning* prec = new PredictionRecordBinning( tax );  // new target record
	{ // initialize combined record
		prec->setQueryIdentifier( identifier );
		prec->setQueryLength( summed_length );
		prec->setQueryFeatureBegin( 1 ); //TODO: range select
		prec->setQueryFeatureEnd( summed_length ); //TODO: range select
	}

	medium_unsigned_int direct_support_thresh = std::max( static_cast< medium_unsigned_int >( min_signal_percentage*summed_support ), min_support );
	medium_unsigned_int direct_support, total_support;
	std::vector< boost::tuple<const TaxonNode*, medium_unsigned_int, medium_unsigned_int, bool> > path;
	const TaxonNode* root_node = taxinter.getRoot();

	//debug
	debug_output << std::endl << "combining " << predictions.size() << " independent predictions for query " << identifier << ", threshold " << direct_support_thresh << " (" << static_cast<uint>( min_signal_percentage*100 ) << " %)" << std::endl;

	for ( typename ContainerT::const_iterator it = predictions.begin(); it != predictions.end(); ++it ) {

		debug_output << static_cast< int >( it->getSupportAt( it->getLowerNode() ) ) << ": ";
		for ( Taxonomy::CPathDownIterator pit( root_node, it->getUpperNode() ); pit != it->getUpperNode(); ++pit ) debug_output << pit->data->annotation->name << ";";
		debug_output << "[";
		for ( Taxonomy::CPathDownIterator pit( it->getUpperNode(), it->getLowerNode() ); pit != it->getLowerNode(); ++pit ) debug_output << pit->data->annotation->name << ";";
		debug_output << it->getLowerNode()->data->annotation->name << "]" << std::endl;
	}

	// set values for root node (TODO: put into initialization loop)
	setPathEndState( tlist );
	getSupport( tlist, direct_support, total_support );

	// walk down each path
	int lower_direct_node_index = -1;
	int running_index = 0;
	while ( ! tlist.empty() ) {
		const TaxonNode* node = &*tlist.front().get<0>();
		if ( direct_support >= direct_support_thresh ) lower_direct_node_index = running_index;
		path.push_back( boost::make_tuple( node, direct_support, total_support, false ) );
		removeIf<3>( tlist ); //remove paths that have ended
		stepDown( tlist ); //forward PathIterators
		++running_index;
		setPathEndState( tlist ); //update bool flags for range end
		path.back().get<3>() = reduceToMajority( tlist ); //set parent branching flag
		getSupport( tlist, direct_support, total_support );
	}

	debug_output << std::endl;

	//debug: output whole path for backtracking
	debug_output << std::setprecision( 3 );
	debug_output << "  L |  direct s. |    total s.| B | name" << std::endl;
	debug_output << "--------------------------------------------" << std::endl;
	for ( std::vector< boost::tuple<const TaxonNode*, medium_unsigned_int, medium_unsigned_int, bool> >::const_iterator it = path.begin(); it != path.end(); ++it ) {
		debug_output << std::fixed << std::setw( 3 ) << static_cast<int>( it->get<0>()->data->root_pathlength ) << " | " << std::setw( 10 ) << static_cast<int>( it->get<1>() ) << " | " << std::setw( 10 ) << static_cast<int>( it->get<2>() ) << " | " << it->get<3>() << " | ";
		if ( it->get<1>() >= direct_support_thresh ) debug_output << "*";
		debug_output	<< it->get<0>()->data->annotation->name << std::endl;
	}

	if ( lower_direct_node_index >= 0 ) { //direct mode
		debug_output << "using direct binning mode..." << std::endl;
		prec->setBinningType( PredictionRecordBinning::direct );

		const TaxonNode* const lower_node = path[ lower_direct_node_index ].get<0>();
		const medium_unsigned_int lower_node_support = path[ lower_direct_node_index ].get<2>(); //return total support (like fallback mode)


		{
			medium_unsigned_int upper_node_support = lower_node_support;
			const TaxonNode* upper_node = lower_node;
			int upper_direct_node_index = lower_direct_node_index;

			for ( int j = lower_direct_node_index; j >= 0; --j ) {
				if ( path[j].get<1>() >= direct_support_thresh ) {
					upper_node_support = path[j].get<2>(); //return total support (like fallback mode)
					upper_node = path[j].get<0>();
					upper_direct_node_index = j;
					if ( path[j].get<3>() ) break;
				}
			}
			prec->setNodeRange( lower_node, lower_node_support, upper_node, upper_node_support );
// 			std::cerr << "node range from " << static_cast<int>( lower_node->data->root_pathlength ) << " to " << static_cast<int>( upper_node->data->root_pathlength ) << std::endl;
			for ( int i = lower_direct_node_index; i > upper_direct_node_index; --i ) prec->setSupportAt( path[i].get<0>(), path[i].get<1>() ); //set support values in between
// 			debug_output << "upper node in direct mode: " << upper_node->data->annotation->name << " (" << upper_node->data->taxid << ")" << std::endl;
// 			debug_output << *prec;
		}

// 		for ( uint j = 0; j < path.size(); ++j ) {
// 			medium_unsigned_int upper_support = path[j].get<1>();
// 			if ( upper_support >= direct_support_thresh ) {
// 				prec->setNodeRange( lower_node, lower_node_support, path[j].get<0>(), upper_support ); //initialize range
// 				for ( uint i = lower_direct_node_index; i > j; --i ) prec->setSupportAt( path[i].get<0>(), path[i].get<1>() ); //set support values in between
//
// 				break;
// 			}
// 		}
		return prec;
	}

	// fallback mode: find first branching point above threshold (majority LCA)
	debug_output << "using fallback binning mode..." << std::endl;
	prec->setBinningType( PredictionRecordBinning::fallback );

	for ( int i = path.size() - 1; i >= 0; --i ) {
		if ( path[i].get<2>() >= direct_support_thresh ) {
			prec->setNodePoint( path[i].get<0>(), path[i].get<2>() );
			return prec;
		}
	}

	prec->setNodePoint( path[0].get<0>(), path[0].get<2>() );
	return prec;
}

#endif //predictionranges_hh_
