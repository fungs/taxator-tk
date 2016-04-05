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

#ifndef taxonfilter_hh_
#define taxonfilter_hh_

#include <string>
#include <fstream>
#include <boost/iterator/iterator_concepts.hpp>
#include <boost/lexical_cast.hpp>
#include "types.hh"
#include "exception.hh"
#include "taxonomyinterface.hh"
#include "constants.hh"



// abstract base class
class TaxonFilter {
	public:
		virtual ~TaxonFilter() {};
		virtual void operator()( std::string& inputstr ) = 0;
		virtual const std::string getInfo() { return description_; };
	private:
		static const std::string description_;
};

const std::string TaxonFilter::description_ = "abstract TaxonFilter";



class NewickTaxonFilter : public TaxonFilter {
	public:
		NewickTaxonFilter( TaxonomyInterface& taxinter, const std::string& outfile, const std::vector< std::string >& rank_names, bool show_names, bool fill_empty_ranks ) :
		maxdepth_(0),
		taxinter_(taxinter),
		outfile_(outfile),
		rank_names_( rank_names.size() ),
		show_names_(show_names),
		fill_empty_ranks_(fill_empty_ranks) {
			newicklists.reserve( taxinter.getMaxDepth() );
			for (small_unsigned_int i = 0; i < rank_names.size(); ++i ) {
				const std::string& rank = taxinter.getRankInternal( rank_names[i] );
				if( rank.empty() ) {
					std::cerr << "Rank '" << rank_names_[i] << "' not found in taxonomy, ignoring..." << std::endl;
					continue;
				}
				ranks_[&rank] = i;
				rank_names_[i] = &rank;
			}
		};
		
		void operator()( std::string& inputstr ) {
			TaxonID taxid = boost::lexical_cast< TaxonID >( inputstr );
      try {
        const TaxonNode* node = taxinter_.getNode( taxid );
        while( node != taxinter_.getRoot() && ! ranks_.count( &(node->data->annotation->rank) ) ) node = node->parent;
        small_unsigned_int& depth = node->data->root_pathlength;
        if ( newicklists.size() <= depth ) {
          newicklists.resize( depth + 1 );
          maxdepth_ = depth;
        }
        newicklists[depth][node]; //standard constructor = std::list<std::string&>>();
      } catch ( TaxonNotFound &e ) {
        if( TaxonID const * taxid = boost::get_error_info<taxid_info>(e) ) std::cerr << "Could not find node with taxid " << *taxid << " in the taxonomy, skipping record..." << std::endl;
        else std::cerr << "Could not find node by its taxid in the taxonomy, skipping record..." << std::endl;
      }
		};
		
		~NewickTaxonFilter() { //writes the actual newick tree
			for (; maxdepth_ > 0; --maxdepth_ ) {
// 				std::cerr << "depth: " << (int)maxdepth_ << std::endl;  // DEBUG
				for ( nodemap::iterator it = newicklists[maxdepth_].begin(); it != newicklists[maxdepth_].end(); ++it ) {
					const TaxonNode* node = it->first;
					nwlist& l = it->second;
					node = newickify( node, l );

					// add code to parent node
					nwlist& pl = newicklists[node->data->root_pathlength][node];
					if ( ! pl.empty() ) pl.push_front( &newick::nodesep );
					pl.splice( pl.begin(), l ); // empty list remains
				}
				newicklists[maxdepth_].clear(); // remove empty lists
			}
			
			// output root as newick tree
			std::ofstream f( outfile_.c_str() );
			assert( newicklists[0].begin()->first == taxinter_.getRoot() );
			nwlist& rootlist = newicklists[0].begin()->second;
			f << newick::nodestart;
			for (nwlist::const_iterator it = rootlist.begin(); it != rootlist.end(); ++it ) f << **it;
			f << newick::nodestop << newick::treestop;
		}

	private:
		typedef std::list<const std::string*> nwlist;
		typedef std::map<const TaxonNode*, nwlist> nodemap;
		
		const TaxonNode* newickify( const TaxonNode* node, nwlist& l ) {
					// put interior node into parenthesis
					if ( ! l.empty() ) {
						l.push_front( &newick::nodestart );
						l.push_back( &newick::nodestop );
					}
					
					// append node name
					if ( show_names_ ) l.push_back( &node->data->annotation->name );
					else l.push_back( &node->data->taxid );  // only if we can ensure TaxID == std::string
					
					// go to next parent
					const TaxonNode* parent = node->parent;
					while( parent != taxinter_.getRoot() && ! ranks_.count( &(parent->data->annotation->rank) ) ) parent = parent->parent;
					
					if ( fill_empty_ranks_ && parent != taxinter_.getRoot() ) {
						small_unsigned_int running_index = ranks_[&(node->data->annotation->rank)];  // exists by definition
						small_unsigned_int parent_index = ranks_[&(parent->data->annotation->rank)];  // exists by definition
						assert( running_index < parent_index );
						while ( parent_index - ++running_index ) {  // insert anonymous node
// 							std::cerr << "inserting dummy node at level " << *rank_names_[running_index] << std::endl;
							l.push_front( &newick::nodestart );
							l.push_back( &newick::nodestop );
						}
					}
					
					return parent;
		}
		
		
		small_unsigned_int maxdepth_;
		const TaxonomyInterface& taxinter_;
		const std::string& outfile_;
		std::vector< const std::string* > rank_names_;
		std::map< const std::string*, small_unsigned_int > ranks_;
		const bool show_names_;
		const bool fill_empty_ranks_;
		std::vector< nodemap > newicklists;
		static const std::string description_;
};

const std::string NewickTaxonFilter::description_ = "NewickTaxonFilter";



#endif // taxonfilter_hh_
