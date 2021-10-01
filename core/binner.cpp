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
#include <unordered_map>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/regex.hpp>
#include "boost/filesystem.hpp"
#include "src/taxontree.hh"
#include "src/ncbidata.hh"
#include "src/constants.hh"
#include "src/predictionrecordbinning.hh"
#include "src/taxonomyinterface.hh"
#include "src/predictionranges.hh"
#include "src/fastnodemap.hh"
#include "src/exception.hh"
#include "src/bioboxes.hh"

using namespace std;

const std::string extractRegex(const std::string& text, const boost::regex& regex) {
  // special case, empty regex equals full globbing"
  if(!regex.size()) return "consensus_sequence";

  boost::cmatch re_results;
  assert( boost::regex_match( text.c_str(), re_results, regex ) );
  assert(  re_results.size() > 1 );
  assert( re_results[1].first != re_results[1].second );  // match empty string is disallowed atm, use empty regex instead
  return std::string(re_results[1].first, re_results[1].second);
}

int main ( int argc, char** argv ) {

    vector< string > ranks, files;
    bool delete_unmarked;
    large_unsigned_int min_support_in_sample( 0 );
    float signal_majority_per_sequence, min_support_in_sample_percentage( 0. );
    string min_support_in_sample_str, log_filename, sample_identifier, glob_identifier_regex;
    large_unsigned_int min_support_per_sequence;
    boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::size_type num_queries_preallocation;

    namespace po = boost::program_options;
    po::options_description visible_options ( "Allowed options" );
    visible_options.add_options()
    ( "help,h", "show help message" )
    ( "citation", "show citation info" )
    ( "advanced-options", "show advanced program options" )
    ( "version,V", "show program version" )
    ( "sample-identifier,n", po::value< std::string >( &sample_identifier)->required(), "unique sample identifier")
    ( "glob-identifier,g", po::value< std::string >( &glob_identifier_regex )->default_value("(.+)"), "grouping regex for substring matching to glob sequence identifiers for consensus binning (for instance to bin amino acid sequences derived from a long nucleotide string or to do consensus genome binning)")
    ( "sequence-min-support,s", po::value< large_unsigned_int >( &min_support_per_sequence )->default_value( 50 ), "minimum number of positions supporting a taxonomic signal for any single sequence. If not reached, a fall-back on a more robust algorthm will be used" )
    ( "signal-majority,j", po::value< float >( &signal_majority_per_sequence )->default_value( .7 ), "minimum combined fraction of support for any single sequence (> 0.5 to be stable)" )
    ( "identity-constrain,i", po::value< vector< string > >(), "minimum required identity for this rank (e.g. -i species:0.8 -i genus:0.7)")
    ( "files,f", po::value< vector< string > >( &files )->multitoken(), "arbitrary number of prediction files (instead of standard input, use \"-\" in list to specify a combination of stdin and files)" )
    ( "logfile,l", po::value< std::string >( &log_filename )->default_value( "binning.log" ), "specify name of file for logging (appending lines)" );

    po::options_description hidden_options("Hidden options");
    hidden_options.add_options()
    ( "ranks,r", po::value< vector< string > >( &ranks )->multitoken(), "set ranks at which to do predictions" )
    ( "sample-min-support,m", po::value< std::string >( &min_support_in_sample_str )->default_value( "0" ), "minimum support in positions (>=1) or fraction of total support (<1) for any taxon" )
    ( "preallocate-num-queries", po::value< boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::size_type >( & num_queries_preallocation )->default_value( 5000 ), "advanced parameter for better memory allocation, set to number of query sequences or similar (no need to be set)" )
    ( "delete-notranks,d", po::value< bool >( &delete_unmarked )->default_value( true ), "delete all nodes that don't have any of the given ranks (make sure that input taxons are at those ranks)" );

    po::options_description all_options;
    all_options.add( visible_options ).add( hidden_options );

    po::variables_map vm;
    po::store ( po::command_line_parser ( argc, argv ).options ( all_options ).run(), vm );

    if ( vm.count ( "help" ) ) {
        cout << visible_options << endl;
        return EXIT_SUCCESS;
    }

    if ( vm.count ( "citation" ) ) {
        cout << citation_note << endl;
        return EXIT_SUCCESS;
    }

    if ( vm.count ( "advanced-options" ) ) {
        cout << hidden_options << endl;
        return EXIT_SUCCESS;
    }

    if ( vm.count ( "version" ) ) {
        cout << program_version << endl;
        return EXIT_SUCCESS;
    }

    po::notify ( vm );  // check required etc.

    if ( ! vm.count ( "ranks" ) ) ranks = default_ranks;

    // test sequence globbing
    const boost::regex globbing_regex(glob_identifier_regex);  // eg "([^_]+)_.*"

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

    try {
        // STEP 0: PARSING INPUT

        // check input files and setup one parser each
        boost::ptr_list< PredictionFileParser< PredictionRecordBinning > > input;
        if ( files.empty() ) {
            input.push_back(new PredictionFileParser< PredictionRecordBinning > ( std::cin, tax.get() ) );
        } else {
            for(auto file_it = files.begin(); file_it != files.end(); ++file_it ) {
                if( *file_it == "-" ) {
                    input.push_back(new PredictionFileParser< PredictionRecordBinning > ( std::cin, tax.get() ) );
                } else {
                    if( boost::filesystem::exists( *file_it ) ) {
                        input.push_back(new PredictionFileParser< PredictionRecordBinning > ( *file_it, tax.get() ) );
                    } else {
                        cerr << "Could not read file \"" << *file_it++ << "\"" << endl;
                        return EXIT_FAILURE;
                    }
                }
            }
        }

        // collect all records in input files and group by glob identifier
        typedef boost::ptr_list< PredictionRecordBinning > RecordGroup;
        std::unordered_map<std::string, RecordGroup> grouped_records;

        {
          std::string a, b;
          std::string* new_name = &a;
          std::string* old_name = &b;
          RecordGroup* records = NULL;
          const RecordGroup empty_record_group;

          for(auto input_stream = input.begin(); input_stream != input.end(); ++input_stream) {
            for(auto rec = input_stream->next(); rec; rec = input_stream->next()) {

              *new_name = extractRegex(rec->getQueryIdentifier(), globbing_regex);

              // debug output
              // std::cerr << "extracted identifier is '" << *new_name << "'" << std::endl;
              // std::cerr << rec->getQueryIdentifier() << " --> " << *new_name << std::endl;
              // std::cerr << "entry is: " << *rec;

              if ( *new_name != *old_name ) { // TODO: make sure that regex matching fails, if no match can be found
                  // insert using emplace
                  records = &(grouped_records.emplace(
                    *new_name,
                    empty_record_group // clone empty list to avoid construction
                    // std::initializer_list<RecordGroup>({})
                  ).first->second);

                  std::swap(old_name, new_name);
              }
              records->push_back(rec); // transfer record ownership
            }
          }
        }

        // STEP 1: RANGE PRUNING
        // in this step the overall sample support for each node is recorded and each
        // range is shrunk such that the remaining nodes have a minimum support (unit is sequence positions)

        //counting support of nodes
        std::cerr << "Analyzing sample composition: ";
        large_unsigned_int minimum_support_found = std::numeric_limits< large_unsigned_int >::max();
        const TaxonNode* const root_node = taxinter.getRoot();
        FastNodeMap< large_unsigned_int > support( taxinter.getMaxDepth() );
        large_unsigned_int& root_support = support[ root_node ];
        for ( auto group_it = grouped_records.begin(); group_it != grouped_records.end(); ++group_it ) {
            RecordGroup& records = group_it->second;
            // const std::string& identifier = group_it->first;

            // std::cerr << "processing glob identifier" << identifier << std::endl;
            for ( RecordGroup::const_iterator prec_it = records.begin(); prec_it != records.end(); ++prec_it ) {
                Taxonomy::PathUpIterator pit = taxinter.traverseUp( prec_it->getLowerNode() );

                // process lowest node
                large_unsigned_int total_node_support = prec_it->getSupportAt( &*pit );
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
        std::cerr << support.size() << " nested taxa with total support of " << support[ root_node ] << " positions" << std::endl;

        // if min_support_in_sample was given as fraction
        if ( min_support_in_sample_percentage ) min_support_in_sample = support[ root_node ]*min_support_in_sample_percentage;

        // shrink ranges from lower end if support is smaller than the minimum required or if it does not comply with user-defined PID per rank.
        std::cerr << "Noise removal: ";
        std::set< const TaxonNode* > pruned_nodes;
        if ( minimum_support_found < min_support_in_sample ) {
            for ( auto group_it = grouped_records.begin(); group_it != grouped_records.end(); ++group_it ) {
                RecordGroup& records = group_it->second;
                for ( RecordGroup::iterator prec_it = records.begin(); prec_it != records.end(); ) {
                    const TaxonNode* lower_node = prec_it->getLowerNode();
                    const TaxonNode* upper_node = prec_it->getUpperNode();

                    Taxonomy::PathUpIterator pit = taxinter.traverseUp( lower_node );
                    while ( pit != upper_node && support[ &*pit ] < min_support_in_sample ) {
                        pruned_nodes.insert( &*pit );
                        ++pit;
                    }

                    if ( pit == upper_node && support[ &*pit ] < min_support_in_sample ) { //remove whole range
                        pruned_nodes.insert( &*pit );
                        prec_it = records.erase( prec_it ); //TODO: mask instead of delete TODO: is this a memory leak?
                        continue;
                    }

                    if ( pit != lower_node ) prec_it->pruneLowerNode( &*pit ); //prune
                    ++prec_it; // because erase already increments
                }
            }
        }
        std::cerr << pruned_nodes.size() << " taxa removed" << std::endl;

        // STEP 2: BINNING and output
        // in this step multiple ranges are combined into a single range by combining
        // evidence for sub-ranges. This algorithm considers only support. Signal
        // strength and interpolation values are ignored. This heuristic seems quite
        // robust

        std::cerr << "Consensus taxonomy assignment: ";
        std::ofstream binning_debug_output( log_filename.c_str() );

        // use Bioboxes output format
        const std::vector<std::tuple<const std::string, const std::string>> custom_header_tags = {std::make_tuple("Version", program_version)};
        const std::vector<std::string> custom_column_tags = {"Support", "Length"};
        std::vector<std::string> extra_cols(2);
        BioboxesBinningFormat binning_output(BioboxesBinningFormat::ColumnTags::taxid, sample_identifier, taxinter.getVersion(), std::cout, "TaxatorTK", custom_header_tags, custom_column_tags);

        // consensus and output
        for ( auto it = grouped_records.begin(); it != grouped_records.end(); ++it ) {
            const std::string& identifier = it->first;
            RecordGroup& records = it->second;
            if( records.empty() ) continue;  // TODO: suppress empty lists
            boost::scoped_ptr< PredictionRecordBinning > prec_sptr;
            const PredictionRecordBinning* prec;
            if ( records.size() > 1 ) { //run combination algo for sequence segments
                prec_sptr.reset( combinePredictionRanges( records, identifier, tax.get(), signal_majority_per_sequence, min_support_per_sequence, binning_debug_output ) );
                prec = prec_sptr.get();
            } else { // pass-through segment prediction for whole sequence
                prec = &records.front();
            }
            // apply user-defined constrain
            if ( prec->getUpperNode() != root_node && ! pid_per_rank.empty() ) {
                const double seqlen = static_cast< double >( prec->getQueryLength() );
                float min_pid = 0.; //enforce consistency when walking down
                map< const string*, float >::const_iterator find_it;
                const TaxonNode* predict_node = root_node;
                const TaxonNode* target_node = prec->getUpperNode();
                const float rank_pid = prec->getSupportAt( target_node )/seqlen;
                Taxonomy::CPathDownIterator pit = taxinter.traverseDown<Taxonomy::CPathDownIterator>( target_node );
                do {
                    pit++;
                    find_it = pid_per_rank.find( &(pit->data->annotation->rank) );
                    if ( find_it != pid_per_rank.end() ) min_pid = max( min_pid, find_it->second );
                    binning_debug_output << "constraint ctrl: " << rank_pid << " >= " << min_pid << " ?" << endl;
                    if ( rank_pid < min_pid ) break;
                    predict_node = &*pit;
                } while ( pit != target_node );
                extra_cols[0] = boost::lexical_cast<std::string>(prec->getSupportAt(predict_node));
                extra_cols[1] = boost::lexical_cast<std::string>(prec->getQueryLength());
                binning_output.writeBodyLine(identifier, predict_node->data->taxid, extra_cols);
            } else {
                extra_cols[0] = boost::lexical_cast<std::string>(prec->getSupportAt(prec->getUpperNode()));
                extra_cols[1] = boost::lexical_cast<std::string>(prec->getQueryLength());
                binning_output.writeBodyLine(identifier, prec->getUpperNode()->data->taxid, extra_cols);
            }
        }
        std::cerr << " done" << std::endl;

        return EXIT_SUCCESS;
    } catch(Exception &e) {
        cerr << "An unrecoverable error occurred." << endl;
        cerr << boost::diagnostic_information(e) << endl;
        return EXIT_FAILURE;
    }
}
