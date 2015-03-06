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
#include <climits> // std numeric limits
#include <stack>
#include <string>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
//#include <boost/algorithm/string.hhp>
#include "boost/filesystem.hpp"
#include "src/taxontree.hh"
#include "src/ncbidata.hh"
#include "src/constants.hh"
#include "src/predictionrecordbinning.hh"
#include "src/taxonomyinterface.hh"
//#include "src/predictionranges.hh"
#include "src/fastnodemap.hh"
//#include "src/taxalist.hh"



using namespace std;

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
			else find_it = supports.insert( std::make_pair< const TaxonNode*, float >( node, it->get<2>()[ it->get<0>()->data->root_pathlength ] ) ).first;

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


class TaxaList{
    public:
        TaxaList(TaxonomyInterface taxint, string whitelist, string blacklist, string realdata): taxinter_(taxint){
        num_taxa = 0;
        total_taxa = 0;
        total_support = 0;
        root_node = taxinter_.getRoot();
        
       

        if(whitelist != ""){
            string line;
            ifstream infile (whitelist.c_str());
            if (infile.is_open()){
                while( getline(infile, line)){
                    string taxa = line.c_str();
                    const TaxonNode* current = taxinter_.getNode(taxa);
                    while(current != root_node){
                        taxa = current->data->taxid;
                        //std::cerr << taxa << "\n";
                        if(whitelist_.find(taxa) == whitelist_.end()){
                            whitelist_.insert(taxa);
                            //taxa_prosperties[taxa] = boost::make_tuple(std::numeric_limits<int>::max(), std::numeric_limits<large_unsigned_int>::max());
                        }
                        current = current->parent;
                    }
                }
            }
        }
        //std::cerr << "whitelist loaded \n";

        if(blacklist != ""){
            string line;
            ifstream infile (whitelist.c_str());
            if (infile.is_open()){
                while( getline(infile, line)){
                    string taxa = line.c_str();
                    if(blacklist_.find(taxa) == blacklist_.end()){
                        blacklist_.insert(taxa);
                    }
                }
            }
        }

        //std::cerr << "blacklist loaded \n";
        
        if(realdata != ""){
            string line;
            ifstream infile (realdata.c_str());
            if (infile.is_open()){
                while( getline(infile, line)){
                    string taxa = line.c_str();
                    const TaxonNode* current = taxinter_.getNode(taxa);
                    while(current != root_node){
                        taxa = current->data->taxid;
                        //std::cerr << taxa << "\n";
                        if(realdata_.find(taxa) == realdata_.end()){
                            realdata_.insert(taxa);
                
                        }
                        current = current->parent;
                    }
                }
            }
            
            //std::set< std::string >::iterator it;
            //for(it = realdata_.begin(); it != realdata_.end(); ++it){
           //     std::cerr << *it << "\n";
            //}
        }
        
        }

        typedef std::set< std::string > TaxaSet;
        typedef TaxaSet::iterator iterator;
        typedef TaxaSet::const_iterator const_iterator;
        iterator begin() {return taxalist.begin();}
        iterator end() {return taxalist.end();}
//
//        void load(std::string filename){
//            //taxafile_name_ = filename;
//            std::string line;
//            std::ifstream infile (filename.c_str());
//            if (infile.is_open()){
//                while( getline(infile, line)){
//                    //int taxa = atoi(line.c_str());
//                    std::string taxa = line.c_str();
//                    if (taxalist.find(taxa) == taxalist.end()){
//                        taxalist.insert(taxa);
//                    }
//                }
//            }
//        }

        void print(){
            std::set< std::string >::iterator it;
            for(it = taxalist.begin(); it != taxalist.end(); ++it){
                std::cerr << "taxa: " << *it << "\n";
            }
        }

        void print_taxa_prosp(string sample_path){
            int in_sample;
            int positive = 0;
            int false_ = 0;
            int false_positive = 0;
            int false_negative = 0;
            int true_positive = 0;
            int true_negative = 0;
            int number_of_taxa = taxa_prosperties.size();
            std::ofstream info_stream((sample_path+string("_info.txt")).c_str());
            std::ofstream list_file((sample_path+string("_list.tsv")).c_str());
            
            info_stream << "total support\t" << total_support << "\n";
            info_stream << "total number of taxa\t" << total_taxa << "\n";
            info_stream << "number of different taxa in sample\t" << number_of_taxa << "\n";
            info_stream << "actual number (after filter)\t" << (number_of_taxa-filtered_taxa) << "\n" 
                        << "number of occuring in the sample(from real data)\t" << realdata_.size() << "\n\n\n";
            
            list_file << "#taxid\tnumber of samples\tsupport\tlength\tival\tsupport%\tmean support\tmean length\tmean ival\trank\tin list\tin sample\tcorrect\n";
               
            for(taxa_list::iterator it = taxa_prosperties.begin(); it != taxa_prosperties.end(); ++it){
                if(realdata_.find(it->first) != realdata_.end()){in_sample = 1;}
                else{in_sample = 0 ;}
                string taxid = it->first, rank;
                int num_taxa, filtered, correct(0);
                large_unsigned_int support, length, dep_support;
                float ival;
                
                boost::tuples::tie(num_taxa, support, length, ival, rank, dep_support ,filtered) = it->second;
                
                if(int(filtered && in_sample)){
                    correct = 1;
                    true_positive += 1;
                }else if(int(!filtered && !in_sample)){
                    correct = 1;
                    true_negative += 1;
                }
                positive += correct;
                if(not correct){
                    if (filtered) false_positive += 1;
                    else false_negative +=1;
                    false_ +=1;
                }
                
                list_file << std::fixed << taxid << "\t" << num_taxa << "\t" << support << "\t" << length << "\t" << ival << "\t"
                          << double(support)/double(length) << "\t" // percentage of support 
                          << double(support)/double(num_taxa) << "\t" // mean support
                          << length/num_taxa << "\t" // mean length
                          << ival/num_taxa << "\t" // mean ival
                          << rank << "\t"//rank
                          << filtered << "\t"
                          << in_sample << "\t"
                          << correct << "\n"
                        ;
            }
            
            info_stream << "positive\t" << positive << "\n";
            info_stream << "true positive\t" << true_positive << "\n";
            info_stream << "true negative\t" << true_negative << "\n";
            info_stream << "false\t" << false_ << "\n";
            info_stream << "false positive\t" << false_positive << "\n";
            info_stream << "false negative\t" << false_negative << "\n";
            
        }

        bool is_in_list(std::string taxid){
            return taxalist.find(taxid) != taxalist.end();
        }

        void add(std::string taxid){
            taxalist.insert(taxid);
            //std::cerr << "add" << taxid << "\n";
        }

        void add(PredictionRecordBinning &prec){

            const TaxonNode* current_node = prec.getUpperNode();
            const TaxonNode* start_node = current_node;

            //for(TaxonomyInterface::CPathUpIterator)

            while(current_node != root_node){
                string taxid = current_node->data->taxid;

                if(blacklist_.find(taxid) != blacklist_.end()){
                    current_node = current_node->parent;
                    continue;
                }

                std::string rank = current_node->data->annotation->rank;
                
                if(ranks.find(rank) == ranks.end()) ranks.insert(rank);
                
                if(taxalist.find(taxid) == taxalist.end()){
                    taxalist.insert(taxid);
                    std::vector< large_unsigned_int > emptylist;
                    supportlist_[taxid] = emptylist;
                    value_list empty_tup;
                    taxa_prosperties[taxid] = empty_tup;
                    taxa_prosperties[taxid].get<0>() = 1;
                    taxa_prosperties[taxid].get<1>() = prec.getSupportAt_full(current_node);
                    taxa_prosperties[taxid].get<2>() = prec.getQueryLength();
                    taxa_prosperties[taxid].get<3>() = prec.getInterpolationValue();
                    taxa_prosperties[taxid].get<4>() = rank;//current_node->data->annotation->rank;
                    taxa_prosperties[taxid].get<5>() = prec.getSupportAt(start_node); //TODO reduce to one variable
                    taxa_prosperties[taxid].get<6>() = 1;
                    
                    //supportlist_[taxid] = emptylist;
                    supportlist_[taxid].push_back(prec.getSupportAt(current_node));
                }
                //else if(taxa_prosperties[taxid].get<0>() != std::numeric_limits<int>::max()){
                else{
                    taxa_prosperties[taxid].get<0>() += 1;
                    taxa_prosperties[taxid].get<1>() += prec.getSupportAt_full(current_node);
                    taxa_prosperties[taxid].get<2>() += prec.getQueryLength();
                    taxa_prosperties[taxid].get<3>() += prec.getInterpolationValue();
                    taxa_prosperties[taxid].get<5>() += prec.getSupportAt(start_node);
                    
                    supportlist_[taxid].push_back(prec.getSupportAt(current_node));
                }

                total_taxa += 1;
                total_support += prec.getSupportAt_full(current_node);
                current_node = current_node->parent;
            }
        }

        void reduce_by_support(float treshold){
            int erased_tax = 0;
            filtered_taxa = 0;
            uint total_support_treshold = total_support*treshold;
            std::set< std::string > to_erase;
            
            //support-filter
            
            for(taxa_list::iterator it = taxa_prosperties.begin(); it != taxa_prosperties.end(); ++it){
                //std::cerr << it->first << " - " << it->second.get<1>() << ": "<< total_support_treshold;
                if( it->second.get<1>() < total_support_treshold && whitelist_.find(it->first) == whitelist_.end()){
                    //std::cerr << "erased\n";
                    erased_tax += it->second.get<0>();
                    to_erase.insert(it->first);
                    it->second.get<6>() = 0;
                    filtered_taxa += 1;
                    //taxalist.erase(it->first);
                    //taxa_prosperties.erase(it);
                }
                //else std::cerr << "\n";
            }
            for(std::set< string >::iterator it = to_erase.begin();it != to_erase.end(); ++it ){
                //taxalist.erase(*it);
                //taxa_prosperties.erase(*it);
            }
            //total_taxa -= erased_tax;
            
            //add parents of all not filtered taxa
            
            for(taxa_list::iterator it = taxa_prosperties.begin(); it != taxa_prosperties.end(); ++it){
                if(it->second.get<6>()){
                    const TaxonNode* current_node = taxinter_.getNode(it->first)->parent;
                    
                    while(current_node != root_node){
                        taxa_prosperties.find(current_node->data->taxid)->second.get<6>() = 1;
                        current_node = current_node->parent;
                    }
                }
            }
            
            filtered_taxa = 0;
            
            for(taxa_list::iterator it = taxa_prosperties.begin(); it != taxa_prosperties.end(); ++it){
                if(!it->second.get<6>()){
                    filtered_taxa += 1;
                }
            }
        }
        
        void reduce_by(string filter_type, float treshold){
            int erased_tax = 0;
            filtered_taxa = 0;
            large_unsigned_int mean_total_support = total_support/total_taxa; 
            std::set< std::string > to_erase;
            int filtered = 0;
            //support-filter
            
            for(taxa_list::iterator it = taxa_prosperties.begin(); it != taxa_prosperties.end(); ++it){
                
                if(filter_type == "number_per_taxa") filtered = it->second.get<0>() < total_taxa*treshold;
                if(filter_type == "support") filtered = it->second.get<1>() < total_support*treshold;
                if(filter_type == "mean_support"){
                    filtered = it->second.get<1>()/it->second.get<0>() < mean_total_support*treshold;
                    //std::cerr << it->second.get<1>()/it->second.get<0>() << " < " << mean_total_support*treshold << " " << filtered << "\n";
                
                }
                if(filter_type == "median_support"){
                    large_unsigned_int median_support = get_median(supportlist_[it->first]);
                    filtered = median_support < mean_total_support*treshold ;
                    //std::cerr << median_support << " < " << mean_total_support*treshold << " " << filtered << "\n";
                
                }
                
                if(filtered && whitelist_.find(it->first) == whitelist_.end()){
                    //std::cerr << "erased\n";    
                    erased_tax += it->second.get<0>();
                    to_erase.insert(it->first);
                    it->second.get<6>() = 0;
                    filtered_taxa += 1;
                    //taxalist.erase(it->first);
                    //taxa_prosperties.erase(it);
                }
                //else std::cerr << "\n";
            }
            for(std::set< string >::iterator it = to_erase.begin();it != to_erase.end(); ++it ){
                //taxalist.erase(*it);
                //taxa_prosperties.erase(*it);
            }
            //total_taxa -= erased_tax;
            
            //add parents of all not filtered taxa
            
            for(taxa_list::iterator it = taxa_prosperties.begin(); it != taxa_prosperties.end(); ++it){
                if(it->second.get<6>()){
                    const TaxonNode* current_node = taxinter_.getNode(it->first)->parent;
                    
                    while(current_node != root_node){
                        taxa_prosperties.find(current_node->data->taxid)->second.get<6>() = 1;
                        current_node = current_node->parent;
                    }
                }
            }
            
            filtered_taxa = 0;
            
            for(taxa_list::iterator it = taxa_prosperties.begin(); it != taxa_prosperties.end(); ++it){
                if(!it->second.get<6>()){
                    filtered_taxa += 1;
                }
            }
        }

        int get_number(){
            return total_taxa;
        }
        
        large_unsigned_int get_median(std::vector<large_unsigned_int> input_list){
            std::sort(input_list.begin(),input_list.end());
            int n = input_list.size();
            if(n%2 == 0){
                return (input_list[n]+input_list[n+1])/2;
            }
            else{
                return input_list[(n+1)/2];
            }
        }

    protected:
        typedef boost::tuple< int,large_unsigned_int, large_unsigned_int, float, std::string, large_unsigned_int, int > value_list;
        typedef std::map<std::string, value_list> taxa_list;
        typedef std::map<std::string, std::vector< large_unsigned_int > > supportlist;
        TaxonomyInterface taxinter_;
        int num_taxa,total_taxa;
        large_unsigned_int total_support, filtered_taxa;
        const TaxonNode* root_node;
        TaxaSet taxalist, whitelist_, blacklist_, realdata_;
        taxa_list taxa_prosperties;
        supportlist supportlist_;
        const std::string taxafile_name_;
        std::set< std::string > ranks;

};

int main ( int argc, char** argv ) {

	vector< string > ranks, files;
	bool delete_unmarked;
	large_unsigned_int min_support_in_sample( 0 );
	float signal_majority_per_sequence, min_support_in_sample_percentage( 0. ), treshold_;
	string min_support_in_sample_str, tax_list_input, log_filename, source_binning, support_mode, blacklist_file, whitelist_file, realdata_file, filter_by_, sample_path_;
	large_unsigned_int min_support_per_sequence;
	boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::size_type num_queries_preallocation;



	namespace po = boost::program_options;
	po::options_description visible_options ( "Allowed options" );
	visible_options.add_options()
	( "help,h", "show help message" )
	( "advanced-options", "show advanced program options" )
	( "sequence-min-support,s", po::value< large_unsigned_int >( &min_support_per_sequence )->default_value( 50 ), "minimum number of positions supporting a taxonomic signal for any single sequence. If not reached, a fall-back on a more robust algorthm will be used" )
	( "signal-majority,j", po::value< float >( &signal_majority_per_sequence )->default_value( .7 ), "minimum combined fraction of support for any single sequence (> 0.5 to be stable)" )
	( "identity-constrain,i", po::value< vector< string > >(), "minimum required identity for this rank (e.g. -i species:0.8 -i genus:0.7)")
	( "files,f", po::value< vector< string > >( &files )->multitoken(), "arbitrary number of prediction files (replaces standard input, use \"-\" to specify a combination of both)" )
	( "logfile,l", po::value< std::string >( &log_filename )->default_value( "binning.log" ), "specify name of file for logging (appending lines)" )
	( "taxalist,t",po::value< std::string >(&tax_list_input)->default_value("none"),"use list for taxa in sample")
	( "source-binning,b",po::value< std::string >(&source_binning)->default_value("none"),"source binning file - used to analyze taxa distribution")
	( "support-mode,e",po::value< std::string >(&support_mode)->default_value("direct"),"support mode - direct or total")
        ( "whitelist",po::value< std::string >(&whitelist_file)->default_value(""),"give a whitelist")
        ( "blacklist",po::value< std::string >(&blacklist_file)->default_value(""),"add blacklist")
        ( "real_data",po::value< std::string >(&realdata_file)->default_value(""),"real data list")
        ( "filter_treshold",po::value< float >(&treshold_)->default_value(.0),"value for treshold")
        ( "filter_by",po::value< string >(&filter_by_)->default_value(""),"which value should be filtered")
        ( "sample_path_name",po::value< string >(&sample_path_)->default_value(""),"path to sample file with name");
        


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
	po::notify ( vm );

	if ( vm.count ( "help" ) ) {
		cout << visible_options << endl;
		return EXIT_SUCCESS;
	}

	if ( vm.count ( "advanced-options" ) ) {
		cout << hidden_options << endl;
		return EXIT_SUCCESS;
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
        TaxaList taxalist (taxinter,whitelist_file,blacklist_file,realdata_file);


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
	// parse taxa-list is given

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

//TODO read files
// generate Taxalist with most abundant taxa
//////    //if( vm.count ( "taxalist" ) ){
//////        //taxalist.load(vm["taxalist"].as<string>());
//////        //taxalist.print();
//////	//}
//
//	if (vm.count( "source-binning" )){
//        std::list< boost::tuple< std::string , int > > binner_taxa_list;
//
//        std::string line;
//        std::ifstream infile (vm["source-binning"].as<string>().c_str());
//        if (infile.is_open()){
//            while( getline(infile, line)){
//                std::vector< std::string > fields;
//                boost::split(fields,line,boost::is_any_of("\t"));
//
//                std::string current_taxa = fields[1];
//                int n;
//                for(std::list< boost::tuple< std::string, int > >::const_iterator it = binner_taxa_list.begin(); it != binner_taxa_list.end(); ++it){
//
//                    if(*it->get<0>() == current_taxa){
//                        *it->get<1>() = *it->get<1>+1);
//                        break;
//                    }
//                    else if(std::next(it) == binner_taxa_list.end()){
//                        binner_taxa_list.push_back(boost::make_tuple(current_taxa,1));
//                        break;
//                    }
//                }
//                //int taxa = atoi(line.c_str());
//                //if (taxalist.find(taxa) == taxalist.end()){
//                //    taxalist.insert(taxa);
//                }
//            }
//        }
//	}

// STEP 1: Find abundant taxa
    std::list< std::string > taxa_list;

    for ( boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::iterator query_it = predictions_per_query.begin(); query_it != predictions_per_query.end(); ++query_it ) {
		for ( boost::ptr_list< PredictionRecordBinning >::iterator prec_it = query_it->begin(); prec_it != query_it->end(); ++prec_it ) {
            //if(std::find(taxa_list.begin(),taxa_list.end(),prec_it->getUpperNode()->data->taxid) == taxa_list.end()){
            //std::cerr << "current taxid: " << prec_it->getUpperNode()->data->taxid << "\n" ;

//            taxalist.add(*prec_it);
//            if (not taxalist.is_in_list(prec_it->getUpperNode()->data->taxid)){
//                //std::cerr << "added";
//                taxalist.add(prec_it->getUpperNode()->data->taxid);
                  taxalist.add(*prec_it);
//                //std::cerr << prec_it->getUpperNode()->data->taxid << "\n";
//            }
		}
    }

    //std::cerr << taxalist.get_number() << "\n";
    //taxalist.reduce_by_support(treshold_);
    taxalist.reduce_by(filter_by_,treshold_);
    //std::cerr << taxalist.get_number() << "\n";
    taxalist.print_taxa_prosp(sample_path_);
    std::exit(0);
        const TaxonNode* root_node = taxinter.getRoot();


	// STEP 2: BINNING
	// in this step multiple ranges are combined into a single range by combining
	// evidence for sub-ranges. This algorithm considers only support. Signal
	// strength and interpolation values are ignored. This heuristic seems quite
	// robust
	// Prediction Record Binning weg !!!!!!!!!
    //taxalist.print();
	//std::cerr << "binning step... ";
	std::ofstream binning_debug_output( log_filename.c_str() );
	for ( boost::ptr_vector< boost::ptr_list< PredictionRecordBinning > >::iterator it = predictions_per_query.begin(); it != predictions_per_query.end(); ++it ) {
		if( it->empty() ) continue;
		boost::scoped_ptr< PredictionRecordBinning > prec_sptr;
		const PredictionRecordBinning* output_prec;
		//if ( it->size() > 1 ) { //run combination algo for sequence segments
            using namespace details;
            TupleRangeCombineList tlist;
            //TaxonomyInterface taxinter( tax.get() );
            boost::ptr_list< PredictionRecordBinning >& predictions = *it;
            PredictionRecordBinning* prec = new PredictionRecordBinning (tax.get());
			{ // copy values
                const PredictionRecordBinning tmp = predictions.front();
                prec->setQueryIdentifier( tmp.getQueryIdentifier() );
                prec->setQueryLength( tmp.getQueryLength() );
                prec->setQueryFeatureBegin( 1 ); //TODO: range select
                prec->setQueryFeatureEnd( tmp.getQueryLength() ); //TODO: range select
            }

            // initialize temporary data structure (list of tuples)
            medium_unsigned_int summed_support = 0;
            for ( boost::ptr_list< PredictionRecordBinning >::iterator itp = predictions.begin(); itp != predictions.end(); ++itp ) {
                const TaxonNode* lower_node = itp->getRtax();
                int i = lower_node->data->root_pathlength;
                medium_unsigned_int support = itp->getSupportAt( i , true );
                //std::cerr << "got support\n";
                summed_support += support;

                const std::size_t depth = i + 1;
                tlist.push_back( boost::make_tuple( taxinter.traverseDown< Taxonomy::CPathDownIterator >( lower_node ), std::vector<medium_unsigned_int>( depth ), std::vector<medium_unsigned_int>( depth ), false ) );

                std::vector< medium_unsigned_int >& direct_support = tlist.back().get<1>(); //TODO: initialize direct_support with taxon_support
                std::vector< medium_unsigned_int >& total_support = tlist.back().get<2>();
                total_support[i] = direct_support[i] = support;

        // 		debug_output << *it;
        // 		debug_output << " summed support at level " << i << " is " << direct_support[i] << " (" << total_support[i] << ")" << std::endl;
                while ( --i >= 0 ) {
                    direct_support[i] = itp->getSupportAt( i , true );
                    total_support[i] = std::max( total_support[i+1], direct_support[i] );
        // 			debug_output << "summed support at level " << i << " is " << direct_support[i] << " (" << total_support[i] << ")" << std::endl;
                }
            }


            medium_unsigned_int direct_support_thresh = std::max( static_cast< medium_unsigned_int >( signal_majority_per_sequence*summed_support ), static_cast< medium_unsigned_int >(min_support_per_sequence) );
            medium_unsigned_int direct_support, total_support;
            std::vector< boost::tuple<const TaxonNode*, medium_unsigned_int, medium_unsigned_int, bool> > path;
            //const TaxonNode* root_node = taxinter.getRoot();

            //
            //medium_unsigned_int total_support_tresh = std::max( )

            // set values for root node (TODO: put into initialization loop)
            setPathEndState( tlist );
            getSupport( tlist, direct_support, total_support );

            // walk down each path
            int lower_direct_node_index = -1;
            int running_index = 0;

            while ( ! tlist.empty() ) {
                const TaxonNode* node = &*tlist.front().get<0>();
                if( vm["support-mode"].as<string>() == "direct" ){
                    if ( direct_support >= direct_support_thresh ) lower_direct_node_index = running_index;
                }
                else{
                    if ( total_support >= direct_support_thresh ) lower_direct_node_index = running_index;
                }
                path.push_back( boost::make_tuple( node, direct_support, total_support, false ) );
                removeIf<3>( tlist ); //remove paths that have ended
                stepDown( tlist ); //forward PathIterators
                ++running_index;
                setPathEndState( tlist ); //update bool flags for range end
                path.back().get<3>() = reduceToMajority( tlist ); //set parent branching flag
                getSupport( tlist, direct_support, total_support );
            }



            bool setp = false;
            //std::cerr << prec->getQueryIdentifier() << "\t" ;
            for ( int i = path.size() - 1; i >= 0; --i ) {
                 //std::cerr << path[i].get<0>()->data->annotation->name << "->" ;
                //int taxid = (atoi(path[i].get<0>()->data->taxid.c_str()));
                if ( taxalist.is_in_list(path[i].get<0>()->data->taxid)) { //path[i].get<2>() >= direct_support_thresh && \\TODO after debug
                    prec->setNodePoint( path[i].get<0>(), path[i].get<2>() );
                    setp = true;
                    break;//return prec;

                }
                else{

                }
            }
            //std::cerr << std::endl;
            if(not setp){
            prec->setNodePoint( path[0].get<0>(), path[0].get<2>() );
            }
			////prec_sptr.reset( combinePredictionRangesRtax( *it, tax.get(), signal_majority_per_sequence, min_support_per_sequence, taxalist, binning_debug_output) );
			output_prec = prec;
		//} else { // pass-through segment prediction for whole sequence
		//	output_prec = &it->front();
// 			output_prec->setBinningType( PredictionRecordBinning::single );
		//}
		// apply user-defined constrain
		if ( output_prec->getUpperNode() != root_node && ! pid_per_rank.empty() ) {
			const double seqlen = static_cast< double >( output_prec->getQueryLength() );
			float min_pid = 0.; //enforce consistency when walking down
			map< const string*, float >::const_iterator find_it;
			const TaxonNode* predict_node = root_node;
			const TaxonNode* target_node = output_prec->getUpperNode();
			const float rank_pid = output_prec->getSupportAt( target_node , true)/seqlen;
			Taxonomy::CPathDownIterator pit = taxinter.traverseDown<Taxonomy::CPathDownIterator>( target_node );
			do {
				pit++;
				find_it = pid_per_rank.find( &(pit->data->annotation->rank) );
				if ( find_it != pid_per_rank.end() ) min_pid = max( min_pid, find_it->second );
				binning_debug_output << "constraint ctrl: " << rank_pid << " >= " << min_pid << " ?" << endl;
				if ( rank_pid < min_pid ) break;
				predict_node = &*pit;
			} while ( pit != target_node );
			std::cout << output_prec->getQueryIdentifier() << tab << predict_node->data->taxid << tab << output_prec->getSupportAt(predict_node, true) << endline;
		} else {
			std::cout << output_prec->getQueryIdentifier() << tab << output_prec->getUpperNode()->data->taxid << tab << output_prec->getSupportAt(output_prec->getUpperNode(), true) << endline;
		}
	}
	//std::cerr << " done" << std::endl;
	return EXIT_SUCCESS;
}

