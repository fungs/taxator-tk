/*
 this program merges different binning results and chooses the deeper predicted
 one

Copyright (C) 2014 Eik Dahms

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


#include <vector>
#include <boost/program_options/cmdline.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/scoped_ptr.hpp>
#include "src/taxontree.hh"
#include "src/ncbidata.hh"
#include "src/bioboxes.hh"
#include "src/taxonomyinterface.hh"
#include "src/utils.hh"

int main ( int argc, char** argv ) {

std::string filenames_input,file_name1,file_name2, use_extra_columns;    
std::list<std::string> filenames, sequence_ids;

namespace po = boost::program_options;
po::options_description desc("Allowed options");
desc.add_options()
    ("help,h", "produce help message")
    ("files,f", po::value<std::string>(&filenames_input)->required(), "filename")
    //("file2,g", po::value<std::string>(&file_name2)->required(), "filename")
    ("extra-columns",po::value<std::string>(&use_extra_columns)->default_value("off"),"To use extra columns \"on\" else \"off\"");

po::variables_map vm;
po::store(po::parse_command_line(argc, argv, desc), vm);

if ( vm.count ( "help" ) ) {
    std::cout << desc << std::endl;
    return EXIT_SUCCESS;
}

po::notify ( vm );  // check required etc.

typedef std::map< std::string, RowValues* > tax_file_map;
std::list<boost::tuple<std::string,tax_file_map>> file_maps;

tokenizeSingleCharDelim(filenames_input, filenames, ",");
for(auto it = filenames.begin(); it != filenames.end(); ++it){
    std::cerr << *it << "\n";
    tax_file_map current_map;
    RowValues* row;
    BioboxesParser bio_parser(*it);
    
    row = bio_parser.getNext();
    while(row){
        if(std::find(sequence_ids.begin(),sequence_ids.end(),row->seqid) == sequence_ids.end()){
            sequence_ids.push_back(row->seqid);
        }
        current_map[row->seqid] = row;
        row = bio_parser.getNext();
    }
    
    file_maps.push_back(boost::tuple<std::string,tax_file_map>(bio_parser.getHeader(),current_map));
}
// create taxonomy
//std::vector< std::string > ranks =  default_ranks;

boost::scoped_ptr< Taxonomy > tax( loadTaxonomyFromEnvironment( &default_ranks ) );
if(!tax) return EXIT_FAILURE;
tax->deleteUnmarkedNodes(); //collapse taxonomy to contain only specified ranks
TaxonomyInterface taxinter (tax.get());

// find classifications in each file for each sequence id and make LCC/LCA

const TaxonNode* outnode;

for(auto id = sequence_ids.begin(); id != sequence_ids.end(); ++id){
    std::list<const TaxonNode*> taxon_nodes;
    for(auto file_iterator = file_maps.begin();file_iterator != file_maps.end(); ++file_iterator){
        if(file_iterator->get<1>().find(*id) != file_iterator->get<1>().end()){
            taxon_nodes.push_back(taxinter.getNode(file_iterator->get<1>()[*id]->taxid));
        }
    }
    outnode = taxinter.getLCC(taxon_nodes);
}


//return;
//BioboxesParser bio_parser_1(file_name1);
//BioboxesParser bio_parser_2(file_name2);
//RowValues* row1;
//RowValues* row2;
//tax_file_map tfmap_1;
//tax_file_map tfmap_2;
////const TaxonNode* outnode;
//
//std::cout << bio_parser_1.getHeader();
//
//while(true){
//    row1 = bio_parser_1.getNext();
//    if(row1){
//        tfmap_1[row1->seqid] = row1;
//    }
//    
//    row2 = bio_parser_2.getNext();
//    if(row2){
//        tfmap_2[row2->seqid] = row2;
//    }
//    if(not row1 and not row2) break;
//}
//
//
//for(tax_file_map::iterator it = tfmap_1.begin(); it != tfmap_1.end(); it++){
//    row1 = it->second;
//    row2 = tfmap_2[it->first];//TODO exception - get - find 
//    if(row1 && row2){
//        assert(row1->seqid == row2->seqid);
//        if(row1->taxid == row2->taxid) outnode = taxinter.getNode(row1->taxid);
//        else outnode = taxinter.getLCC(row1->taxid,row2->taxid);
//       
////        const TaxonNode* node1 = taxinter.getNode(row1->taxid);
////        const TaxonNode* node2 = taxinter.getNode(row2->taxid);
////        if(node1->data->root_pathlength > node2->data->root_pathlength){
////            outnode = node1;
////        }
////        else{
////            outnode = node2;
////        }
//        
//        std::cout << row1->seqid ;
//        
//        if(use_extra_columns == "on"){
//            if(row1->taxid == row2->taxid){
//                std::cout << "\t" <<row1->taxid;
//                for(std::vector<std::string>::iterator it = row1->extra_cols.begin()+2;it != row1->extra_cols.end();it++){
//                    std::cout << "\t" << *it;
//                }
//            }
//            else if(outnode->data->taxid == row1->taxid){
//                std::cout << "\t" <<row1->taxid;
//                for(std::vector<std::string>::iterator it = row1->extra_cols.begin()+2;it != row1->extra_cols.end();it++){
//                    std::cout << "\t" << *it;
//                }
//            }
//            else if (outnode->data->taxid == row2->taxid){// new
//                std::cout << "\t" << row2->taxid;
//                for(std::vector<std::string>::iterator it = row2->extra_cols.begin()+2;it != row2->extra_cols.end();it++){
//                    std::cout << "\t" << *it;
//                }
//            }
//            else{
//                std::cout << "\t" << outnode->data->taxid;
//                for(std::vector<std::string>::iterator it = row2->extra_cols.begin()+2;it != row2->extra_cols.end();it++){
//                    std::cout << "\t" << *it;
//                }
//            }
//        }
//        
//        else std::cout << "\t" << outnode->data->taxid;
//        
//        std::cout << "\n";
//    }
//    //else break;
//};

}
