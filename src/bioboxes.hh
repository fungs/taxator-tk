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

#ifndef bioboxes_hh_
#define bioboxes_hh_

#include <iostream>
#include <ostream>
#include <vector>
#include <tuple>
#include <fstream>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/algorithm/string.hpp>


class BioboxesBinningFormat{  // implements Bioboxes.org binning format 0.9
public:
    enum class ColumnTags {
        taxid,
        binid,
        taxid_binid
    };
    
    BioboxesBinningFormat(
        const ColumnTags cols,
        const std::string& sampleid,
        const std::string& taxonomyid = "",
        std::ostream& ostr = std::cout,
        const std::string& custom_tag_prefix = std::string(),
        const std::vector<std::tuple<const std::string, const std::string>> custom_header_tags = std::vector<std::tuple<const std::string, const std::string>>(),  // TODO: easier syntax?
        const std::vector<std::string> custom_column_tags = std::vector<std::string>()  // TODO: easier syntax?
        );
    
    ~BioboxesBinningFormat();
    
    void writeBodyLine(
        const std::string& sequenceid,
        const std::string& singleid
    );
    
    void writeBodyLine(
        const std::string& sequenceid,
        const std::string& singleid,
        const std::vector<std::string>& columns_custom
    );

    void writeBodyLine(
        const std::string& sequenceid,
        const std::string& binid,
        const std::string& taxid
    );
    
    void writeBodyLine(
        const std::string& sequenceid,
        const std::string& binid,
        const std::string& taxid,
        const std::vector<std::string>& columns_custom
    );
    
private:
    void writeHeader(const std::string& sampleid, const std::string& taxonomyid);
    
    void writeHeaderCustom(const std::vector<std::tuple<const std::string, const std::string>>& custom_header_tags);
    
    void writeHeaderColumnTags(const ColumnTags cols);
    
    void writeHeaderColumnTagsCustom(const std::vector<std::string>& custom_column_tags);
    
    std::ostream& ostr_;
    const ColumnTags cols_;
    const std::string custom_tag_prefix_;
    const std::string format_version_ = "0.9.1";
};


struct RowValues {
        std::string seqid;
        std::string taxid;
        std::string binid;
        std::string taxid_binid;
        std::vector<std::string> extra_cols;
        std::string line;
    };


class BioboxesParser{
public:
    enum class ColumnTags {
        taxid,
        binid,
        taxid_binid
    };
    
    
    //typedef std::vector<std::string> row_values;
    //row_values current_vals;
    
    
    BioboxesParser(std::string filename) : filehandle(filename), handle(filehandle){
        //TODO make input more Variable // read header exactly
        std::string line;
        while(std::getline( handle, line )){
            header += line + "\n";
            if(boost::starts_with(line,"@@")){
                std::vector<std::string> SplitVec;
                boost::split( SplitVec, line, boost::is_any_of("\t"), boost::token_compress_on );
                num_columns = SplitVec.size();
                break;
            }
        }
    };
    
    inline void destroyRow( const RowValues* row ) const {
        delete row;
    }
    
    bool ignoreLine(std::string line){
        return line.empty();
    }
    
    RowValues* getNext(){
        std::string line;
        RowValues* row = new RowValues();
        while( std::getline( handle, line ) ) {    
            if(ignoreLine( line )) continue;
            std::vector<std::string> SplitVec; // #2: Search for tokens
            boost::split( SplitVec, line, boost::is_any_of("\t"), boost::token_compress_on ); 
            row->seqid = SplitVec[0];
            row->taxid = SplitVec[1];
            row->extra_cols = SplitVec;
            row->line = line;
            return row;
            std::cout << line << "\n";
        }
        destroyRow(row);
        return NULL;
    };
    
    std::string getHeader(){
        return header;
    };
        
private:
    std::ifstream filehandle;
    std::istream& handle;
    std::string header = "";
    int num_columns;
    //std::ostream& ostr_;
    const ColumnTags cols_;
    const std::string custom_tag_prefix_;
    const std::string format_version_ = "0.9.1";
    
};

#endif // bioboxes_hh_
