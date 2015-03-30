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
        const std::vector<std::tuple<std::string, std::string>> custom_header_tags = std::vector<std::tuple<std::string, std::string>>(),  // TODO: easier syntax?
        const std::vector<std::string> custom_column_tags = std::vector<std::string>()  // TODO: easier syntax?
        );
    
    void writeBodyLine(
        const std::string& sequenceid,
        const std::string& singleid
    );

    void writeBodyLine(
        const std::string& sequenceid,
        const std::string& binid,
        const std::string& taxid
    );
    
private:
    void writeHeader(const std::string& sampleid, const std::string& taxonomyid);
    
    void writeHeaderCustom(const std::vector<std::tuple<std::string, std::string>>& custom_header_tags);
    
    void writeHeaderColumnTags(const ColumnTags cols);
    
    void writeHeaderColumnTagsCustom(const std::vector<std::string>& custom_column_tags);
    
    std::ostream& ostr_;
    const ColumnTags cols_;
    const std::string format_version_ = "0.9.0";
    const std::string custom_tag_prefix_ = "TaxatorTK";
};

#endif // bioboxes_hh_
