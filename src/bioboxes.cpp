#include "bioboxes.hh"
#include "constants.hh"

BioboxesBinningFormat::BioboxesBinningFormat(
    const BioboxesBinningFormat::ColumnTags cols,  // TODO: make template argument
    const std::string& sampleid,
    const std::string& taxonomyid,
    std::ostream& ostr,
    const std::string& custom_tag_prefix,
    const std::vector<std::tuple<const std::string, const std::string>> custom_header_tags,
    const std::vector<std::string> custom_column_tags)
        : ostr_(ostr), cols_(cols), custom_tag_prefix_(custom_tag_prefix)
{
    // TODO: check values with regexp
    
    ostr_ << "# This is the bioboxes.org binning output format at" << endline
          << "# https://github.com/bioboxes/rfc/tree/master/data-format"
          << endline << endline;
          
    writeHeader (sampleid, taxonomyid);
    writeHeaderCustom(custom_header_tags);
    ostr_ << endline;
    
    writeHeaderColumnTags(cols);
    writeHeaderColumnTagsCustom(custom_column_tags);
    ostr_ << endline << std::flush;
}


BioboxesBinningFormat::~BioboxesBinningFormat()
{
    ostr_ << std::flush;
}


void BioboxesBinningFormat::writeHeader(const std::string& sampleid, const std::string& taxonomyid)
{
    // output header tags
    ostr_ << "@Version:" << format_version_ << endline;
    ostr_ << "@SampleID:" << sampleid << endline;
    if(!taxonomyid.empty()) ostr_ << "@TaxonomyID:" << taxonomyid << endline;
}


void BioboxesBinningFormat::writeHeaderCustom(const std::vector<std::tuple<const std::string, const std::string>>& custom_header_tags)
{
    for(auto it = custom_header_tags.begin(); it != custom_header_tags.end(); ++it) {
        ostr_ << "@_" << custom_tag_prefix_ << '_' << std::get<0>(*it) << ':' << std::get<1>(*it) << endline;
    }
}


void BioboxesBinningFormat::writeHeaderColumnTags(const BioboxesBinningFormat::ColumnTags cols)
{
    ostr_ << "@@SequenceID";
    if(cols == ColumnTags::taxid) ostr_ << tab << "TaxID";
    else if(cols == ColumnTags::binid) ostr_ << tab << "BinID";
    else if(cols == ColumnTags::taxid_binid) ostr_ << tab << "TaxID" << tab << "BinID";
}


void BioboxesBinningFormat::writeHeaderColumnTagsCustom(const std::vector< std::string >& custom_column_tags)
{
    for(auto it = custom_column_tags.begin(); it != custom_column_tags.end(); ++it) {
        ostr_ << tab << '_' << custom_tag_prefix_ << '_' << *it;
    }
}


void BioboxesBinningFormat::writeBodyLine(const std::string& sequenceid, const std::string& singleid)
{
    ostr_ << sequenceid << tab << singleid << endline;
}


void BioboxesBinningFormat::writeBodyLine(const std::string& sequenceid, const std::string& singleid, const std::vector< std::string >& columns_custom)
{
    ostr_ << sequenceid << tab << singleid;
    for(auto it = columns_custom.begin(); it != columns_custom.end(); ++it) {
        ostr_ << tab << *it;
    }
    ostr_ << endline;
}



void BioboxesBinningFormat::writeBodyLine(const std::string& sequenceid, const std::string& binid, const std::string& taxid)
{
    ostr_ << sequenceid << tab << binid << tab << taxid << endline;
}


void BioboxesBinningFormat::writeBodyLine(const std::string& sequenceid, const std::string& binid, const std::string& taxid, const std::vector< std::string >& columns_custom)
{
    ostr_ << sequenceid << tab << binid << tab << taxid;
    for(auto it = columns_custom.begin(); it != columns_custom.end(); ++it) {
        ostr_ << tab << *it;
    }
    ostr_ << endline;
}
