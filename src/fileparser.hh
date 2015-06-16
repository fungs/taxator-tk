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

#ifndef fileparser_hh_
#define fileparser_hh_

#include "exception.hh"
#include "utils.hh"


template< typename FactoryType >
class FileParser {
public:
    typedef typename FactoryType::value_type RecordType;
    
    FileParser( const std::string& filename, FactoryType& factory ) : filehandle_(filename.c_str()),
                                                                      handle_(filehandle_),
                                                                      factory_(factory) {
        feed();
    }
    
    FileParser( std::istream& strm, FactoryType& factory ) : handle_(strm),
                                                             factory_(factory) {
        feed();
    }

    RecordType* next() {
        try {
            RecordType* ret = factory_.create(line_);
            feed();
            return ret;
        }
        catch (Exception &e) {
            e << line_info{line_num_};
            BOOST_THROW_EXCEPTION(e);
        }
        return NULL;  // should never be reached
    }

    inline void destroy( const RecordType* rec ) const { factory_.destroy(rec); }
    inline bool eof() { return eof_; }

private:
    void feed() {
        while (std::getline(handle_, line_)) {
            ++line_num_;
            if (!ignoreLine(line_)) return;
        }
        eof_ = true;
    }
    
    std::ifstream filehandle_;
    std::istream& handle_;
    std::string line_;
    FactoryType& factory_;
    
    unsigned int line_num_ = 0;
    bool eof_ = false;
};


template< typename InType, typename FactoryType >  // TODO: not yet working!!
auto make_file_parser(InType& in, FactoryType& fac) -> FileParser<FactoryType> {
    return FileParser<FactoryType>(in, fac);
}

#endif  // fileparser_hh_
