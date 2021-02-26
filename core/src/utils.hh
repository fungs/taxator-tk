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

#ifndef utils_hh_
#define utils_hh_

#include "constants.hh"
#include "types.hh"
#include <boost/lexical_cast.hpp>
#include <boost/tuple/tuple.hpp>
#include <iostream>
#include <map>
#include <string>
#include <list>
#include <fstream>
#include <limits>
#include <assert.h>



//TODO: clean up template functions


inline bool emptyLine( const std::string& line ) {
  return line.empty();
}


inline bool ignoreLine( const std::string& line ) {
  return !emptyLine(line) && line[0] == default_comment_symbol;
}


inline bool maskedLine( const std::string& line ) {
	return ! line.empty() && line[0] == default_mask_symbol;
}



// easy tokenizer without having to use boost
template < class ContainerT >
void tokenizeSingleCharDelim(const std::string& str, ContainerT& tokens, const std::string& delimiters = " ", int fieldnum = 0, const bool trimempty = false) {
  const unsigned int stringlength = str.size();
  if (! fieldnum) fieldnum = stringlength;  // in case not provided or 0 given
  std::string::size_type pos, lastpos = 0;
  while ( fieldnum && lastpos < stringlength ) {
    pos = str.find_first_of(delimiters, lastpos);
    if ( pos == std::string::npos ) {
      pos = str.length();
      if ( pos != lastpos || !trimempty ) tokens.push_back( typename ContainerT::value_type( str.data() + lastpos, (typename ContainerT::value_type::size_type)pos - lastpos ) );
      lastpos = pos;
			break;
    }
    if ( pos != lastpos || !trimempty ) {
      tokens.push_back( typename ContainerT::value_type( str.data() + lastpos, (typename ContainerT::value_type::size_type)pos - lastpos ) );
      --fieldnum;
    }
    lastpos = pos + 1;
  }
  tokens.push_back( typename ContainerT::value_type( str.data() + lastpos, (typename ContainerT::value_type::size_type)stringlength - lastpos ) );  //append rest
}



template < class ContainerT >
void tokenizeMultiCharDelim(const std::string& str, ContainerT& tokens, const std::string& delimiters = " ", int fieldnum = 0, const bool trimempty = false) {
  const unsigned int stringlength = str.size();

  if( ! fieldnum ) { // in case not provided or 0 given
    fieldnum = stringlength;
  }

  const int delimsize( delimiters.size() );

  std::string::size_type pos, lastpos = 0;
  while( fieldnum && lastpos < stringlength ) {
    pos = str.find( delimiters, lastpos );
    if( pos == std::string::npos ) {
      pos = str.length();

      if( pos != lastpos || !trimempty ) {
        tokens.push_back( typename ContainerT::value_type( str.data() + lastpos, (typename ContainerT::value_type::size_type)pos - lastpos ) );
      }
      lastpos = pos;
      return;
    }
    if( pos != lastpos || !trimempty ) {
      tokens.push_back( typename ContainerT::value_type( str.data() + lastpos, (typename ContainerT::value_type::size_type)pos - lastpos ) );
      --fieldnum;
    }
    lastpos = pos + delimsize;
  }
  tokens.push_back( typename ContainerT::value_type( str.data() + lastpos, (typename ContainerT::value_type::size_type)stringlength - lastpos ) ); //append rest
}



template< typename KeyT, typename ValueT >
void loadMapFromFile( const std::string& filename, std::map< KeyT, ValueT >& map_fill, const std::string& SEP = "\t" ) {
	std::string line;
	std::list< std::string > fields;
	std::list< std::string >::iterator field_it;
	std::ifstream file_handle( filename.c_str() );

	if( ! file_handle.good() ) {
		std::cerr << "something wrong with file " << filename << std::endl;
		file_handle.close();
		return;
	}

	KeyT key;
	ValueT value;
	while( std::getline( file_handle, line ) ) {
		if( ! (emptyLine(line) || ignoreLine(line)) ) {
			tokenizeSingleCharDelim( line, fields, SEP, 2 );
			field_it = fields.begin();
			try {
				key = boost::lexical_cast< KeyT >( *field_it++ );
				value = boost::lexical_cast< ValueT >( *field_it );
				map_fill.insert( std::make_pair( key, value ) );
			} catch( boost::bad_lexical_cast& ) {
				std::cerr << "loadMapFromFile(): could not parse this line: '" << line << "'" << std::endl;
			}
			fields.clear();
		}
	}
	file_handle.close();
}



template< typename ContainerT >
typename ContainerT::iterator firstUnmaskedIter( ContainerT& cont ) {
	typename ContainerT::iterator rec_it = cont.begin();
	while ( rec_it != cont.end() && (*rec_it)->isFiltered() ) { ++rec_it; };
	return rec_it;
}


template< typename TType, unsigned int position >
struct compareTupleFirstLT {
	bool operator() ( const TType& a, const TType& b ) {
		return boost::get<position>( a ) < boost::get<position>( b );
	}
};



// provides memory efficient storage of query sequences via reference counting
template< typename UIntType = small_unsigned_int > //default maximum 256 references per unique string!
class ReferencedStringStore {
	public:
		ReferencedStringStore() : last_( id2counter_.end() ), max_capacity_( std::numeric_limits<UIntType>::max() ) {}

		const std::string& add( const std::string& id ) {
			if ( &id != &last_->first ) {
				last_ = id2counter_.find( id );
				if ( last_ == id2counter_.end() ) last_ = id2counter_.insert( std::make_pair( id, 0 ) ).first;
			}

			assert( last_->second < max_capacity_ );
			last_->second += 1;
			return last_->first;
		}

		bool remove( const std::string& id ) {
			if ( &id != &last_->first ) {
				typename std::map< const std::string, UIntType >::iterator it = id2counter_.find( id );
				if ( it == id2counter_.end() || &id != &it->first ) return false;
				last_ = it;
			}

			if ( last_->second > 0 ) {
				last_->second -= 1;
			} else {
				id2counter_.erase( last_ );
				last_ = id2counter_.end();
			}
			return true;
		}

	private:
		typename std::map< const std::string, UIntType > id2counter_; //max. 256 references
		typename std::map< const std::string, UIntType >::iterator last_; //save some time
		const UIntType max_capacity_;
};

// ----------------------------------------------------------------------------
// Function write()
// ----------------------------------------------------------------------------

/*!
 * @fn Align#write
 * @deprecated Old-style I/O.
 * @brief Writing of Gaps to Streams in human-readable format.
 *
 * @signature void write(stream, align);
 *
 * @param[in,out] stream The Stream to write to.
 * @param[in]     align  The Align object to write out.
 */

// simplified seqan alignment print/write without without breaking it
template <typename TFile, typename TSource, typename TSpec>
inline void
_write(TFile & target,
      Align<TSource, TSpec> const & source)
{
    typedef Align<TSource, TSpec> const TAlign;
    typedef typename Row<TAlign>::Type TRow;
    typedef typename Position<typename Rows<TAlign>::Type>::Type TRowsPosition;
    typedef typename Position<TAlign>::Type TPosition;

    TRowsPosition row_count = length(rows(source));
    TPosition begin_ = 0;
    TPosition end_ = std::min(length(row(source, 0)), length(row(source, 1)));

    unsigned int baseCount = 0;
    unsigned int leftSpace = 6;
    while (begin_ < end_)
    {
        unsigned int windowSize_ = 50;
        if ((begin_ + windowSize_) > end_)
            windowSize_ = end_ - begin_;

        // Print header line
        char buffer[20];
        int len = snprintf(buffer, 20, "%7u", (unsigned)baseCount);
        write(target, buffer, len);
        baseCount += windowSize_;
        writeValue(target, ' ');
        for (TPosition i = 1; i <= windowSize_; ++i)
        {
            if ((i % 10) == 0)
                writeValue(target, ':');
            else if ((i % 5) == 0)
                writeValue(target, '.');
            else
                writeValue(target, ' ');
        }
        writeValue(target, ' ');
        writeValue(target, '\n');

        // Print sequences
        for (TRowsPosition i = 0; i < 2 * row_count - 1; ++i)
        {
            for (unsigned int j = 0; j < leftSpace + 2; ++j)
                writeValue(target, ' ');
            if ((i % 2) == 0)
            {
                TRow & row_ = row(source, i / 2);
                typedef typename Iterator<typename Row<TAlign>::Type const, Standard>::Type TIter;
                TIter begin1_ = iter(row_, begin_);
                TIter end1_ = iter(row_, begin_ + windowSize_);
                for (; begin1_ != end1_; ++begin1_)
                {
                    if (isGap(begin1_))
                        writeValue(target, gapValue<char>());
                    else
                        writeValue(target, getValue(begin1_));
                }
            }
            else
            {
                for (unsigned int j = 0; j < windowSize_; ++j)
                {
                    if ((!isGap(row(source, (i - 1) / 2), begin_ + j)) &&
                        (!isGap(row(source, (i + 1) / 2), begin_ + j)) &&
                        (row(source, (i - 1) / 2)[begin_ + j] == row(source, (i + 1) / 2)[begin_ + j]))
                    {
                        writeValue(target, '|');
                    }
                    else
                    {
                        writeValue(target, ' ');
                    }
                }
            }
            writeValue(target, '\n');
        }
        writeValue(target, '\n');
        begin_ += 50;
    }
    writeValue(target, '\n');
}

#endif // utils_hh_
