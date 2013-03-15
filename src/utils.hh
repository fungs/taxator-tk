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



//TODO: clean up template functions


inline bool ignoreLine( const std::string& line ) {
	return line.empty() || line[0] == default_comment_symbol;
}



inline bool maskedLine( const std::string& line ) {
	return ! line.empty() && line[0] == default_mask_symbol;
}



// easy tokenizer without having to use boost
template < class ContainerT >
void tokenizeSingleCharDelim(const std::string& str, ContainerT& tokens, const std::string& delimiters = " ", int fieldnum = 0, const bool trimempty = false) {

	const unsigned int stringlength = str.size();

	if( ! fieldnum ) { // in case not provided or 0 given
		fieldnum = stringlength;
	}

	std::string::size_type pos, lastpos = 0;
	while( fieldnum && lastpos < stringlength ) {
		pos = str.find_first_of(delimiters, lastpos);
		if( pos == std::string::npos ) {
			pos = str.length();

			if( pos != lastpos || !trimempty ) {
				tokens.push_back( typename ContainerT::value_type( str.data() + lastpos, (typename ContainerT::value_type::size_type)pos - lastpos ) );
			}
			lastpos = pos;
			break;
		} else {
			if( pos != lastpos || !trimempty ) {
				tokens.push_back( typename ContainerT::value_type( str.data() + lastpos, (typename ContainerT::value_type::size_type)pos - lastpos ) );
				--fieldnum;
			}
		}
		lastpos = pos + 1;
	}
 	tokens.push_back( typename ContainerT::value_type( str.data() + lastpos, (typename ContainerT::value_type::size_type)stringlength - lastpos ) ); //append rest
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
			break;
		} else {
			if( pos != lastpos || !trimempty ) {
				tokens.push_back( typename ContainerT::value_type( str.data() + lastpos, (typename ContainerT::value_type::size_type)pos - lastpos ) );
				--fieldnum;
			}
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
		if( ! ignoreLine( line ) ) {
			tokenizeSingleCharDelim( line, fields, SEP, 2 );
			field_it = fields.begin();
			try {
				key = boost::lexical_cast< KeyT >( *field_it++ );
				value = boost::lexical_cast< ValueT >( *field_it );
				map_fill.insert( std::make_pair( key, value ) );
			} catch( boost::bad_lexical_cast e ) {
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

#endif // utils_hh_
