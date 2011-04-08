#ifndef utils_hh_
#define utils_hh_

#include <iostream>
#include <map>
#include <string>
#include <list>
#include <fstream>
#include <boost/lexical_cast.hpp>

//TODO: clean up template functions


inline bool ignoreLine( const std::string& line ) {
	return line.empty() || line[0] == '#';
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
};



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
};



template< typename KeyT, typename ValueT >
void loadMapFromFile( const std::string& filename, std::map< KeyT, ValueT >& map_fill, const std::string& SEP = "\t" ) {
	std::string line;
	std::list< std::string > fields;
	std::list< std::string >::iterator field_it;
	std::ifstream file_handle( filename.c_str() );

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
};

#endif // utils_hh_
