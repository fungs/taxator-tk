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

#ifndef sequencestorage_hh_
#define sequencestorage_hh_

#include <seqan/sequence.h>
#include <seqan/file.h>
#include <boost/progress.hpp>
#include <set>
#include <string>
#include "ncbidata.hh"
#include <assert.h>


// This currently works with standard and packed strings
template < typename StorageStringType = seqan::Dna5String, typename CopyStringType = seqan::Dna5String, typename Format = seqan::Fasta >
class RandomSeqStorRO {
	public:
		RandomSeqStorRO ( const std::string& filename ) : format_( Format() ) {
			std::cerr << "analyzing file '" << filename << "'... ";
			seqan::MultiSeqFile db_sequences;
			if ( seqan::open( db_sequences.concat, filename.c_str(), seqan::OPEN_RDONLY ) ) {
				seqan::split( db_sequences, format_ );
				std::cerr << "done" << std::endl;
				large_unsigned_int num_records = seqan::length( db_sequences );
				std::cerr << "importing sequences from '" << filename << "' (total=" << num_records << ")" << std::endl;
				{
					boost::progress_display eta( num_records - 1, std::cerr ); //progress bar
					for( large_unsigned_int i = 0; i < num_records; ++i ) {
						StorageStringType seq;
						seqan::assignSeq( seq, db_sequences[i], format_ );
						std::string id;
						seqan::assignSeqId( id, db_sequences[i], format_ );
						id2pos_[ id ] = seqan::assignValueById( data_, seq );
						++eta;
					}
				}
				std::cerr << std::endl;
			} else {
				std::cerr << "Could not open input FASTA file \"" << filename << "\"" << std::endl;
			}
		}

		RandomSeqStorRO ( const std::string& filename, const std::set< std::string >& whitelist ) : format_( Format() ) {
			std::cerr << "analyzing file '" << filename << "'... ";
			seqan::MultiSeqFile db_sequences;
			if ( seqan::open( db_sequences.concat, filename.c_str(), seqan::OPEN_RDONLY ) ) {
				seqan::split( db_sequences, format_ );
				std::cerr << "done" << std::endl;
				large_unsigned_int num_records = seqan::length( db_sequences );
				large_unsigned_int effective_num_records = std::min< large_unsigned_int >( num_records, whitelist.size() );
				std::cerr << "importing sequences from '" << filename << "' (total=" << effective_num_records << ")" << std::endl;
				{
					boost::progress_display eta( effective_num_records - 1, std::cerr ); //progress bar
					for( large_unsigned_int i = 0; i < num_records; ++i ) {
						StorageStringType seq;
						seqan::assignSeq( seq, db_sequences[i], format_ );
						std::string id;
						seqan::assignSeqId( id, db_sequences[i], format_ );
						
						if ( whitelist.count( id ) ) {
							id2pos_[ id ] = seqan::assignValueById( data_, seq );
							++eta;
						}
					}
					assert( seqan::length( data_ ) <= effective_num_records );
				}
				std::cerr << std::endl;
			} else {
				std::cerr << "Could not open input FASTA file \"" << filename << "\"" << std::endl;
			}
		}
		
		const StorageStringType& getSequence ( const std::string& id ) const {
			std::map< std::string, large_unsigned_int >::const_iterator find_it = id2pos_.find( id );
			if( find_it != id2pos_.end() ) {
				return seqan::value( data_, find_it->second );
			}
			std::cerr << "Could not find sequence with id " << id << " in RandomSeqStorRO" << std::endl;
			return empty_string_;
		};
		
		const CopyStringType getSequence ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
			const StorageStringType& db_seq = getSequence ( id );
			if ( &db_seq != &empty_string_ ) {
				stop = std::min< large_unsigned_int >( stop, seqan::length( db_seq ) );
				CopyStringType seq = seqan::infix ( db_seq, start - 1, stop );
				assert( seqan::length( seq ) == (stop - start + 1) );
				return seq;
			}
			return empty_string_;
		};
		
		const CopyStringType getSequenceComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
			const StorageStringType& db_seq = getSequence ( id );
			if ( &db_seq != &empty_string_ ) {
				stop = std::min< large_unsigned_int >( stop, seqan::length( db_seq ) );
				CopyStringType seq = seqan::ModifiedString< seqan::ModifiedString< CopyStringType, seqan::ModView< seqan::FunctorComplement< seqan::Dna > > >, seqan::ModReverse> ( seqan::infix ( db_seq, start - 1, stop ) );
				assert( seqan::length( seq ) == (stop - start + 1) );
				return seq;
			}
			return empty_string_;
		};
		
		const CopyStringType getSequenceAuto ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
			if ( start < stop ) {
				return getSequence( id, start, stop );
			}
			return getSequenceComplement( id, stop, start );
		};

	protected:
		seqan::StringSet< StorageStringType > data_;
		std::map< std::string, large_unsigned_int > id2pos_; //hash_map aka unordered_map would be more apt
		const StorageStringType empty_string_;
		Format format_;
};



template< typename StringType, bool skip = true, typename Format = seqan::Fasta >
class SequentialSeqStorRO {
	public:
		SequentialSeqStorRO ( const std::string& filename ) : strm_( filename.c_str() ), format_( Format() ) {};
		
		const StringType& getSequence ( const std::string& id ) {
			
			if ( id == last_id_ ) return last_entry_;
			
			seqan::readMeta( strm_, last_id_, format_ );
			
			if( skip ) {
				while ( last_id_ != id && ! strm_.eof() ) {
					std::cerr << "skipping sequence id: " << last_id_ << std::endl;
					seqan::read( strm_, last_entry_, format_ ); //eat it (is there no better way in seqan?)
					seqan::readMeta( strm_, last_id_, format_ );
				};
			}
			
			seqan::read( strm_, last_entry_, format_ );
			assert( last_id_ == id );
			
			return last_entry_;
		};

		const StringType getSequence ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) {
			return seqan::infix ( getSequence ( id ), start - 1, stop );
		};
		
		const StringType getSequenceComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) {
			return seqan::ModifiedString< seqan::ModifiedString< StringType, seqan::ModView< seqan::FunctorComplement< seqan::Dna > > >, seqan::ModReverse> ( seqan::infix ( getSequence ( id ), start, stop ) );
		};
		
		const StringType getSequenceAuto ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) {
			if ( start < stop ) {
				return seqan::infix ( getSequence ( id ), start, stop );
			}
			return seqan::ModifiedString< seqan::ModifiedString< StringType, seqan::ModView< seqan::FunctorComplement< seqan::Dna > > >, seqan::ModReverse> ( seqan::infix ( getSequence ( id ), stop, start + 1 ) );
		};
		
	protected:
		std::ifstream strm_;
		Format format_;
		StringType last_entry_;
		std::string last_id_;
};



void populateIdentSet( std::set< std::string >& whitelist, const std::string& filename ) {
	std::ifstream flatfile( filename.c_str() );
	std::string line;
	while( std::getline( flatfile, line ) ) {
		whitelist.insert( line ); //supposing that newline characters are stripped (UNIX)
	}
	flatfile.close();
}



#endif // sequencestorage_hh_
