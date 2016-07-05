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
#include <seqan/seq_io.h>
#include <boost/progress.hpp>
#include <boost/concept_check.hpp>
#include <boost/filesystem.hpp>
#include <set>
#include <string>
#include "ncbidata.hh"
#include <assert.h>
#include "exception.hh"


// This currently works with standard and packed strings
template <typename WorkingStringType>
class RandomSeqStoreROInterface {
public:
    virtual const WorkingStringType getSequence ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const = 0;
    virtual const WorkingStringType getSequenceReverseComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const = 0;
    virtual ~RandomSeqStoreROInterface() {};
    
    const WorkingStringType getSequenceAuto ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
      if ( start < stop ) return getSequence( id, start, stop );
      return getSequenceReverseComplement( id, stop, start );
    };
};


template < typename StorageStringType = seqan::Dna5String, typename WorkingStringType = seqan::Dna5String, typename Format = seqan::Fasta >
class RandomInmemorySeqStoreRO : public RandomSeqStoreROInterface<WorkingStringType> {
public:
    RandomInmemorySeqStoreRO ( const std::string& filename ) : format_( Format() ) {
        
        if( ! boost::filesystem::exists( filename ) ) BOOST_THROW_EXCEPTION(FileNotFound{} << file_info{filename});
        
        std::cerr << "Loading '" << filename;
        seqan::SeqFileIn db_sequences(filename.c_str());

        seqan::readRecords(ids_, seqs_, db_sequences);
        large_unsigned_int num_records = seqan::length( ids_ );
        
        std::cerr  << "' (total=" << num_records << ")" << std::endl;
        
        auto id = seqan::begin(ids_);
        
        for(large_unsigned_int i = 0; i < num_records; ++i){
            id2pos_[ *id ] = i;
            seqan::goNext(id);
        }
        
        std::cerr << "fasta-file loaded" << std::endl;
    }

    RandomInmemorySeqStoreRO ( const std::string& filename, const std::set< std::string >& whitelist ) : format_( Format() ) {
        
        if( ! boost::filesystem::exists( filename ) ) BOOST_THROW_EXCEPTION(FileNotFound{} << file_info{filename});
        
        std::cerr << "Loading '" << filename;
        seqan::SeqFileIn db_sequences(filename.c_str());
        
        seqan::readRecords(ids_, seqs_, db_sequences);
        large_unsigned_int num_records = seqan::length( ids_ );
        
        auto id = seqan::begin(ids_);
        auto seq = seqan::begin(seqs_);
        
        //only seqs, make index of ids
        
        for(large_unsigned_int i = 0; i < num_records; ++i){
            if ( whitelist.count( std::string::c_str(*id) ) ) {
                id2pos_[ *id ] = seqan::assignValueById( data_, *seq );
            }
            id++;
            seq++;
        }
    }
    
    const StorageStringType& getSequence ( const std::string& id ) const {
        seqan::CharString id_ss = id;
        std::map< seqan::CharString, large_unsigned_int >::const_iterator find_it = id2pos_.find( id_ss );
        if( find_it == id2pos_.end() ) BOOST_THROW_EXCEPTION(SequenceNotFound {} << seqid_info{id});
        return(seqan::value( seqs_, find_it->second ));
        
    };

    const WorkingStringType getSequence ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
        const StorageStringType& db_seq = getSequence ( id );
        stop = std::min< large_unsigned_int >( stop, seqan::length( db_seq ) );
        if( start > seqan::length( db_seq ) ) BOOST_THROW_EXCEPTION(SequenceRangeError{} << general_info{"invalid position"} << seqid_info{id} << position_info{start});
        WorkingStringType seq = seqan::infix ( db_seq, start - 1, stop );
        assert( seqan::length( seq ) == (stop - start + 1) );
        return seq;
    };

    const WorkingStringType getSequenceReverseComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
        const StorageStringType& db_seq = getSequence ( id );
//        stop = std::min< large_unsigned_int >( stop, seqan::length( db_seq ) );
//        if( start > seqan::length( db_seq ) ) BOOST_THROW_EXCEPTION(SequenceRangeError{} << general_info{"invalid position"} << seqid_info{id} << position_info{start});
//        WorkingStringType cst = seqan::infix ( db_seq, start - 1, stop );
//        seqan::ModifiedString< seqan::ModifiedString< WorkingStringType, seqan::ModView< seqan::FunctorComplement< seqan::Dna > > >, seqan::ModReverse> seq( cst );
//        assert( seqan::length( seq ) == (stop - start + 1) );
        return db_seq;
    };
    
    
protected:
    seqan::StringSet< StorageStringType > data_;
    std::map< seqan::CharString, large_unsigned_int > id2pos_;
    const StorageStringType empty_string_;
    Format format_;
    seqan::StringSet<seqan::CharString> ids_;
    seqan::StringSet<WorkingStringType> seqs_;
};


// AA Template Specialisation --------------------------------------------------------------------------------

template < typename StorageStringType, typename Format >
class RandomInmemorySeqStoreRO<StorageStringType, seqan::AminoAcid, Format> : public RandomSeqStoreROInterface<seqan::AminoAcid> {
public:
    RandomInmemorySeqStoreRO ( const std::string& filename ) : format_( Format() ) {
        
        if( ! boost::filesystem::exists( filename ) ) BOOST_THROW_EXCEPTION(FileNotFound{} << file_info{filename});
        
        std::cerr << "Loading '" << filename;
        seqan::SeqFileIn db_sequences(filename.c_str());

        seqan::readRecords(ids_, seqs_, db_sequences);
        large_unsigned_int num_records = seqan::length( ids_ );
        
        std::cerr  << "' (total=" << num_records << ")" << std::endl;
        
        auto id = seqan::begin(ids_);
        
        for(large_unsigned_int i = 0; i < num_records; ++i){
            id2pos_[ *id ] = i;
            seqan::goNext(id);
        }
        
        std::cerr << "fasta-file loaded" << std::endl;
    }

    RandomInmemorySeqStoreRO ( const std::string& filename, const std::set< std::string >& whitelist ) : format_( Format() ) {
        
        if( ! boost::filesystem::exists( filename ) ) BOOST_THROW_EXCEPTION(FileNotFound{} << file_info{filename});
        
        std::cerr << "Loading '" << filename;
        seqan::SeqFileIn db_sequences(filename.c_str());
        
        seqan::readRecords(ids_, seqs_, db_sequences);
        large_unsigned_int num_records = seqan::length( ids_ );
        
        auto id = seqan::begin(ids_);
        auto seq = seqan::begin(seqs_);
        
        //only seqs, make index of ids
        
        for(large_unsigned_int i = 0; i < num_records; ++i){
            if ( whitelist.count( std::string::c_str(*id) ) ) {
                id2pos_[ *id ] = seqan::assignValueById( data_, *seq );
            }
            id++;
            seq++;
        }
    }
    
    const StorageStringType& getSequence ( const std::string& id ) const {
        seqan::CharString id_ss = id;
        std::map< seqan::CharString, large_unsigned_int >::const_iterator find_it = id2pos_.find( id_ss );
        if( find_it == id2pos_.end() ) BOOST_THROW_EXCEPTION(SequenceNotFound {} << seqid_info{id});
        return(seqan::value( seqs_, find_it->second ));
        
    };

    const seqan::AminoAcid getSequence ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
        const StorageStringType& db_seq = getSequence ( id );
        stop = std::min< large_unsigned_int >( stop, seqan::length( db_seq ) );
        if( start > seqan::length( db_seq ) ) BOOST_THROW_EXCEPTION(SequenceRangeError{} << general_info{"invalid position"} << seqid_info{id} << position_info{start});
        seqan::AminoAcid seq = seqan::infix ( db_seq, start - 1, stop );
        assert( seqan::length( seq ) == (stop - start + 1) );
        return seq;
    };

    const seqan::AminoAcid getSequenceReverseComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
        const StorageStringType& db_seq = getSequence ( id );
//        stop = std::min< large_unsigned_int >( stop, seqan::length( db_seq ) );
//        if( start > seqan::length( db_seq ) ) BOOST_THROW_EXCEPTION(SequenceRangeError{} << general_info{"invalid position"} << seqid_info{id} << position_info{start});
//        seqan::AminoAcid cst = seqan::infix ( db_seq, start - 1, stop );
//        seqan::ModifiedString< seqan::ModifiedString< seqan::AminoAcid, seqan::ModView< seqan::FunctorComplement< seqan::Dna > > >, seqan::ModReverse> seq( cst );
//        assert( seqan::length( seq ) == (stop - start + 1) );
        return db_seq;
    };
    
    
protected:
    seqan::StringSet< StorageStringType > data_;
    std::map< seqan::CharString, large_unsigned_int > id2pos_;
    const StorageStringType empty_string_;
    Format format_;
    seqan::StringSet<seqan::CharString> ids_;
    seqan::StringSet<seqan::AminoAcid> seqs_;
};


template< typename StringType, bool skip = true, typename Format = seqan::Fasta >
class SequentialSeqStoreRO : public RandomSeqStoreROInterface<StringType> {
public:
    SequentialSeqStoreRO ( const std::string& filename ) : strm_( filename.c_str() ), format_( Format() ) {};

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
        
        if(last_id_ != id) BOOST_THROW_EXCEPTION(SequenceNotFound {} << seqid_info{id});

        return last_entry_;
    };

    const StringType getSequence ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) {
        return seqan::infix ( getSequence ( id ), start - 1, stop );
    };

    const StringType getSequenceReverseComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) {
        return seqan::ModifiedString< seqan::ModifiedString< StringType, seqan::ModView< seqan::FunctorComplement< seqan::Dna > > >, seqan::ModReverse> ( seqan::infix ( getSequence ( id ), start, stop ) );
    };

protected:
    std::ifstream strm_;
    Format format_;
    StringType last_entry_;
    std::string last_id_;
};

// AA Template Specialisation --------------------------------------------------------------------------------

template< bool skip , typename Format >
class SequentialSeqStoreRO< seqan::AminoAcid, skip, Format> : public RandomSeqStoreROInterface<seqan::AminoAcid> {
public:
    SequentialSeqStoreRO ( const std::string& filename ) : strm_( filename.c_str() ), format_( Format() ) {};

    const seqan::AminoAcid& getSequence ( const std::string& id ) {

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
        
        if(last_id_ != id) BOOST_THROW_EXCEPTION(SequenceNotFound {} << seqid_info{id});

        return last_entry_;
    };

    const seqan::AminoAcid getSequence ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) {
        return seqan::infix ( getSequence ( id ), start - 1, stop );
    };

    const seqan::AminoAcid getSequenceReverseComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) {
        //return seqan::ModifiedString< seqan::ModifiedString< seqan::AminoAcid, seqan::ModView< seqan::FunctorComplement< seqan::Dna > > >, seqan::ModReverse> ( seqan::infix ( getSequence ( id ), start, stop ) );
        return seqan::infix ( getSequence ( id ), start, stop );
    };

protected:
    std::ifstream strm_;
    Format format_;
    seqan::AminoAcid last_entry_;
    std::string last_id_;
};


template< typename StringType >
class RandomIndexedSeqstoreRO : public RandomSeqStoreROInterface<StringType> {
public:
    RandomIndexedSeqstoreRO( const std::string& fasta_filename, const std::string& index_filename ) : index_filename_( index_filename ), write_on_exit_( false ) {
        if ( ! boost::filesystem::exists( index_filename ) )  {
            if ( seqan::build( index_, fasta_filename.c_str() ) ) { //TODO: propagate error
                BOOST_THROW_EXCEPTION(GeneralError{} << general_info{"could not build fasta index"} << file_info{index_filename});
                return;
            } else write_on_exit_ = true;
        } else if ( ! seqan::open( index_, fasta_filename.c_str(), index_filename.c_str() ) ) {
            BOOST_THROW_EXCEPTION(FileError{} << file_info{index_filename});
            return;
        }

        //make a thread-safe lookup for identifiers, broken in SEQAN as of version 1.4.1
        unsigned int idx = 0;
        for (auto it = seqan::begin(index_.seqNameStore); !seqan::atEnd(it); seqan::goNext(it)) {
            refid2position_[*it] = idx++;
        }
    }

    const StringType getSequence ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
        assert( start <= stop );
        unsigned int seq_num;
        StringType seq;
        
        /*if ( ! seqan::getIdByName( index_, id.c_str(), seq_num ) ) {
        	std::cerr << "Sequence " << id << " not found in sequence file." << std::endl; //TODO. propagate error
        	return seq;
        }*/
        std::map<seqan::CharString, unsigned int>::const_iterator it = refid2position_.find( id.c_str() );
        if( it != refid2position_.end() ) seq_num = it->second;
        else {
            BOOST_THROW_EXCEPTION(SequenceNotFound {} << seqid_info{id});
        }
        
        stop = std::min< large_unsigned_int >( stop, seqan::sequenceLength( index_, seq_num) );
        seqan::readRegion( seq, index_, seq_num, start - 1, stop );
        assert( seqan::length( seq ) == (stop - start + 1) );
        return seq;
    }
   
    const StringType getSequenceReverseComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
        assert( start <= stop );
        StringType seq = getSequence( id , start, stop );
        seqan::reverseComplement( seq );
        return seq;
    }
 
//    const StringType getSequenceReverseComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
//        assert( start <= stop );
//        StringType seq = getSequence( id , start, stop );
//        //seqan::reverseComplement( seq );
//        return seq;
//    }
    
//    template<typename seqan::String<seqan::Dna5>>
//    const StringType getSequenceReverseComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
//        assert( start <= stop );
//        StringType seq = getSequence( id , start, stop );
//        seqan::reverseComplement( seq );
//        return seq;
//    }
//    
//    template<typename seqan::String<seqan::AminoAcid>>
//    const StringType getSequenceReverseComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
//        assert( start <= stop );
//        StringType seq = getSequence( id , start, stop );
//        //seqan::reverseComplement( seq );
//        return seq;
//    }
    
    
    ~RandomIndexedSeqstoreRO() {
        if ( write_on_exit_ && ! boost::filesystem::exists( index_filename_ ) )
            if( seqan::save( index_, index_filename_.c_str() ) ) BOOST_THROW_EXCEPTION(FileError{} << file_info{index_filename_});
    }

protected:
    const std::string index_filename_;
    seqan::FaiIndex index_;
    bool write_on_exit_;
    std::map<seqan::CharString, unsigned int> refid2position_;
};

// AA Template Specialisation --------------------------------------------------------------------------------

template<> class RandomIndexedSeqstoreRO <seqan::String<seqan::AminoAcid>> : public RandomSeqStoreROInterface<seqan::String<seqan::AminoAcid>>{
public:
        RandomIndexedSeqstoreRO( const std::string& fasta_filename, const std::string& index_filename ) : index_filename_( index_filename ), write_on_exit_( false ) {
        if ( ! boost::filesystem::exists( index_filename ) )  {
            if ( seqan::build( index_, fasta_filename.c_str() ) ) { //TODO: propagate error
                BOOST_THROW_EXCEPTION(GeneralError{} << general_info{"could not build fasta index"} << file_info{index_filename});
                return;
            } else write_on_exit_ = true;
        } else if ( ! seqan::open( index_, fasta_filename.c_str(), index_filename.c_str() ) ) {
            BOOST_THROW_EXCEPTION(FileError{} << file_info{index_filename});
            return;
        }

        //make a thread-safe lookup for identifiers, broken in SEQAN as of version 1.4.1
        unsigned int idx = 0;
        for (auto it = seqan::begin(index_.seqNameStore); !seqan::atEnd(it); seqan::goNext(it)) {
            refid2position_[*it] = idx++;
        }
        }

        const seqan::String<seqan::AminoAcid> getSequence ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
        assert( start <= stop );
        unsigned int seq_num;
        seqan::String<seqan::AminoAcid> seq;
        
        /*if ( ! seqan::getIdByName( index_, id.c_str(), seq_num ) ) {
        	std::cerr << "Sequence " << id << " not found in sequence file." << std::endl; //TODO. propagate error
        	return seq;
        }*/
        std::map<seqan::CharString, unsigned int>::const_iterator it = refid2position_.find( id.c_str() );
        if( it != refid2position_.end() ) seq_num = it->second;
        else {
            BOOST_THROW_EXCEPTION(SequenceNotFound {} << seqid_info{id});
        }
        
        stop = std::min< large_unsigned_int >( stop, seqan::sequenceLength( index_, seq_num) );
        seqan::readRegion( seq, index_, seq_num, start - 1, stop );
        assert( seqan::length( seq ) == (stop - start + 1) );
        return seq;
        }
    
    
    const seqan::String<seqan::AminoAcid> getSequenceReverseComplement ( const std::string& id, large_unsigned_int start, large_unsigned_int stop ) const {
        assert( start <= stop );
        seqan::String<seqan::AminoAcid> seq = getSequence( id , start, stop );
        //seqan::reverseComplement( seq );
        return seq;
    }
    
    protected:
    const std::string index_filename_;
    seqan::FaiIndex index_;
    bool write_on_exit_;
    std::map<seqan::CharString, unsigned int> refid2position_;

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
