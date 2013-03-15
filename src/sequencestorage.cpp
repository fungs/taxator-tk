#include <assert.h>
#include "sequencestorage.hh"
#include "ncbidata.hh"


RandomMMSeqStorRO::RandomMMSeqStorRO( const std::string& filename ) {
	seqan::MultiSeqFile db_sequences;
	seqan::open( db_sequences.concat, filename.c_str(), seqan::OPEN_RDONLY );
	seqan::split( db_sequences, seqan::Fasta() );
	medium_unsigned_int num_records = seqan::length( db_sequences );
	std::cerr << "indexing DB sequence positions (total=" << num_records << ")..."; 
	for( medium_unsigned_int i = 0; i < num_records; ++i ) {
		StringType seq;
		std::cerr << "assigning sequence (pos=" << i << ")... ";
		seqan::assignSeq( seq, db_sequences[i], fasta_format_ );
		std::cerr << " done" << std::endl;
		
		std::string id;
		seqan::assignSeqId( id, db_sequences[i], fasta_format_ );
		std::cerr << "id is: " << id << std::endl;
		id2pos_[ extractFastaCommentField( id, "gi" ) ] = seqan::assignValueById( data_, seq );
		std::cerr << "seqan assign value by id successful" << std::endl;
	}
	std::cerr << " done" << std::endl;
}



const RandomMMSeqStorRO::StringType& RandomMMSeqStorRO::getSequence ( const std::string& id ) const {
	std::map< std::string, large_unsigned_int >::const_iterator find_it = id2pos_.find( id );
	if( find_it != id2pos_.end() ) {
		return seqan::value( data_, find_it->second );
// 		seqan::assignSeq( seq, data_[find_it->second], fasta_format_ );
	}
	return StringType(); //TODO: exception instead, invalid reference
}



// const seqan::Dna5String RandomSeqStorRO::getSequence ( const std::string& id, unsigned int start, unsigned int stop ) const {
// 
// }



SequentialSeqStorRO::SequentialSeqStorRO ( const std::string& filename ) {

}




const SequentialSeqStorRO::StringType& SequentialSeqStorRO::getSequence ( const std::string& id ) const {
	return StringType(); //TODO: exception instead, invalid reference
}



// const seqan::Dna5String SequentialSeqStorRO::getSequence ( const std::string& id, unsigned int start, unsigned int stop ) const {
// 	return seqan::Dna5String();
// }
