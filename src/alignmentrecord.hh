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

#ifndef alignmentrecord_hh_
#define alignmentrecord_hh_

#include <fstream>
#include <string>
#include <vector>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include "types.hh"
#include "utils.hh"
#include "ncbidata.hh"
#include "accessconv.hh"
#include "taxonomyinterface.hh"



class AlignmentRecord {
	public:
		~AlignmentRecord() {};
		inline const std::string& getQueryIdentifier() const { return query_identifier_; };
		inline large_unsigned_int getQueryStart() { return query_start_; };
		inline large_unsigned_int getQueryStop() { return query_stop_; };
		inline large_unsigned_int getQueryLength() { return query_length_; };
		inline const std::string& getReferenceIdentifier() const { return reference_identifier_; };
		inline large_unsigned_int getReferenceStart() { return reference_start_; };
		inline large_unsigned_int getReferenceStop() { return reference_stop_; };
		inline float getScore() const { return score_; };
		inline double getEValue() const { return evalue_; };
		inline large_unsigned_int getIdentities() { return identities_; };
		inline large_unsigned_int getAlignmentLength() { return alignment_length_; };
		inline const std::string& getAlignmentCode() { return alignment_code_; };
		inline bool isFiltered() const { return blacklist_this_; };
		inline float getPID() const { return identities_/float( std::max( query_length_, alignment_length_ ) ); };
		
// 		inline void setReferenceIdentifier( const std::string& ident ) { reference_identifier_ = ident; };
// 		inline void setReferenceStart( unsigned int pos ) { reference_start_ = pos; };
// 		inline void setReferenceStop( unsigned int pos ) { reference_stop_ = pos; };
// 		inline void setQueryIdentifier( const std::string& ident ) { query_identifier_ = ident; };
// 		inline void setQueryStart( unsigned int pos ) { query_start_ = pos; };
// 		inline void setQueryStop( unsigned int pos ) { query_stop_ = pos; };
// 		inline void setPID( float pid ) { pid_ = pid; };
// 		inline void setScore( float score ) { score_ = score; };
// 		inline void setEValue( double evalue ) { evalue_ = evalue; };
		
		inline void filterOut() { blacklist_this_ = true; };

		bool parse( const std::string& line ) {
			if ( line.size() > 1 ) {
				std::vector< std::string > fields;
				if ( line[0] == '*' ) {
					blacklist_this_ = true;
					tokenizeSingleCharDelim( line.substr( 1 ), fields, default_field_separator, 12, false );
				} else {
					blacklist_this_ = false;
					tokenizeSingleCharDelim( line, fields, default_field_separator, 12, false );
				}
				return parse( fields );
			}
			return false;
		}
		
		virtual bool parse( const std::vector< std::string >& fields ) {
			if ( fields.size() >= 12 ) {
				try {
					query_start_ = boost::lexical_cast< large_unsigned_int >( fields[1] );
					query_stop_ = boost::lexical_cast< large_unsigned_int >( fields[2] );

					if( query_start_ > query_stop_ ) {
						std::cerr << "reverse query positions are not allowed (only reference positions can be swapped to indicate the reverse complement, adjust your input file format)" << std::endl;
						return false;
					}
					
					query_length_ = boost::lexical_cast< large_unsigned_int >( fields[3] );
					
					reference_start_ = boost::lexical_cast< large_unsigned_int >( fields[5] );
					reference_stop_ = boost::lexical_cast< large_unsigned_int >( fields[6] );

				} catch ( boost::bad_lexical_cast e ) {
					std::cerr << "could not parse position number or query length" << std::endl;
					return false;
				}

				try {
					score_ = boost::lexical_cast< float >( fields[7] );
				} catch( boost::bad_lexical_cast e ) {
					std::cerr << "could not parse score" << std::endl;
					return false;
				}
				
				try {
					evalue_ = boost::lexical_cast< double >( fields[8] );
				} catch( boost::bad_lexical_cast e ) {
					std::cerr << "could not parse E-value" << std::endl;
					return false;
				}
				
				try {
					identities_ = boost::lexical_cast< large_unsigned_int >( fields[9] );
				} catch( boost::bad_lexical_cast e ) {
					std::cerr << "could not identities" << std::endl;
					return false;
				}
				
				try {
					alignment_length_ = boost::lexical_cast< large_unsigned_int >( fields[10] );
				} catch( boost::bad_lexical_cast e ) {
					std::cerr << "could not parse alignment length" << std::endl;
					return false;
				}
				
				alignment_code_ = fields[11];
				
				// easy things that cannot go wrong (I know: what can go wrong, will go wrong)
				query_identifier_ = fields[0];
				reference_identifier_ = fields[4];
				
				return true;
			} else {
				std::cerr << "could not parse alignment because input line has too few entries" << std::endl;
			}
			return false;
		}
		
		void print( std::ostream& strm = std::cout ) const {
			if ( blacklist_this_ ) {
				strm << '*';
			}
			
			strm << query_identifier_ << default_field_separator
			     << query_start_ << default_field_separator
			     << query_stop_ << default_field_separator
			     << query_length_ << default_field_separator
			     << reference_identifier_ << default_field_separator
			     << reference_start_ << default_field_separator
			     << reference_stop_ << default_field_separator
			     << score_ << default_field_separator
			     << evalue_ << default_field_separator
			     << identities_ << default_field_separator
			     << alignment_length_ << default_field_separator
			     << alignment_code_ << default_field_separator
			     << endline;
		}
		
	private:
		std::string reference_identifier_;
		std::string query_identifier_;
		large_unsigned_int query_start_;
		large_unsigned_int query_stop_;
		large_unsigned_int query_length_;
		large_unsigned_int reference_start_;
		large_unsigned_int reference_stop_;
		float score_;
		double evalue_;
		large_unsigned_int identities_;
		large_unsigned_int alignment_length_;
		std::string alignment_code_;
		bool blacklist_this_;
};



//overload ostream operator for class AlignmentRecord ->print()
std::ostream& operator<<( std::ostream& strm, const AlignmentRecord& rec );

//overload istream operator for class AlignmentRecord ->parse()
std::istream& operator>>( std::istream& strm, AlignmentRecord& rec );




class AlignmentRecordTaxonomy : public AlignmentRecord {
	public:
		AlignmentRecordTaxonomy( StrIDConverter& converter, const Taxonomy* tax ) : acc2taxid_( converter ), taxinter( tax ) {}
		
		bool parse( const std::vector< std::string >& fields ) {
			if( this->AlignmentRecord::parse( fields ) ) {
				TaxonID taxid = acc2taxid_[ getReferenceIdentifier() ];				
				reference_node_ = taxinter.getNode( taxid );
				if( reference_node_ ) {
					return true;
				} else {
					std::cerr << "Could not find node with taxonomic id " << taxid << " in taxonomy" << std::endl;
				}
			}
			return false;
		}
		
		inline const TaxonNode* getReferenceNode() const { return reference_node_; }
	
	private:
		const TaxonNode* reference_node_;
// 		StrIDConverterFlatfileMemory& acc2taxid_;
		StrIDConverter& acc2taxid_;
		TaxonomyInterface taxinter;
};



template< typename T >
class AlignmentRecordFactory {
	typedef T AlignmentRecordType;
}; //TODO: add virtual create() and destroy functions



template<> //specialization for AlignmentRecord
class AlignmentRecordFactory< AlignmentRecord > {
	public:
		AlignmentRecordFactory() {}

		AlignmentRecord* create( const std::string& line ) {
			AlignmentRecord* rec = new AlignmentRecord;
			if ( rec->parse( line ) ) {
				return rec;
			}
			delete rec;
			return NULL;
		}
		
		void destroy( const AlignmentRecord* rec ) {
			delete rec;
		}
};



template<> //specialization for AlignmentRecordTaxonomy
class AlignmentRecordFactory< AlignmentRecordTaxonomy > {
	public:
		AlignmentRecordFactory( StrIDConverter& acc2taxid, const Taxonomy* tax ) : acc2taxid_( acc2taxid ), tax_( tax ) {}

		AlignmentRecordTaxonomy* create( const std::string& line ) {			
			
			AlignmentRecordTaxonomy* rec = new AlignmentRecordTaxonomy( acc2taxid_, tax_ );
			if ( rec->AlignmentRecord::parse( line ) ) { //is this really necessary?
				return rec;
			}
			delete rec;
			return NULL;
		}
		
		void destroy( const AlignmentRecordTaxonomy* rec ) {
			delete rec;
		}
		
	private:
		StrIDConverter& acc2taxid_;
		const Taxonomy* tax_;
};



template< typename RecordType = AlignmentRecord >
class AlignmentFileParser {
	public:
		typedef AlignmentRecordFactory< RecordType > FactoryType;
		
		AlignmentFileParser( const std::string& filename, FactoryType& factory ) : filehandle_( filename.c_str() ), handle_( filehandle_ ), factory_( factory ), line_num_( 0 ) {}
		
		AlignmentFileParser( std::istream& strm, FactoryType& factory ) : handle_( strm ), factory_( factory ), line_num_( 0 ) {}
		
		RecordType* next() {
// 			std::cerr << "calling parser.next()..." << std::endl;
			while (	std::getline( handle_, line_ ) ) {
				++line_num_;
				if ( ignoreLine( line_ ) ) continue;
				RecordType* rec = factory_.create( line_ );
				if ( rec ) {
					return rec;
				}
				std::cerr << "there was an error parsing line " << line_num_ << " in alignments, skipping..." << std::endl;
			}
			//std::cout << "UNEXPECTED INPUT STREAM ERROR WITH ALIGNMENTS FILE, RETURNING NULL POINTER OR EOF (check out .eof() function" << std::endl;
			return NULL;
		}
		
		inline void destroy( const RecordType* rec ) const { factory_.destroy( rec ); }
		inline bool eof() { return handle_.eof(); }
		
	private:
		std::ifstream filehandle_;
		std::istream& handle_;
		std::string line_;
		FactoryType& factory_;
		unsigned int line_num_;
};



template< typename RecordType >
class RecordSetGenerator {
	public:
		typedef AlignmentFileParser< RecordType > ParserType;
		
		RecordSetGenerator( ParserType& parser ) : parser_( parser ), last_record_( parser.next() ) {
			last_query_id_ = last_record_ ? &( last_record_->getQueryIdentifier() ) : NULL;
		}

		template< typename ContainerT >
		void getNext( ContainerT& recordset ) {
			if( last_record_ ) {
				RecordType* record = last_record_;
				const std::string& query_id = *last_query_id_;

				do {
					if( query_id == record->getQueryIdentifier() ) { //still the same query
						recordset.push_back( record );
					} else {
						last_query_id_ = &(record->getQueryIdentifier());
						last_record_ = record;
						break;
					}
					record = parser_.next();
				} while( record );
				last_record_ = record;
			}
		}

		bool notEmpty() { return last_record_; };

	private:
		ParserType& parser_;
		RecordType* last_record_;
		const std::string* last_query_id_;
};



template< typename ContainerT1, typename ContainerT2 >
void records2Nodes( const ContainerT1& recordset, const TaxonomyInterface& taxinter, StrIDConverter& acc2taxid, ContainerT2& refnodes ) {
	typename ContainerT1::const_iterator it = recordset.begin();
	const TaxonNode* node;
	while( it != recordset.end() ) {
		if( !(*it)->isFiltered() ) {
			TaxonID taxid = acc2taxid[ (*it)->getReferenceIdentifier() ];
			node = taxinter.getNode( taxid );
			if( node ) { //silently ignore non-mapping :)
				refnodes.push_back( node );
			}
		}
		++it;
	}
}



template< typename ContainerT1, typename ContainerT2 >
void records2Nodes( const ContainerT1& recordset, ContainerT2& refnodes ) {
	typename ContainerT1::const_iterator record_it = recordset.begin();
	const TaxonNode* node;
	while( record_it != recordset.end() ) {
		if( !(*record_it)->isFiltered() ) {
			node = (*record_it)->getReferenceNode();
			if( node ) {
				refnodes.push_back( node );
			}
		}
		++record_it;
	}
}



template< typename ContainerT >
void deleteRecords( ContainerT& recordset ) {
	while( recordset.size() ) {
		delete recordset.back();
		recordset.pop_back();
	}
}



template< typename ContainerT, typename QueueLikeContainer >
void separateAlignmentsByRange( ContainerT& recordset, QueueLikeContainer& workload ) {
	
	if ( recordset.empty() ) return;
	
	typedef typename ContainerT::value_type AlignmentRecordTypePtr; //expect stdcontainer
	
	// walk over original recordset and determine split point(s)
	std::size_t i = 0;
	std::vector< boost::tuple< large_unsigned_int, large_unsigned_int, AlignmentRecordTypePtr > > ranges( recordset.size() ); //temporary space
	for( typename ContainerT::const_iterator it = recordset.begin(); it != recordset.end(); ++it ) {
		ranges[i++] = boost::make_tuple( (*it)->getQueryStart(), (*it)->getQueryStop(), *it );
	}
// 	recordset.clear(); //remove pointers just to be sure
	
	// sort vector (in increasing order)
	std::sort( ranges.begin(), ranges.end() ); //TODO: sort by start is enough (maye use ordered map)
	
	// push into queue as separate sets to be treated independently by prediction algorithm
	ContainerT rset;
	
	large_unsigned_int start = boost::get<0>( ranges[0] );
	large_unsigned_int stop = boost::get<1>( ranges[0] );
	large_unsigned_int rstart = start;
	large_unsigned_int rstop = stop;
	rset.push_back( boost::get<2>( ranges[0] ) );
	
	for ( i = 1; i < ranges.size(); ++i ) {
		
		start = boost::get<0>( ranges[i] );
		stop = boost::get<1>( ranges[i] );
		
		if ( start > rstop ) { //split point detected
// 			std::cerr << "range from " << rstart << " to " << rstop << " with " << rset.size() << " alignments" << std::endl;
			workload.push( rset ); //copy pointer list to working queue
			rset.clear();
			rstop = stop;
			rstart = start;
		} else {
			rstop = std::max( rstop, stop );
		}
		
		rset.push_back( boost::get<2>( ranges[i] ) );
	}
	
	if ( ! rset.empty() ) {
// 		std::cerr << "range from " << rstart << " to " << rstop << " with " << rset.size() << " alignments" << std::endl;
		workload.push( rset );
	}
}

#endif // alignmentrecord_hh_
