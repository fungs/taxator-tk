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

#ifndef predictionrecord_hh_
#define predictionrecord_hh_

#include <string>
#include <fstream>
#include <iostream>
#include "types.hh"
#include "constants.hh"
#include "taxontree.hh"
#include "taxonomyinterface.hh"
#include "types.hh"
#include "utils.hh"


class PredictionRecordBase { //TODO: rename to something like feature
	public:
		PredictionRecordBase( const Taxonomy* tax ) : query_identifier_ ( NULL ), query_length_( 0 ), lower_node_( NULL ), upper_node_( NULL ), interpolation_value_( -1. ), signal_strength_( 0. ), taxinter_( tax ) {};
		virtual ~PredictionRecordBase() {}

		void initialize( const std::string& query_identifier, large_unsigned_int query_length ) { initialize( query_identifier, query_length, 0, query_length ); }
		
		void initialize( const std::string& query_identifier, large_unsigned_int query_length, large_unsigned_int feature_begin, large_unsigned_int feature_end ) {
			setQueryIdentifier( query_identifier );
			setQueryLength( query_length );
			setQueryFeatureBegin( feature_begin );
			setQueryFeatureEnd( feature_end );
		}
		
		//pure getters
		const std::string& getQueryIdentifier() const { return *query_identifier_; }
		large_unsigned_int getQueryLength() const { return query_length_; }
		large_unsigned_int getQueryFeatureBegin() const { return query_feature_begin_; }
		large_unsigned_int getQueryFeatureEnd() const { return query_feature_end_; }
		large_unsigned_int getQueryFeatureWidth() const { return query_feature_end_ - query_feature_begin_ + 1; }
		
		medium_unsigned_int getSupportAt( const TaxonNode* node ) const { return getSupportAt( node->data->root_pathlength ); } //TODO: range check
		
		medium_unsigned_int getSupportAt( small_unsigned_int depth ) const {
// 			std::cerr << "taxon support size is: " << taxon_support_.size() << std::endl;
// 			std::cerr << "upper node depth is: " << static_cast< int >( upper_node_->data->root_pathlength ) << std::endl;
// 			std::cerr << "lower node depth is: " << static_cast< int >( lower_node_->data->root_pathlength ) << std::endl;
// 			std::cerr << "lower node is: " << lower_node_->data->taxid << std::endl;
			int index = depth - upper_node_->data->root_pathlength;
			if ( index >= 0 && index < static_cast< const int >( taxon_support_.size() ) ) return taxon_support_.at( index );
			return 0;
		} //TODO: at->[]
		
		float getInterpolationValue() const { return interpolation_value_; }
		float getSignalStrength() const { return signal_strength_; }
		const TaxonNode* getUpperNode() const { return upper_node_; }
		const TaxonNode* getLowerNode() const { return lower_node_; }
		
		//pure setters
		virtual void setQueryIdentifier( const std::string& id ) = 0;
		void setQueryLength( large_unsigned_int i ) { query_length_ = i; }
		void setQueryFeatureBegin( large_unsigned_int i ) { query_feature_begin_ = i; }
		void setQueryFeatureEnd( large_unsigned_int i ) { query_feature_end_ = i; }
		void setInterpolationValue( float f ) { interpolation_value_ = f; }
		void setSignalStrength( float f ) { signal_strength_ = f; }
		
		void setNodeRange( const TaxonNode* lower_node, medium_unsigned_int lower_node_support, const TaxonNode* upper_node, medium_unsigned_int upper_node_support ) {
			lower_node_ = lower_node;
			upper_node_ = upper_node;
			taxon_support_ = std::vector< medium_unsigned_int >( lower_node->data->root_pathlength - upper_node->data->root_pathlength + 1, upper_node_support );
			taxon_support_.back() = lower_node_support;
		}
		
		void setNodeRange( const TaxonNode* lower_node, const TaxonNode* upper_node ) {
			medium_unsigned_int support = getQueryFeatureWidth();
			setNodeRange( lower_node, support, upper_node, support );
		}
		
		void setNodePoint( const TaxonNode* node, medium_unsigned_int support ) { setNodeRange( node, support, node, support ); } //alias for point estimate
		
		void setNodePoint( const TaxonNode* node ) { //alias for point estimate
			medium_unsigned_int support = getQueryFeatureWidth();
			setNodeRange( node, support, node, support );
		}
		
		void pruneLowerNode( const TaxonNode* node ) {
			assert( taxinter_.isParentOf( node, lower_node_ ) && node->data->root_pathlength >= upper_node_->data->root_pathlength ); //TODO: debug mode
			taxon_support_.resize( node->data->root_pathlength - upper_node_->data->root_pathlength + 1 );
			lower_node_ = node;
		}
		
		void setSupportAt( const TaxonNode* node, medium_unsigned_int support ) { setSupportAt( node->data->root_pathlength, support ); }
		void setSupportAt( small_unsigned_int depth, medium_unsigned_int support ) { taxon_support_.at( depth - upper_node_->data->root_pathlength ) = support; } //TODO: at->[]
		
		//de-serialization
		virtual bool parse( const std::string& line ) { //read GFF3-style
// 			std::cerr << "parsing" << tab << line << std::endl;
		
			if ( line.empty() ) return false;
			
			std::vector< std::string > fields;
			tokenizeSingleCharDelim( line, fields, tab_as_str, 9, false );
			
			if ( fields.size() < 9 ) {
				std::cerr << "could not parse alignment because input line has too few entries" << std::endl;
				return false;
			}
			
			if ( fields[1] != "taxator-tk" ) return false;
			
			try {
				setQueryFeatureBegin( boost::lexical_cast< large_unsigned_int >( fields[3] ) );
				setQueryFeatureEnd( boost::lexical_cast< large_unsigned_int >( fields[4] ) );

				if ( query_feature_begin_ > query_feature_end_ ) {
					std::cerr << "reverse feature positions are not allowed in taxator-tk GFF3 file" << std::endl;
					return false;
				}
				
			

			} catch ( boost::bad_lexical_cast e ) {
				std::cerr << "could not parse feature position number" << std::endl;
				return false;
			}

			try {
				setSignalStrength( boost::lexical_cast< float >( fields[5] ) );
			} catch( boost::bad_lexical_cast e ) {
				std::cerr << "could not parse signal strength (score field) in input" << std::endl;
				return false;
			}
			
// 			std::cerr << "basic fields are set, trying to parse column 9" << std::endl;
			
			{ //parse variable field (column 9)
				//split on ';'
				std::vector< std::string > key_value_fields;
				std::vector< std::string > key_value;
				tokenizeSingleCharDelim( fields[8], key_value_fields, ";", 0, true );
				for ( std::vector< std::string >::const_iterator it = key_value_fields.begin(); ! it->empty(); ++it ) {
					tokenizeSingleCharDelim( *it, key_value, "=", 2, false );
// 					std::cerr << "parsing key \"" << key_value[0] << "\" with value \"" << key_value[1] << "\"" << std::endl;
					if ( ! parseKeyValue( key_value[0], key_value[1] ) ) return false; //set values
					key_value.clear();
				}
			}
			
			// easy things that cannot go wrong (I know: what can go wrong, will go wrong)
			setQueryIdentifier( fields[0] );
			
// 			std::cerr << "parsing column 9 finished" << std::endl;
			return true;
		}

		
		//serialization
		virtual void print( std::ostream& strm = std::cout ) const { //write GFF3-style
			printColumns1to8( strm );
			printFeatureSeqLen( strm );
			strm << ';';
			printFeatureTax( strm );
			if ( interpolation_value_ >= 0 ) {
				strm << ';';
				printFeatureIVal( strm );
			}
			strm << endline;
		}
		
	protected:
		const std::string* query_identifier_;
		large_unsigned_int query_length_;
		large_unsigned_int query_feature_begin_;
		large_unsigned_int query_feature_end_;
		const TaxonNode* lower_node_;
		const TaxonNode* upper_node_; //can be removed with knowledge of lower_node_ and taxon_support_.size()
		float interpolation_value_;
		float signal_strength_;
		TaxonomyInterface taxinter_;
		std::vector< medium_unsigned_int > taxon_support_; //internal encoding of support, TODO: change to small_unsigned_int?
		
		void printColumns1to8( std::ostream& strm ) const {
			strm << *query_identifier_ << tab << "taxator-tk" << tab << "sequence_feature" << tab << query_feature_begin_ << tab << query_feature_end_ << tab << signal_strength_ << tab << '.' << tab << '.' << tab;
		}
		
		void printFeatureSeqLen( std::ostream& strm ) const { strm << "seqlen=" << query_length_; }
		void printFeatureIVal( std::ostream& strm ) const { strm << "ival=" << interpolation_value_; }
		void printFeatureTax( std::ostream& strm ) const {
			assert( lower_node_ );
			assert( upper_node_ );
			assert( ! taxon_support_.empty() );
			
			strm << "tax=";
			medium_unsigned_int last_support = 0;
			Taxonomy::PathUpIterator pit( lower_node_ );
			unsigned int i = taxon_support_.size() - 1;
			while ( pit != upper_node_ ) {
				if ( taxon_support_[i] != last_support ) {
					strm << pit->data->taxid << ':' << taxon_support_[i] << '-';
					last_support = taxon_support_[i];
				}
				--i;
				++pit;
			}
			strm << pit->data->taxid << ':' << taxon_support_[i];
		}
		
		virtual bool parseKeyValue( const std::string& key, const std::string& value ) {
			try {
				if ( key == "seqlen" ) {
					setQueryLength( boost::lexical_cast< large_unsigned_int >( value ) );
					return true;
				}
				if ( key == "ival" ) {
					setInterpolationValue( boost::lexical_cast< float >( value ) );
					return true;
				}
				if ( key == "tax" ) {
					taxon_support_.clear();
					std::vector< std::string > taxpath;
					std::vector< std::string > taxid_support;
					TaxonID taxid;
					medium_unsigned_int support;
					
					tokenizeSingleCharDelim( value, taxpath, "-", 0, false );
					
					std::vector< std::string >::const_iterator it = taxpath.begin();
					tokenizeSingleCharDelim( *it, taxid_support, ":", 2, false );
					taxid = boost::lexical_cast< TaxonID >( taxid_support[0] );
					if ( taxid_support[1].empty() ) support = getQueryFeatureWidth();
					else support = boost::lexical_cast< medium_unsigned_int >( taxid_support[1] );
					const TaxonNode* last_node = taxinter_.getNode( taxid );
					lower_node_ = last_node;
					std::list< medium_unsigned_int > tmp_taxon_support;
					
					while ( ! (++it)->empty() ) { //last field is empty by definition, if well formed
						taxid_support.clear();
						tokenizeSingleCharDelim( *it, taxid_support, ":", 2, false );
						taxid = boost::lexical_cast< TaxonID >( taxid_support[0] );
						const TaxonNode* node = taxinter_.getNode( taxid );
						
						//sanity check path
						if ( ! taxinter_.isParentOf( node, last_node ) ) {
							std::cerr << "incorrect taxonomic path, " << node->data->taxid << " (" << last_node->data->annotation->name << ") is not parent of " << last_node->data->taxid << " (" << node->data->annotation->name << ")" << std::endl;
							return false;
						}
						for ( Taxonomy::PathUpIterator pit( last_node ); pit != node; ++pit ) {
// 							std::cerr << "pushing support value for node: " << pit->data->taxid << std::endl;
							tmp_taxon_support.push_front( support );
						}
						
						if ( ! taxid_support[1].empty() ) support = boost::lexical_cast< medium_unsigned_int >( taxid_support[1] );
						last_node = node;
					}
					tmp_taxon_support.push_front( support );
					upper_node_ = last_node;
					
					// assign to taxon_support_
					taxon_support_.reserve( tmp_taxon_support.size() );
					std::copy( tmp_taxon_support.begin(), tmp_taxon_support.end(), std::back_inserter( taxon_support_ ) );
					
					assert( lower_node_->data->root_pathlength - upper_node_->data->root_pathlength + 1 == taxon_support_.size() );
					return true;
				}
			} catch( boost::bad_lexical_cast e ) {
				std::cerr << "could not parse value of attribute " << key << std::endl;
				return false;
			}
			return true;
		}
};



class PredictionRecordSaveMem : public PredictionRecordBase {
	public:
		PredictionRecordSaveMem( const Taxonomy* tax, ReferencedStringStore<>& qid_store ) : PredictionRecordBase( tax ), qid_store_( qid_store ) {};
		~PredictionRecordSaveMem() { qid_store_.remove( *query_identifier_ ); }
		void setQueryIdentifier( const std::string& id ) {
			if ( query_identifier_ ) qid_store_.remove( *query_identifier_ );
			query_identifier_ = &qid_store_.add( id );
		}
		
	private:
		ReferencedStringStore<>& qid_store_;
};



class PredictionRecord : public PredictionRecordBase {
	public:
		PredictionRecord ( const Taxonomy* tax ) : PredictionRecordBase( tax ) {}
		~PredictionRecord() {	if( query_identifier_ ) delete query_identifier_; }

		void setQueryIdentifier( const std::string& id ) {
			if ( query_identifier_ ) delete query_identifier_;
			query_identifier_ = new const std::string( id );
		}
};



template< class PredictionRecordType >
class PredictionFileParser {
	public:
		PredictionFileParser( const std::string& filename, const Taxonomy* tax ) : filehandle( filename.c_str() ), handle( filehandle ), tax_( tax ) {};
		PredictionFileParser( std::istream& strm, const Taxonomy* tax ) : handle( strm ), tax_( tax ) {};
		void destroyRecord( const PredictionRecordType* rec ) const { delete rec; }
		bool eof() const { return handle.eof(); };
		
		PredictionRecordType* next() {
			PredictionRecordType* rec = new PredictionRecordType( tax_ );
			std::string line;
			while( std::getline( handle, line ) ) {
				if( rec->parse( line ) ) return rec;
			}
			destroyRecord( rec );
			return NULL;
		}

	protected:
		std::ifstream filehandle;
		std::istream& handle;
		const Taxonomy* tax_;
};



std::ostream& operator<<( std::ostream& strm, const PredictionRecordBase& prec );


std::istream& operator>>( std::istream& strm, PredictionRecordBase& prec );


#endif // predictionrecord_hh_
