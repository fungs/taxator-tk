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
#include "exception.hh"
#include "fileparser.hh"



class AlignmentRecord {
public:
    virtual ~AlignmentRecord() {};
    inline const std::string& getQueryIdentifier() const {
        return query_identifier_;
    };
    inline large_unsigned_int getQueryStart() const {
        return query_start_;
    };
    inline large_unsigned_int getQueryStop() const {
        return query_stop_;
    };
    inline large_unsigned_int getQueryLength() const {
        return query_length_;
    };
    inline const std::string& getReferenceIdentifier() const {
        return reference_identifier_;
    };
    inline large_unsigned_int getReferenceStart() const {
        return reference_start_;
    };
    inline large_unsigned_int getReferenceStop() const {
        return reference_stop_;
    };
    inline float getScore() const {
        return score_;
    };
    inline double getEValue() const {
        return evalue_;
    };
    inline large_unsigned_int getIdentities() const {
        return identities_;
    };
    inline large_unsigned_int getAlignmentLength() const {
        return alignment_length_;
    };
    inline const std::string& getAlignmentCode() const {
        return alignment_code_;
    };
    inline bool isFiltered() const {
        return blacklist_this_;
    };
    inline float getPID() const {
        return identities_/float( std::max( query_length_, alignment_length_ ) );
    };

    inline void filterOut() {
        blacklist_this_ = true;
    };
    
    inline bool operator<(const AlignmentRecord& other) const {
        if (this->getScore() < other.getScore()) return true;
        if (this->getScore() > other.getScore()) return false;
        return this->getIdentities() < other.getIdentities();
    }

    void parse( const std::string& line ) {
        if (line.size() <= 1) BOOST_THROW_EXCEPTION(ParsingError {} << general_info {"alignment line too short"});
        std::vector< std::string > fields;
        if ( line[0] == '*' ) {
            blacklist_this_ = true;
            tokenizeSingleCharDelim( line.substr( 1 ), fields, default_field_separator, 12, false );
        } else {
            blacklist_this_ = false;
            tokenizeSingleCharDelim( line, fields, default_field_separator, 12, false );
        }
        parse( fields );
    }

    virtual void parse( const std::vector< std::string >& fields ) {
        if ( fields.size() >= 11 ) {
            try {
                query_start_ = boost::lexical_cast< large_unsigned_int >( fields[1] );
                query_stop_ = boost::lexical_cast< large_unsigned_int >( fields[2] );

                if( query_start_ > query_stop_ ) BOOST_THROW_EXCEPTION(ParsingError {} << general_info {"reverse query positions not allowed (only reference positions can be swapped to indicate the reverse complement, adjust input"});

                query_length_ = boost::lexical_cast< large_unsigned_int >( fields[3] );

                reference_start_ = boost::lexical_cast< large_unsigned_int >( fields[5] );
                reference_stop_ = boost::lexical_cast< large_unsigned_int >( fields[6] );

            } catch(boost::bad_lexical_cast&) {
                BOOST_THROW_EXCEPTION(ParsingError {} << general_info {"bad position number or query length"});
            }

            try {
                score_ = boost::lexical_cast< float >( fields[7] );
            } catch(boost::bad_lexical_cast&) {
                BOOST_THROW_EXCEPTION(ParsingError {} << general_info {"bad score"});
            }

            try {
                evalue_ = boost::lexical_cast< double >( fields[8] );
            } catch(boost::bad_lexical_cast&) {
                BOOST_THROW_EXCEPTION(ParsingError {} << general_info {"bad E-value"});
            }

            try {
                identities_ = boost::lexical_cast< large_unsigned_int >( fields[9] );
            } catch(boost::bad_lexical_cast&) {
                BOOST_THROW_EXCEPTION(ParsingError {} << general_info {"bad identity value"});
            }

            try {
                alignment_length_ = boost::lexical_cast< large_unsigned_int >( fields[10] );
            } catch(boost::bad_lexical_cast&) {
                BOOST_THROW_EXCEPTION(ParsingError {} << general_info {"bad alignment length"});
            }

	    if(fields.size() >= 12){
            	alignment_code_ = fields[11];
	    }

            // easy things that cannot go wrong
            query_identifier_ = fields[0];
            reference_identifier_ = fields[4];

        } else BOOST_THROW_EXCEPTION(ParsingError {} << general_info {"bad number of fields in alignment line"});
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

    void parse( const std::vector< std::string >& fields ) {
        this->AlignmentRecord::parse( fields );

        TaxonID taxid;
        try {
            taxid = acc2taxid_[getReferenceIdentifier()];
        }
        catch(Exception &e) {
            BOOST_THROW_EXCEPTION(e << general_info {"bad taxon mapping for alignment reference sequence"});
        }

        try {
            reference_node_ = taxinter.getNode( taxid );
        }
        catch(Exception &e) {
            BOOST_THROW_EXCEPTION(e << general_info {"bad alignment reference taxon"});
        }
    }

    inline const TaxonNode* getReferenceNode() const {
        return reference_node_;
    }

private:
    const TaxonNode* reference_node_;
    StrIDConverter& acc2taxid_;
    TaxonomyInterface taxinter;
};



template< typename T >
class AlignmentRecordFactory {
    typedef T value_type;
}; //TODO: add virtual create() and destroy functions



template<>
class AlignmentRecordFactory< AlignmentRecord > {
public:
    typedef AlignmentRecord value_type;
    
    AlignmentRecordFactory() {}

    AlignmentRecord* create(const std::string& line) {
        AlignmentRecord* rec = new AlignmentRecord;
        try {
            rec->parse(line);
        } catch (Exception &e) {  // prevent memory leak
            destroy(rec);
            BOOST_THROW_EXCEPTION(e);
        }
        return rec;
    }

    inline void destroy( const AlignmentRecord* rec ) { delete rec; }
};



template<>
class AlignmentRecordFactory< AlignmentRecordTaxonomy > {
public:
    typedef AlignmentRecordTaxonomy value_type;
    
    AlignmentRecordFactory( StrIDConverter& acc2taxid, const Taxonomy* tax ) : acc2taxid_( acc2taxid ), tax_( tax ) {}
    
    AlignmentRecordTaxonomy* create( const std::string& line ) {
        AlignmentRecordTaxonomy* rec = new AlignmentRecordTaxonomy( acc2taxid_, tax_ );
        try {
            rec->AlignmentRecord::parse( line );
        } catch (Exception &e) {  // prevent memory leak
            destroy(rec);
            BOOST_THROW_EXCEPTION(e);
        }
        return rec;
    }

private:
    inline void destroy( const AlignmentRecordTaxonomy* rec ) { delete rec; }
    StrIDConverter& acc2taxid_;
    const Taxonomy* tax_;
};



template< typename RecordType, typename RecordSetType >
class RecordSetGenerator {
public:
    virtual ~RecordSetGenerator() {};
    virtual void getNext( RecordSetType& rset ) = 0;
    virtual bool notEmpty () {
        return record_ ;
    };

private:
    RecordType* record_;
};


// template<typename RecordType, typename RecordSetType, bool split_alignments>
// class RecordSetGeneratorUnsorted : public RecordSetGenerator< RecordType, RecordSetType > {
// public:
//     typedef FileParser< RecordType > ParserType;
// 
//     RecordSetGeneratorUnsorted(ParserType& parser) : parser_( parser ), record_( parser.next()) , tmpindex_(0)  {
//         last_query_id_ = record_ ? &( record_->getQueryIdentifier() ) : NULL;
//     };
// 
//     typedef typename RecordSetType::value_type AlignmentRecordTypePtr;
//     std::vector< boost::tuple< large_unsigned_int, large_unsigned_int, AlignmentRecordTypePtr > > ranges;
// 
//     void getNext(RecordSetType& rset) {
//         if(!split_alignments) {
//             if(record_) {  // always true unless called on empty input
//                 const std::string& query_id = *last_query_id_;
//                 rset.push_back(record_);
// 
//                 while(true) {
//                     if(parser_.eof()) {
//                         record_ = NULL;
//                         break;
//                     }
// 
//                     RecordType* record = parser_.next();
// 
//                     if( query_id == record->getQueryIdentifier() ) { //still the same query
//                         rset.push_back( record );
//                     } else {
//                         last_query_id_ = &(record->getQueryIdentifier());
//                         record_ = record;
//                         break;
//                     }
//                 };
//             } else BOOST_THROW_EXCEPTION(EOFError {} << general_info {"alignment recordset parser: trying to read from empty input"});
//         }
// 
//         else {  // split alignments
//             if(ranges.empty()) {  // read new query
//                 if(record_) {  // always true unless called on empty input
//                     const std::string& query_id = *last_query_id_;
//                     ranges.push_back(boost::make_tuple(record_->getQueryStart(), record_->getQueryStop(), record_));
// //                     rset.push_back(record_);
// 
//                     while(true) {
//                         if(parser_.eof()) {
//                             record_ = NULL;
//                             break;
//                         }
// 
//                         RecordType* record = parser_.next();
// 
//                         if( query_id == record->getQueryIdentifier() ) { //still the same query
//                             ranges.push_back(boost::make_tuple(record->getQueryStart(), record->getQueryStop(),record));
//                         } else {
//                             last_query_id_ = &(record->getQueryIdentifier());
//                             record_ = record;
//                             break;
//                         }
//                     };
//                     std::sort( ranges.begin(), ranges.end() ); //TODO: sort by start is enough (maybe use ordered map)
//                 } else BOOST_THROW_EXCEPTION(EOFError {} << general_info {"alignment recordset parser: trying to read from empty input"});
//             }
// 
//             // push into queue as separate sets to be treated independently by prediction algorithm
//             large_unsigned_int start = boost::get<0>( ranges[tmpindex_] );
//             large_unsigned_int stop = boost::get<1>( ranges[tmpindex_] );
//             large_unsigned_int rstop = stop;
//             rset.push_back( boost::get<2>( ranges[tmpindex_] ) );
//             for (std::size_t i = tmpindex_+1; i < ranges.size(); ++i) {
// 
//                 start = boost::get<0>(ranges[i]);
//                 stop = boost::get<1>(ranges[i]);
// 
//                 if (start > rstop) { //split point detected
//                     tmpindex_ = i;
//                     return;
//                 } else {
//                     rstop = std::max(rstop, stop);
//                 }
//                 rset.push_back(boost::get<2>(ranges[i]));
//             }
//             ranges.clear();
//             tmpindex_ = 0;
//         }
//     };
// 
//     bool notEmpty() {
//         if(!split_alignments) return record_ ;
//         else return (record_ || (ranges.size() > tmpindex_));
//     };
// 
// private:
//     ParserType& parser_;
//     RecordType* record_;
//     const std::string* last_query_id_;
//     large_unsigned_int tmpindex_;
// 
// };


template<typename RecordType, typename RecordSetType, bool split_alignments>
class RecordSetGeneratorUnsorted;


template<typename RecordType, typename RecordSetType>
class RecordSetGeneratorUnsorted<RecordType, RecordSetType, true> : public RecordSetGenerator< RecordType, RecordSetType > {
public:
    typedef FileParser< AlignmentRecordFactory< RecordType> > ParserType;
    RecordSetGeneratorUnsorted(ParserType& parser);
    void getNext(RecordSetType& rset);
    bool notEmpty();

private:
    ParserType& parser_;
    RecordType* record_;
    const std::string* last_query_id_;
    typedef typename RecordSetType::value_type AlignmentRecordTypePtr;
    std::vector< boost::tuple< large_unsigned_int, large_unsigned_int, AlignmentRecordTypePtr > > ranges;
    
    large_unsigned_int tmpindex_ = 0;

};


// specialization which splits the alignments
template< typename RecordType, typename RecordSetType >
RecordSetGeneratorUnsorted<RecordType, RecordSetType, true>::RecordSetGeneratorUnsorted(ParserType& parser) : parser_(parser) {
        if (parser_.eof()) {
            record_ = NULL;
            last_query_id_ = NULL;
        }
        else {
            record_ = parser_.next();
            last_query_id_ = &(record_->getQueryIdentifier());
        }
}


template< typename RecordType, typename RecordSetType >
bool RecordSetGeneratorUnsorted<RecordType, RecordSetType, true>::notEmpty() {
        return (record_ || (ranges.size() > tmpindex_));
}


template<typename RecordType, typename RecordSetType>
void RecordSetGeneratorUnsorted<RecordType, RecordSetType, true>::getNext(RecordSetType& rset) {
    if(ranges.empty()) {  // read new query
        if(record_) {  // always true unless called on empty input
            const std::string& query_id = *last_query_id_;
            ranges.push_back(boost::make_tuple(record_->getQueryStart(), record_->getQueryStop(), record_));

            while(true) {
                if(parser_.eof()) {
                    record_ = NULL;
                    break;
                }

                RecordType* record = parser_.next();

                if( query_id == record->getQueryIdentifier() ) { //still the same query
                    ranges.push_back(boost::make_tuple(record->getQueryStart(), record->getQueryStop(),record));
                } else {
                    last_query_id_ = &(record->getQueryIdentifier());
                    record_ = record;
                    break;
                }
            };
            std::sort( ranges.begin(), ranges.end() ); //TODO: sort by start is enough (maybe use ordered map)
        } else BOOST_THROW_EXCEPTION(EOFError {} << general_info {"alignment recordset parser: read from empty input"});
    }

    // push into queue as separate sets to be treated independently by prediction algorithm
    large_unsigned_int start = boost::get<0>( ranges[tmpindex_] );
    large_unsigned_int stop = boost::get<1>( ranges[tmpindex_] );
    large_unsigned_int rstop = stop;
    rset.push_back( boost::get<2>( ranges[tmpindex_] ) );
    for (std::size_t i = tmpindex_+1; i < ranges.size(); ++i) {

        start = boost::get<0>(ranges[i]);
        stop = boost::get<1>(ranges[i]);

        if (start > rstop) { //split point detected
            tmpindex_ = i;
            return;
        } else {
            rstop = std::max(rstop, stop);
        }
        rset.push_back(boost::get<2>(ranges[i]));
    }
    ranges.clear();
    tmpindex_ = 0;
}


// specialization which doesn't split the alignments
template<typename RecordType, typename RecordSetType>
class RecordSetGeneratorUnsorted<RecordType, RecordSetType, false> : public RecordSetGenerator< RecordType, RecordSetType > {
public:
    typedef FileParser< AlignmentRecordFactory< RecordType> > ParserType;
    RecordSetGeneratorUnsorted(ParserType& parser);
    void getNext(RecordSetType& rset);
    bool notEmpty();

private:
    ParserType& parser_;
    RecordType* record_;
    const std::string* last_query_id_;

};


template< typename RecordType, typename RecordSetType >
RecordSetGeneratorUnsorted<RecordType, RecordSetType, false>::RecordSetGeneratorUnsorted(ParserType& parser) : parser_(parser) {
    if (parser_.eof()) {
        record_ = NULL;
        last_query_id_ = NULL;
    }
    else {
        record_ = parser_.next();
        last_query_id_ = &(record_->getQueryIdentifier());
    }
}

    
template< typename RecordType, typename RecordSetType >
bool RecordSetGeneratorUnsorted<RecordType, RecordSetType, false>::notEmpty() {
        return record_ ;
}


template< typename RecordType, typename RecordSetType >
void RecordSetGeneratorUnsorted< RecordType, RecordSetType, false >::getNext(RecordSetType& rset) {
    if(record_) {  // always true unless called on empty input
        const std::string& query_id = *last_query_id_;
        rset.push_back(record_);

        while(true) {
            if(parser_.eof()) {
                record_ = NULL;
                break;
            }

            RecordType* record = parser_.next();

            if( query_id == record->getQueryIdentifier() ) { //still the same query
                rset.push_back( record );
            } else {
                last_query_id_ = &(record->getQueryIdentifier());
                record_ = record;
                break;
            }
        };
    } else BOOST_THROW_EXCEPTION(EOFError {} << general_info {"alignment recordset parser: read from empty input"});
}


template< typename RecordType, typename RecordSetType, bool split_alignments >
class RecordSetGeneratorSorted : public RecordSetGenerator< RecordType, RecordSetType > {
public:
    typedef FileParser< AlignmentRecordFactory< RecordType> > ParserType;

    RecordSetGeneratorSorted( ParserType& parser ) : parser_(parser) {
        if (parser_.eof()) {
            record_ = NULL;
            last_query_id_ = NULL;
        }
        else {
            record_ = parser_.next();
            last_query_id_ = &(record_->getQueryIdentifier());
            rstop_ = record_->getQueryStop();
        }
    }


    void getNext( RecordSetType& rset ) {
        typedef typename RecordSetType::value_type AlignmentRecordTypePtr;

        while( record_ ) {  // TODO: rework
            AlignmentRecordTypePtr record = record_;
            const std::string& query_id = *last_query_id_;

            large_unsigned_int start = record->getQueryStart();
            large_unsigned_int stop = record->getQueryStop();

            if(query_id == record->getQueryIdentifier()) {
                if (split_alignments) {
                    if (start > rstop_) {  // split
                        last_query_id_ = &(record->getQueryIdentifier());  //TODO: check if necessary!
                        rstop_ = stop;
                        return;
                    }
                    else {  // no split
                        rstop_ = std::max( record->getQueryStop() , rstop_ );
                        last_query_id_ = &(record->getQueryIdentifier());
                        rset.push_back(record);
                        if (parser_.eof()) record_ = NULL;
                        else record_ = parser_.next();
                    }
                }
                else rset.push_back(record);
            }
            else {  // new recordset
                last_query_id_ = &(record->getQueryIdentifier());
                rstop_ = stop;
                return;
            }
        }
    }

    bool notEmpty() {
        return record_;
    };

private:
    ParserType& parser_;
    RecordType* record_;
    const std::string* last_query_id_;
    large_unsigned_int rstop_;
};


template< typename ContainerT1, typename ContainerT2 >
void records2Nodes( const ContainerT1& recordset, const TaxonomyInterface& taxinter, StrIDConverter& acc2taxid, ContainerT2& refnodes ) {
    typename ContainerT1::const_iterator it = recordset.begin();
    const TaxonNode* node;
    while( it != recordset.end() ) {
        if( !(*it)->isFiltered() ) {
            TaxonID taxid = acc2taxid[ (*it)->getReferenceIdentifier() ];
            refnodes.insert(refnodes.end(), taxinter.getNode( taxid ) );
        }
        ++it;
    }
}



template< typename ContainerT1, typename ContainerT2 >
void records2Nodes( const ContainerT1& recordset, ContainerT2& refnodes ) {
    typename ContainerT1::const_iterator record_it = recordset.begin();
    while( record_it != recordset.end() ) {
        if( !(*record_it)->isFiltered() ) refnodes.insert(refnodes.end(), (*record_it)->getReferenceNode() );
        ++record_it;
    }
}



template< typename ContainerT >
void deleteRecords( ContainerT& recordset ) {
    for(typename ContainerT::const_iterator rec_it = recordset.begin(); rec_it != recordset.end(); ++rec_it) delete *rec_it;
    recordset.clear();
}


#endif // alignmentrecord_hh_

