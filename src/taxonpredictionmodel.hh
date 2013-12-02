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

#ifndef taxonpredictionmodel_hh_
#define taxonpredictionmodel_hh_

#include <cmath>
#include <boost/tuple/tuple.hpp>
#include <ostream>
#include "types.hh"
#include "alignmentrecord.hh"
#include "taxonomyinterface.hh"
#include "alignmentsfilter.hh"
#include "predictionrecord.hh"
#include "constants.hh"



template< typename ContainerT >
class TaxonPredictionModel {
	public:
		TaxonPredictionModel( const Taxonomy* tax ) : taxinter_( tax ), root_( taxinter_.getRoot() ) {};
		virtual void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) = 0;
	protected:
		void initPredictionRecord( ContainerT& recordset, PredictionRecord& prec ) {
			prec.initialize( recordset.front()->getQueryIdentifier(), recordset.front()->getQueryLength() );
		};
		const TaxonNode* unclassified() { return root_; };
		void setUnclassified( PredictionRecord& prec ) {
			prec.setNodePoint( root_, 0 );
		}
		TaxonomyInterface taxinter_;
		const TaxonNode* root_;
};



//only for debugging!
template< typename ContainerT >
class DummyPredictionModel : public TaxonPredictionModel< ContainerT > { //TODO: include winscore
	public:
		DummyPredictionModel( const Taxonomy* tax ) : TaxonPredictionModel< ContainerT >( tax ) {}

		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
};



template< typename ContainerT >
class LCASimplePredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		LCASimplePredictionModel( const Taxonomy* tax ) : TaxonPredictionModel< ContainerT >( tax ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			std::list< const TaxonNode* > refnodes;
			records2Nodes( recordset, refnodes );
			if( ! refnodes.empty() ) {
				prec.setNodePoint( this->taxinter_.getLCA( refnodes ) );
				return;
			}
			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
};



template< typename ContainerT >
class MeganLCAPredictionModel : public TaxonPredictionModel< ContainerT > { //TODO: include winscore
	public:
		MeganLCAPredictionModel( const Taxonomy* tax, bool iuc, const float toppercent, const float minscore, const int minsupport, const double maxevalue, const float winscore = 0.0) : TaxonPredictionModel< ContainerT >( tax ),
			msp( minsupport ), wsc( winscore ), ms_me_toppercent( minscore, maxevalue, toppercent ), ignore_unclassified( iuc ) {};
		
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			ms_me_toppercent.filter( recordset );
			std::list< const TaxonNode* > refnodes;
			
			typename ContainerT::iterator rec_it = firstUnmaskedIter( recordset );
			if ( rec_it != recordset.end() ) {
				large_unsigned_int qrstart = (*rec_it)->getQueryStart();
				large_unsigned_int qrstop = (*rec_it)->getQueryStop();
// 				large_unsigned_int qlength = (*rec_it)->getQueryLength();
				
				// determine range of query used
				large_unsigned_int n = 0;
				if( qrstart > qrstop ) std::swap( qrstart, qrstop ); //TODO: is this allowed? should we through an error/exception

				refnodes.push_back( (*rec_it)->getReferenceNode() );
				++n;
				++rec_it;
				
				for ( ;rec_it != recordset.end(); ++rec_it ) { 
					if ( ! (*rec_it)->isFiltered() ) {
						if( (*rec_it)->getQueryStart() <= (*rec_it)->getQueryStop() ) { //TODO: is this allowed? should we through an error/exception
							qrstart = std::min( (*rec_it)->getQueryStart(), qrstart );
							qrstop = std::max( (*rec_it)->getQueryStop(), qrstop );
						} else {
							qrstart = std::min( (*rec_it)->getQueryStop(), qrstart );
							qrstop = std::max( (*rec_it)->getQueryStart(), qrstop );
						}
						refnodes.push_back( (*rec_it)->getReferenceNode() );
						++n;
					}
				}
				
				// set region coverage on query
				prec.setQueryFeatureBegin( qrstart );
				prec.setQueryFeatureEnd( qrstop );
			}
			
			if( ignore_unclassified ) { //TODO: avoid extra pass
				std::list< const TaxonNode* >::iterator it = refnodes.begin();
				while( it != refnodes.end() ) {
					if( (*it)->data->is_unclassified ) {
						it = refnodes.erase( it );
					} else {
						++it;
					}
				}
			}

			if( refnodes.size() >= msp ) {
				prec.setNodePoint( this->taxinter_.getLCA( refnodes ) );
				return;
			}

			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
	private:
		unsigned int msp;
		float wsc;
		MinScoreMaxEvalueTopPercentFilter< ContainerT > ms_me_toppercent;
		bool ignore_unclassified;
};



template< typename ContainerT > //makes use of incomplete knowledge
class ICMeganLCAPredictionModel : public TaxonPredictionModel< ContainerT > { //TODO: include winscore
	public:
		ICMeganLCAPredictionModel( const Taxonomy* tax, const float toppercent, const float minscore, const int minsupport, const double maxevalue, const float winscore = 0.0) : TaxonPredictionModel< ContainerT >( tax ),
		msp( minsupport ), wsc( winscore ), ms_me_toppercent( minscore, maxevalue, toppercent ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			ms_me_toppercent.filter( recordset );
			std::list< const TaxonNode* > refnodes;
			records2Nodes( recordset, refnodes );
			if( refnodes.size() >= msp ) {
				prec.setNodePoint( this->taxinter_.getICLCA( refnodes ) );
				return;
			}

			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
	private:
		unsigned int msp;
		float wsc;
		MinScoreMaxEvalueTopPercentFilter< ContainerT > ms_me_toppercent;
};



template< typename ContainerT >
class NBestLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
	public:
		NBestLCAPredictionModel( const Taxonomy* tax, const int n = 1) : TaxonPredictionModel< ContainerT >( tax ), findnbest( n ) {};
		void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
			this->initPredictionRecord( recordset, prec);
			findnbest.filter( recordset );
			std::list< const TaxonNode* > refnodes;
			records2Nodes( recordset, refnodes );
			if( ! refnodes.empty() ) {
				prec.setNodePoint( this->taxinter_.getLCA( refnodes ) );
				return;
			}
			TaxonPredictionModel< ContainerT >::setUnclassified( prec );
		};
	private:
		NumBestBitscoreFilter< ContainerT > findnbest;
};


// the following classes are experimental

class PredictionsParser {
	public:
		PredictionsParser( std::istream& istr ) : source( istr ), not_empty( true ) {
			parseLine();
		};

		PredictionsParser( const std::string& filename ) : filehandle( filename.c_str() ), source( filehandle ), not_empty( true ) {
			parseLine();
		};

		boost::tuple< std::string, unsigned int > getNext() {
			boost::tuple< std::string, unsigned int> ret_type = boost::make_tuple( qid, taxid_pred );
			parseLine();
			return ret_type;
		};

		bool notEmpty() {
			return not_empty;
		};

	private:
		void parseLine() {
			std::list< std::string > fields;
			std::list< std::string >::iterator field_it;

			while( std::getline( source, line ) ) {
				if( ignoreLine( line ) ) {
					std::cerr << "PredictionsParser: Ignoring comment line" << std::endl;
					continue;
				}

				tokenizeSingleCharDelim( line, fields, default_field_separator, 2 );
				field_it = fields.begin();

				try {
					qid = *field_it++;
					taxid_pred = boost::lexical_cast< unsigned int >( *field_it );
					return;
				} catch ( boost::bad_lexical_cast e ) {
					std::cerr << "PredictionsParser: Could not parse line '" << line << '\'' << std::endl;
					fields.clear();
				}
			}
			not_empty = false;
		}

    std::ifstream filehandle;
		std::istream& source;
		std::string qid;
		std::string line;
		unsigned int taxid_pred;
		bool not_empty;
};



#endif // taxonpredictionmodel_hh_
