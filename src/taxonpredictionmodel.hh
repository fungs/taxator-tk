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
      prec.setBestReferenceTaxon(root_);
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



template< typename ContainerT, bool treat_unclassified=false >
class LCASimplePredictionModel : public TaxonPredictionModel< ContainerT > {
  public:
    LCASimplePredictionModel( const Taxonomy* tax ) : TaxonPredictionModel< ContainerT >( tax ) {};
    void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
      this->initPredictionRecord( recordset, prec);
      std::set< const TaxonNode* > refnodes;
      std::set< const TaxonNode* > refnodes_best;
      
      auto rec_it = firstUnmaskedIter( recordset );
      if ( rec_it != recordset.end() ) {
        large_unsigned_int qrstart = (*rec_it)->getQueryStart();
        large_unsigned_int qrstop = (*rec_it)->getQueryStop();
        float maxscore = (*rec_it)->getScore();
        
        { // determine range of query used
          large_unsigned_int n = 0;
          if( qrstart > qrstop ) std::swap( qrstart, qrstop ); //TODO: Should we through an error/exception instead?
          
          refnodes.insert( (*rec_it)->getReferenceNode() );
          ++n;
          ++rec_it;
          
          for ( ;rec_it != recordset.end(); ++rec_it ) { 
            if ( ! (*rec_it)->isFiltered() ) {
              if( (*rec_it)->getQueryStart() <= (*rec_it)->getQueryStop() ) { //TODO: Should we through an error/exception instead?
                qrstart = std::min( (*rec_it)->getQueryStart(), qrstart );
                qrstop = std::max( (*rec_it)->getQueryStop(), qrstop );
              } else {
                qrstart = std::min( (*rec_it)->getQueryStop(), qrstart );
                qrstop = std::max( (*rec_it)->getQueryStart(), qrstop );
              }
              if( (*rec_it)->getScore() > maxscore) maxscore = (*rec_it)->getScore();
              refnodes.insert( (*rec_it)->getReferenceNode() );
              ++n;
            }
          }
          prec.setQueryFeatureBegin( qrstart );
          prec.setQueryFeatureEnd( qrstop );
        }
        
        // find best reference taxon
        for (rec_it = recordset.begin() ;rec_it != recordset.end(); ++rec_it ) {
          if( !(*rec_it)->isFiltered() && (*rec_it)->getScore() == maxscore) refnodes_best.insert((*rec_it)->getReferenceNode());
        }
        
        if(treat_unclassified) {
          const TaxonNode* node = this->taxinter_.getLCC( refnodes );
          prec.setNodePoint( node );
          if(refnodes.size() != refnodes_best.size()) prec.setBestReferenceTaxon(this->taxinter_.getLCC(refnodes_best));
          else prec.setBestReferenceTaxon(node);
        } else {
          const TaxonNode* node = this->taxinter_.getLCA( refnodes );
          prec.setNodePoint( node );
          if(refnodes.size() != refnodes_best.size()) prec.setBestReferenceTaxon(this->taxinter_.getLCA(refnodes_best));
          else prec.setBestReferenceTaxon(node);
        }
        return;
      }
      
      TaxonPredictionModel< ContainerT >::setUnclassified( prec );
    }
};



template< typename ContainerT, bool treat_unclassified=false >
class MeganLCAPredictionModel : public TaxonPredictionModel< ContainerT > { //TODO: include winscore
  public:
    MeganLCAPredictionModel(const Taxonomy* tax, bool ignore_unclassified=false, const float toppercent=1.0, const float minscore=0.0, const int minsupport=1, const double maxevalue=std::numeric_limits<double>::max(), const float winscore=0.0) :
      TaxonPredictionModel< ContainerT >(tax),
      lca_simple_model_(tax),
      ms_me_toppercent_(minscore, maxevalue, toppercent),
      minsupport_(minsupport), winscore_(winscore), ignore_unclassified_(ignore_unclassified) {};
    
    void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
      ms_me_toppercent_.filter(recordset);
      if( ignore_unclassified_ ) remove_unclassified.filter(recordset);
      if( ms_me_toppercent_.getSupport() >= minsupport_ ) {
        lca_simple_model_.predict(recordset, prec, logsink);
        return;
      }

      this->initPredictionRecord( recordset, prec);
      TaxonPredictionModel< ContainerT >::setUnclassified( prec );
    };
  private:
    LCASimplePredictionModel<ContainerT, treat_unclassified> lca_simple_model_;
    RemoveUnclassifiedFilter<ContainerT> remove_unclassified;
    MinScoreMaxEvalueTopPercentFilter< ContainerT > ms_me_toppercent_;
    const unsigned int minsupport_;
    const float winscore_;
    const bool ignore_unclassified_;
};



// template< typename ContainerT > //makes use of incomplete knowledge
// class ICMeganLCAPredictionModel : public TaxonPredictionModel< ContainerT > { //TODO: include winscore
//   public:
//     ICMeganLCAPredictionModel(const Taxonomy* tax, bool ignore_unclassified=false, const float toppercent=1.0, const float minscore=0.0, const int minsupport=1, const double maxevalue=std::numeric_limits<double>::max(), const float winscore=0.0) :
//     TaxonPredictionModel< ContainerT >(tax),
//     minsupport_(minsupport), winscore_(winscore), ms_me_toppercent_(minscore, maxevalue, toppercent), ignore_unclassified_(ignore_unclassified) {};
//     
//     void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
//       this->initPredictionRecord( recordset, prec);
//       ms_me_toppercent_.filter( recordset );
//       std::set< const TaxonNode* > refnodes;
//       std::set< const TaxonNode* > refnodes_best;
//       
//       typename ContainerT::iterator rec_it = firstUnmaskedIter( recordset );
//       if ( rec_it != recordset.end() ) {
//         large_unsigned_int qrstart = (*rec_it)->getQueryStart();
//         large_unsigned_int qrstop = (*rec_it)->getQueryStop();
//         float maxscore = (*rec_it)->getScore();
//         
//         // determine range of query used
//         large_unsigned_int n = 0;
//         if( qrstart > qrstop ) std::swap( qrstart, qrstop ); //TODO: Should we through an error/exception instead?
//         
//         refnodes.insert( (*rec_it)->getReferenceNode() );
//         ++n;
//         ++rec_it;
//         
//         for ( ;rec_it != recordset.end(); ++rec_it ) { 
//           if ( ! (*rec_it)->isFiltered() ) {
//             if( (*rec_it)->getQueryStart() <= (*rec_it)->getQueryStop() ) { //TODO: Should we through an error/exception instead?
//               qrstart = std::min( (*rec_it)->getQueryStart(), qrstart );
//               qrstop = std::max( (*rec_it)->getQueryStop(), qrstop );
//             } else {
//               qrstart = std::min( (*rec_it)->getQueryStop(), qrstart );
//               qrstop = std::max( (*rec_it)->getQueryStart(), qrstop );
//             }
//             if( (*rec_it)->getScore() > maxscore) maxscore = (*rec_it)->getScore();
//             refnodes.insert( (*rec_it)->getReferenceNode() );
//             ++n;
//           }
//         }
//         
//         for (rec_it = recordset.begin() ;rec_it != recordset.end(); ++rec_it ) {
//           if( !(*rec_it)->isFiltered() && (*rec_it)->getScore() == maxscore) refnodes_best.insert((*rec_it)->getReferenceNode());
//         }
//         prec.setBestReferenceTaxon(this->taxinter_.getICLCA(refnodes_best));
//         
//         // set region coverage on query
//         prec.setQueryFeatureBegin( qrstart );
//         prec.setQueryFeatureEnd( qrstop );
//       }
//       
//       if( ignore_unclassified_ ) { //TODO: avoid extra pass
//         for( auto it = refnodes.begin(); it != refnodes.end(); ) {
//           if( (*it)->data->is_unclassified ) it = refnodes.erase( it );
//           else ++it;
//         }
//       }
//       
//       if( refnodes.size() >= minsupport_ ) {
//         prec.setNodePoint( this->taxinter_.getICLCA( refnodes ) );
//         return;
//       }
//       
//       TaxonPredictionModel< ContainerT >::setUnclassified( prec );
//     };
//   private:
//     const unsigned int minsupport_;
//     const float winscore_;
//     MinScoreMaxEvalueTopPercentFilter< ContainerT > ms_me_toppercent_;
//     const bool ignore_unclassified_;
//   };



template< typename ContainerT, bool treat_unclassified=false >
class NBestLCAPredictionModel : public TaxonPredictionModel< ContainerT > {
  public:
    NBestLCAPredictionModel( const Taxonomy* tax, const int n = 1) :
    TaxonPredictionModel< ContainerT >( tax ),
    lca_simple_model_(tax),
    findnbest( n ) {};
    
    void predict( ContainerT& recordset, PredictionRecord& prec, std::ostream& logsink ) {
      this->initPredictionRecord( recordset, prec);
      findnbest.filter( recordset );
      lca_simple_model_.predict(recordset, prec, logsink);
    };
  private:
    LCASimplePredictionModel<ContainerT, treat_unclassified> lca_simple_model_;
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
        } catch ( boost::bad_lexical_cast& ) {
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
