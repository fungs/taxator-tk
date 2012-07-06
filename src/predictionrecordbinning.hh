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

#ifndef predictionrecordbinning_hh_
#define predictionrecordbinning_hh_

#include "predictionrecord.hh"

class PredictionRecordBinning : public PredictionRecordBase {
	public:
		enum BinningType { none, single, direct, fallback };
		
		PredictionRecordBinning( const Taxonomy* tax ) : PredictionRecordBase( tax ), binning_type_( none ) {}
		~PredictionRecordBinning() {	if( query_identifier_ ) delete query_identifier_; }
		
		void setQueryIdentifier( const std::string& id ) {
			if ( query_identifier_ ) delete query_identifier_;
			query_identifier_ = new const std::string( id );
		}
		
		//serialization
		virtual void print( std::ostream& strm = std::cout ) const { //write GFF3-style
			printColumns1to8( strm );
			printFeatureSeqLen( strm );
			strm << ';';
			printFeatureTax( strm );
			
			switch( binning_type_ ) {
				case single:
					if ( interpolation_value_ >= 0 ) {
						strm << ';';
						printFeatureIVal( strm );
					}
					strm << ";binning=single";
					break;
				case direct:
					strm << ";binning=direct";
					break;
				case fallback:
					strm << ";binning=fallback";
					break;
				default:
					break;
			}
			
			strm << endline;
		}
		
		void setBinningType( BinningType t ) { binning_type_ = t; };
		
	private:
		BinningType binning_type_;
};

#endif // predictionrecordbinning_hh_
