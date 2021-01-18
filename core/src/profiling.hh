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

#include <string>
#include <ctime>
#include "types.hh"

#ifndef profiling_hh_
#define profiling_hh_

class StopWatchCPUTime {
	public:
		StopWatchCPUTime( const std::string& info ) : info_( info ), stopped_(true), sum_( 0 ), counter_( 0 ), conversion_to_milliseconds_( CLOCKS_PER_SEC/1000 ) {}
		
// 		~StopWatchCPUTime() {
// 			std::cerr << info_ << " took: " << static_cast< int >( sum_ ) << '/' << sum_/static_cast< double >( counter_ ) << " (total/avg. in milliseconds)" << std::endl;
// 		}
		
		void start() {
            if(stopped_) {
                timestamp_ = clock();
                stopped_ = false;
            }
		}
		
		void stop() {
			if(!stopped_) {
                sum_ += (std::clock() - timestamp_)/conversion_to_milliseconds_;
                ++counter_;
                stopped_ = true;
            }
		}
		
		large_unsigned_int read() {
            if(stopped_) return sum_;
            return (std::clock() - timestamp_)/conversion_to_milliseconds_;
        }
		
	private:
		const std::string info_;
        bool stopped_;
        clock_t timestamp_;
		large_unsigned_int sum_;
		large_unsigned_int counter_;
		const unsigned int conversion_to_milliseconds_;
};

#endif //profiling_hh_