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

#ifndef concurrentoutstream_hh_
#define concurrentoutstream_hh_

#include <boost/thread.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include<ostream>
#include<sstream>

// this implementation uses stringstreams for buffering and is meant to be
// used until there is a decent logging mechanism in a stable boost release

class ConcurrentOutStream {
	public:
		ConcurrentOutStream( std::ostream& os, const uint threads, const uint buffer_size ) :
			os_( os ),
			max_buffer_size_( buffer_size ),
 			buffers_( threads ) //TODO: ensure exact buffer vector length
			{
				for ( uint i=0; i<threads; ++i ) buffers_.push_back( new std::ostringstream ); //because streams are not copyable
			};
		
		~ConcurrentOutStream() {
			for ( uint i=0; i<buffers_.size(); ++i ) flushSerial( i ); //no locking required
		}
		
		std::ostream& operator()( const uint channel ) { return buffers_[channel]; }
		
		void flush( const uint channel ) {
			if ( buffers_[channel].str().size() < max_buffer_size_ ) tryFlush( channel );
			else forceFlush( channel );
		}
		
		const uint channels() { return buffers_.size(); };
	
	protected:
		void tryFlush( const uint channel ) { // write if ostream not busy
			if ( mutex_.try_lock() ) {
// 				os_ << "normal write " << buffers_[channel].str().size() << std::endl;
				flushSerial( channel );
				mutex_.unlock();
			}
		}
		
		void forceFlush( const uint channel ) {
			if ( ! buffers_[channel].str().empty() ) {
// 				os_ << "forced write " << buffers_[channel].str().size() << std::endl;
				boost::mutex::scoped_lock( mutex_ );
				flushSerial( channel );
			}
		}
		
		void flushSerial( const uint channel ) {
			os_ << buffers_[channel].str();
			buffers_[channel].str("");
		}
		
		std::ostream& os_;
		const uint max_buffer_size_;
		boost::ptr_vector< std::ostringstream > buffers_;
		boost::mutex mutex_;
};

#endif // concurrentoutstream_hh_
