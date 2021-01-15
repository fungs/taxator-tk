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

#ifndef ostreamwrapperts_hh_
#define ostreamwrapperts_hh_

#include <boost/thread/mutex.hpp>
#include <boost/thread/locks.hpp>
#include <boost/thread/condition.hpp>
#include<ostream>

class OStreamWrapperTS : public std::ostream {
	public:
		OStreamWrapperTS( std::ostream& strm ) : strm_( strm ) {}
		
		OStreamWrapperTS& operator<<( std::istream& strm ) {
			strm_ << strm;
			return *this;
		}
		
		boost::mutex& getMutex() { return m_mutex_; }
		
	private:
		std::ostream& strm_;
		boost::mutex m_mutex_;
};

#endif //ostreamwrapperts_hh_
