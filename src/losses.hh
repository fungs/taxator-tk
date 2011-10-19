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

#ifndef losses_hh_
#define losses_hh_

#include <string>
#include <vector>
#include <cmath>
#include <map>
#include <boost/tuple/tuple.hpp>

class BaseLoss {
	public:
		BaseLoss( unsigned int size ) : values_( size, 0.0 ) {};
		virtual void add( const int a, const int b, const int c) = 0;

		virtual void add( const boost::tuple< int, int, int >& bt ) {
			int a, b, c;
			boost::tie( a, b, c ) = bt;
			add( a, b, c );
		};

		virtual std::vector< float > get() = 0;

		virtual int support() = 0;

		virtual void clear() {
			for( std::vector< float >::iterator it = values_.begin(); it != values_.end(); ++it ) {
				*it = 0.0;
			}
		};

		virtual unsigned int valNum() {
			return values_.size() + 1;
		};

	protected:
		std::vector< float > values_;
};



class BaseLossFactory {
	public:
		virtual BaseLoss* operator()() = 0;
};



template< class T >
class LossFactory : public BaseLossFactory {
	public:
		virtual BaseLoss* operator()() {
			return new T;
		}
};



class LossMap {
	public:
		typedef std::map< std::string, BaseLoss* > MapType;
		typedef MapType::iterator iterator;

		LossMap( BaseLossFactory& fac ) : factory_( fac ) {};
		BaseLoss& operator[]( const std::string& id ) {
			iterator it = store_.find( id );
			if( it != store_.end() ) {
				return *it->second;
			}
			return *store_.insert( std::make_pair( id, factory_() ) ).first->second;
		}

		BaseLoss& front() {
			return *store_.begin()->second;
		}

		iterator begin() {
			return store_.begin();
		}

		iterator end() {
			return store_.end();
		}

		~LossMap() {
			for( iterator it = store_.begin(); it != store_.end(); ++it ) {
				delete it->second;
			}
		}

	private:
		MapType store_;
		BaseLossFactory& factory_;
};



class L1Loss : public BaseLoss {
	public:
		L1Loss( const float al = 0.5) : BaseLoss( 3 ), count( 0 ), alpha( al ), oneminusalpha( 1.0 - al ) {};

		virtual void add( const int a, const int b, const int c ) {
			values_[0] += a;
			values_[1] += b;
			values_[2] += c;
			++count;
		};

		virtual void clear() {
			count = 0;
			BaseLoss::clear();
		};

		virtual std::vector< float > get() {
			std::vector< float > ret( 4 );
			float tmp_a = values_[0] / float( count );
			float tmp_b = values_[1] / float( count );
			float tmp_c = values_[2] / float( count );
			ret[0] = tmp_a;
			ret[1] = tmp_b;
			ret[2] = tmp_c;
			ret[3] = 2.0 * (alpha*tmp_a + oneminusalpha*tmp_b);
			return ret;
		};

		virtual int support() {
			return count;
		}

	protected:
	  unsigned int count;
		const float alpha;
		const float oneminusalpha;
};

class L2Loss : public BaseLoss {
	public:
		L2Loss( const float al = 0.5) : BaseLoss( 3 ), count( 0 ), alpha( al ), oneminusalpha( 1.0 - al ) {};

		virtual void add( const int a, const int b, const int c ) {
			values_[0] += a*a;
			values_[1] += b*b;
			values_[2] += c*c;
			++count;
		};

		virtual void clear() {
			count = 0;
			BaseLoss::clear();
		};

		virtual std::vector< float > get() {
			std::vector< float > ret( 4 );
			float tmp_a = values_[0] / float( count );
			float tmp_b = values_[1] / float( count );
			float tmp_c = values_[2] / float( count );
			ret[0] = std::sqrt( tmp_a );
			ret[1] = std::sqrt( tmp_b );
			ret[2] = std::sqrt( tmp_c );
			ret[3] = std::sqrt( 2.0 * (alpha*tmp_a + oneminusalpha*tmp_b) );
			return ret;
		};

		virtual int support() {
			return count;
		}

	protected:
	  unsigned int count;
		const float alpha;
		const float oneminusalpha;
};

#endif //losses_hh_
