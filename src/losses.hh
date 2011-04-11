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
		BaseLoss( unsigned int size ) : values( size, 0.0 ) {};
		virtual void add( const int a, const int b, const int c) = 0;

		virtual void add( const boost::tuple< int, int, int >& bt ) {
			int a, b, c;
			boost::tie( a, b, c ) = bt;
			add( a, b, c );
		};

		virtual std::vector< float > get() = 0;

		virtual int support() = 0;

		virtual void clear() {
			for( std::vector< float >::iterator it = values.begin(); it != values.end(); ++it ) {
				*it = 0.0;
			}
		};

		virtual unsigned int valNum() {
			return values.size();
		};

	protected:
		std::vector< float > values;
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

		LossMap( BaseLossFactory* fac ) : factory( *fac ) {};
		BaseLoss& operator[]( const std::string& id ) {
			iterator it = store.find( id );
			if( it != store.end() ) {
				return *it->second;
			}
			return *store.insert( std::make_pair( id, factory() ) ).first->second;
		}

		BaseLoss& front() {
			return *store.begin()->second;
		}

		iterator begin() {
			return store.begin();
		}

		iterator end() {
			return store.end();
		}

		~LossMap() {
			//std::pair< std::string,BaseLoss* > p;
			for( iterator it = store.begin(); it != store.end(); ++it ) {
				delete it->second;
			}

		}

	private:
		MapType store; //destroys ojects by itself
		BaseLossFactory& factory;
};



class L1Loss : public BaseLoss {
	public:
		L1Loss( const float al = 0.5) : BaseLoss( 3 ), count( 0 ), alpha( al ), oneminusalpha( 1.0 - al ) {};

		virtual void add( const int a, const int b, const int c ) {
			values[0] += a;
			values[1] += c;
			++count;
		};

		virtual void clear() {
			count = 0;
			BaseLoss::clear();
		};

		virtual std::vector< float > get() {
			std::vector< float > ret( values.size() );
			float tmp_a = values[0] / float( count );
			float tmp_c = values[1] / float( count );
			ret[0] = tmp_a;
			ret[1] = tmp_c;
			ret[2] = 2.0 * (alpha*tmp_a + oneminusalpha*tmp_c);
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
			values[0] += a*a;
			values[1] += c*c;
			++count;
		};

		virtual void clear() {
			count = 0;
			BaseLoss::clear();
		};

		virtual std::vector< float > get() {
			std::vector< float > ret( values.size() );
			float tmp_a = values[0] / float( count );
			float tmp_c = values[1] / float( count );
			ret[0] = std::sqrt( tmp_a );
			ret[1] = std::sqrt( tmp_c );
			ret[2] = std::sqrt( 2.0 * (alpha*tmp_a + oneminusalpha*tmp_c) );
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
