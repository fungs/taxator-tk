// Comparison of bounded buffers based on different containers.
// Copyright (c) 2003-2008 Jan Gaspar
// Use, modification, and distribution is subject to the Boost Software
// License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)

// slight changes to original version!

#ifndef boundedbuffer_hh_
#define boundedbuffer_hh_

#include <boost/circular_buffer.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>
#include <boost/thread/thread.hpp>
#include <boost/progress.hpp>
#include <boost/bind.hpp>

template <class T>
class BoundedBuffer {
	public:

		typedef boost::circular_buffer<T> container_type;
		typedef typename container_type::size_type size_type;
		typedef typename container_type::value_type value_type;

		explicit BoundedBuffer(size_type capacity) : m_unread_( 0 ), m_container_( capacity ) {}

		void push(const value_type& item) {
			boost::mutex::scoped_lock lock( m_mutex_ );
			m_not_full_.wait(lock, boost::bind(&BoundedBuffer<value_type>::is_not_full, this));
			m_container_.push_front(item);
			++m_unread_;
			lock.unlock();
			m_not_empty_.notify_one();
// 			std::cerr << "push()" << std::endl;
		}

		value_type pop() {
			boost::mutex::scoped_lock lock( m_mutex_ );
			m_not_empty_.wait(lock, boost::bind(&BoundedBuffer<value_type>::is_not_empty, this));
			value_type retobj = m_container_[ --m_unread_ ];
			lock.unlock();
			m_not_full_.notify_one();
			if ( empty() ) empty_.notify_all();
// 			std::cerr << "pop()" << std::endl;
			return retobj;
		}
		
		void waitUntilEmpty() {
			if ( empty() ) return;
			boost::unique_lock< boost::mutex > lock( m_mutex_ );
			empty_.wait( lock );
		}
		
		size_type size() { return m_container_.size(); }
		bool empty() { return ! m_unread_; }
		
	private:
		BoundedBuffer(const BoundedBuffer&);              // Disabled copy constructor
		BoundedBuffer& operator = (const BoundedBuffer&); // Disabled assign operator

		bool is_not_empty() const { return m_unread_; }
		bool is_not_full() const { return m_unread_ < m_container_.capacity(); }

		size_type m_unread_;
		container_type m_container_;
		boost::mutex m_mutex_;
		boost::condition m_not_empty_;
		boost::condition m_not_full_;
		boost::condition empty_;
};

#endif //boundedbuffer_hh_