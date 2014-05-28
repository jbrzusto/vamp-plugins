/* -*- mode: c++ -*- */
/*
    Copyright 2014 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _SORTED_FIFO_H
#define _SORTED_FIFO_H

#include <boost/circular_buffer.hpp>
#include <set>

// SortedFifo: maintain a fifo of size at most n in sorted order.
// Adding to a full SortedFifo removes the oldest element.

// Storage: O(n) * sizeof (TYPE)
// Complexity: process() is O(log(n)); other operations constant.

template < typename TYPE >
class SortedFifo {
 public:
  SortedFifo (size_t n = 1) :
  buf(n),
    ordbuf()
  {
  };

  bool push_back(TYPE d) { // add a value; return true if the
    // fifo is full (whether or not it already was)

    // if the container is full, drop the first element
    // from the order buffer; (dropping from the circular buffer
    // is automatic)
    if (buf.full()) {
      // get the oldest element, and remove one copy
      // from ordbuf; we don't just do
      // ordbuf.erase(buf[0]) because that would erase
      // *all* elements with the same value, not just one

      ordbuf.erase(ordbuf.find(buf[0]));
    }
    // insert the new element in the circular buffer
    buf.push_back(d);

    // and the ordered multiset
    ordbuf.insert(d);

    return buf.full();
  };

  TYPE & back() {
    return buf.back();
  };

  TYPE & operator[](unsigned i) {
    return buf[i];
  };

  bool full() {
    return buf.full();
  };

  typename std::multiset < TYPE > :: reverse_iterator max() {
    // iterator from largest to smallest element
    return ordbuf.rbegin();
  };

  typename std::multiset < TYPE > :: iterator min() {
    // iterator from smallest to largest element
    return ordbuf.begin();
  };

 protected:
  boost::circular_buffer < TYPE > buf; // elements in FIFO order
  std::multiset < TYPE > ordbuf; // elements in magnitude order
};

#endif //  _SORTED_FIFO_H
