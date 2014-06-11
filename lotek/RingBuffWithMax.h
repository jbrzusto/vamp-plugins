/* -*- mode: c++ -*- */
/*
    Copyright 2014 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _RING_BUFF_WITH_MAX_H
#define _RING_BUFF_WITH_MAX_H

#include <boost/circular_buffer.hpp>

// RingBuffWithMax: a ring buffer that maintains the location of its (latest) max
// element.

// Storage: O(n) * sizeof (TYPE)

template < typename TYPE >
class RingBuffWithMax : public boost::circular_buffer < TYPE > {

private:
  typedef boost::circular_buffer < TYPE > super;

public:
  RingBuffWithMax (size_t n = 1) :
    boost::circular_buffer < TYPE > (n),
  maxloc(-1)
  {
  };

  bool push_back(TYPE d) { // add a value; return true if the
    // fifo is full (whether or not it already was)
    bool recalc_needed = false;

    if (maxloc < 0 || d >= max()) {
      maxloc = super::size() - 1;
    } else if (maxloc > 0) {
      --maxloc;
    } else {
      recalc_needed = true;
      maxloc = super::size() - 1; // assume new element is max
    }
    super::push_back( d );
    if (recalc_needed) {
      for (int i = super::size() - 1; i >= 0; --i) {
        if ((*this)[i] > max())
          maxloc = i;
      }
    }
    
    return this->full();
  };

  TYPE max() {
    if (maxloc >= 0)
      return (*this)[maxloc];
    throw std::runtime_error("no max available");
  };
      
 protected:
  int maxloc; // location of max element; -1 means not known
};

#endif //  _RING_BUFF_WITH_MAX_H
