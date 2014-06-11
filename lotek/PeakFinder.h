/* -*- mode: c++ -*- */
/*
    Copyright 2011, 2014 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _PEAK_FINDER_H
#define _PEAK_FINDER_H

#include "RingBuffWithMax.h"

// PeakFinder: given a window of radius n (i.e. size 2n+1), search an
//             input sequence for locations when a local max occurs.
//             A local max occurs for sample j when input[j] >
//             input[k] for all k in [j - n, ..., j - 1, j + 1, ..., j
//             + n] Processing a new element from the input sequence
//             returns a boolean flag indicating whether a local max
//             occurred n elements earlier (which is the soonest we
//             can tell).  
//
// Storage: O(n) * sizeof (TYPE)
// Complexity: process() is O(log(n)); other operations constant.

template < typename TYPE >
class PeakFinder {
 public:
  PeakFinder (size_t n = 1) : 
    n(n),
    buf(2 * n + 1)
  {
  };

  operator TYPE () { // the value of peak; only valid if the most recent call to process returned true
    return buf[n];
  };

  bool operator () (const TYPE & d) { // process a value from the data stream
    if (buf.push_back(d)) {
      if (buf[n] == buf.max())
        return true;
    }
    return false;
  };

  TYPE & operator[](unsigned i) {
    return buf[i];
  };

  bool full() {
    return buf.full();
  };

 protected:
  int n;
  RingBuffWithMax < TYPE > buf; // elements in magnitude order
};

#endif //  _PEAK_FINDER_H
