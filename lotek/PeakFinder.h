/* -*- mode: c++ -*- */
/*
    Copyright 2011, 2014 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _PEAK_FINDER_H
#define _PEAK_FINDER_H

#include "SortedFifo.h"

// PeakFinder: given a window of radius n (i.e. size 2n+1), search an
//             input sequence for locations when a local max or min occurs.
//             A local max occurs for sample j when input[j] >
//             input[k] for all k in [j - n, ..., j - 1, j + 1, ..., j
//             + n] Processing a new element from the input sequence
//             returns a boolean flag indicating whether a local max
//             occurred n elements earlier (which is the soonest we
//             can tell).  Users can specify whether the finder is to
//             look for a local max, a local min, or either.
//
// Storage: O(n) * sizeof (TYPE)
// Complexity: process() is O(log(n)); other operations constant.

template < typename TYPE >
class PeakFinder {
 public:
  PeakFinder (size_t n = 1, bool find_max = true, bool find_min = false) :
    n(n),
    buf(2 * n + 1),
    find_max(find_max),
    find_min(find_min)
  {
  };

  operator TYPE () { // the value of peak; only valid if the most recent call to process returned true
    return buf[n];
  };

  bool operator () (const TYPE & d) { // process a value from the data stream
    if (buf.push_back(d)) {
      TYPE mid = buf[n];
      if (find_max) {
        // if looking for a max, see whether this is one
        auto max = buf.max();
        if (mid == *max && mid != *++max) {
          was_max = true;
          return true;
        }
      } 
      if (find_min) {
        // if looking for a min, see whether this is one
        auto min = buf.min();
        if (mid == *min && mid != *++min) {
          was_max = false;
          return true;
        }
      }
    }
    return false;
  };

  TYPE & operator[](unsigned i) {
    return buf[i];
  };

  bool full() {
    return buf.full();
  };

  operator bool() {
    return was_max;
  };

 protected:
  int n;
  SortedFifo < TYPE > buf; // elements in magnitude order
  bool find_max; // if true, find a local maximum
  bool find_min; // if true, find a local minimum
  bool was_max; // if true, last peak was a max; otherwise, last peak was a min
};

#endif //  _PEAK_FINDER_H
