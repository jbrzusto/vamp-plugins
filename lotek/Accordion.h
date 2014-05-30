/* -*- mode: c++ -*- */
/*

  Accordion.h: map an input sequence x[i] to x[j[i]] such that the sequence j[i]
  forms an 'accordion fold' of specified width and overlap over the natural numbers.
  e.g. if width = 6 and overlap = 4, the sequence

     x[1], x[2], x[3], ...

  maps to

 <----------- width = 6 ---------->   

 x[1], x[2], x[3], x[4], x[5], x[6],  

             <--- overlap = 4 ---->                  
             <----------- width = 6 ----------> 
             x[3], x[4], x[5], x[6], x[7], x[8],

                         <--- overlap = 4 ---->   
                         <----------- width = 6 ------------->   
                         x[5], x[6], x[7], x[8], x[9], x[10], ...

  This can be used for overlapping DFT frames.

  Samples in each consecutive window of given width are available as
  two <pointer, length> pairs after the first 'width' samples, and
  thereafter every width - overlap samples.

  Copyright 2014 John Brzustowski
  Licence: GPL version 2 or later

*/

#ifndef _ACCORDION_H
#define _ACCORDION_H

#include <boost/circular_buffer.hpp>

template < typename TYPE > 

class Accordion {

 public:

  typedef typename boost::circular_buffer < TYPE > ::const_array_range array_range_t;

  Accordion (size_t width, size_t overlap) :
    width(width), 
    overlap(overlap),
    buff(width),
    i(width)
    {
    };

  void clear() {
    i = width;
    buff.clear();
  };

  bool operator () (TYPE d) {
    // process a sample and return true each time the appropriate 'width'
    // samples are available

    buff.push_back(d);
    if (--i == 0) {
      i = width - overlap;
      return true;
    }
    return false;
  };

  array_range_t array_one () {
    // return the first linear array of samples in the latest window
    // only valid immediately after operator() returns true
    return buff.array_one();
  };

  array_range_t array_two () {
    // return the second linear array of samples in the latest window
    // only valid immediately after operator() returns true
    return buff.array_two();
  };

  // the next pair of functions return the two segments corresponding
  // to any 'valid' samples in the buffer (i.e. samples which will 
  // eventually be output, before being overwritten)
  // (Remaining samples have been output the correct number of times
  // already, and will be overwritten.)

  array_range_t partial_array_one () {
    // return the first linear array of 'valid' samples in the latest window
    auto a1 = buff.array_one();
    auto a2 = buff.array_two();
    size_t n = width - i;
    
    if (n > a2.second) {
      size_t m = n - a2.second;
      // some valid samples are in first linear array
      return array_range_t (a1.first + a1.second - m, m);
    } else {
      // the whole array is valid, but there are more in array_two
      return array_range_t (0, 0);
    }
  };

  array_range_t partial_array_two () {
    // return the second linear array of 'valid' samples in the latest window
    auto a2 = buff.array_two();
    size_t n = width - i;
    if (n == 0 || a2.second == 0)
      return array_range_t (0, 0);

    if (a2.second > n) {
      size_t m = a2.second - n;
      // some samples in this linear array are not valid
      return array_range_t (a2.first + m, a2.second - m);
    } else {
      // this whole linear array is valid
      return a2;
    }
  };

 protected:
  size_t width;
  size_t overlap;
  boost::circular_buffer < TYPE > buff;
  size_t i;
};

#endif //  _ACCORDION_H
