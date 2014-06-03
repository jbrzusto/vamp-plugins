/* -*- mode: c++ -*- */
/*
    Copyright 2014 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _FIXED_PULSE_DETECTOR_H
#define _FIXED_PULSE_DETECTOR_H

#include "EdgeDetector.h"

// FixedPulseDetector: find pulses of fixed width in a data stream.
// a pulse is a block of 'width' consecutive samples whose mean is
// significantly (in magnitude or probability) higher than the mean
// of the immediately adjacent blocks of the same width.

template < typename TYPE >

class FixedPulseDetector {
 public:

  typedef typename boost::circular_buffer < TYPE > ::const_array_range array_range_t;
  
  FixedPulseDetector(size_t width=0, TYPE min_diff=0, double max_prob=0) :
    width(min_width),
    min_diff(min_diff),
    max_prob(max_prob),
    buff(3 * min_width), 
    ej(min_width, min_diff, max_prob),
    awaiting_fall(false),
    last_location(0),
    last_width(0)
      {
      };
  
  bool operator() (TYPE d) { 
    // process a value from the data stream; return true if a pulse was detected

    // add to ring buffer
    buff.push_back(d);

    // if waiting for a falling edge, bump up the length
    if (awaiting_fall)
      ++samples_since_rise;

    // check for an edge
    if (ej(d)) {
      bool found_rise = ej.rising();
      if (found_rise) {
          // found a rising edge; cancel any started pulse and start a new one
        awaiting_fall = true;
        samples_since_rise = 0;
        left_bg = ej.left();
      } else if (awaiting_fall) {
        // we found the falling edge we were waiting for
        awaiting_fall = false;
        if (samples_since_rise <= max_width && samples_since_rise >= min_width) {
          last_width = samples_since_rise;
          last_location = last_width + 2 * min_width - 1;
          right_bg = ej.right();

          // set up (pointer, length) pairs to 1 or 2 segments in
          // the ring buffer containing the samples for this pulse
          
          array_range_t a1 = buff.array_one();
          size_t first = buff.size() - last_location - 1;
          size_t len_in_first;
          if (first < a1.second) {
            // beginning of pulse is in ring buffer's first array
            len_in_first = std::min(a1.second - first, last_width);
            array1 = array_range_t (a1.first + first, len_in_first);
          } else {
            // no part of pulse is in ring buffer's first array
            len_in_first = 0;
            array1 = array_range_t (0, 0);
          }
          size_t len_in_second = last_width - len_in_first;
          if (len_in_second > 0) {
            // there are some pulse samples in the ring buffer's second array
            array_range_t a2 = buff.array_two();
            if (len_in_second > a2.second)
              std::cerr << "oops!  mis-aligned frame" << std::endl;
            array2 = array_range_t(a2.first, len_in_second);
          } else {
            array2 = array_range_t(0, 0);
          }
          return true;
        }
      } else {
        // we found a falling edge but had not seen a rising edge, so ignore it
      }
    }
    return false;
  };

  size_t location () {
    // how many samples back from most recently processed sample is first sample in
    // detected pulse?  
    // only valid immediately after a call to operator() returns true
    return last_location;
  };

  size_t width () {
    // how many samples are in the pulse?
    return last_width;
  };

  TYPE bg () {
    // background value; mean of means in two background windows on either
    // side of pulse.  Valid only immediately after a call to operator() returns true.
    return (left_bg + right_bg) / 2.0;
  };

  typename boost::circular_buffer < TYPE > :: const_array_range pulse_array1 () {
    return array1;
  };

  typename boost::circular_buffer < TYPE > :: const_array_range pulse_array2 () {
    return array2;
  };

 protected:
  size_t min_width;  // minimum number of samples in pulse
  size_t max_width; // maximum number of samples in pulse
  TYPE min_diff; // minimum difference between averages of min_width samples on either side of pulse edge
  double max_prob; // maximum probability of difference between averages of min_width samples on either side of pulse edge
  boost::circular_buffer < TYPE > buff; // buffer of samples large enough to contain longest pulse when it has just been detected
  EdgeDetector < TYPE > ej; // edge detector 
  bool awaiting_fall; // true iff a rising edge has been detected but no falling edge has been seen since and max_width samples have not elapsed since detection
  size_t samples_since_rise; // number of samples elapsed since a rising edge was detected
  size_t last_location;
  size_t last_width;
  TYPE left_bg; // value of mean in left background window
  TYPE right_bg; // value of mean in right background window
  array_range_t array1;
  array_range_t array2;
};

#endif //  _FIXED_PULSE_DETECTOR_H
