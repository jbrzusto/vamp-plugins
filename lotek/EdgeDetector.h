/* -*- mode: c++ -*- */
/*
    Copyright 2014 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _EDGE_DETECTOR_H
#define _EDGE_DETECTOR_H

#include <boost/circular_buffer.hpp>
#include "MovingAverager.h"
#include "DynamicQuantile.h"
#include "PeakFinder.h"

using namespace boost::accumulators;

// EdgeDetector:  find edges in a data stream.  User species the window
// size, minimum difference, and probability threshold.  An edge is detected
// where the difference between two adjacent half-windows is larger, in
// absolute value, than the minimum difference, or lower, in probability
// than the probability threshold, and where the difference is a 
// local maximum over twice the window size.

typedef accumulator_set < float, stats < tag::p_square_quantile > > quant_est;

template < typename TYPE >

class EdgeDetector {
 public:
  
 EdgeDetector(size_t win_size, TYPE min_diff, double max_prob) :
  win_size(win_size),
    min_diff(min_diff),
    max_prob(max_prob),
    ma_count(win_size),
    ma_count_from(win_size),
    win(win_size),
    ma_buff(win_size + 1),
    dq(1.0 - max_prob, 10000),
    pk (win_size, true, true)
      {
      };
  
  bool operator() (TYPE d) { 
    // process a value from the data stream; return true if an edge was detected
    win.process(d);

    if (win.have_average()) {
      // add this average to our buffer
      double ma_right = win.get_average();
      ma_buff.push_back(ma_right);

      if (ma_buff.full()) {
        // we have enough samples to have already buffered the average
        // for the left window, so subtract from the right window's average
        double ma_left = ma_buff[0];
        double diff = ma_right - ma_left;
        double posdiff = fabs(diff);

        if (--ma_count == 0) {
          // add this difference to the probability distribution
          // we don't do this every frame to reduce processing requirements;
          // in any case, the difference changes little between adjacent frames
          
          dq(posdiff);
          
          // reset the counter
          ma_count = ma_count_from;
        }
        
        // send the difference through the peak finder
        if (pk.process(diff)) {
          // a local min/max in the between-window difference was found.
          // check whether it is large enough (in magnitude or probability)
          // to count as an edge
          if (posdiff >= dq || posdiff >= min_diff) {
            // record the quantile
            last_quant = dq;
            // record the average from the left window
            left_ave = ma_left;
            // note whether this was a rising or falling edge
            was_rising = diff > 0;
            return true;
          }
        }
      }
    }
    return false;
  };

  TYPE diff () {
    // return the mean in the right window minus the mean in the left window
    // only valid immediately after a call to operator() returns true
    return right() - left();
  };

  bool full() {
    // return true if we've seen sufficient samples to find an edge
    return pk.full();
  };

  TYPE left () {
    // return the mean in the left window
    // only valid immediately after a call to operator() returns true
    return left_ave;
  };

  TYPE right () {
    // return the mean in the right window
    return ma_buff[win_size];
  };

  bool rising() {
    // return true if the last edge was rising; otherwise false
    return was_rising;
  };

  double quantile() {
    // return quantile value when last edge was detected
    return last_quant;
  };

 protected:
  size_t win_size;
  TYPE min_diff;
  double max_prob;
  size_t ma_count;
  size_t ma_count_from;
  double left_ave; // moving average of left window for last peak
  bool was_rising; // direction of last detected edge
  double last_quant; // quantile at last detected edge

  MovingAverager < TYPE, double > win; // moving average of most recent win_size samples
  boost::circular_buffer < TYPE > ma_buff; // buffer of recent moving averages, so we need only compute the average of the right window, and lookup the average of the left window, which was computed originally when those samples were in the right window
  DynamicQuantile < double > dq; // maintain dynamic quantile for absolute values of window difference
  PeakFinder < double > pk; // find peaks in the difference between right and left windows
  
};

#endif //  _EDGE_DETECTOR_H
