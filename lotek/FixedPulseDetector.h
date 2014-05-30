/* -*- mode: c++ -*- */
/*

// FixedPulseDetector: find pulses of fixed width in a data stream.
// a pulse is a window of 'width' consecutive samples whose mean is
// significantly (in magnitude or probability) higher than the mean
// of the immediately adjacent windows of the same width.
// Moreover, the pulse must have the most significant such difference
// among adjacent locations to 'width' samples in each direction.

    Copyright 2014 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _FIXED_PULSE_DETECTOR_H
#define _FIXED_PULSE_DETECTOR_H

#include <boost/circular_buffer.hpp>
#include "MovingAverager.h"
#include "DynamicQuantile.h"
#include "PeakFinder.h"

using namespace boost::accumulators;

typedef accumulator_set < float, stats < tag::p_square_quantile > > quant_est;

template < typename TYPE >

class FixedPulseDetector {
public:
  
  FixedPulseDetector(size_t width, TYPE min_SNR, double max_prob) :
    width(width),
    min_SNR(min_SNR),
    max_prob(max_prob),
    min_needed(2 * width + 1),
    ma_count(width),
    ma_count_from(width),
    win(width),
    ma_buff(3 * width + 1),
    dq(1.0 - max_prob, 10000000),
    pk (width, true, false) // only look for local maxima
  {
  };
  
  bool operator() (const TYPE & d) { 
    // process a value from the data stream; return true if a pulse was detected
    if (win(d)) {
      // add this average to our buffer
      double ma_right = win;
      ma_buff.push_back(ma_right);
      if (ma_buff.size() >= min_needed) {
        // we have enough moving window averages buffered so that the left
        // background window is the first ma in the buffer, and the middle
        // (signal) window is also there.

        double ma_left = ma_buff[ma_buff.size() - 2 * width - 1];
        double ma_sig = ma_buff[ma_buff.size() - width - 1];
        double diff = ma_sig - (ma_right + ma_left) / 2.0;

        if (diff > 0) {
          // this might be a pulse: the signal window has a larger
          // mean than the background windows on both sides

          if (--ma_count == 0) {
            // add this difference to the probability distribution
            // we don't do this every frame to reduce processing requirements;
            // in any case, the difference changes little between adjacent frames
            
            dq(diff);
            
            // reset the counter
            ma_count = ma_count_from;
          }
        }
        
        // send the difference through the peak finder
        if (pk(diff)) {
          // a local max in the between-window difference was found.
          // Check whether it is large enough (in magnitude or
          // probability) to count as a pulse

          was_unlikely = ((double) pk >= dq);
          was_big = SNR() >= min_SNR;

          if (was_big || was_unlikely) {
            // record the quantile
            last_quant = dq;
            return true;
          }
        }
      }
    }
    return false;
  };

  TYPE signal () {
    // return the mean in the central (signal) window
    // only valid immediately after a call to operator() returns true
    return ma_buff[width];
  };

  TYPE bkgd () {
    // return the mean of the left and right background windows
    // only valid immediately after a call to operator() returns true
    return (ma_buff[0] + ma_buff[2 * width]) / 2.0;
  };

  double SNR () {
    // return the mean of the left and right background windows
    // only valid immediately after a call to operator() returns true
    double bkg = bkgd();
    if (bkg <= 0)
      bkg = 2.5118864E-10; // -96 dB unlikely case, but we protect against it anyway
    double sig = signal();
    if (sig < bkg)
      return 0.0;
    return (sig - bkg) / bkg; // or should we just do sig / bkg?
  };

  bool big () {
    // return true if the pulse had a mean difference above the minimum
    // only valid immediately after operator() returns true
    return was_big;
  };

  bool unlikely () {
    // return true if the pulse had a mean difference with lower probability
    // than max_prob
    return was_unlikely;
  };

  size_t location () {
    // how many samples back from most recently processed sample is first sample in
    // detected pulse?  
    // only valid immediately after a call to operator() returns true
    return 3 * width - 1;
  };

  double quantile() {
    // return quantile value when last edge was detected
    return last_quant;
  };

protected:
  size_t width;
  TYPE min_SNR;
  double max_prob;
  size_t min_needed;
  size_t ma_count;
  size_t ma_count_from;
  double last_quant; // quantile at last detected edge
  bool was_big; // true if difference between signal and bkgd was signficant
  bool was_unlikely; // true if the difference between signal and bkgd was low probability
 
  MovingAverager < TYPE, double > win; // moving average of most recent width samples
  boost::circular_buffer < TYPE > ma_buff; // buffer of recent moving averages, so we need only compute the average of the right window, and lookup the average of the left window, which was computed originally when those samples were in the right window
  DynamicQuantile < double > dq; // maintain dynamic quantile for absolute values of window difference
  PeakFinder < double > pk; // find peaks in the difference between right and left windows
  
};

#endif //  _FIXED_PULSE_DETECTOR_H
