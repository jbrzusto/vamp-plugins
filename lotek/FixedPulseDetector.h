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
#include "PeakFinder.h"

template < typename TYPE >

class FixedPulseDetector {
public:
  
  FixedPulseDetector(size_t width, TYPE min_SNR, double min_z) :
    width(width),
    z_scale(sqrt((2.0 * width - 1) / (2.0 * width))),
    min_SNR(min_SNR),
    min_z(min_z / z_scale),  // pre-scale to simplify SD calculation lator
    min_needed(2 * width + 1),
    ma_count(width),
    ma_count_from(width),
    win(width),
    win2(width),
    ma_buff(3 * width + 1),
    ma_buff2(3 * width + 1),
    pk (width, true, false) // only look for local maxima
  {
  };
  
  bool operator() (const TYPE & d) { 
    // process a value from the data stream; return true if a pulse was detected
    win(d);
    TYPE d2 = d * d;
    if (win2(d2)) {
      // add moving averages of samples and squares to the buffer
      double ma_right = win;
      ma_buff.push_back(ma_right);
      double ma_right2 = win2;
      ma_buff2.push_back(ma_right2);

      if (ma_buff.size() >= min_needed) {
        // we have enough moving window averages buffered so that the left
        // background window is the first ma in the buffer, and the middle
        // (signal) window is also there.

        double ma_left = ma_buff[ma_buff.size() - 2 * width - 1];
        double ma_sig = ma_buff[ma_buff.size() - width - 1];
        double diff = ma_sig - (ma_right + ma_left) / 2.0;

        // send the difference through the peak finder
        if (pk(diff)) {
          // a local max in the between-window difference was found.
          // Its signal and noise window values are not those calculated
          // immediately above, but were buffered some time ago.

          // Check whether it is large enough (in magnitude or
          // probability) to count as a pulse

          was_big = SNR() >= min_SNR;

          // do the more expensive Z-score calculation
          // the formula for bksd below would need to be multiplied
          // by sqrt((N-1)/N) to be correct, but we have
          // already divided min_z by that factor so that
          // the threshold comparison is correct.
          
          double pkbkgd = bkgd();
          double bksd = sqrt(bkgd2() - pkbkgd*pkbkgd); // * sqrt((2*width-1)/(2*width))
          last_z = (signal() - pkbkgd) / bksd;

          was_unlikely = last_z >= min_z;

          return was_unlikely || was_big;
        }
      }
    }
    return false;
  };

  TYPE signal () {
    // return the mean of the values in the central (signal) window
    // only valid immediately after a call to operator() returns true
    return ma_buff[width];
  };

  TYPE bkgd () {
    // return the mean of the values in the left and right background windows
    // only valid immediately after a call to operator() returns true
    return (ma_buff[0] + ma_buff[2 * width]) / 2.0;
  };

  TYPE bkgd2 () {
    // return the mean of the squares of values in left and right background windows
    // only valid immediately after a call to operator() returns true
    return (ma_buff2[0] + ma_buff2[2 * width]) / 2.0;
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
    // return true if the pulse had z-score larger than the min_z; only valid
    // immediately after operator() returns true; and only valid if big()
    // returned false
    return was_unlikely;
  };

  size_t location () {
    // how many samples back from most recently processed sample is first sample in
    // detected pulse?  
    // only valid immediately after a call to operator() returns true
    return 3 * width - 1;
  };

  double Z() {
    // return z score for the most recently detected pulse
    return last_z * z_scale;
  };

protected:
  size_t width;
  double z_scale; // sqrt((2 * width) / (2 * width - 1))
  TYPE min_SNR;  // minimum signal to noise ratio in linear units
  double min_z;
  size_t min_needed;
  size_t ma_count;
  size_t ma_count_from;
  double last_z; // quantile at last detected edge
  bool was_big; // true if difference between signal and bkgd was signficant
  bool was_unlikely; // true if the difference between signal and bkgd was low probability
 
  MovingAverager < TYPE, double > win; // moving average of most recent width samples
  MovingAverager < TYPE, double > win2; // moving average of squares of most recent width samples
  boost::circular_buffer < TYPE > ma_buff; // buffer of recent moving averages of samples values, so we need only compute the average of the right window, and lookup the average of the left window, which was computed originally when those samples were in the right window
  boost::circular_buffer < TYPE > ma_buff2; // buffer of recent moving averages of samples squares
  PeakFinder < double > pk; // find peaks in the difference between right and left windows
  
};

#endif //  _FIXED_PULSE_DETECTOR_H
