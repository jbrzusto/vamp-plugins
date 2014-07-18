/* -*- mode: c++ -*- */
/*

// FixedPulseDetector: find pulses of fixed width in a data stream.
// a pulse is a window of 'width' consecutive samples whose mean is
// significantly (in magnitude or probability) higher than the mean
// of the immediately adjacent windows of the same width.
// Moreover, the pulse must have the most significant such difference
// among adjacent locations to 'width' samples in each direction.

// The signal and background windows are allowed to be separated,
// to prevent partial leakage of energy from the pulse into the
// background due to imprecision in locating the pulse in the input
// domain (e.g. because the input domain is frequency vs. time
// with a stepping, overlapping FFT window).

    Copyright 2014 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _FIXED_PULSE_DETECTOR_H
#define _FIXED_PULSE_DETECTOR_H

#include <boost/circular_buffer.hpp>
#include "MovingSum.h"
#include "PeakFinder.h"
#include "vamp-plugins-common.h"

template < typename TYPE >

class FixedPulseDetector {
public:
  
  FixedPulseDetector(size_t width, size_t bkgd_width, size_t avoid_width, TYPE min_SNR, double min_z, double max_noise_for_Z) :
    width(width),
    bkgd_width(bkgd_width),
    avoid_width(avoid_width),
    left_offset(width + bkgd_width + 2 * avoid_width + 1),
    sig_offset(bkgd_width + avoid_width + 1),
    z_scale(sqrt(2.0 * bkgd_width - 1.0)),
    min_SNR(min_SNR),
    min_z(min_z / z_scale),  // pre-scale to simplify SD calculation later
    max_noise_for_Z(max_noise_for_Z),
    min_needed(width + 2 * avoid_width + bkgd_width + 1),
    signal_scale(1.0 / width),
    bkgd_scale(0.5 / width),
    winma(width),
    bkgdma(bkgd_width),
    bkgd2ma(bkgd_width),
    buf_size(2 * bkgd_width + 2 * avoid_width + width + 1),
    ms_win_buf(buf_size),
    ms_bkgd_buf(buf_size),
    ms_bkgd2_buf(buf_size),
    pk (2*width - 1)
  {
  };
  
  bool operator() (const TYPE & d) { 
    // process a value from the data stream; return true if a pulse was detected
    winma(d);
    bkgdma(d);
    TYPE d2 = d * d;
    if (bkgd2ma(d2)) {
      // add moving averages of samples and squares to the buffer
      ms_win_buf.push_back(winma);
      double ms_right = bkgdma;
      ms_bkgd_buf.push_back(ms_right);
      double ms_right2 = bkgd2ma;
      ms_bkgd2_buf.push_back(ms_right2);

      if (ms_win_buf.size() >= min_needed) {
        // we have enough moving window averages buffered so that the left
        // background window is the first ma in the buffer, and the middle
        // (signal) window is also there.

        double ms_left = ms_bkgd_buf[ms_bkgd_buf.size() - left_offset ];
        double ms_sig = ms_win_buf[ms_win_buf.size() - sig_offset];
        //        double diff = ms_sig - (ms_right + ms_left) / 2.0;
        // detect only left edge:
        double diff = ms_sig - (ms_left);

        // send the difference through the peak finder
        if (pk(diff)) {
          // a local max in the between-window difference was found.
          // Its signal and noise window values are not those calculated
          // immediately above, but were buffered some time ago.

          // Check whether it is large enough (in magnitude or
          // probability) to count as a pulse

          was_big = SNR() >= min_SNR;

          double pkbkgd = bkgd();
          if (pkbkgd <= max_noise_for_Z) {
            //
            // low noise environment, so look for weaker pulses
            // by using the 'Z score' (signal - noise) / (SE noise)
            //
            // the formula for bkse below would need to be divided by
            // by sqrt(N-1) to be correct (and so last_z should be
            // multiplied by it) but we have already divided min_z by
            // that factor so that the threshold comparison is
            // correct.
          
            double bkse = sqrt((bkgd2() - pkbkgd * pkbkgd));
            last_z = (signal() - pkbkgd) / bkse;
            was_unlikely = last_z >= min_z;
          } else {
            last_z = 0;
            was_unlikely = false;
          }
          return was_unlikely || was_big;
        }
      }
    }
    return false;
  };

  TYPE sigdif () {
    return signal() - bkgd();
  };

  TYPE signal () {
    // return the mean of the values in the central (signal) window
    // only valid immediately after a call to operator() returns true
    return ms_win_buf[width + avoid_width] * signal_scale;
  };

  TYPE bkgd () {
    // return the mean of the values in the left and right background windows
    // only valid immediately after a call to operator() returns true
    return ms_bkgd_buf[0] * bkgd_scale;
  };

  TYPE bkgd2 () {
    // return the mean of the squares of values in left and right background windows
    // only valid immediately after a call to operator() returns true
    return (ms_bkgd2_buf[0] + ms_bkgd2_buf[width + bkgd_width + 2 * avoid_width]) * bkgd_scale;
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
    return 2 * bkgd_width + width + avoid_width - 1;
 };

  double Z() {
    // return z score for the most recently detected pulse
    return last_z * z_scale;
  };

protected:
  size_t width;
  size_t bkgd_width; // size of bkgd window on each side of pulse
  size_t avoid_width; // size of zone separating pulse window from background window on each side
  size_t left_offset; // how far back from the end of buffer is the left bkgd window value
  size_t sig_offset; // how far back from the end of the buffer is the signal value
  double z_scale; // multiple to convert simple z-score formula to correct formula!
  TYPE min_SNR;  // minimum signal to noise ratio in linear units
  double min_z;
  double max_noise_for_Z; // maximum noise at which a pulse can be accepted on Z score alone (i.e. ignoring SNR)
  size_t min_needed;
  double last_z; // quantile at last detected edge
  bool was_big; // true if difference between signal and bkgd was signficant
  bool was_unlikely; // true if the difference between signal and bkgd was low probability

  float signal_scale; // factor to multiply signal window sum to get mean
  float bkgd_scale; // factor to multiply bkgd window sum to get mean
 
  MovingSum < TYPE > winma; // moving average of most recent width samples
  MovingSum < TYPE > bkgdma; // moving average of most recent bkgd samples
  MovingSum < TYPE > bkgd2ma; // moving average of squares of most recent bkgd samples

  size_t buf_size; // length of buffers we need to retain for signal and noise moving averages
  boost::circular_buffer < TYPE > ms_win_buf; // buffer of recent moving averages of samples values, so we need only compute the average of the right window, and lookup the average of the left window, which was computed originally when those samples were in the right window
  boost::circular_buffer < TYPE > ms_bkgd_buf; // buffer of recent moving averages of samples in bkgd window
  boost::circular_buffer < TYPE > ms_bkgd2_buf; // buffer of recent moving averages of samples squares
  PeakFinder < double > pk; // find peaks in the difference between right and left windows
  
};

#endif //  _FIXED_PULSE_DETECTOR_H
