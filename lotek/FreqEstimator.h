/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    FreqEstimator.h - estimate main frequency of a pulse from raw audio data.  FIXME: hardcoded S16_LE

    Copyright 2014 John Brzustowski

    License: GPL v 2.0 or later.

*/

#ifndef _FREQ_ESTIMATOR_H_
#define _FREQ_ESTIMATOR_H_

#include <fftw3.h>
#include <stdint.h>
#include <vector>
#include <map>


// whenever an estimator is requested, we see whether we already
// have a plan for that number of samples and channels, so that we 
// don't need to reallocate storage.  This is NOT thread safe!

struct cachedFFTWPlan {
  fftwf_complex *m_input; // for 1 channel: samples stored as floats; for 2 channels: floating point samples, as interleaved I/C pairs
  fftwf_complex *m_output; // DFT output from pulse samples
  fftwf_plan m_plan; // FFT plan for I/Q sample
};

class FreqEstimator {

 public:

  // constructor.  max_frames is the largest number of frames of data
  // allowed (typically made large to allow zero-padding for better
  // frequency resolution) channels is the number of channels, 1 or 2

  FreqEstimator (unsigned max_frames, unsigned channels=2);

  ~FreqEstimator();

  // estimate the frequency with maximum power, and (if pwr_out != NULL)
  // return an estimate of its power in *pwr_out
  // the estimate is in cycles per sample, so the caller has to 
  // multiply by the sampling rate to get frequency in Hz.
  float get (int16_t *samples, unsigned frames, float *pwr_out = 0);

 protected:

  unsigned m_max_frames;
  unsigned m_channels;
  unsigned m_fft_bins;

  unsigned m_first_zero_frame; // index of first frame in m_ff
  float m_freq_scale;  // amount to multiply a bin index to get cycles per sample
  static const char * m_fftw_wisdom_filename;
  static bool m_fftw_wisdom_loaded;

  fftwf_complex *m_input; // for 1 channel: samples stored as floats; for 2 channels: floating point samples, as interleaved I/C pairs
  fftwf_complex *m_output; // DFT output from pulse samples
  fftwf_plan m_plan; // FFT plan for I/Q sample

  static void generateWindowingCoefficients(int N, std::vector < float > &window, float &win_sum, float &win_sumsq);

  static std::map < std::pair < int, short >, cachedFFTWPlan > m_cached_plans;
};


#endif // _FREQ_ESTIMATOR_H_
