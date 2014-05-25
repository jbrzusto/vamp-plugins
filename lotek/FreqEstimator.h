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

    FreqEstimator (int max_frames, int channels=2);

    ~FreqEstimator();

    // estimate the frequency with maximum power, and (if pwr_out != NULL)
    // return an estimate of its power in *pwr_out
    // the estimate is in cycles per sample, so the caller has to 
    // multiply by the sampling rate to get frequency in Hz.
    float get (int16_t *samples, int frames, float *pwr_out = 0);

    void generateWindowingCoefficients();

protected:

    int m_max_frames;
    int m_channels;
    int m_fft_bins;

    int m_first_zero_frame; // index of first frame in m_ff
    float m_freq_scale;  // amount to multiply a bin index to get cycles per sample
    static const char * m_fftw_wisdom_filename;
    static bool m_fftw_wisdom_loaded;

    fftwf_complex *m_input; // for 1 channel: samples stored as floats; for 2 channels: floating point samples, as interleaved I/C pairs
    fftwf_complex *m_output; // DFT output from pulse samples
    fftwf_plan m_plan; // FFT plan for I/Q sample

    bool m_have_window;
    std::vector < float > m_window;
    double m_win_sum;
    double m_win_sumsq;
    double m_power_scale;

    float binRatio(float d); // ratio of spectral power in bin x+d to true power for the function exp(2 * pi * %i * x * t)
    float twoBinRatio(float d); // ratio of spectral power in bin x+d to bin x+d+1 for the function exp(2 * pi * %i * x * t)
    float estimateBinOffset(int bin); // given 'bin' has max power, estimate the offset to the true frequency using adjacent bins, and the ratio of bin power to true power
    std::vector < float > m_ratio2_to_offset; // lookup table that inverts twoBinRatio function
    std::vector < float > m_ratio3_to_offset; // lookup table that inverts twoBinRatio function

    double m_ratioMapScale2;
    double m_ratioMapScale3;

    void calculateRatioMaps(float prec = 0.00001); // set up ratio to offset lookup tables

    static std::map < std::pair < int, short >, cachedFFTWPlan > m_cached_plans;
};


#endif // _FREQ_ESTIMATOR_H_
