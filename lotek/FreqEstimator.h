/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  FreqEstimator.h - estimate main frequency of a fixed-width pulse from complex I/Q data.

  Copyright 2014 John Brzustowski

  License: GPL v 2.0 or later.

*/

#ifndef _FREQ_ESTIMATOR_H_
#define _FREQ_ESTIMATOR_H_

#include <fftw3.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <complex>

class FreqEstimator {

public:

    // constructor. 

    // frames it the width of the pulse, in samples

    FreqEstimator (int frames);

    ~FreqEstimator();

    // estimate the frequency with maximum power
    // the estimate is in cycles per sample, so the caller has to 
    // multiply by the sampling rate to get frequency in Hz.

    float get (std::complex < float > * seg1, int n1, std::complex < float > * seg2, int n2);

protected:

    int m_frames;

    std::complex < float > *m_input; // for 1 channel: samples stored as floats; for 2 channels: floating point samples, as interleaved I/C pairs
    std::complex < float > *m_output; // DFT output from pulse samples
    fftwf_plan m_plan; // FFT plan for I/Q sample

    std::vector < float > m_window;

    double m_win_sum;
    double m_win_sumsq;

    void generateWindowingCoefficients();

    double estimateBinOffset(int bin); // given 'bin' has max power, estimate the offset to the true frequency using adjacent bins, and the ratio of bin power to true power
    double cubicMaximize(double y0, double y1, double y2, double y3); // estimate location of maximum (using cubic curve) given function values at x = 0, 1, 2, 3
};


#endif // _FREQ_ESTIMATOR_H_
