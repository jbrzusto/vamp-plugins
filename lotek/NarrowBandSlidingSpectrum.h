/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  NarrowBandSlidingSpectrum.h - calculate a small portion of the
  Fourier spectrum of an complex (I/Q) stream.  No windowing, and
  an estimate of power in each bin is available after each new sample.
  Uses direct accumulation, rather than libfftw.  Also maintains
  a 'total' bin, which is simply the sum of power in all bins.

  Note: we use std::complex < float > as the hardcoded sample type.

  Copyright 2014 John Brzustowski

  License: GPL v 2.0 or later.

*/

#ifndef _NARROW_BAND_SLIDING_SPECTRUM_H_
#define _NARROW_BAND_SLIDING_SPECTRUM_H_

#include <stdint.h>
#include <vector>
#include <map>
#include <complex>
#include <boost/circular_buffer.hpp>
#include "KahanMovingAverager.h"
#include "vamp-plugins-common.h"

// FIXME: maybe use this when working: #include "MovingAveragerWithRecalc.h"

// FIXME: hardcoded float sample type

class NarrowBandSlidingSpectrum {

public:

    // constructor.
    // win_size: number of samples over which DFT estimate is computed
    //           an estimate of DFT for each bin and each sample is computed; win_size
    //           gives the number of consecutive samples used in each estimate.
    //           This is the same as using an overlap of n-1 in the usual size-n FFT.
    //
    // freq_low: frequency of lowest bin, in cycles per sample
    // freq_step: frequency increment between bins, in cycles per sample
    // num_bins: total number of DFT bins with frequencies:
    //           freq_low, freq_low + freq_step, freq_low + 2 * freq_step, ..., freq_low + (num_bins - 1) * freq_step

    NarrowBandSlidingSpectrum (int win_size, float freq_low, float freq_step, int num_bins) :
        win_size (win_size),
        freq_low (freq_low),
        freq_step (freq_step),
        num_bins (num_bins),
        phase (0)
    {
        for (int i = 0; i < num_bins; ++i) {
            //            ma_bin.push_back (MovingAveragerWithRecalc < std::complex < float > > (win_size));
            ma_bin.push_back (KahanMovingAverager < std::complex < float >, std::complex < double > > (win_size));
        }
        ma_power = KahanMovingAverager < float, double > (win_size);

        generateUnits();

    };

    ~NarrowBandSlidingSpectrum() {
        for (int i = 0; i < num_bins; ++i)
            delete [] units[i];
        delete [] units;
    };

    void generateUnits() {
        units = new std::complex < float > * [num_bins];
        for (int i = 0; i < num_bins; ++i) {
            units[i] = new std::complex < float > [win_size];

            for (int j = 0; j < win_size; ++j)
                // Note: do exp calculation in double
                units[i][j] = exp (std::complex < double > (0, - 2.0 * M_PI * (freq_low + i * (double) freq_step) * j / win_size));
        }
    };
        
    bool operator () (const std::complex < float >  &d) {
        // process a sample; return true if estimates of the spectral components are available
        float total = 0.0;

        for (int i = 0; i < num_bins; ++i) {
            // multiply by complex unit for bin i, given phase
            std::complex < float > bin_component = d * units[i][phase];
            ma_bin[i](bin_component);
            total += Pwr(bin_component);
        }
        // add total component; this is used to find the pulse

        bool have_est = ma_power(total);
        if (++phase == win_size)
            phase = 0;

        return have_est;
    };

    operator float () {
        // return the total power across all bins
        // only valid if operator() has returned true
        return ma_power.sum();
    };

    std::complex < double > operator [] (unsigned i) {
        // return the component from the i'th bin.
        // only valid if operator() has returned true
        return ma_bin[i].sum();
    };

protected:

    int win_size;
    float freq_low;
    float freq_step;
    int num_bins;

    std::vector < KahanMovingAverager < std::complex < float >, std::complex < double> > > ma_bin; // moving average of component in each bin
    KahanMovingAverager < float, double > ma_power; // moving average of power across all bins

    std::complex < float > ** units;  // elementary functions for required frequencies

    int phase; // phase for current sample; selects which element of each unit function is being used

};

#endif // _NARROW_BAND_SLIDING_SPECTRUM_H_
