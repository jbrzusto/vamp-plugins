/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  SpectralPulseFinder.h - find fixed-width pulses in the power
  spectrum of a stream of complex values.  Pulses are distinguished by
  having SNR exceeding a threshold, or S minus N having high
  Z score (using sd(N) as the denominator)

  Copyright 2014 John Brzustowski

  License: GPL v 2.0 or later.

*/

#ifndef _SPECTRAL_PULSE_FINDER_H_
#define _SPECTRAL_PULSE_FINDER_H_

#include <stdint.h>
#include <vector>
#include <map>
#include <complex>

#include "vamp-plugins-common.h"
#include "SlidingSpectrum.h"
#include "FixedPulseDetector.h"

// FIXME: hardcoded float sample type
// NOTE: bin numbers are relative to the padded fft.

class SpectralPulseFinder {

public:

    typedef std::vector < int > :: const_iterator bin_iter_t;

    // constructor.
    // fft_size: number of samples per FFT
    // pad: number of zeroes to pad FFT with (for higher precision peak frequency estimate)
    // overlap: number of (non-pad) samples to overlap between consecutive fft
    // minbin, maxbin: indices of FFT bins in which to search for pulses; negative values
    //     mean bins above win_size / 2; e.g. -1 is bin win_size - 1.
    //     so the range minbin = -3, maxbin = 3 is permitted

    SpectralPulseFinder (int width,
                         int bkgd,
                         int win_size, 
                         int pad, 
                         int overlap,
                         int minbin,
                         int maxbin,
                         double min_SNR_dB,
                         double min_z,
                         double max_noise_for_Z) : 
        width (width),
        bkgd (bkgd),
        win_size (win_size),
        pad (pad), 
        overlap (overlap),
        minbin (minbin),
        maxbin (maxbin),
        numbins (maxbin - minbin + 1),
        min_SNR_dB (min_SNR_dB),
        min_z (min_z),
        max_noise_for_Z (max_noise_for_Z),
        ss(win_size, pad, overlap),
        pd()
    {
        // Note: pulse width is specified in the time domain, but pulses
        // are detected in the sliding spectrum.  In the latter, the 'width'
        // of a pulse being sought corresponds to the number of FFTs
        // that will cover the samples in the time-domain pulse, taking into
        // account window overlap.  Padding does not need to be considered here.

        width_in_spectrum = ceil(1.0 + width / (double) (win_size - overlap));

        bkgd_in_spectrum = ceil(1.0 + bkgd / (double) (win_size - overlap));

        for (int i = 0; i < numbins; ++i) 
            pd.push_back(FixedPulseDetector < float > ( width_in_spectrum,
                                                        bkgd_in_spectrum,
                                                        exp10(min_SNR_dB / 10.0), // min_SNR in linear units
                                                        min_z,
                                                        max_noise_for_Z
                                                        ));
    };

    ~SpectralPulseFinder() {
    };

    int ith_bin(int i) {
        // return the non-zero index of the ith bin of interest
        i += minbin;
        if (i >= 0) 
            return i;
        else 
            return i + win_size * pad;
    }
            
    double Power(const std::complex < float >  &d) {
        return d.real()*d.real() + d.imag() * d.imag();
    };

    bool operator () (const std::complex < float >  &d) {
        // process a sample; return true if a pulse has been found

        bestbin = -1;
        float bestsig = -1;

        if (ss(d)) {
            // a new fourier spectrum is available; push power values for each
            // bin of interest into the appropriate pulse detector


            for (int i = 0; i < numbins; ++i) {
                int ii = ith_bin(i);
                if ( pd[i](std::abs(Power(ss[ii])))) {
                    if (bestbin == -1 || pd[i].signal() > bestsig) {
                        bestbin = i;
                        bestsig = pd[i].signal();
                    }
                }
            }
        }
        return bestbin >= 0;
    };

    int pulsebin () {
        // return the index (relative to minbin) of the bin with the pulse
        // only valid immediately after operator() returns true
        return bestbin;
    };

    double signal (unsigned i) {
        // return the signal of the pulse in the i'th (relative to minbin) FFT bin
        return pd[i].signal();
    };

    double noise (unsigned i) {
        // return the noise of the pulse in the i'th FFT bin
        return pd[i].bkgd();
    };

    double SNR (unsigned i) {
        // return the SNR for the pulse in the ith FFT bin
        return pd[i].SNR();
    };

    double Z (unsigned i) {
        // return the Z score for the pulse in the ith FFT bin
        return pd[i].Z();
    };

    size_t location () {
        // how many samples back from most recently processed sample is first sample in
        // detected pulse?  
        // only valid immediately after a call to operator() returns true
        return (width_in_spectrum + 2 * bkgd_in_spectrum) * (win_size - overlap);
    };
       
    protected:

    size_t width; // samples per pulse
    size_t bkgd; // samples in the (half) background window
    int win_size; // samples per FFT (not including zero padding)
    int pad;     // number of zero samples to pad with (for tighter freq. estimation)
    int overlap; // overlap between consecutive FFT windows, in samples
    int minbin; // min bin we look for pulses in
    int maxbin; // max bin we look for pulses in
    int numbins; // number of bins we're looking in

    double min_SNR_dB; // minimum signal to noise ratio for a pulse
    double min_z; // maximum probability for a pulse
    double max_noise_for_Z; // maximum noise level at which z score is effective

    int width_in_spectrum; // width of pulse, in numbers of FFTs
    int bkgd_in_spectrum; // width of background, in numbers of FFTs

    SlidingSpectrum ss;

    std::vector < FixedPulseDetector < float > > pd; // vector of pulse detectors, one per bin in the range minbin..maxbin
    
    int bestbin;
};

#endif // _SPECTRAL_PULSE_FINDER_H_
