/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  SpectralPulseFinder.h - find fixed-width pulses in the power
  spectrum of a stream of complex values.  Pulses are distinguished by
  having SNR exceeding a threshold, or S minus N having high
  Z score (using SE(N) as the denominator)

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
#include "VectorRingBuffer.h"

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
        numbins (pad * win_size),
        overlap (overlap),
        step (win_size - overlap),
        avoid (win_size / step),
        minbin (minbin),
        maxbin (maxbin),
        numwatchbin (maxbin - minbin + 1),
        min_SNR_dB (min_SNR_dB),
        min_z (min_z),
        max_noise_for_Z (max_noise_for_Z),
        width_in_spectrum(std::max(1.0, floor((width - win_size) / (double) step))),
        bkgd_in_spectrum (std::max(1.0, ceil(1.0 + bkgd / (double) step))),
        hist_size (2 * bkgd_in_spectrum + avoid + width_in_spectrum + 1),
        ss(win_size, (pad - 1) * win_size, overlap),
        pd()
    {
    // Note: pulse width is specified in the time domain, but pulses
        // are detected in the sliding spectrum.  In the latter, the 'width'
        // of a pulse being sought corresponds to the number of FFTs
        // that will cover the samples in the time-domain pulse, taking into
        // account window overlap.  Padding does not need to be considered here.

        for (int i = 0; i < numbins; ++i) 
            pd.push_back(FixedPulseDetector < float > ( width_in_spectrum,
                                                        bkgd_in_spectrum,
                                                        avoid,
                                                        exp10(min_SNR_dB / 10.0), // min_SNR in linear units
                                                        min_z,
                                                        max_noise_for_Z * (win_size) // scale for FFT window size
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
            return i + numbins;
    }
            
    bool operator () (const std::complex < float >  &d) {
        // process a sample; return true if a pulse has been found

        bestbin = -1;
        float bestsig = -1;

        if (ss(d)) {
            // a new fourier spectrum is available; calculate power
            // values for bins of interest, buffer them, get best power, see
            // if a peak has been found in best power

            std::complex < float > * fft_out = ss;

            for (int i = 0; i < numwatchbin; ++i) {
                int ii = ith_bin(i);
                float p = Pwr(fft_out[ii]);
                if (pd[i](p)) {
                    float sigdif = pd[i].sigdif();
                    if (sigdif > bestsig || bestbin < 0) {
                        bestbin = i;
                        bestsig = sigdif;
                    }
                }
            };
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
        return (hist_size - 1) * step;
    };
       
    protected:

    size_t width; // samples per pulse
    size_t bkgd; // samples in the (half) background window
    int win_size; // samples per FFT (not including zero padding)
    int pad;     // number of zero samples to pad with (for tighter freq. estimation)
    int numbins; // real number of FFT bins (i.e. including padding)
    int overlap; // overlap between consecutive FFT windows, in samples
    int step; // win_size - overlap; number of new samples per FFT
    int avoid; // number of FFT steps to skip on either side of signal window to avoid leakage into background window
    int minbin; // min bin we look for pulses in
    int maxbin; // max bin we look for pulses in
    int numwatchbin; // number of bins we're looking in

    double min_SNR_dB; // minimum signal to noise ratio for a pulse
    double min_z; // maximum probability for a pulse
    double max_noise_for_Z; // maximum noise level at which z score is effective

    int width_in_spectrum; // width of pulse, in numbers of FFTs
    int bkgd_in_spectrum; // width of background, in numbers of FFTs

    int hist_size; // number of bin power time slots we need to keep to
    // use bin power once pulse has been detected

    SlidingSpectrum ss;

    std::vector < FixedPulseDetector < float > > pd; // pulse detectors, to find pulse in values of largest power across bins

    int bestbin;

};

#endif // _SPECTRAL_PULSE_FINDER_H_
