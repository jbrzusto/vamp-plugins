/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  NarrowBandNarrowBandSpectralPulseFinder.h - find fixed-width pulses in the power
  spectrum of a stream of complex values.  Pulses are distinguished by
  having SNR exceeding a threshold, or S minus N having high
  Z score (using sd(N) as the denominator)

  This version doesn't use FFT but rather does direct calculation of 
  a limited set of FFT bins using floating point accumulators.

  Copyright 2014 John Brzustowski

  License: GPL v 2.0 or later.

*/

#ifndef _NARROW_BAND_SPECTRAL_PULSE_FINDER_H_
#define _NARROW_BAND_SPECTRAL_PULSE_FINDER_H_

#include <stdint.h>
#include <vector>
#include <map>
#include <complex>

#include "vamp-plugins-common.h"
#include "NarrowBandSlidingSpectrum.h"
#include "PeakFinder.h"

// FIXME: hardcoded float sample type
// NOTE: bin numbers are relative to the padded fft.

class NarrowBandSpectralPulseFinder {

public:

    typedef std::vector < int > :: const_iterator bin_iter_t;

    // constructor.
    // width: pulse width, in samples
    // bin_low: lowest DFT bin of interest (i.e. cycles per pulse width)
    // bin_high: highest DFT bin of interest (i.e. cycles per pulse width)
    // min_SNR_dB: minimum signal to noise for a pulse
    // min_Z: minimum z-score for a pulse
    // max_noise_for_Z: max noise at which min_Z criterion is valid

    NarrowBandSpectralPulseFinder (int width,
                                   int bin_low,
                                   int bin_high,
                                   double min_SNR_dB,
                                   double min_Z,
                                   double max_noise_for_Z) : 
        width (width),
        bin_low (bin_low),
        bin_high (bin_high),
        num_bins (bin_high - bin_low + 1),
        min_SNR undB(min_SNR_dB),
        min_Z (min_Z),
        max_noise_for_Z (max_noise_for_Z),
        nbss(width, bin_low, 1, num_bins),
        total_pwr_buff (2 * width + 1),
        pk(width),
        sig_scale (1.0 / (width * width)),
        bkgd_scale (0.5 / width)
    {
        for (int i = 0; i < num_bins; ++i)
            bin_pwr_buff.push_back (boost::circular_buffer < float > (3 * width + 1));
    };

    ~NarrowBandSpectralPulseFinder() {
    };

    bool operator () (const std::complex < float >  &d) {
        // process a sample; return true if a pulse has been found

        if (nbss(d)) { // process this sample
            // bin estimates are available
            
            // buffer running power estimates for each bin
            for (int i = 0; i < num_bins; ++i) 
                bin_pwr_buff[i].push_back(Pwr(nbss[i]));

            total_pwr_buff.push_back(nbss);
            if (total_pwr_buff.full()) {
                // we have enough estimates of total power in bins
                // to compare central window to mean(left, right)
                float diff = total_pwr_buff[width] - (total_pwr_buff[0] + total_pwr_buff[2 * width]) / 2.0;
                if (pk(diff)) {
                    return process_pulse ();
                }
            }
        } 
        return false;
    };

    float signal () {
        // return the power of the pulse (best bin only)
        // only valid immediately after operator() returns true
        return signal(best_bin);
    };

    float signal (int i) {
        // return the power of the pulse (bin i)
        // only valid immediately after operator() returns true
        return bin_pwr_buff[i][width] * sig_scale;
    };

    float noise () {
        // return the bkgd noise of the pulse (best bin only)
        // only valid immediately after operator() returns true
        return noise(best_bin);
    };

    float noise (int i) {
        // return the bkgd noise of the pulse (bin i)
        // only valid immediately after operator() returns true
        float n = (bin_pwr_buff[i][0] + bin_pwr_buff[i][2 * width]) * bkgd_scale;
        if (n == 0)
            n = 2.5118864E-10; // -96 dB unlikely case, but we protect against it anyway
        return n;
    };

    double SNR () {
        // return the SNR for the pulse
        // only valid immediately after operator() returns true
        float n = noise();
        return (signal() - n) / n;
    };

    // FIXME: add Z functionality
    // double Z (unsigned i) {
    //     // return the Z score for the pulse in the ith FFT bin
    //     return pd.Z();
    // };

    size_t location () {
        // how many samples back from most recently processed sample is first sample in
        // detected pulse?  
        // only valid immediately after a call to operator() returns true
        return (5 * width + 1) / 2;
    };
       
    float freq () {
        // return the estimated frequency of the pulse, as a bin number

        // use a cubic smoother to get the finer offset from the highest
        // power bin and its neighbours.  The 4th neighbour is on
        // the same side of the highest power bin as the higher of the
        // two immediate neighbours.

        if (num_bins < 4)
            return best_bin;

        int first_bin = best_bin - 1;
        if (first_bin <= 0)
            // Note: '<=' so that the case first_bin==0 isn't available in
            // the third clause below
            first_bin = 0;
        else if (first_bin > num_bins - 4)
            first_bin = num_bins - 4;
        else if (smoothed_signal(first_bin) > smoothed_signal(first_bin + 2))
            ;//-- first_bin;

        float offset = cubicMaximize(smoothed_signal(first_bin),
                                     smoothed_signal(first_bin + 1),
                                     smoothed_signal(first_bin + 2),
                                     smoothed_signal(first_bin + 3));
        if (fabs(offset) < 3)
            return bin_low + first_bin + offset;
        else
            // smoothed offset is too large
            return best_bin;
    };

protected:
    float smoothed_signal (int i) {
        // return the power of the pulse in the i'th bin, averaged
        // across 3 timesteps: pulse centre, and +/- 1/4 pulse width;
        // this seems to behave similarly as explicitly windowing the
        // samples before doing an FFT (which we're not doing)

        float ss = 0.0;
        for (int j = -1; j <= 1; ++j) 
            ss += bin_pwr_buff[i][width + j * (signed) width / 4];
        return ss * sig_scale;
    };

    bool process_pulse () {
        // we have a candidate pulse - the difference between power in the 'signal' and
        // bkgd windows is at a local max.
        // return true if it should be accepted as a pulse

        // check whether the SNR criterion is met
        // we'll be accepting this pulse.
        // find bin with highest power

        best_bin = 0;
        best_bin_power = bin_pwr_buff[0][width];
        for (int i = 1; i < num_bins; ++i) {
            if (bin_pwr_buff[i][width] > best_bin_power) {
                best_bin_power = bin_pwr_buff[i][width];
                best_bin = i;
            }
        }

        if (SNR() < min_SNR)
            return false;

        return true;
    };

    static double cubicMaximize(double y0, double y1, double y2, double y3)
    {
        // Find coefficients of cubic

        double a, b, c;

        a = y0 / -6.0 + y1 / 2.0 - y2 / 2.0 + y3 / 6.0;

        if (a == 0.0)
            return double(-1000); // error

        b = y0 - 5.0 * y1 / 2.0 + 2.0 * y2 - y3 / 2.0;
        c = -11.0 * y0 / 6.0 + 3.0 * y1 - 3.0 * y2 / 2.0 + y3 / 3.0;

        // Take derivative

        double da, db, dc;

        da = 3 * a;
        db = 2 * b;
        dc = c;

        // Find zeroes of derivative using quadratic equation

        double discriminant = db * db - 4 * da * dc;
        if (discriminant < 0.0) {
            if (discriminant < -1.0)
                return double(-1000);              // error
            else
                discriminant = 0.0;
        }              
    
        double x1 = (-db + sqrt(discriminant)) / (2 * da);
        double x2 = (-db - sqrt(discriminant)) / (2 * da);

        // The one which corresponds to a local _maximum_ in the
        // cubic is the one we want - the one with a negative
        // second derivative

        double dda = 2 * da;
        double ddb = db;

        if (dda * x1 + ddb < 0)
            {
                return x1;
            }
        else
            {
                return x2;
            }
    };
        
    protected:

    size_t width; // samples per pulse
    int bin_low; // min bin we look for pulses in
    int bin_high; // max bin we look for pulses in
    int num_bins; // number of bins we're looking in

    double min_SNR; // minimum signal to noise ratio for a pulse
    double min_Z; // maximum probability for a pulse
    double max_noise_for_Z; // maximum noise level at which z score is effective

    NarrowBandSlidingSpectrum nbss;

    std::vector < boost::circular_buffer < float > > bin_pwr_buff; // buffer of power for bins, for retrieving values at pulse location

    boost::circular_buffer < float > total_pwr_buff;  // buffer of total power across bins

    PeakFinder < double > pk; // find peaks in the difference between power at time step i and mean power at time steps (i - width), (i + width)
    
    // parameters related to most recently found pulse
    int best_bin; // bin of highest power at peak
    float best_bin_power; // power in this best bin

    float sig_scale; // multiplier to get signal power 
    float bkgd_scale; // multipler to get noise power

};

#endif // _NARROW_BAND_SPECTRAL_PULSE_FINDER_H_
