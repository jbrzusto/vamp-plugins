/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  SlidingSpectrum.h - calculate the Fourier spectrum of an complex
  (I/Q) stream, with specified window size and overlap.

  Note: we use std::complex < float > as the hardcoded sample type,
  and use reinterpret_cast < fftwf_complex * > to interface with
  fftwf_..() functions.

  Copyright 2014 John Brzustowski

  License: GPL v 2.0 or later.

*/

#ifndef _SLIDING_SPECTRUM_H_
#define _SLIDING_SPECTRUM_H_

#include <fftw3.h>
#include <stdint.h>
#include <vector>
#include <map>
#include <complex>
#include <boost/circular_buffer.hpp>

#include "Accordion.h"

// FIXME: hardcoded float sample type

class SlidingSpectrum {

public:

    typedef enum {WINDOW_RECTANGULAR = 0, WINDOW_HAMMING = 1} WinType;

    // constructor.
    // fft_size: number of samples per FFT
    // pad: number of zeroes to pad FFT with (for higher precision peak frequency estimate)
    // overlap: number of (non-pad) samples to overlap between consecutive fft
    // wintype: type of windowing to apply to samples

    SlidingSpectrum (int win_size, int pad, int overlap, WinType win_type = WINDOW_HAMMING) :
        win_size (win_size),
        pad (pad), 
        overlap (overlap),
        win_type (win_type),
        fft_bins (win_size + pad),
        acc (win_size, overlap),
        have_window (false)
    {
        if (! fftw_wisdom_loaded) {
            // silently fail if wisdom cannot be found
            FILE *f = fopen(fftw_wisdom_filename.c_str(), "r");
            if (f) {
                (void) fftwf_import_wisdom_from_file(f);
                fclose(f);
                fftw_wisdom_loaded = true;
            }
        }

        // allocate input and output arrays with optimal alignment

        input  = reinterpret_cast < std::complex < float > * > (fftwf_malloc(fft_bins * sizeof(fftw_complex)));

        // zero padding portion of input array 
        for (int i = win_size; i < fft_bins; ++i)
                 input[i] = 0.0;
             
        output  = reinterpret_cast < std::complex < float > * > (fftwf_malloc(fft_bins * sizeof(fftw_complex)));

        plan  = fftwf_plan_dft_1d (fft_bins, 
                                  reinterpret_cast < fftwf_complex *> (input),
                                  reinterpret_cast < fftwf_complex *> (output),
                                  -1, FFTW_PATIENT);

        if (win_type == WINDOW_HAMMING)
            generateHammingWindow();
    };

    ~SlidingSpectrum() {
        // silently fail if we can't export wisdom
        if (! fftw_wisdom_loaded) {
            FILE *f = fopen(fftw_wisdom_filename.c_str(), "wb");
            if (f) {
                (void) fftwf_export_wisdom_to_file(f);
                fclose(f);
            }
        }
        fftwf_destroy_plan(plan);
        fftwf_free(input);
        fftwf_free(output);
    };

    bool operator () (std::complex < float >  &d) {
        // process a sample; return true if new FFT output is available
        if (acc (d)) {

            // window the samples and do the fft
            auto a = acc.array_one();

            // loop to process both linear array segments
            int i;
            for (int jj = 0; jj < 2; ++jj) {
                if (a.first != 0 && a.second > 0) {
                    if (have_window) {
                        for (size_t ii = 0; ii < a.second; ++ii, ++i) {
                            input[i] = a.first[ii] * window[i];
                        }
                    } else {
                        for (size_t ii = 0; ii < a.second; ++ii, ++i) {
                            input[i] = a.first[ii];
                        }
                    }
                }
                // get the 2nd linear segment and repeat (once only)
                a = acc.array_two();
            }
            // do FFT
            fftwf_execute(plan);

            return true;
        }
        return false;
    };

    size_t size() {
        return fft_bins;
    };

    std::complex < float > & operator [] (unsigned i) {
        // return the i'th bin of the most recent FFT output
        // only valid if operator() has returned true
        return output[i];
    };

    static std::string fftw_wisdom_filename;

protected:

    int win_size;
    int pad;
    int overlap;
    WinType win_type;
    int fft_bins;

    Accordion < std::complex < float > > acc; // to generate window overlaps

    static bool fftw_wisdom_loaded;

    std::complex < float > *input; // for 2 channels: floating point samples, as interleaved I/C pairs
    std::complex < float > *output; // DFT output from pulse samples
    fftwf_plan plan; // FFT plan for I/Q sample

    bool have_window;
    std::vector < float > window;
    double win_sum;
    double win_sumsq;

    void generateHammingWindow()
    {
        // generate a Hamming window of size N in float array window, 
        // Return sums and sums of squares of window coefficients in win_sum and win_sumsq.

        window = std::vector < float > (win_size);
        win_sum = win_sumsq = 0.0;
        for (int i = 0; i < win_size; ++i) {
            window[i] = 0.54 - 0.46 * cosf(2 * M_PI * i / (win_size - 1));
            win_sum += window[i];
            win_sumsq += window[i] * window[i];
        }
        have_window = true;
    };

};

std::string SlidingSpectrum::fftw_wisdom_filename = "fftwf3_wisdom.dat";
bool SlidingSpectrum::fftw_wisdom_loaded = false;

#endif // _SLIDING_SPECTRUM_H_
