/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    FreqEstimator.h - estimate main frequency of a pulse from stereo
    raw audio data as an I/Q pair.
    
    Assumes input data have had DC offset corrected with a slow-moving
    average, with a decay constant much larger than the pulse length.
    That way, the 0 bin of the DFT contains energy from the pulse,
    rather than the DC offset.  This is required in case the pulse
    frequency is less than the frequency of the 1-st non-DC DFT bin.
    
    Copyright 2014 John Brzustowski

    License: GPL v 2.0 or later.

*/

#include "FreqEstimator.h"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <string.h>

#ifdef MINGW
#define fftw_free(X) fftwf_free(X)
#endif

FreqEstimator::FreqEstimator (int frames) :
    m_frames(frames)
{
    if (! m_fftw_wisdom_loaded) {
        // silently fail if wisdom cannot be found
        FILE *f = fopen(m_fftw_wisdom_filename, "r");
        if (f) {
            (void) fftwf_import_wisdom_from_file(f);
            fclose(f);
            m_fftw_wisdom_loaded = true;
        }
    }

    m_input =  reinterpret_cast < std::complex < float > *> (fftwf_malloc(m_frames * sizeof(fftw_complex)));
    memset(m_input, 0, sizeof(std::complex < float > ) * m_frames);

    m_output = reinterpret_cast < std::complex < float > *> (fftwf_malloc(m_frames * sizeof(fftw_complex)));
    m_plan = fftwf_plan_dft_1d(m_frames,
                               reinterpret_cast < fftwf_complex * > (m_input), 
                               reinterpret_cast < fftwf_complex * > (m_output)
                               , -1, FFTW_PATIENT);

    // silently fail if we can't export wisdom
    FILE *f = fopen(m_fftw_wisdom_filename, "wb");
    if (f) {
        (void) fftwf_export_wisdom_to_file(f);
        fclose(f);
    }
    generateWindowingCoefficients();
}

FreqEstimator::~FreqEstimator()
{
    fftwf_destroy_plan(m_plan);
    fftwf_free(m_input);
    fftwf_free(m_output);
}


#define Pwr(x) (x.real() * x.real() + x.imag() * x.imag())

// allow for incoming frames to be in two separate contiguous segments, as they might come from
// a ring buffer.  Must have n1 + n2 = frames

float FreqEstimator::get (std::complex < float > * seg1, int n1, std::complex < float > * seg2, int n2)
{
    int i = 0;

    for (int j = 0; j < n1; ++j, ++i) {
        m_input[i] = seg1[j] * m_window[i];
    }
    for (int j = 0; j < n2; ++j, ++i) {
        m_input[i] = seg2[j] * m_window[i];
    }
    // do FFT
    fftwf_execute(m_plan);

    // find the max power bin
            
    int bin_high = m_frames;
            
    float max_power = 0.0;
    int max_bin = -1;
    for (int j = 0; j < bin_high; ++j) {
        float pwr = Pwr (m_output[j]);
        if (pwr > max_power) {
            max_power = pwr;
            max_bin = j;
        }
    }

    float bin_offset = estimateBinOffset(max_bin);
 
    // win + bin_offset is the estimate of the (possibly non-integer) bin with the highest
    // power.  We convert that to an estimate in cycles per sample, keeping in mind that
    // the data were zero-padded out to m_frames

    float bin_est = max_bin + bin_offset;

    if (bin_est > m_frames / 2.0)
        // frequencies above Nyquist are considered negative frequencies
        bin_est = - (m_frames - bin_est);

    return bin_est;
}

void
FreqEstimator::generateWindowingCoefficients()
{
    // generate a Hamming window of size N in float array window, 
    // Return sums and sums of squares of window coefficients in win_sum and win_sumsq.

    m_window = std::vector < float > (m_frames);
    m_win_sum = m_win_sumsq = 0.0;
    for (int i = 0; i < m_frames; ++i) {
        m_window[i] = 0.54 - 0.46 * cosf(2 * M_PI * i / (m_frames - 1));
        m_win_sum += m_window[i];
        m_win_sumsq += m_window[i] * m_window[i];
    }
}

double
FreqEstimator::cubicMaximize(double y0, double y1, double y2, double y3)
{
   // Find coefficients of cubic

   double a, b, c;

   a = y0 / -6.0 + y1 / 2.0 - y2 / 2.0 + y3 / 6.0;

   if (a == 0.0)
       return double(-1); // error

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

#define WRAP_BIN(x) (((x) + m_frames) % m_frames)
double 
FreqEstimator::estimateBinOffset(int max_bin)
{

    // use a cubic estimator to find the peak frequency estimate using nearby bins

    float pwr[4];
    for (int i = -1; i < 3; ++i)
        pwr[i+1] = Pwr(m_output[WRAP_BIN(max_bin + i)]);

    // get the estimate of the peak beat frequency (in bin units)
    return cubicMaximize(pwr[0], pwr[1], pwr[2], pwr[3]) - 1;
}

const char * FreqEstimator::m_fftw_wisdom_filename = "./fftw_wisdom.dat";
bool FreqEstimator::m_fftw_wisdom_loaded = false;
