/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    FreqEstimator.h - estimate main frequency of a pulse from raw audio data.  FIXME: hardcoded S16_LE
    
    Copyright 2014 John Brzustowski

    License: GPL v 2.0 or later.

*/

#include "FreqEstimator.h"

#include <cmath>
#include <complex>
#include <stdexcept>
#include <iostream>

#ifdef MINGW
#define fftw_free(X) fftwf_free(X)
#endif

FreqEstimator::FreqEstimator (unsigned max_frames, unsigned channels) :
    m_max_frames(max_frames),
    m_channels(channels),
    m_first_zero_frame(max_frames), // indicate input float array needs zeroing
    m_freq_scale(1.0 / max_frames)
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

    switch (channels) {
    case 1:
        m_input =  (fftwf_complex *) fftwf_malloc(m_max_frames * sizeof(float));
        m_fft_bins = (m_max_frames + 1) / 2;
        m_output = (fftwf_complex *) fftwf_malloc( m_fft_bins * sizeof(fftw_complex));
        m_plan = fftwf_plan_dft_r2c_1d(m_max_frames, (float *) m_input, m_output, FFTW_PATIENT);
        break;
    case 2:
        m_input =  (fftwf_complex *) fftwf_malloc(m_max_frames * sizeof(fftw_complex));
        m_fft_bins = m_max_frames;
        m_output = (fftwf_complex *) fftwf_malloc(m_fft_bins * sizeof(fftw_complex));
        m_plan = fftwf_plan_dft_1d(m_max_frames, m_input, m_output, -1, FFTW_PATIENT);
        break;
    default:
        throw std::runtime_error("FreqEstimator: channels must be 1 or 2");
    }
}

FreqEstimator::~FreqEstimator()
{
    // silently fail if we can't export wisdom
    FILE *f = fopen(m_fftw_wisdom_filename, "wb");
    if (f) {
        (void) fftwf_export_wisdom_to_file(f);
        fclose(f);
    }
    fftwf_destroy_plan(m_plan);
    fftwf_free(m_input);
    fftwf_free(m_output);
}


#define Pwr(x) (((float *)(&x))[0] * ((float *)(&x))[0] + ((float *)(&x))[1] * ((float *)(&x))[1])

float FreqEstimator::get (int16_t *samples, unsigned frames, float *pwr_out)
{
    unsigned i, j = m_channels * std::min(frames, m_max_frames);
    float *input = (float *) m_input;

    // convert samples from S16_LE to float

    for (i = 0; i < j; ++i)
        input[i] = samples[i];
    
    // pad with zeroes, if necessary

    for(j = m_channels * m_first_zero_frame; i < j; ++i)
        input[i] = 0.0;

    // note start of zero block in input buffer, so we don't have
    // to re-zero it next time

    m_first_zero_frame = frames;

    // do FFT
    fftwf_execute(m_plan);

    // find the max power bin
            
    unsigned bin_high = m_fft_bins;
            
    float max_power = 0.0;
    int max_bin = -1;
    for (j = 1; j < bin_high; ++j) {  // skip DC bin (j = 0)
        float pwr = Pwr (m_output[j]);
        if (pwr > max_power) {
            max_power = pwr;
            max_bin = j;
        }
    }
            
    // use a quadratic estimator to find the peak frequency estimate using
    // 3 bins that include the max power bin.  These will be numbered
    // win - 1, win, win+1, and win + 2. We avoid the DC bin, and
    // don't overstep the array

    // FIXME: what if max_bin == 1 or max_bin == m_fft_bin - 1 ?
    int win = max_bin;

    if (win == 1)
        ++ win;
    else if (win + 1 >= m_fft_bins)
        -- win;

    // get the power values in the window of 3 bins, for peak estimation
    // apparently, this works better in log units

    float pwr[3] = {log10f(Pwr(m_output[win - 1])),
                    log10f(Pwr(m_output[win + 0])),
                    log10f(Pwr(m_output[win + 1]))};

    // this is the quadratic estimate of peak location given the three values
    // computed above.

    float bin_offset = (pwr[0] - pwr[2]) / (2.0* (pwr[2] - 2.0 * pwr[1] + pwr[0]));

    // interpolate the quadratic estimate.
    if (pwr_out)
        *pwr_out = exp10f(pwr[1] - (pwr[0]-pwr[2]) * (pwr[0] - pwr[2]) / (8.0 * (pwr[2]-2.0 * pwr[1]+pwr[0]))) / frames;
 
    // win + bin_offset is the estimate of the (possibly non-integer) bin with the highest
    // power.  We convert that to an estimate in cycles per sample, keeping in mind that
    // the data were zero-padded out to m_max_frames

    if (win + bin_offset > m_fft_bins / 2)
        return - (m_fft_bins - win - bin_offset) * m_freq_scale;

    return (win + bin_offset) * m_freq_scale;
}

void
FreqEstimator::generateWindowingCoefficients(int N, std::vector < float > &window, float &win_sum, float &win_sumsq)
{
    // generate a Hamming window of size N in float array window, 
    // Return sums and sums of squares of window coefficients in win_sum and win_sumsq.

    window = std::vector < float > (N);
    win_sum = win_sumsq = 0.0;
    for (int i = 0; i < N; ++i) {
        window[i] = 0.54 - 0.46 * cosf(2 * M_PI * i / (N - 1));
        win_sum += window[i];
        win_sumsq += window[i] * window[i];
    }
}

const char * FreqEstimator::m_fftw_wisdom_filename = "./fftw_wisdom.dat";
bool FreqEstimator::m_fftw_wisdom_loaded = false;

std::map < std::pair < int, short >, cachedFFTWPlan > 
FreqEstimator::m_cached_plans;
