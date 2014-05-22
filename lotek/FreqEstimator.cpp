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

FreqEstimator::FreqEstimator (int max_frames, int channels) :
    m_max_frames(max_frames),
    m_channels(channels),
    m_first_zero_frame(max_frames), // indicate input float array needs zeroing
    m_freq_scale(1.0 / max_frames),
    m_have_window(false),
    m_power_scale(1.0 / sqrtf(max_frames))
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
    //    generateWindowingCoefficients();
    calculateRatioMaps();
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

float FreqEstimator::get (int16_t *samples, int frames, float *pwr_out)
{
    int i, j = m_channels * std::min(frames, m_max_frames);
    float *input = (float *) m_input;

    // convert samples from S16_LE to float

    if (m_have_window) {
        for (i = 0; i < j; ++i)
            input[i] = samples[i] * m_window[i / 2];
    } else {
        for (i = 0; i < j; ++i)
            input[i] = samples[i];
    }        
    
    // pad with zeroes, if necessary

    for(j = m_channels * m_first_zero_frame; i < j; ++i)
        input[i] = 0.0;

    // note start of zero block in input buffer, so we don't have
    // to re-zero it next time

    m_first_zero_frame = frames;

    // do FFT
    fftwf_execute(m_plan);

    // find the max power bin
            
    int bin_high = m_fft_bins;
            
    float max_power = 0.0;
    int max_bin = -1;
    for (j = 1; j < bin_high; ++j) {  // skip DC bin (j = 0)
        float pwr = Pwr (m_output[j]);
        if (pwr > max_power) {
            max_power = pwr;
            max_bin = j;
        }
    }
            
    float prevPwr, nextPwr;

    if (max_bin > 1) {
        prevPwr = Pwr(m_output[max_bin - 1]);
    } else {
        prevPwr = 0;
    }
    if (max_bin >= bin_high) {
        nextPwr = 0;
    } else {
        nextPwr = Pwr(m_output[max_bin + 1]);
    }


    float bin_offset = estimateBinOffset(prevPwr, max_power, nextPwr);

    // sum the energy in the three bins, then scale
    if (pwr_out)
        *pwr_out = 10 * log10f((prevPwr + max_power + nextPwr) * m_power_scale / 32767 / 32767);
 
    // win + bin_offset is the estimate of the (possibly non-integer) bin with the highest
    // power.  We convert that to an estimate in cycles per sample, keeping in mind that
    // the data were zero-padded out to m_max_frames

    float bin_est = max_bin + bin_offset;

    if (bin_est > m_fft_bins / 2)
        // frequencies above 1/2 Nyquist are considered negative frequencies
        bin_est = - (m_fft_bins - bin_est);

    return bin_est * m_freq_scale;
}

void
FreqEstimator::generateWindowingCoefficients()
{
    // generate a Hamming window of size N in float array window, 
    // Return sums and sums of squares of window coefficients in win_sum and win_sumsq.

    m_window = std::vector < float > (m_max_frames);
    m_win_sum = m_win_sumsq = 0.0;
    for (int i = 0; i < m_max_frames; ++i) {
        m_window[i] = 0.54 - 0.46 * cosf(2 * M_PI * i / (m_max_frames - 1));
        m_win_sum += m_window[i];
        m_win_sumsq += m_window[i] * m_window[i];
    }
    m_have_window = true;
    m_power_scale = 1.0 / (sqrtf(m_max_frames) * (m_win_sum * m_win_sum / 2));
}

float
FreqEstimator::binRatio(float d)
{
    // ratio of power in adjacent bins for a pure tone of given offset
    // frequency (in units of bins) from being an exact bin.

    // Consider the DFT of exp(2 \pi \i f i / n) on the sequence i = 0, 1, ..., n - 1

    // if f = m + d where m is an integer in 2..n -2 and abs(d) <= 0.5 then 
    // this function returns Pow (m - 1) ^2 / Pow ( m ) ^2. 
    // Pow(m) is the power in bin m.
    // If -0.5 <= d <= 0, the function returns the ratio [max(Pow(m-1), Pow(m+1)) / Pow(m)]^2
    // If 0 <= d <= 0.5, the function returns the ratio [min(Pow(m-1), Pow(m+1)) / Pow(m)]^2
    // In particular, the function returns 0 for d = 0 and 1 for d = -0.5

    double t1, t2, t3, t31, t32, t4, t5, t51, t6, t7, t72, t8, t82;
    if (fabsf(d) < 1.0e-5)
        return 0.0;

    t1 = 1.0 / m_max_frames;
    t2 = 2.0 * M_PI * d * t1;
    t3 = cos(t2);
    t31 = t3 - 1.0;
    t32 = sin(t2);
    t4 = 2.0 * M_PI * (d + 1.0) * t1;
    t5 = cos(t4);
    t51 = t5 - 1.0;
    t6 = 2.0 * M_PI * d;
    t7 = cos(t6);
    t72 = sin(t6);
    t8 = t72 * t72 + t7 * t7 - 2.0 * t7 + 1.0;
    t82 = sin(t4);
    return t31 * t31 * t8 * (t82 * t82 +t5 * t5 - 2.0 * t5 + 1.0) / (t51 * t51 * t8 * (t32 * t32 + t3 * t3 - 2.0 * t3 + 1.0));

};


float 
FreqEstimator:: estimateBinOffset(float pwrPrev, float pwrMax, float pwrNext)
{
    // pwrMax: the power in the bin where a local maximum occurs
    // pwrPrev: the power in the preceding bin, if that is not the DC bin, else -1
    // pwrNext: the power in the following bin, if that exists, else -1
    // at most one of pwrPrev and pwrNext can be -1; other parameters must be non-negative

    if (pwrPrev <= 0.0)
        return m_ratio2_to_offset [pwrNext / pwrMax * m_ratioMapScale2];
    if (pwrNext <= 0.0)
        return m_ratio2_to_offset [pwrPrev / pwrMax * m_ratioMapScale2];
    
    // in the general case, we get ratios on both sides, with each giving an offset.
    // we average these.

    float off1, off2;
    if (pwrPrev >= pwrNext) { 
        // prev is 2nd highest bin; offset is positive
        off1 = m_ratio2_to_offset [pwrPrev / pwrMax * m_ratioMapScale2];
        off2 = m_ratio3_to_offset [pwrNext / pwrMax * m_ratioMapScale3];
        return - (off1 + off2) / 2;
    } else {  
        // next is 2nd highest bin; offset is positive
        off1 = m_ratio2_to_offset [pwrNext / pwrMax * m_ratioMapScale2];
        off2 = m_ratio3_to_offset [pwrPrev / pwrMax * m_ratioMapScale3];
        return (off1 + off2) / 2;
    }
};

void
FreqEstimator::calculateRatioMaps(float prec) 
{
    // invert the binRatio function to a given precision, on each branch
    // m_ratio2_to_offset will map [0,1] -> [0, 0.5] using the negative side of binRatio()
    // m_ratio3_to_offset will map [0,1] -> [0, 0.3334095 ] using the positive side of binRatio()

    // we use bisection on this monotonic function to search for inverse values

    int n1 = binRatio(-0.5) / prec;
    int n2 = binRatio( 0.5) / prec;
    
    m_ratioMapScale2 = n1 / binRatio(-0.5);
    m_ratioMapScale3 = n2 / binRatio(0.5);

    m_ratio2_to_offset = std::vector < float > (n1);
    m_ratio3_to_offset = std::vector < float > (n2);

    int i;
    float y = 0;
    for (i = 0; i < n1; ++i) {
        float ylo = y, yhi = 0.5;
        for (;;) {
            float diff = binRatio(-y) - i / m_ratioMapScale2;
            if (fabsf(diff) < prec / 10)
                break;
            if (diff < 0.0) {
                ylo = y;
                y = (yhi + y) / 2;
            } else {
                yhi = y;
                y = (ylo + y) / 2;
            }
        }
        m_ratio2_to_offset[i] = y;
    }

    y = 0;
    for (i = 0; i < n2; ++i) {
        float ylo = y, yhi = 0.5;
        for (;;) {
            float diff = binRatio(y) - i / m_ratioMapScale3;
            if (fabsf(diff) < prec / 10)
                break;
            if (diff < 0.0) {
                ylo = y;
                y = (yhi + y) / 2;
            } else {
                yhi = y;
                y = (ylo + y) / 2;
            }
        }
        m_ratio3_to_offset[i] = y;
    }

};
    
std::vector < float > m_ratio2_to_offset; // lookup table converting ratio (of 2nd bin to max) to offset
std::vector < float > m_ratio3_to_offset; // lookup table converting ratio (of 3rd bin to max) to offset

const char * FreqEstimator::m_fftw_wisdom_filename = "./fftw_wisdom.dat";
bool FreqEstimator::m_fftw_wisdom_loaded = false;

std::map < std::pair < int, short >, cachedFFTWPlan > 
FreqEstimator::m_cached_plans;
