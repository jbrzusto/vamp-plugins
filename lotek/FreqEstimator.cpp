/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    FreqEstimator.h - estimate main frequency of a pulse from raw audio data.
    Assumes input data have had DC offset corrected with a slow-moving average,
    with a decay constant much larger than the pulse length.  That way, 
    the 0 bin of the DFT contains energy from the pulse, rather than the 
    DC offset.  This is required in case the pulse frequency is less than
    the frequency of the 1-st non-DC DFT bin.
    
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

FreqEstimator::FreqEstimator (int max_frames, int channels, bool standard) :
    m_max_frames(max_frames),
    m_channels(channels),
    m_standard(standard),
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
    //    calculateRatioMaps();
    calculateBalRatioMap();
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


#define Pwr(x) (((x)[0]*(x)[0] + (x)[1]*(x)[1]) / m_max_frames)

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
    for (j = 0; j < bin_high; ++j) {  // skip DC bin (j = 0)
        float pwr = m_output[j][0] = Pwr (m_output[j]);
        if (pwr > max_power) {
            max_power = pwr;
            max_bin = j;
        }
    }

    float bin_offset;
    if (m_standard) {
        bin_offset = standardEstimateBinOffset(max_bin);
    } else {
        bin_offset = estimateBinOffset(max_bin);
    }

    // if requested, estimate the power at the true frequency
    if (pwr_out)
        *pwr_out = 10 * log10f(max_power * binRatio(bin_offset) * m_power_scale / 32767 / 32767);
 
    // win + bin_offset is the estimate of the (possibly non-integer) bin with the highest
    // power.  We convert that to an estimate in cycles per sample, keeping in mind that
    // the data were zero-padded out to m_max_frames

    float bin_est = max_bin + bin_offset;

    if (bin_est > m_fft_bins)
        // frequencies above Nyquist are considered negative frequencies
        bin_est = - (m_fft_bins - bin_est);

    return bin_est * m_freq_scale;
}

void
FreqEstimator::generateWindowingCoefficients()
{
    // generate a Hamming window of size N in float array window, 
    // Return sums and sums of squares of window coefficients in win_sum and win_sumsq.
    // This is mainly for testing, to demonstrate that windowing hinders
    // peak frequency estimation.

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

double
FreqEstimator::balancedBinRatio(double d)
{
    // (Pwr(x+d)-Pwr(x+d+1)) / (Pwr(x+d)-Pwr(x+d-1)) for the function

    //    return (powerAtOffset(d) - powerAtOffset(d - 1)) / (powerAtOffset(d) - powerAtOffset(d + 1));

    // more stable?

    // double ad = 2 * M_PI * d;
    // // double adp1 = 2 * M_PI * (d + 1);
    // // double adm1 = 2 * M_PI * (d - 1);
    // double adp1 = ad;
    // double adm1 = ad;

    // double adn = ad / m_max_frames;
    // double adp1n = adp1 / m_max_frames;
    // double adm1n = adm1 / m_max_frames;

    // return (1.0 - ((cos(adn) - 1.0) / (cos(adm1n) - 1.0)) * ((cos(adm1) - 1.0) / (cos(ad) - 1.0))) / 
    //     (1.0 - ((cos(adn) - 1.0) / (cos(adp1n) - 1.0)) * ((cos(adp1) - 1.0) / (cos(ad) - 1.0)));

    // (Pwr(d+1)-Pwr(d-1))/Pwr(d)
    if (fabs(d) < 1e-6)
        return 0.0;
    return (powerAtOffset(d+1) - powerAtOffset(d-1)) / powerAtOffset(d);
}

float
FreqEstimator::binRatio(float d)
{
    // ratio of true power at x to power in DFT bin at x + d for the function
    // exp (2 * pi * %i% * x)
    //
    // courtesy of maxima:
    //  f(d);
    //
    //                              n - 1   2 %i %pi d i
    //                              ====    ------------
    //                              \            n
    // (%o123)                       >    %e
    //                              /
    //                              ====
    //                              i = 0
    // (%i125) f90(optimize(factor(simplify_sum(cabs(f(0)/f(d))))));
    // block([%1,%2,%3,%4],%1:2*%pi*d/n,%2:cos(%1),%3:2*%pi*d,%4:cos(%3)&
    // &,2*abs(%2-1)*abs(n)/sqrt((sin(%3)**2+%4**2-2*%4+1)*(sin(%1)**2+%2&
    // &**2-2*%2+1)))
    
    float t1, t2, t3, t4, t5, t6;
    t1 = 2 * M_PI * d;
    t3 = t1 / m_max_frames;
    t2 = cos(t3);
    t4 = cos(t1);
    t5 = sin(t1);
    t6 = sin(t3);
    return sqrt((fabs(t2 - 1.0) * m_max_frames) * 2.0 / ((t5 * t5 + t4 * t4 - 2.0 * t4 + 1.0) * (t6 * t6 + t2 * t2 - 2.0 * t2 + 1.0)));
}

float
FreqEstimator::twoBinRatio(float d)
{
    // ratio of power in DFT bin x+d to power in bin x+d+1 for the function
    // exp (2 * pi * %i% * x)
    //
    // should be equal to binRatio(d) / binRatio (d + 1), but calculated more efficiently.
    //
    // courtesy of maxima:
    //
    // (%i127) f(d);
    //                              n - 1   2 %i %pi d i
    //                              ====    ------------
    //                              \            n
    // (%o127)                       >    %e
    //                              /
    //                              ====
    //                              i = 0
    // (%i128) f90(optimize(factor(simplify_sum(cabs(f(d)/f(d+1))))));
    // block([%1,%2,%3,%4,%5,%6,%7,%8],%1:1/n,%2:2*%pi*d*%1,%3:cos(%2),%&
    // &4:2*%pi*(d+1)*%1,%5:cos(%4),%6:2*%pi*d,%7:cos(%6),%8:sin(%6)**2+%&
    // &7**2-2*%7+1,abs(%5-1)*sqrt(%8*(sin(%2)**2+%3**2-2*%3+1))/(abs(%3-&
    // &1)*sqrt(%8*(sin(%4)**2+%5**2-2*%5+1))))
    //
    //
    // Details:
    //
    // Consider the DFT of exp(2 \pi \i f i / n) on the sequence i = 0, 1, ..., n - 1
    //
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


double 
FreqEstimator:: estimateBinOffset(int bin)
{
    // given 'bin' has max power, estimate the offset to the true
    // frequency using power in adjacent bins.

    // Requires that the power has been stored in the real component
    // of m_output.

    //  FIXME: allow for a constant noise value
    // added to each bin.

    // We compare the ratios of power of each adjacent bin to the central
    // one to determine which branch of the curve to apply to each side
    // i.e. to determine whether the true frequency is above or below the
    // exact frequency corresponding to bin.

    // Then we apply the appropriate map to each ratio, and average the resulting
    // estimated offsets.  FIXME: we can detect uniform noise by noticing when
    // these offsets don't match well, then subtract noise until they do.

    // In the case of the endpoint bins 0, and -1, we wrap around, and
    // we use the so-called DC bin on the assumption the original signal
    // has been (slowly) corrected for DC bias, leaving any low frequency
    // energy from the pulse in bin 0.

    double pwrMax = m_output[bin][0];

    int prevBin = bin - 1;
    if (prevBin == -1)
        prevBin = m_max_frames - 1;
    int nextBin = bin +1;
    if (nextBin == m_max_frames)
        nextBin = 0;

    double ratio = (m_output[nextBin][0] - m_output[prevBin][0]) / pwrMax;

    double off;
    if (ratio > 0) {
        off = - m_balanced_ratio_to_offset[ ratio * m_balRatioMapScale];
    } else {
        off = m_balanced_ratio_to_offset[ - ratio * m_balRatioMapScale];
    };
    return off;
};

void
FreqEstimator::calculateRatioMaps(float prec) 
{
    // invert the twoBinRatio function to a given precision, on each branch
    // m_ratio2_to_offset will map [0,1] -> [0, 0.5] using the negative side of twoBinRatio()
    // m_ratio3_to_offset will map [0,1] -> [0, 0.3334095 ] using the positive side of twoBinRatio()

    // we use bisection on this monotonic function to search for inverse values

    // FIXME: redo assuming linear interpolation is used on resulting values

    int n1 = twoBinRatio(-0.5) / prec;
    int n2 = twoBinRatio( 0.5) / prec;
    
    m_ratioMapScale2 = n1 / twoBinRatio(-0.5);
    m_ratioMapScale3 = n2 / twoBinRatio(0.5);

    m_ratio2_to_offset = std::vector < float > (n1);
    m_ratio3_to_offset = std::vector < float > (n2);

    int i;
    double y = 0;
    for (i = 0; i < n1; ++i) {
        double ylo = y, yhi = std::min(y + 0.1, (double) 0.5);
        for (;;) {
            double diff = twoBinRatio(-y) - i / m_ratioMapScale2;
            if (fabs(diff) < prec / 2)
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
        double ylo = y, yhi = std::min(y + 0.1, (double) 0.5);
        for (;;) {
            double diff = twoBinRatio(y) - i / m_ratioMapScale3;
            if (fabs(diff) < prec / 2)
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

double 
FreqEstimator::powerAtOffset(double d) {
    return ((cos( (double) 2.0 * M_PI * d) - 1) /
            (cos( (double) 2.0 * M_PI * d / m_max_frames) - 1));
};

void
FreqEstimator::calculateBalRatioMap(double prec)
{
    // invert the balancedBinRatio function to a given precision
    // this will map 0..1 to 0.5 .. 0.0
    // FIXME: redo assuming linear interpolation is used on resulting values
    
    // int n = 1.0 / prec;

    // m_balRatioMapScale = n;

    // double scale = 1.0 / (double) n;

    // m_balanced_ratio_to_offset = std::vector < double > (n);

    // int i;
    // double y = 0.5;
    // for (i = 0; i < n; ++i) {
    //     double yhi = 0.5, ylo = 0.0;
    //     for (;;) {
    //         double diff = balancedBinRatio(y) - i * scale;
    //         if (fabs(diff) < prec/3 && fabs(yhi - ylo) < prec/3)
    //             break;
    //         if (diff >= 0.0) {
    //             ylo = y;
    //             y = (yhi + y) / 2;
    //         } else {
    //             yhi = y;
    //             y = (ylo + y) / 2;
    //         }
    //     }
    //     m_balanced_ratio_to_offset[i] = y;
    //     std::cerr << i * m_balRatioMapScale << ',' << y << std::endl;
    // }


    int n = 1.0 / prec;

    double ymax = 0.51362681;
    m_balRatioMapScale = n;

    double scale = 1.0 / (double) n;

    m_balanced_ratio_to_offset = std::vector < double > (n);

    int i;
    double y = 0.0;
    for (i = 0; i < n; ++i) {
        double yhi = ymax, ylo = 0.0;
        for (;;) {
            double diff = balancedBinRatio(-y) - i * scale;
            if (fabs(diff) < prec/3 && fabs(yhi - ylo) < prec/3)
                break;
            if (diff <= 0.0) {
                ylo = y;
                y = (yhi + y) / 2;
            } else {
                yhi = y;
                y = (ylo + y) / 2;
            }
        }
        m_balanced_ratio_to_offset[i] = -y;
    }

};

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

double 
FreqEstimator::standardEstimateBinOffset(int max_bin)
{

    // use a cubic estimator to find the peak frequency estimate using nearby bins
    
    int bin_low = std::max(0, std::min(m_max_frames / 2 - 4, max_bin - 1));  // avoid the DC bin
                    
    if (bin_low + 3 <= m_max_frames / 2) {
        // get the estimate of the peak beat frequency (in bin units)
        return bin_low - max_bin + cubicMaximize(m_output[bin_low][0], m_output[bin_low+1][0], m_output[bin_low+2][0], m_output[bin_low+3][0]);
    } else {
        return 0;
    }
}
    
std::vector < float > m_ratio2_to_offset; // lookup table converting ratio (of 2nd bin to max) to offset
std::vector < float > m_ratio3_to_offset; // lookup table converting ratio (of 3rd bin to max) to offset
std::vector < float > m_balanced_ratio_to_offset;

const char * FreqEstimator::m_fftw_wisdom_filename = "./fftw_wisdom.dat";
bool FreqEstimator::m_fftw_wisdom_loaded = false;

std::map < std::pair < int, short >, cachedFFTWPlan > 
FreqEstimator::m_cached_plans;
