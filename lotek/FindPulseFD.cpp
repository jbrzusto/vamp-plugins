/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
  Vamp

  An API for audio analysis and feature extraction plugins.

  Centre for Digital Music, Queen Mary, University of London.
  Copyright 2006 Chris Cannam.

  VAMP license:

  Permission is hereby granted, free of charge, to any person
  obtaining a copy of this software and associated documentation
  files (the "Software"), to deal in the Software without
  restriction, including without limitation the rights to use, copy,
  modify, merge, publish, distribute, sublicense, and/or sell copies
  of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:
    
  The above copyright notice and this permission notice shall be
  included in all copies or substantial portions of the Software.
    
  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
  ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
  CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
  WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    
  Except as contained in this notice, the names of the Centre for
  Digital Music; Queen Mary, University of London; and Chris Cannam
  shall not be used in advertising or otherwise to promote the sale,
  use or other dealings in this Software without prior written
  authorization.
    
  FindPulseFD.cpp - find pulses from Lotek tags - frequency domain
  Copyright 2012 John Brzustowski

  License: GPL v 2.0 or later.  This is required in order to use fftw.

*/

#include "FindPulseFD.h"

using std::stringstream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;

const char * FindPulseFD::fftw_wisdom_filename = "./fftw_wisdom.dat";

// from Audacity 2.0.1's src/FreqWindow.cpp:

float 
FindPulseFD::cubicMaximize(float y0, float y1, float y2, float y3)
{
    // Find coefficients of cubic

    float a, b, c;

    a = y0 / -6.0 + y1 / 2.0 - y2 / 2.0 + y3 / 6.0;

    if (a == 0.0)
        return float(-1); // error

    b = y0 - 5.0 * y1 / 2.0 + 2.0 * y2 - y3 / 2.0;
    c = -11.0 * y0 / 6.0 + 3.0 * y1 - 3.0 * y2 / 2.0 + y3 / 3.0;

    // Take derivative

    float da, db, dc;

    da = 3 * a;
    db = 2 * b;
    dc = c;

    // Find zeroes of derivative using quadratic equation

    float discriminant = db * db - 4 * da * dc;
    if (discriminant < 0.0) {
        if (discriminant < -1.0)
            return float(-1000);              // error
        else
            discriminant = 0.0;
    }              

    float x1 = (-db + sqrt(discriminant)) / (2 * da);
    float x2 = (-db - sqrt(discriminant)) / (2 * da);

    // The one which corresponds to a local _maximum_ in the
    // cubic is the one we want - the one with a negative
    // second derivative

    float dda = 2 * da;
    float ddb = db;

    if (dda * x1 + ddb < 0)
        {
            return x1;
        }
    else
        {
            return x2;
        }
}

// from Audacity 2.0.1's src/FreqWindow.cpp:

float 
FindPulseFD::cubicInterpolate(float y0, float y1, float y2, float y3, float x)
{
    float a, b, c, d;

    a = y0 / -6.0 + y1 / 2.0 - y2 / 2.0 + y3 / 6.0;
    b = y0 - 5.0 * y1 / 2.0 + 2.0 * y2 - y3 / 2.0;
    c = -11.0 * y0 / 6.0 + 3.0 * y1 - 3.0 * y2 / 2.0 + y3 / 3.0;
    d = y0;

    float xx = x * x;
    float xxx = xx * x;

    return (a * xxx + b * xx + c * x + d);
};

void
FindPulseFD::generateWindowingCoefficients(int N, std::vector < float > &window, float &win_sum, float &win_sumsq)
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

FindPulseFD::FindPulseFD(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_stepSize(0),
    m_blockSize(0),
    m_plen(m_default_plen),
    m_min_pulse_SNR(exp10(m_default_min_pulse_SNR_dB / 10.0)),
    m_fft_win_size (m_default_fft_win_size),
    m_noise_win_size (m_default_noise_win_size),
    m_min_pulse_sep (m_default_min_pulse_sep),
    m_min_freq (m_default_min_freq),
    m_max_freq (m_default_max_freq),
    m_have_fft_plan(false)
{
    // silently fail if wisdom cannot be found
    FILE *f = fopen(fftw_wisdom_filename, "r");
    if (f) {
        (void) fftwf_import_wisdom_from_file(f);
        fclose(f);
    }
}

FindPulseFD::~FindPulseFD()
{
    if (m_have_fft_plan) {
        for (int i=0; i < 2; ++i) {
            fftwf_destroy_plan(m_plan[i]);
            fftw_free(m_windowed[i]); 
        }
        fftw_free(m_fft);
        m_have_fft_plan = false;
    }
    // silently fail if we can't export wisdom
    FILE *f = fopen(fftw_wisdom_filename, "wb");
    if (f) {
        (void) fftwf_export_wisdom_to_file(f);
        fclose(f);
    }
}

string
FindPulseFD::getIdentifier() const
{
    return "findpulseFD";
}

string
FindPulseFD::getName() const
{
    return "Find Pulses in Frequency Domain";
}

string
FindPulseFD::getDescription() const
{
    return "Find pulses (e.g. from Lotek telemetry tags)";
}

string
FindPulseFD::getMaker() const
{
    return "flightcalls.org jbrzusto@fastmail.fm";
}

int
FindPulseFD::getPluginVersion() const
{
    return 1;
}

string
FindPulseFD::getCopyright() const
{
    return "GPL version 2 or later";
}

bool
FindPulseFD::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_stepSize = stepSize;
    m_blockSize = blockSize;

    m_plen_in_samples = (m_plen / 1000.0) * m_inputSampleRate;

    m_pf_size = roundf(m_plen_in_samples / (float) (m_fft_win_size));

    int num_bins = m_fft_win_size;

    m_last_timestamp = std::vector < Vamp::RealTime > (num_bins, Vamp::RealTime(-1, 0));

    m_freq_bin_pulse_finder = std::vector < PulseFinder < double > > (num_bins, PulseFinder < double > (131072, m_pf_size, m_pf_size * m_noise_win_size, m_pf_size * m_min_pulse_sep));

    // allocate time-domain sample buffers for each channel which are large enough to contain
    // the samples for a pulse when it has been detected.  Because pulses are not
    // detected until the sliding window determines them to be maximal in the
    // frequency domain, we need to keep a rather long window.


    // allocate fft output buffer for pulse finding
    m_fft = (fftwf_complex *) fftwf_malloc(m_fft_win_size * sizeof(fftwf_complex));

    // allocate input window and generate the plan
    for (unsigned i=0; i < 2; ++i) {
        m_windowed[i] = (fftwf_complex *) fftwf_malloc(m_fft_win_size * sizeof(fftwf_complex));
        m_plan[i] = fftwf_plan_dft_1d(m_fft_win_size, m_windowed[i], m_fft, -1, FFTW_PATIENT);
    }

    // cap frequency limits at Nyquist
    if (m_min_freq > m_inputSampleRate / 2000)
        m_min_freq = m_inputSampleRate / 2000;
    if (m_max_freq > m_inputSampleRate / 2000)
        m_max_freq = m_inputSampleRate / 2000;
    
    m_first_freq_bin = floorf(m_min_freq * 1000.0 / (m_inputSampleRate / m_fft_win_size) - 0.5);
    m_last_freq_bin =  ceilf(m_max_freq * 1000.0 / (m_inputSampleRate / m_fft_win_size) + 0.5);
    if (m_last_freq_bin < m_first_freq_bin)
        m_last_freq_bin = m_first_freq_bin;

    m_num_windowed_samples[0] = 0;

    // ugly hack to simplify accumulation in odd-phase sliding window
    m_num_windowed_samples[1] = m_fft_win_size / 2;

    // get windowing coefficients for sliding frequency bin window
    generateWindowingCoefficients(m_fft_win_size, m_window, m_win_s1, m_win_s2);

    // get windowing coefficients for pulse-sized window; we don't estimate power
    // from this, so don't need the moments

    float ignore1, ignore2;
    generateWindowingCoefficients(m_plen_in_samples, m_pulse_window, ignore1, ignore2);

    m_probe_scale = m_channels * (m_win_s1 * m_win_s1 / 2);

    m_dcma[0] = MovingAverager < float, float > (5 * m_fft_win_size);
    if (m_channels == 2)
        m_dcma[1] = MovingAverager < float, float > (5 * m_fft_win_size);

    m_odd_phase_window_is_bogus = true;

    return true;
}

void
FindPulseFD::reset()
{
}

FindPulseFD::ParameterList
FindPulseFD::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "plen";
    d.name = "Pulse Length (unit: milliseconds)";
    d.description = "Duration of a transmitted pulse in milliseconds";
    d.unit = "milliseconds";
    d.minValue = 0.1;
    d.maxValue = 500;
    d.defaultValue = FindPulseFD::m_default_plen;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "minsnr";
    d.name = "Minimum Pulse SNR (unit: dB)";
    d.description = "Minimum pulse signal-to-noise ratio";
    d.unit = "dB";
    d.minValue = 0;
    d.maxValue = 96;
    d.defaultValue = FindPulseFD::m_default_min_pulse_SNR_dB;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "fftsize";
    d.name = "Size of FFT window (unit: samples)";
    d.description = "The number of samples in each window for which an FFT is computed.  Windows are non-overlapping.  A window of 96 samples for 96kHz sampling means one FFT calculated each millisecond";
    d.unit = "samples";
    d.minValue = 10;
    d.maxValue = 1000;
    d.defaultValue = FindPulseFD::m_default_fft_win_size;
    d.isQuantized = true;
    d.quantizeStep = 1;
    list.push_back(d);

    d.identifier = "noisesize";
    d.name = "Size of noise window (unit: pulse length)";
    d.description = "Size of window on each side of signal that is used to estimate noise.  In multiples of signal pulse length.";
    d.unit = "pulses";
    d.minValue = 1;
    d.maxValue = 100;
    d.defaultValue = FindPulseFD::m_default_noise_win_size;
    d.isQuantized = true;
    d.quantizeStep = 1;
    list.push_back(d);

    d.identifier = "pulsesep";
    d.name = "Minimum separation of pulses (unit: pulse length)";
    d.description = "Minimum separation between a pulse and adjacent pulses in order to be detected, in units of pulse length.";
    d.unit = "pulses";
    d.minValue = 1;
    d.maxValue = 100;
    d.defaultValue = FindPulseFD::m_default_min_pulse_sep;
    d.isQuantized = true;
    d.quantizeStep = 1;
    list.push_back(d);

    d.identifier = "minfreq";
    d.name = "Minimum Tag Offset Frequency (unit: kHz)";
    d.description = "Minimum frequency by which tag differs from receiver, in kHz";
    d.unit = "kHz";
    d.minValue = - m_inputSampleRate / 2000;
    d.maxValue = m_inputSampleRate / 2000;
    d.defaultValue = FindPulseFD::m_default_min_freq;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxfreq";
    d.name = "Maximum Tag Offset Frequency (unit: kHz)";
    d.description = "Maximum frequency by which tag differs from receiver, in kHz";
    d.unit = "kHz";
    d.minValue = - m_inputSampleRate / 2000;
    d.maxValue = m_inputSampleRate / 2000;
    d.defaultValue = FindPulseFD::m_default_max_freq;
    d.isQuantized = false;
    list.push_back(d);

    return list;
}

float
FindPulseFD::getParameter(string id) const
{
    if (id == "plen") {
        return m_plen;
    } else if (id == "minsnr") {
        return 10 * log10(m_min_pulse_SNR);
    } else if (id == "fftsize") {
        return m_fft_win_size;
    } else if (id == "noisesize") {
        return m_noise_win_size;
    } else if (id == "pulsesep") {
        return m_min_pulse_sep;
    } else if (id == "minfreq") {
        return m_min_freq;
    } else if (id == "maxfreq") {
        return m_max_freq;
    }
    return 0.f;
}

void
FindPulseFD::setParameter(string id, float value)
{
    if (id == "plen") {
        FindPulseFD::m_default_plen = m_plen = value;
    } else if (id == "minsnr") {
        FindPulseFD::m_default_min_pulse_SNR_dB =  value;
        m_min_pulse_SNR = exp10(value / 10.0);
    } else if (id == "fftsize") {
        FindPulseFD::m_default_fft_win_size = m_fft_win_size = value;
    } else if (id == "noisesize") {
        FindPulseFD::m_default_noise_win_size = m_noise_win_size = value;
    } else if (id == "pulsesep") {
        FindPulseFD::m_default_min_pulse_sep = m_min_pulse_sep = value;
    } else if (id == "minfreq") {
        FindPulseFD::m_default_min_freq = m_min_freq = value;
    } else if (id == "maxfreq") {
        FindPulseFD::m_default_max_freq = m_max_freq = value;
    }
}


FindPulseFD::OutputList
FindPulseFD::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor zc;

    zc.identifier = "pulses";
    zc.name = "Pulses";
    zc.description = "The locations and features of pulses";
    zc.unit = "";
    zc.hasFixedBinCount = true;
    zc.binCount = 0;
    zc.sampleType = OutputDescriptor::VariableSampleRate;
    zc.sampleRate = m_inputSampleRate;
    list.push_back(zc);

    return list;
}

FindPulseFD::FeatureSet
FindPulseFD::process(const float *const *inputBuffers,
                     Vamp::RealTime timestamp)
{
    FeatureSet returnFeatures;

    if (m_stepSize == 0) {
	cerr << "ERROR: FindPulseFD::process: "
	     << "FindPulseFD has not been initialised"
	     << endl;
	return returnFeatures;
    }

    for (unsigned int i=0; i < m_blockSize; ++i) {
        for (unsigned short ch = 0; ch < m_channels; ++ch) {
            m_dcma[ch].process(inputBuffers[ch][i]);
        }

        if (! m_dcma[0].have_average())
            continue;
        
        // append each weighted sample to each window
        for (unsigned short w = 0; w < 2; ++w) {

            float avg1 = m_dcma[0].get_average();
            float avg2 = m_dcma[1].get_average();
            // to omit DC cancelling, do: float avg1 = 0.0, avg2 = 0.0;
            m_windowed[w][m_num_windowed_samples[0]][0] = (inputBuffers[0][i] - avg1) * m_window[m_num_windowed_samples[0]];
            m_windowed[w][m_num_windowed_samples[0]][1] = (inputBuffers[1][i] - avg2) * m_window[m_num_windowed_samples[0]];
        }

        for (unsigned short w=0; w < 2; ++w) {
            ++m_num_windowed_samples[w];

            if (m_num_windowed_samples[w] == m_fft_win_size) {
                m_num_windowed_samples[w] = 0;
                if (w == 1 && m_odd_phase_window_is_bogus) {
                    // discard the first odd phase window; it only had half the samples
                    m_odd_phase_window_is_bogus = false;
                    break;
                }

                fftwf_execute(m_plan[w]);

                // process power from each fft bin, and look for the pulse with largest SNR

                int best;
                float highest_probe_signal = 0;
                bool have_pulse = false;
                float bin_probe;

                for (int j = m_first_freq_bin; j <= m_last_freq_bin; ++j) {
                    int jj = j;
                    if (jj < 0)
                        jj = m_fft_win_size + jj;

                    // for each bin, process power in that bin through the pulse finder
                    
                    float pwr = m_fft[jj][0] * m_fft[jj][0] + m_fft[jj][1] * m_fft[jj][1];

                    if (m_freq_bin_pulse_finder[jj].process(pwr) 
                        && m_freq_bin_pulse_finder[jj].pulse_SNR() >= m_min_pulse_SNR
                        && (bin_probe = m_freq_bin_pulse_finder[jj].pulse_signal()) >= highest_probe_signal ) {
                        // found a peak in this bin (a while back; current power value confirms it)
                        highest_probe_signal = bin_probe;
                        best = j;
                        have_pulse = true;
                    }
                }            

                if (have_pulse) {
                    // estimate the frequency more finely, then see whether it's within the appropriate range

                    float pwr[4];  // for power in adjacent bins
                    if (best == m_first_freq_bin)
                        ++best;
                    while (best >= m_last_freq_bin - 1)
                        -- best;

                    float bin_est;
                    if (best > m_first_freq_bin) {
                        for (int i = -1; i < 3; ++i) {
                            int jj = best + i;
                            if (jj < 0)
                                jj = jj + m_fft_win_size;
                            pwr[i + 1] = m_freq_bin_pulse_finder[jj].pulse_signal();
                        }
                        
                        bin_est = best + cubicMaximize(pwr[0], pwr[1], pwr[2], pwr[3]);
                        if (bin_est > m_fft_win_size / 2.0)
                            bin_est -= m_fft_win_size;
                    } else {
                        bin_est = best;
                    }

                    if (fabs(bin_est - best) < 1.0 ) {
                        float freq = bin_est * ((float) m_inputSampleRate / m_fft_win_size) / 1000;
                        if (freq >= m_min_freq && freq <= m_max_freq) {
                            
                            // dump the feature
                            Feature feature;
                            feature.hasTimestamp = true;
                            feature.hasDuration = false;
                    
                            // The pulse timestamp is taken to be the centre of the fft window
                    
                            feature.timestamp = timestamp +
                                Vamp::RealTime::frame2RealTime((signed int) i - m_fft_win_size * (1 + (m_noise_win_size + m_min_pulse_sep) * m_pf_size + m_pf_size / 2) / 2.0 +2, (size_t) m_inputSampleRate);                    

                            std::stringstream ss;
                            ss.precision(5);

                            if (best < 0)
                                best += m_fft_win_size;

                            ss << " freq: " << freq
                               << " kHz; SNR: " << 10 * log10(m_freq_bin_pulse_finder[best].pulse_SNR())
                               << " dB; sig: " << 10 * log10((m_freq_bin_pulse_finder[best].pulse_signal() - m_freq_bin_pulse_finder[best].pulse_noise()) / m_fft_win_size)
                               << " dB; noise: " << 10 * log10(m_freq_bin_pulse_finder[best].pulse_noise() / m_fft_win_size)
                               << " dB;";
                            
                            ss.precision(2);
                            if (m_last_timestamp[best].sec >= 0) {
                                Vamp::RealTime gap =  feature.timestamp - m_last_timestamp[best];
                                if (gap.sec < 1) {
                                    ss.precision(0);
                                    ss << "; Gap: " << gap.msec() << " ms";
                                } else {
                                    ss.precision(1);
                                    ss << "; Gap: " << gap.sec + (double) gap.msec()/1000 << " s";
                                }
                            }
                            m_last_timestamp[best] = feature.timestamp;
                            feature.label = ss.str();
                            returnFeatures[0].push_back(feature);
                        }
                    }
                }
            }
        }
    }
    return returnFeatures;
}

FindPulseFD::FeatureSet
FindPulseFD::getRemainingFeatures()
{
    return FeatureSet();
}

float FindPulseFD::m_default_plen = 2.5; // milliseconds
float FindPulseFD::m_default_min_pulse_SNR_dB = 5; // dB
int FindPulseFD::m_default_fft_win_size = 128; // 0.5 milliseconds @ 48kHz
int FindPulseFD::m_default_noise_win_size = 5; // pulse lengths
int FindPulseFD::m_default_min_pulse_sep = 1; //pulse lengths
float FindPulseFD::m_default_min_freq = -2.0; // 2 kHz
float FindPulseFD::m_default_max_freq = 8.0; // 24 kHz
