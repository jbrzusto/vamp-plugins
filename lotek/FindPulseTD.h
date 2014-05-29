/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

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

#ifndef _FIND_PULSETD_PLUGIN_H_
#define _FIND_PULSETD_PLUGIN_H_

#include <boost/circular_buffer.hpp>
#include "vamp-sdk/Plugin.h"
#include "ProbPulseFinder.h"
#include "MovingAverager.h"
#include <complex>
#include <fftw3.h>
#include <boost/circular_buffer.hpp>
#include <cmath>
#include <sstream>

#ifdef MINGW
#define exp10(X) powf(10, X)
#define fftw_free(X) fftwf_free(X)
#endif

/**
 * Look for pulses from tags - Time Domain version
*/

class FindPulseTD : public Vamp::Plugin
{
public:
    static float cubicMaximize(float y0, float y1, float y2, float y3);
    static float cubicInterpolate(float y0, float y1, float y2, float y3, float x);
    static void generateWindowingCoefficients(int N, std::vector < float > &window, float &win_sum, float &win_sumsq);

    FindPulseTD(float inputSampleRate);
    virtual ~FindPulseTD();

    bool initialise(size_t channels, size_t stepSize, size_t blockSize);
    void reset();

    InputDomain getInputDomain() const { return TimeDomain; }
    size_t getMinChannelCount() const {return 1;}
    size_t getMaxChannelCount() const {return 2;}
    size_t getPreferredStepSize() const {return 1024;}
    size_t getPreferredBlockSize() const {return 1024;}
    std::string getIdentifier() const;
    std::string getName() const;
    std::string getDescription() const;
    std::string getMaker() const;
    int getPluginVersion() const;
    std::string getCopyright() const;

    ParameterList getParameterDescriptors() const;
    float getParameter(std::string id) const;
    void setParameter(std::string id, float value);

    OutputList getOutputDescriptors() const;

    FeatureSet process(const float *const *inputBuffers,
                       Vamp::RealTime timestamp);

    FeatureSet getRemainingFeatures();
    

protected:
    size_t m_channels;
    size_t m_stepSize;
    size_t m_blockSize;
    
    // parameters
    float m_plen;        // length of pulse we're trying to detect, in ms
    int m_plen_in_samples; // length of pulse, measured in samples
    float m_min_pulse_z; // minimum pulse sig-noise z score
    int m_noise_win_size; // size of noise window on each side of pulse, in multiples of pulse length
    int m_min_pulse_sep; // minimum separation between pulses, in multiples of pulse length
    float m_min_freq;  // only accept pulses from bins whose centre frequency is at least this
    float m_max_freq;  // only accept pulses from bins whose centre frequency is at most this

    // parameter defaults
    static float m_default_plen;
    static float m_default_min_pulse_z;
    static int   m_default_noise_win_size;
    static int   m_default_min_pulse_sep;
    static float m_default_min_freq;
    static float m_default_max_freq;

    // internal registers
    float m_probe_scale; // divisor to convert raw probe value to power
    float m_min_probe; // scaled value of m_min_pulse_power_dB
    int m_pf_size; // size of peak finder moving average window (in units of fft windows)
    Vamp::RealTime  m_last_timestamp; // timestamp of previous pulse

    // the following members are used to calculate a finer estimate of dfreq once a pulse has been found
    std::vector < float > m_pulse_window; // windowing function for FFT on pulse candidates

    float *m_windowed_fine[2]; // windowed samples
    boost::circular_buffer < float > m_sample_buf[2]; // ring buffer of time domain samples from each channel
    fftwf_plan m_plan_fine[2]; // FFT plans for pulse samples on each channel
    fftwf_complex *m_fft_fine[2]; // DFT output from pulse samples

    int m_num_samples;  // number of samples processed

    ProbPulseFinder < double > m_pulse_finder;

    MovingAverager < float, float > m_dcma[2]; // moving averager for removing DC on each channel

    static const char * fftw_wisdom_filename;
};


#endif // _FIND_PULSETD_PLUGIN_H_
