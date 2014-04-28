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

    StatPulse.cpp - find low-probability pulses with width in a specified range.

    Copyright 2012 John Brzustowski

    License: GPL v 2.0 or later.  This is required in order to use fftw.

*/

#ifndef _STATPULSE_PLUGIN_H_
#define _STATPULSE_PLUGIN_H_

#include <boost/circular_buffer.hpp>
#include "vamp-sdk/Plugin.h"
#include "PulseFinder.h"
#include "FreqEstimator.h"
#include <complex>
#include <fftw3.h>
#include <boost/circular_buffer.hpp>
#include <cmath>
#include <sstream>

#ifdef MINGW
#define exp10(X) powf(10, X)
#define fftw_free(X) fftwf_free(X)
#endif

class StatPulse : public Vamp::Plugin
{
public:

    StatPulse(float inputSampleRate);
    virtual ~StatPulse();

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

    float m_min_plen; // minimum pulse length, in milliseconds
    float m_max_plen; // maximum pulse length, in milliseconds
    float m_bkgd_len; // length of background window, in milliseconds
    float m_max_pulse_prob; // maximum P-value for a pulse to be accepted
    float m_min_freq;  // only accept pulses whose estimated offset frequency is at least this (can be negative)
    float m_max_freq;  // only accept pulses whose estimated offset frequency is at most this (can be negative)
    float m_accept_raw_samples; // defaults to 0; set to 1.0 by host when it can provide raw samples

    // parameter defaults
    static float m_default_min_plen;
    static float m_default_max_plen;
    static float m_default_bkgd_len;
    static float m_default_max_pulse_prob;
    static float m_default_min_freq;
    static float m_default_max_freq;

    // internal registers
    int m_plen_samples;  // maximum pulse length, in samples
    int m_bkgd_samples;  // length of background window, in samples

    float m_probe_scale; // divisor to convert raw probe value to power
    float m_min_probe; // scaled value of m_min_pulse_power_dB
    int m_pf_size; // size of peak finder moving average window, in samples
    Vamp::RealTime  m_last_timestamp; // timestamp of previous pulse

    // the following members are used to calculate a finer estimate of dfreq once a pulse has been found
    std::vector < float > m_pulse_window; // windowing function for FFT on pulse candidates

    boost::circular_buffer < int16_t > m_sample_buf; // ring buffer of time domain samples from each channel, interleaved as I/C complex pair; fixme: hardcoded S16_LE

    int m_num_samples;  // number of samples processed

    PulseFinder < double > m_pulse_finder;

    FreqEstimator freq_est;
    static const char * fftw_wisdom_filename;
};


#endif // _STATPULSE_PLUGIN_H_
