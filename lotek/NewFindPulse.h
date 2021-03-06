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

    NewFindPulse.h - probabilistic time domain pulse finder

    Copyright 2014 John Brzustowski

    License: GPL v 2.0 or later.  This is required in order to use fftw.

*/

#ifndef _NEW_FIND_PULSE_PLUGIN_H_
#define _NEW_FIND_PULSE_PLUGIN_H_

#include "vamp-sdk/Plugin.h"
#include "PulseDetector.h"
#include <complex>
#include <cmath>
#include <sstream>

#ifdef MINGW
#define exp10(X) powf(10, X)
#define fftw_free(X) fftwf_free(X)
#endif

/**
 * Look for pulses from Lotek tags - Frequency Domain version
*/

class NewFindPulse : public Vamp::Plugin
{
public:

    NewFindPulse(float inputSampleRate);
    virtual ~NewFindPulse();

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
    float m_min_plen;        // minimum length of pulse we're trying to detect, in ms
    float m_max_plen;        // minimum length of pulse we're trying to detect, in ms

    double m_max_pulse_prob; // probability of pulse must be at most this to be accepted
    double m_min_power_diff; // minimum power difference at both pulse edges in order to be accepted (linear units)

    float m_min_freq;  // only accept pulses whose offset frequency is at least this (kHz)
    float m_max_freq;  // only accept pulses from bins whose offset frequency is at most this (kHz)

    // parameter defaults
    static float m_default_min_plen;
    static float m_default_max_plen;
    static double m_default_max_pulse_prob;
    static double m_default_min_power_diff;
    static float m_default_min_freq;
    static float m_default_max_freq;

    // internal registers
    int m_min_plen_samples;
    int m_max_plen_samples;

    // sample buffer, with interleaved floats from two channels, from which fft estimates will be calculated
    boost::circular_buffer < float > m_sample_buf; 

    MovingAverager < float, double > m_dcma[2]; // moving averager for removing DC on each channel

    // pulse detector
    PulseDetector < float > m_pd;
    
};


#endif // _NEW_FIND_PULSE_PLUGIN_H_
