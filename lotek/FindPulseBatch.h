/* -*- mode:c++; c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp

    An API for audio analysis and feature extraction plugins.

    Centre for Digital Music, Queen Mary, University of London.
    Copyright 2006 Chris Cannam.

    FindPulseBatch.h 
    Copyright 2011 John Brzustowski

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

    This plugin applies the same algorithm as FindPulse, but generates
    outputs in number form suitable for use by Sonic Annotator, rather
    than by Audacity.

*/

#ifndef _FIND_PULSE_BATCH_PLUGIN_H_
#define _FIND_PULSE_BATCH_PLUGIN_H_

#include <boost/circular_buffer.hpp>
#include "vamp-sdk/Plugin.h"
#include "MovingAverager.h"
#include <complex>
/**
 * Look for pulses from Lotek tags
*/

class FindPulseBatch : public Vamp::Plugin
{
public:
    FindPulseBatch(float inputSampleRate);
    virtual ~FindPulseBatch();

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
    
    // paramters
    float m_min_duration;
    float m_max_duration;
    float m_min_power;
    float m_max_freq_rsd;
    int m_pulse_power_win_size;
    int m_bkgd_power_win_size;

    // parameter defaults
    static float m_default_min_duration;
    static float m_default_max_duration;
    static float m_default_min_power;
    static float m_default_max_freq_rsd;
    static int m_default_pulse_power_win_size;
    static int m_default_bkgd_power_win_size;

    // internal registers
    float m_last_phase;
    bool m_in_pulse;
    int m_duration;
    float m_total_power;
    float m_power_ratio;
    int m_max_samples;
    MovingAverager < float, float > m_pulse_power_ma;
    MovingAverager < float, float > m_bkgd_power_ma;
    MovingAverager < std::complex < float >, std::complex < float > > m_pulse_phasor_ma;

    // constants
    static const float NO_LAST_PHASE = -10;
};


#endif // _FIND_PULSE_BATCH_PLUGIN_H_
