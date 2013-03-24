/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp

    An API for audio analysis and feature extraction plugins.

    Centre for Digital Music, Queen Mary, University of London.
    Copyright 2006 Chris Cannam.

    AveragePower.cpp
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

#include "AveragePower.h"

using std::stringstream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;

#include <cmath>
#include <sstream>

AveragePower::AveragePower(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_stepSize(0),
    m_blockSize(0),
    m_power_win_size(AveragePower::m_default_power_win_size),
    m_decim_rate(AveragePower::m_default_decim_rate)
{
}

AveragePower::~AveragePower()
{
}

string
AveragePower::getIdentifier() const
{
    return "averagepower";
}

string
AveragePower::getName() const
{
    return "Power Moving Average";
}

string
AveragePower::getDescription() const
{
    return "Output the moving average of power";
}

string
AveragePower::getMaker() const
{
    return "flightcalls.org jbrzusto@fastmail.fm";
}

int
AveragePower::getPluginVersion() const
{
    return 1;
}

string
AveragePower::getCopyright() const
{
    return "GPL version 2 or later";
}

bool
AveragePower::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_stepSize = stepSize;
    m_blockSize = blockSize;
    m_power_ma = MovingAverager < float, float > ((size_t) m_power_win_size);
    m_decim_counter = m_decim_rate;
    return true;
}

void
AveragePower::reset()
{
}

AveragePower::ParameterList
AveragePower::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "powwinsize";
    d.name = "Size of power-averaging window";
    d.description = "Number of samples in moving-average window for power";
    d.unit = "samples";
    d.minValue = 1;
    d.maxValue = 1024;
    d.defaultValue = AveragePower::m_default_power_win_size;
    d.isQuantized = true;
    d.quantizeStep = 1;
    list.push_back(d);

    d.identifier = "decimrate";
    d.name = "Decimation rate";
    d.description = "Rate at which power average is reported: 1 = each sample, 2 = every 2 samples, etc.";
    d.unit = "samples";
    d.minValue = 1;
    d.maxValue = 1024;
    d.defaultValue = AveragePower::m_default_decim_rate;
    d.isQuantized = true;
    d.quantizeStep = 1;
    list.push_back(d);

    return list;
}

float
AveragePower::getParameter(string id) const
{
    if (id == "powwinsize") {
        return m_power_win_size;
    } else if (id == "decimrate") {
        return m_decim_rate;
    }
    return 0.f;
}

void
AveragePower::setParameter(string id, float value)
{
    if (id == "powwinsize") {
        AveragePower::m_default_power_win_size = m_power_win_size = value;
    } else if (id == "decimrate") {
        AveragePower::m_default_decim_rate = m_decim_counter = m_decim_rate = value;
    }
}


AveragePower::OutputList
AveragePower::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor zc;

    zc.identifier = "averagepower";
    zc.name = "Power Moving Average ";
    zc.description = "The moving average of power over the given sample window size, decimated";
    zc.unit = "";
    zc.hasFixedBinCount = true;
    zc.binCount = 1;
    zc.sampleType = OutputDescriptor::FixedSampleRate;
    zc.sampleRate = m_inputSampleRate / m_decim_rate;
    list.push_back(zc);

    return list;
}

AveragePower::FeatureSet
AveragePower::process(const float *const *inputBuffers,
                      Vamp::RealTime timestamp)
{
    FeatureSet returnFeatures;

    if (m_stepSize == 0 || m_decim_counter == 0) {
	cerr << "ERROR: AveragePower::process: "
	     << "AveragePower has not been initialised"
	     << endl;
	return returnFeatures;
    }

    for (unsigned int i=0; i < m_blockSize; ++i) {
        float sample_power = inputBuffers[0][i] * inputBuffers[0][i];
        if (m_channels > 1) 
            sample_power +=  inputBuffers[1][i] * inputBuffers[1][i];
        m_power_ma.process(sample_power);
        if (--m_decim_counter == 0) {
            m_decim_counter = m_decim_rate;
            Feature feature;
            feature.hasTimestamp = false;
            feature.hasDuration = false;
            feature.values.clear();
            feature.values.push_back(10 * log10f(m_power_ma.get_average()));
            returnFeatures[0].push_back(feature);
        }
    }
    return returnFeatures;
}

AveragePower::FeatureSet
AveragePower::getRemainingFeatures()
{
    return FeatureSet();
}

int AveragePower::m_default_power_win_size = 24; // = 0.25ms at 96kHz
int AveragePower::m_default_decim_rate = 24;     // = 0.25ms at 96kHz
