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
    
  NewFindPulse.cpp - probabilistic time domain pulse finder

  Copyright 2014 John Brzustowski

  License: GPL v 2.0 or later.  This is required in order to use fftw.

*/

#include "NewFindPulse.h"

using std::stringstream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;

NewFindPulse::NewFindPulse(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_stepSize(0),
    m_blockSize(0),
    m_min_plen(m_default_min_plen),
    m_max_plen(m_default_max_plen),
    m_max_pulse_prob(m_default_max_pulse_prob),
    m_min_power_diff(m_default_min_power_diff),
    m_min_freq (m_default_min_freq),
    m_max_freq (m_default_max_freq)
{
}

NewFindPulse::~NewFindPulse()
{
}

string
NewFindPulse::getIdentifier() const
{
    return "NewFindPulse";
}

string
NewFindPulse::getName() const
{
    return "Probabilistic Time Domain Pulse Finder";
}

string
NewFindPulse::getDescription() const
{
    return "Find pulses (e.g. from Lotek telemetry tags)";
}

string
NewFindPulse::getMaker() const
{
    return "sensorgnome.org jbrzusto@fastmail.fm";
}

int
NewFindPulse::getPluginVersion() const
{
    return 1;
}

string
NewFindPulse::getCopyright() const
{
    return "GPL version 2 or later";
}

bool
NewFindPulse::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_stepSize = stepSize;
    m_blockSize = blockSize;

    m_min_plen_samples = (m_min_plen / 1000.0) * m_inputSampleRate;
    m_max_plen_samples = (m_max_plen / 1000.0) * m_inputSampleRate;

    // cap frequency limits at Nyquist
    if (m_min_freq > m_inputSampleRate / 2000)
        m_min_freq = m_inputSampleRate / 2000;
    if (m_max_freq > m_inputSampleRate / 2000)
        m_max_freq = m_inputSampleRate / 2000;
    
    m_dcma[0] = MovingAverager < float, double > (100 * m_min_plen_samples);
    if (m_channels == 2)
        m_dcma[1] = MovingAverager < float, double > (100 * m_min_plen_samples);

    m_pd = PulseDetector < float > (m_min_plen_samples, m_max_plen_samples, m_min_power_diff * 32767 * 32767, m_max_pulse_prob);

    m_sample_buf = boost::circular_buffer < float > (2 * (m_max_plen_samples + (m_min_plen_samples * 2 + 1)));

    return true;
}

void
NewFindPulse::reset()
{
}

NewFindPulse::ParameterList
NewFindPulse::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "minplen";
    d.name = "Minimum Pulse Length (unit: milliseconds)";
    d.description = "Minimum duration of a transmitted pulse in milliseconds";
    d.unit = "milliseconds";
    d.minValue = 0.1;
    d.maxValue = 500;
    d.defaultValue = NewFindPulse::m_default_min_plen;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxplen";
    d.name = "Maximum Pulse Length (unit: milliseconds)";
    d.description = "Maximum duration of a transmitted pulse in milliseconds";
    d.unit = "milliseconds";
    d.maxValue = 0.1;
    d.maxValue = 500;
    d.defaultValue = NewFindPulse::m_default_max_plen;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxprob";
    d.name = "Maximum Pulse Probability (threshold; e.g. 0.00001)";
    d.description = "Maximum NULL hypothesis probability of a pulse in order to be accepted";
    d.unit = "(none)";
    d.minValue = 0;
    d.maxValue = 0.05;
    d.defaultValue = NewFindPulse::m_default_max_pulse_prob;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "minpowdiff";
    d.name = "Minimum Power Difference for Pulse";
    d.description = "Minimum power rise/fall (in linear units) at pulse edges";
    d.unit = "linear power units";
    d.minValue = 0;
    d.maxValue = 1;
    d.defaultValue = NewFindPulse::m_default_min_power_diff;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "minfreq";
    d.name = "Minimum Tag Offset Frequency (unit: kHz)";
    d.description = "Minimum frequency by which tag differs from receiver, in kHz";
    d.unit = "kHz";
    d.minValue = - m_inputSampleRate / 2000;
    d.maxValue = m_inputSampleRate / 2000;
    d.defaultValue = NewFindPulse::m_default_min_freq;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxfreq";
    d.name = "Maximum Tag Offset Frequency (unit: kHz)";
    d.description = "Maximum frequency by which tag differs from receiver, in kHz";
    d.unit = "kHz";
    d.minValue = - m_inputSampleRate / 2000;
    d.maxValue = m_inputSampleRate / 2000;
    d.defaultValue = NewFindPulse::m_default_max_freq;
    d.isQuantized = false;
    list.push_back(d);

    return list;
}

float
NewFindPulse::getParameter(string id) const
{
    if (id == "minplen") {
        return m_min_plen;
    } else if (id == "maxplen") {
        return m_max_plen;
    } else if (id == "maxprob") {
        return m_max_pulse_prob;
    } else if (id == "minpowdiff") {
        return m_min_power_diff;
    } else if (id == "minfreq") {
        return m_min_freq;
    } else if (id == "maxfreq") {
        return m_max_freq;
    }
    return 0.f;
}

void
NewFindPulse::setParameter(string id, float value)
{
    if (id == "minplen") {
        NewFindPulse::m_default_min_plen = m_min_plen = value;
    } else if (id == "maxplen") {
        NewFindPulse::m_default_max_plen = m_max_plen = value;
    } else if (id == "maxprob") {
        NewFindPulse::m_default_max_pulse_prob = m_max_pulse_prob = value;
    } else if (id == "mindiff") {
        NewFindPulse::m_default_min_power_diff = m_min_power_diff = value;
    } else if (id == "minfreq") {
        NewFindPulse::m_default_min_freq = m_min_freq = value;
    } else if (id == "maxfreq") {
        NewFindPulse::m_default_max_freq = m_max_freq = value;
    }
}


NewFindPulse::OutputList
NewFindPulse::getOutputDescriptors() const
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

NewFindPulse::FeatureSet
NewFindPulse::process(const float *const *inputBuffers,
                     Vamp::RealTime timestamp)
{
    FeatureSet returnFeatures;

    if (m_stepSize == 0) {
	cerr << "ERROR: NewFindPulse::process: "
	     << "NewFindPulse has not been initialised"
	     << endl;
	return returnFeatures;
    }

    for (unsigned int i=0; i < m_blockSize; ++i) {
        // get DC-corrected power values

        float pwr = 0.0;
        for (unsigned short ch = 0; ch < m_channels; ++ch) {
            // update DC offset
            m_dcma[ch].process(inputBuffers[ch][i]);
            if (m_dcma[ch].have_average()) {
                // correct DC offset and calculate contribution to power
                float val = inputBuffers[ch][i] - m_dcma[ch].get_average();
                m_sample_buf.push_back(val);
                pwr += val * val;
            }
        }

        if (! m_dcma[0].have_average())
            continue;
        
        // process power through pulse detector

        if (m_pd(pwr)) {
            // we detected a pulse figure out power etc.
            
            // dump the feature
            Feature feature;
            feature.hasTimestamp = true;
            feature.hasDuration = false;
            
            // The pulse timestamp is taken to be the centre of the pulse window

            float offset = m_pd.location() - m_pd.width() / 2.0;
            feature.timestamp = timestamp + Vamp::RealTime::frame2RealTime(i - offset, (size_t) m_inputSampleRate);
        
            // get pulse power by summing values in pulse; FIXME: need to get BKGD levels
            // to subract from this; PulseDetector needs to get these from Edge Detector
            // ALSO, PulseDetector should sum power in pulse to get average, keep low-side
            // power average for both edges.

            auto a1 = m_pd.pulse_array1();
            auto a2 = m_pd.pulse_array2();

            float pulsepwr = 0.0;

            if (a1.first && a1.second) {
                for (unsigned i=0; i < a1.second; ++i)
                    pulsepwr += a1.first[i];
            }
            if (a2.first && a2.second) {
                for (unsigned i=0; i < a2.second; ++i)
                    pulsepwr += a2.first[i];
            }
            std::stringstream ss;
            ss.precision(5);
 
            ss << "pwr: " <<  10 * log10(pulsepwr / m_pd.width()) << " dB"
               << "; bgkd: " <<  10 * log10(m_pd.bg()) << " dB"
               << ";  dur: " << m_pd.width() / (double) m_inputSampleRate * 1000.0 << " msec";
            
            feature.label = ss.str();
            returnFeatures[0].push_back(feature);
        }
    }
    return returnFeatures;
}

NewFindPulse::FeatureSet
NewFindPulse::getRemainingFeatures()
{
    return FeatureSet();
}

float NewFindPulse::m_default_min_plen = 2.2; // milliseconds
float NewFindPulse::m_default_max_plen = 2.8; // milliseconds
double NewFindPulse::m_default_min_power_diff = 1e-2; // linear power units
double NewFindPulse::m_default_max_pulse_prob = 1e-4; // probability
float NewFindPulse::m_default_min_freq = -4.0; // -4 kHz
float NewFindPulse::m_default_max_freq =  4.0; // +4 kHz
