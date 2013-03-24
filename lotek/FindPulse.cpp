/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Vamp

    An API for audio analysis and feature extraction plugins.

    Centre for Digital Music, Queen Mary, University of London.
    Copyright 2006 Chris Cannam.

    FindPulse.cpp
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
*/

#include "FindPulse.h"

using std::stringstream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;

#include <cmath>
#include <sstream>

FindPulse::FindPulse(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_stepSize(0),
    m_blockSize(0),
    m_min_duration(FindPulse::m_default_min_duration),
    m_max_duration(FindPulse::m_default_max_duration),
    m_min_power(FindPulse::m_default_min_power),
    m_max_freq_rsd(FindPulse::m_default_max_freq_rsd),
    m_pulse_power_win_size(FindPulse::m_default_pulse_power_win_size),
    m_bkgd_power_win_size(FindPulse::m_default_bkgd_power_win_size)
{
}

FindPulse::~FindPulse()
{
}

string
FindPulse::getIdentifier() const
{
    return "findpulse";
}

string
FindPulse::getName() const
{
    return "Find Pulses";
}

string
FindPulse::getDescription() const
{
    return "Find pulses (e.g. from Lotek telemetry tags)";
}

string
FindPulse::getMaker() const
{
    return "flightcalls.org jbrzusto@fastmail.fm";
}

int
FindPulse::getPluginVersion() const
{
    return 1;
}

string
FindPulse::getCopyright() const
{
    return "GPL version 2 or later";
}

bool
FindPulse::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_stepSize = stepSize;
    m_blockSize = blockSize;
    m_duration = 0;
    m_total_power = 0;
    m_in_pulse = false;
    m_last_timestamp = Vamp::RealTime();
    m_pulse_power_ma = MovingAverager < float, float > ((size_t) m_pulse_power_win_size);
    m_bkgd_power_ma = MovingAverager < float, float > ((size_t) m_bkgd_power_win_size);
    m_power_ratio = powf(10, m_min_power / 10);
    m_max_samples = m_max_duration / 1000.0 * m_inputSampleRate;
    m_pulse_phasor_ma = MovingAverager < std::complex < float > , std::complex < float > > ((size_t) m_max_samples + 1);

    return true;
}

void
FindPulse::reset()
{
    m_duration = 0;
    m_total_power = 0;
    m_in_pulse = false;
    m_last_timestamp = Vamp::RealTime();
}

FindPulse::ParameterList
FindPulse::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "minduration";
    d.name = "Minimum pulse duration";
    d.description = "Minimum duration of a pulse in milliseconds; shorter pulses are discarded";
    d.unit = "milliseconds";
    d.minValue = 0.1;
    d.maxValue = 50;
    d.defaultValue = FindPulse::m_default_min_duration;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxduration";
    d.name = "Maximum pulse duration";
    d.description = "Maximum duration of a pulse in milliseconds; longer pulses are discarded";
    d.unit = "milliseconds";
    d.minValue = 0.1;
    d.maxValue = 50;
    d.defaultValue = FindPulse::m_default_max_duration;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "minpower";
    d.name = "Power threshold (db above background)";
    d.description = "a pulse consists of consecutive samples with power at least this much above background, in dB";
    d.unit = "dB";
    d.minValue = 0.001;
    d.maxValue = 100;
    d.defaultValue = FindPulse::m_default_min_power;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxfreqrsd";
    d.name = "Max frequency relative SD";
    d.description = "Maximum relative standard deviation in frequency for a pulse to qualify as real";
    d.unit = "";
    d.minValue = 0;
    d.maxValue = 1;
    d.defaultValue = FindPulse::m_default_max_freq_rsd;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "pulsepowwinsize";
    d.name = "Size of power-averaging window for pulse";
    d.description = "Number of samples in moving-average window for pulse power";
    d.unit = "samples";
    d.minValue = 1;
    d.maxValue = 1024;
    d.defaultValue = FindPulse::m_default_pulse_power_win_size;
    d.isQuantized = true;
    d.quantizeStep = 1;
    list.push_back(d);

    d.identifier = "bkgdpowwinsize";
    d.name = "Size of power-averaging window for background";
    d.description = "Number of samples in moving-average window for background power";
    d.unit = "samples";
    d.minValue = 1;
    d.maxValue = 4096;
    d.defaultValue = FindPulse::m_default_bkgd_power_win_size;
    d.isQuantized = true;
    d.quantizeStep = 1;
    list.push_back(d);

    return list;
}

float
FindPulse::getParameter(string id) const
{
    if (id == "minduration") {
        return m_min_duration;
    } else if (id == "maxduration") {
        return m_max_duration;
    } else if (id == "minpower") {
        return m_min_power;
    } else if (id == "maxfreqrsd") {
        return m_max_freq_rsd;
    } else if (id == "pulsepowwinsize") {
        return m_pulse_power_win_size;
    } else if (id == "bkgdpowwinsize") {
        return m_bkgd_power_win_size;
    }
    return 0.f;
}

void
FindPulse::setParameter(string id, float value)
{
    if (id == "minduration") {
        FindPulse::m_default_min_duration = m_min_duration = value;
    } else if (id == "maxduration") {
        FindPulse::m_default_max_duration = m_max_duration = value;
    } else if (id == "minpower") {
        FindPulse::m_default_min_power = m_min_power = value;
    } else if (id == "maxfreqrsd") {
        FindPulse::m_default_max_freq_rsd = m_max_freq_rsd = value;
    } else if (id == "pulsepowwinsize") {
        FindPulse::m_default_pulse_power_win_size = m_pulse_power_win_size = value;
    } else if (id == "bkgdpowwinsize") {
        FindPulse::m_default_bkgd_power_win_size = m_bkgd_power_win_size = value;
    }
}


FindPulse::OutputList
FindPulse::getOutputDescriptors() const
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

FindPulse::FeatureSet
FindPulse::process(const float *const *inputBuffers,
                      Vamp::RealTime timestamp)
{
    FeatureSet returnFeatures;

    if (m_stepSize == 0) {
	cerr << "ERROR: FindPulse::process: "
	     << "FindPulse has not been initialised"
	     << endl;
	return returnFeatures;
    }

    for (unsigned int i=0; i < m_blockSize; ++i) {
        float sample_power = inputBuffers[0][i] * inputBuffers[0][i];
        if (m_channels > 1) 
            sample_power +=  inputBuffers[1][i] * inputBuffers[1][i];
        m_pulse_power_ma.process(sample_power);
        if (! m_in_pulse)
            m_bkgd_power_ma.process(sample_power);
        if (m_pulse_power_ma.have_average() && m_bkgd_power_ma.have_average()) {
            float power = m_pulse_power_ma.get_average();
            float power_thresh = m_bkgd_power_ma.get_average() * m_power_ratio;
            if (m_in_pulse && (power < power_thresh || m_duration > m_max_samples)) {
                bool pulse_valid = false;
                m_in_pulse = false;
                // we've finished a pulse; check whether it meets filtering criteria
                float duration_msec = m_duration / m_inputSampleRate * 1000;
                if (duration_msec >= m_min_duration && duration_msec <= m_max_duration && m_duration > 1) {
                    float rsd = sqrtf(1.0 - abs(m_pulse_phasor_ma.get_buffer_average()));
                    if (rsd < m_max_freq_rsd) {
                        pulse_valid = true;
                        Feature feature;
                        feature.hasTimestamp = true;
                        feature.hasDuration = true;
                        feature.duration = Vamp::RealTime::fromSeconds(duration_msec / 1000.0);
                        feature.timestamp = timestamp +
                            Vamp::RealTime::frame2RealTime((signed int) i - m_duration / 2, (size_t)m_inputSampleRate);
                        std::stringstream ss;
                        ss.precision(3);
                        ss << "dFreq: " <<
                            arg(m_pulse_phasor_ma.get_buffer_average()) / (2 * M_PI) * m_inputSampleRate / 1000 << " kHz; SNR: " <<
                            10 * log10f(m_total_power / (m_bkgd_power_ma.get_average() * m_duration)) << " dB; Dur: " <<
                            duration_msec << " ms; RSD: ";
                        ss.precision(2);
                        ss << rsd;
                        if (m_last_timestamp.sec > 0) {
                            Vamp::RealTime gap =  feature.timestamp - m_last_timestamp;
                            if (gap.sec < 1) {
                                ss.precision(0);
                                ss << "; Gap: " << gap.msec() << " ms";
                            } else {
                                ss.precision(1);
                                ss << "; Gap: " << gap.sec + (double) gap.msec()/1000 << " s";
                            }
                        }
                        m_last_timestamp = feature.timestamp;
                        feature.label = ss.str();
                        returnFeatures[0].push_back(feature);
                    }
                }
                if (! pulse_valid) {
                        // add this invalid pulse's power to background
                    float average_power = m_total_power / m_duration;
                    for (int i=m_duration; i > 0; --i) 
                        m_bkgd_power_ma.process(average_power);
                }
            } else if (power >= power_thresh && power_thresh > 0) {
                if (! m_in_pulse) {
                    m_total_power = 0;
                    m_duration = 0;
                    m_last_phase = NO_LAST_PHASE;
                    m_in_pulse = true;
                    m_pulse_phasor_ma.clear();
                }
                float phase = (m_channels > 1) ? atan2f(inputBuffers[1][i], inputBuffers[0][i]) : 0;
                if (m_last_phase != NO_LAST_PHASE) {
                    float d_phase = phase - m_last_phase;
                    m_pulse_phasor_ma.process(std::complex < float > (cosf(d_phase), sinf(d_phase)));
                }
                m_last_phase = phase;
                m_total_power += sample_power;
                ++m_duration;
            }
        }
    }
    return returnFeatures;
}

FindPulse::FeatureSet
FindPulse::getRemainingFeatures()
{
    return FeatureSet();
}

float FindPulse::m_default_min_duration = 1;
float FindPulse::m_default_max_duration = 5;
float FindPulse::m_default_min_power = 3;
float FindPulse::m_default_max_freq_rsd = 0.7;
int FindPulse::m_default_pulse_power_win_size = 24;
int FindPulse::m_default_bkgd_power_win_size = 1024;
