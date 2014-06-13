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
    
  FindNarrowBandPulsePlugin.cpp - probabilistic time domain pulse finder

  Copyright 2014 John Brzustowski

  License: GPL v 2.0 or later.

*/

#include "FindNarrowBandPulsePlugin.h"

using std::stringstream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;

FindNarrowBandPulsePlugin::FindNarrowBandPulsePlugin(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_stepSize(0),
    m_blockSize(0),
    m_plen(m_default_plen),
    m_min_freq (m_default_min_freq),
    m_max_freq (m_default_max_freq),
    m_min_SNR_dB(m_default_min_SNR_dB),
    m_min_Z(m_default_min_Z),
    m_max_noise_for_Z (undB(m_default_max_noise_for_Z_dB)),
    m_batch_host (false),
    m_nbpf(0)
{
    
}

FindNarrowBandPulsePlugin::~FindNarrowBandPulsePlugin()
{
    if (m_nbpf)
        delete m_nbpf;
}

string
FindNarrowBandPulsePlugin::getIdentifier() const
{

    return "findnbpulse";
}

string
FindNarrowBandPulsePlugin::getName() const
{
    return "Narrow Band Pulse Finder (SNR or Z)";
}

string
FindNarrowBandPulsePlugin::getDescription() const
{
    return "Find pulses in a narrow portion of the spectrum with minimum SNR or Z-score";
}

string
FindNarrowBandPulsePlugin::getMaker() const
{
    return "sensorgnome.org  jbrzusto@fastmail.fm";
}

int
FindNarrowBandPulsePlugin::getPluginVersion() const
{
    return 1;
}

string
FindNarrowBandPulsePlugin::getCopyright() const
{
    return "(C) 2014 John Brzustowski; Licence: GPL version 2 or later";
}

bool
FindNarrowBandPulsePlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_stepSize = stepSize;
    m_blockSize = blockSize;

    m_plen_samples = (m_plen / 1000.0) * m_inputSampleRate;

    m_freq_scale = m_inputSampleRate / (float) m_plen_samples / 1000.0; // DFT bin width, in kHz, given DFT of size m_plen_samples and sampling rate
    m_min_bin = floor(m_min_freq / m_freq_scale);
    m_max_bin = ceil(m_max_freq / m_freq_scale);

    m_nbpf = new NarrowBandSpectralPulseFinder (m_plen_samples, m_min_bin, m_max_bin, m_min_SNR_dB, m_min_Z, m_max_noise_for_Z);

    return true;
}

void
FindNarrowBandPulsePlugin::reset()
{
}

FindNarrowBandPulsePlugin::ParameterList
FindNarrowBandPulsePlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "plen";
    d.name = "Pulse Length (unit: milliseconds)";
    d.description = "Duration of a transmitted pulse in milliseconds";
    d.unit = "milliseconds";
    d.minValue = 0.1;
    d.maxValue = 500;
    d.defaultValue = FindNarrowBandPulsePlugin::m_default_plen;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "minfreq";
    d.name = "Lowest Tag Offset Frequency (unit: kHz)";
    d.description = "Lowest allowed frequency by which tag differs from receiver, in kHz";
    d.unit = "kHz";
    d.minValue = - m_inputSampleRate / 2000;
    d.maxValue = m_inputSampleRate / 2000;
    d.defaultValue = FindNarrowBandPulsePlugin::m_default_min_freq;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxfreq";
    d.name = "Highest Tag Offset Frequency (unit: kHz)";
    d.description = "Highest allowed frequency by which tag differs from receiver, in kHz";
    d.unit = "kHz";
    d.minValue = - m_inputSampleRate / 2000;
    d.maxValue = m_inputSampleRate / 2000;
    d.defaultValue = FindNarrowBandPulsePlugin::m_default_max_freq;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "minsnr";
    d.name = "Minimum Signal to Noise Ratio";
    d.description = "Minimum ratio of signal (with bkgd subtracted) to bkgd, in dB";
    d.unit = "dB";
    d.minValue = -30;
    d.maxValue = 100;
    d.defaultValue = FindNarrowBandPulsePlugin::m_default_min_SNR_dB;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "minz";
    d.name = "Minimum Z Score";
    d.description = "Minimum Z score of pulse vs background; only in effect at low noise";
    d.unit = "(none)";
    d.minValue = 0;
    d.maxValue = 1000;
    d.defaultValue = FindNarrowBandPulsePlugin::m_default_min_Z;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxnoisez";
    d.name = "Maximum Noise For Using Z-Score (unit: dB)";
    d.description = "Noise Level (dB) below which Z score can be used instead of SNR";
    d.unit = "(none)";
    d.minValue = -100;
    d.maxValue = -30;
    d.defaultValue = FindNarrowBandPulsePlugin::m_default_max_noise_for_Z_dB;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "__batch_host__";
    d.name = "IGNORE: set automatically by audacity or vamp-alsa-host";
    d.description = "set to 1 when host needs batch output";
    d.unit = "";
    d.minValue = 0;
    d.maxValue = 1;
    d.defaultValue = FindNarrowBandPulsePlugin::m_batch_host;
    d.isQuantized = true;
    d.quantizeStep = 1.0;
    list.push_back(d);

    return list;
}

float
FindNarrowBandPulsePlugin::getParameter(string id) const
{
    if (id == "plen") {
        return m_plen;
    } else if (id == "minsnr") {
        return m_min_SNR_dB;
    } else if (id == "maxnoisez") {
        return dB(m_max_noise_for_Z);
    } else if (id == "minz") {
        return m_min_Z;
    } else if (id == "minfreq") {
        return m_min_freq;
    } else if (id == "maxfreq") {
        return m_max_freq;
    } else if (id == "__batch_host__") {
        return m_batch_host;
    } else {
        throw std::runtime_error("invalid parameter name");
    }
}

void
FindNarrowBandPulsePlugin::setParameter(string id, float value)
{
    if (id == "plen") {
        FindNarrowBandPulsePlugin::m_default_plen = m_plen = value;
    } else if (id == "minfreq") {
        FindNarrowBandPulsePlugin::m_default_min_freq = m_min_freq = value;
    } else if (id == "maxfreq") {
        FindNarrowBandPulsePlugin::m_default_max_freq = m_max_freq = value;
    } else if (id == "minsnr") {
        FindNarrowBandPulsePlugin::m_default_min_SNR_dB = m_min_SNR_dB = value;
    } else if (id == "maxnoisez") {
        FindNarrowBandPulsePlugin::m_default_max_noise_for_Z_dB = value;
        m_max_noise_for_Z = undB(value);
    } else if (id == "minz") {
        FindNarrowBandPulsePlugin::m_default_min_Z = m_min_Z = value;
    } else if (id == "maxnoiseforz") {
        FindNarrowBandPulsePlugin::m_default_max_noise_for_Z_dB = m_max_noise_for_Z = value;
    } else if (id == "__batch_host__") {
        // kludge: parameter that affects whether this plugin
        // produces output for a batch-style host (e.g. vamp-alsa-host)
        // or for display in a GUI-style host (e.g. audacity)
        // The default value for m_batch_host is false, so it will stay
        // thus unless a host is aware of this parameter and sets it to
        // true.
        FindNarrowBandPulsePlugin::m_batch_host = value;
    } else {
        throw std::runtime_error("invalid parameter name");
    }
}


FindNarrowBandPulsePlugin::OutputList
FindNarrowBandPulsePlugin::getOutputDescriptors() const
{
    OutputList list;

    OutputDescriptor zc;

    zc.identifier = "pulses";
    zc.name = "Pulses";
    zc.description = "The locations and features of pulses";
    zc.unit = "";
    zc.hasFixedBinCount = true;
    zc.binCount = m_batch_host ? 3 : 0;
    zc.sampleType = OutputDescriptor::VariableSampleRate;
    zc.sampleRate = m_inputSampleRate;
    list.push_back(zc);

    return list;
}

FindNarrowBandPulsePlugin::FeatureSet
FindNarrowBandPulsePlugin::process(const float *const *inputBuffers,
                     Vamp::RealTime timestamp)
{
    FeatureSet returnFeatures;

    if (m_stepSize == 0) {
	cerr << "ERROR: FindNarrowBandPulsePlugin::process: "
	     << "FindNarrowBandPulsePlugin has not been initialised"
	     << endl;
	return returnFeatures;
    }

    for (unsigned i=0; i < m_blockSize; ++i) {
        // grab sample as complex I/Q pair
        
        float I = inputBuffers[0][i];
        float Q = inputBuffers[1][i];

        std::complex < float > sample (I, Q);

        // send it to the pulse finder

        if ((*m_nbpf) (sample)) {
            // found a pulse

            // how many samples back was the centre of the pulse?
            // this is only as precise as the fft step (fft_size - overlap)

            int centre = (int) i - (int) m_nbpf->location();
            Vamp::RealTime ts = timestamp + Vamp::RealTime::frame2RealTime(centre, (size_t) m_inputSampleRate);
 
            // dump the feature
            Feature feature;
            feature.hasTimestamp = true;
            feature.hasDuration = false;
            
            // The pulse timestamp is taken to be the centre of the pulse window
            
            feature.timestamp = ts;
     
            float sig = dB (m_nbpf->signal());
            float noise = dB (m_nbpf->noise());
            float freq = m_nbpf->freq() * m_freq_scale;

            if (m_batch_host) {
                feature.values.push_back(freq);
                feature.values.push_back(sig);
                feature.values.push_back(noise);
            } else {
                std::stringstream ss;
                ss.precision(5);
                    
                ss << "freq: " << freq << " (kHz)"
                   << "; pwr: " << sig << " dB"
                   << "; bgkd: " << noise  << " dB";
            
                feature.label = ss.str();
            }
            returnFeatures[0].push_back(feature);
        }
    }
    return returnFeatures;
}

FindNarrowBandPulsePlugin::FeatureSet
FindNarrowBandPulsePlugin::getRemainingFeatures()
{
    return FeatureSet();
}


float FindNarrowBandPulsePlugin::m_default_plen = 2.5; // milliseconds
float FindNarrowBandPulsePlugin::m_default_min_freq = 2.0; // 2 kHz
float FindNarrowBandPulsePlugin::m_default_max_freq = 6.0; // 6 kHz
double FindNarrowBandPulsePlugin::m_default_min_SNR_dB = 10; // minimum SNR 
double FindNarrowBandPulsePlugin::m_default_min_Z = 150; // z-score
double FindNarrowBandPulsePlugin::m_default_max_noise_for_Z_dB = -70; // z-score
