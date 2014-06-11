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
    
  SNRZFindPulse.cpp - probabilistic time domain pulse finder

  Copyright 2014 John Brzustowski

  License: GPL v 2.0 or later.  This is required in order to use fftw.

*/

#include "FindPulseTDSNRZ.h"

using std::stringstream;
using std::string;
using std::vector;
using std::cerr;
using std::endl;

FindPulseTDSNRZ::FindPulseTDSNRZ(float inputSampleRate) :
    Plugin(inputSampleRate),
    m_stepSize(0),
    m_blockSize(0),
    m_plen(m_default_plen),
    m_bkgd(m_default_bkgd),
    m_min_Z(m_default_min_Z),
    m_min_SNR_dB(m_default_min_SNR_dB),
    m_max_noise_for_Z (undB(m_default_max_noise_for_Z_dB)),
    m_min_freq (m_default_min_freq),
    m_max_freq (m_default_max_freq),
    m_batch_host (false),
    m_fpd(0),
    m_fest(0),
    m_dc_offset_I(0.0),
    m_dc_offset_Q(0.0)
{
}

FindPulseTDSNRZ::~FindPulseTDSNRZ()
{
    if (m_fpd)
        delete m_fpd;
    if (m_fest)
        delete m_fest;
}

string
FindPulseTDSNRZ::getIdentifier() const
{
    // note: to allow use by legacy deployment files,
    // we give this the same name as was used by the very
    // different previous incarnation into early 2014

    return "findpulsetdsnrz";
}

string
FindPulseTDSNRZ::getName() const
{
    return "Time Domain Pulse Finder (SNR or Z)";
}

string
FindPulseTDSNRZ::getDescription() const
{
    return "Find pulses in the time domain with minimum SNR or Z-score";
}

string
FindPulseTDSNRZ::getMaker() const
{
    return "sensorgnome.org  jbrzusto@fastmail.fm";
}

int
FindPulseTDSNRZ::getPluginVersion() const
{
    return 1;
}

string
FindPulseTDSNRZ::getCopyright() const
{
    return "GPL version 2 or later";
}

bool
FindPulseTDSNRZ::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

    m_channels = channels;
    m_stepSize = stepSize;
    m_blockSize = blockSize;

    m_plen_samples = (m_plen / 1000.0) * m_inputSampleRate;
    m_bkgd_samples = (m_bkgd / 1000.0) * m_inputSampleRate;

    // cap frequency limits at Nyquist
    if (m_min_freq > m_inputSampleRate / 2000)
        m_min_freq = m_inputSampleRate / 2000;
    if (m_max_freq > m_inputSampleRate / 2000)
        m_max_freq = m_inputSampleRate / 2000;
    
    m_bin_step = m_inputSampleRate / (1000.0 * m_num_bins);
    m_min_bin = floor(m_min_freq / m_bin_step) ;
    m_max_bin = ceil(m_max_freq / m_bin_step) ;
    
    m_num_seek_bins = m_max_bin - m_min_bin + 1;

    m_fpd = new FixedPulseDetector < float > (m_plen_samples, m_bkgd_samples, undB(m_min_SNR_dB), m_min_Z, undB(m_max_noise_for_Z));

    m_samp_buff = boost::circular_buffer < std::complex < float > > (m_plen_samples + 2 * m_bkgd_samples);

    m_fest = new FreqEstimator (m_plen_samples * DEFAULT_FFT_PADDING);

    return true;
}

void
FindPulseTDSNRZ::reset()
{
}

FindPulseTDSNRZ::ParameterList
FindPulseTDSNRZ::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "plen";
    d.name = "Pulse Length (unit: milliseconds)";
    d.description = "Duration of a transmitted pulse in milliseconds";
    d.unit = "milliseconds";
    d.minValue = 0.1;
    d.maxValue = 500;
    d.defaultValue = FindPulseTDSNRZ::m_default_plen;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "bkgd";
    d.name = "Background Window Length (single-sided; unit: milliseconds)";
    d.description = "Duration of the background window on each side of the pulse";
    d.unit = "milliseconds";
    d.minValue = 0.1;
    d.maxValue = 500;
    d.defaultValue = FindPulseTDSNRZ::m_default_bkgd;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "minsnr";
    d.name = "Minimum Signal to Noise Ratio";
    d.description = "Minimum ratio of signal (with bkgd subtracted) to bkgd, in dB";
    d.unit = "dB";
    d.minValue = -40;
    d.maxValue = 100;
    d.defaultValue = FindPulseTDSNRZ::m_default_min_SNR_dB;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxnoisez";
    d.name = "Maximum Noise For Using Z-Score (unit: dB)";
    d.description = "Noise Level (dB) below which Z score can be used instead of SNR";
    d.unit = "(none)";
    d.minValue = -100;
    d.maxValue = -30;
    d.defaultValue = FindPulseTDSNRZ::m_default_max_noise_for_Z_dB;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "minz";
    d.name = "Minimum Z Score";
    d.description = "Minimum Z score of pulse vs background; only in effect at low noise";
    d.unit = "(none)";
    d.minValue = 0;
    d.maxValue = 1000;
    d.defaultValue = FindPulseTDSNRZ::m_default_min_Z;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "minfreq";
    d.name = "Minimum Tag Offset Frequency (unit: kHz)";
    d.description = "Minimum frequency by which tag differs from receiver, in kHz";
    d.unit = "kHz";
    d.minValue = - m_inputSampleRate / 2000;
    d.maxValue = m_inputSampleRate / 2000;
    d.defaultValue = FindPulseTDSNRZ::m_default_min_freq;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "maxfreq";
    d.name = "Maximum Tag Offset Frequency (unit: kHz)";
    d.description = "Maximum frequency by which tag differs from receiver, in kHz";
    d.unit = "kHz";
    d.minValue = - m_inputSampleRate / 2000;
    d.maxValue = m_inputSampleRate / 2000;
    d.defaultValue = FindPulseTDSNRZ::m_default_max_freq;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "__batch_host__";
    d.name = "IGNORE: set automatically by audacity or vamp-alsa-host";
    d.description = "set to 1 when host needs batch output";
    d.unit = "";
    d.minValue = 0;
    d.maxValue = 1;
    d.defaultValue = FindPulseTDSNRZ::m_batch_host;
    d.isQuantized = true;
    d.quantizeStep = 1.0;
    list.push_back(d);

    return list;
}

float
FindPulseTDSNRZ::getParameter(string id) const
{
    if (id == "plen") {
        return m_plen;
    } else if (id == "bkgd") {
        return m_bkgd;
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
FindPulseTDSNRZ::setParameter(string id, float value)
{
    if (id == "plen") {
        FindPulseTDSNRZ::m_default_plen = m_plen = value;
    } else if (id == "bkgd") {
        FindPulseTDSNRZ::m_default_bkgd = m_bkgd = value;
    } else if (id == "minsnr") {
        FindPulseTDSNRZ::m_default_min_SNR_dB = m_min_SNR_dB = value;
    } else if (id == "maxnoisez") {
        FindPulseTDSNRZ::m_default_max_noise_for_Z_dB = value;
        m_max_noise_for_Z = undB(value);
    } else if (id == "minz") {
        FindPulseTDSNRZ::m_default_min_Z = m_min_Z = value;
    } else if (id == "minfreq") {
        FindPulseTDSNRZ::m_default_min_freq = m_min_freq = value;
    } else if (id == "maxfreq") {
        FindPulseTDSNRZ::m_default_max_freq = m_max_freq = value;
    } else if (id == "maxnoiseforz") {
        FindPulseTDSNRZ::m_default_max_noise_for_Z_dB = m_max_noise_for_Z = value;
    } else if (id == "__batch_host__") {
        // kludge: parameter that affects whether this plugin
        // produces output for a batch-style host (e.g. vamp-alsa-host)
        // or for display in a GUI-style host (e.g. audacity)
        // The default value for m_batch_host is false, so it will stay
        // thus unless a host is aware of this parameter and sets it to
        // true.
        FindPulseTDSNRZ::m_batch_host = value;
    } else {
        throw std::runtime_error("invalid parameter name");
    }
}


FindPulseTDSNRZ::OutputList
FindPulseTDSNRZ::getOutputDescriptors() const
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

FindPulseTDSNRZ::FeatureSet
FindPulseTDSNRZ::process(const float *const *inputBuffers,
                     Vamp::RealTime timestamp)
{
    FeatureSet returnFeatures;

    if (m_stepSize == 0) {
	cerr << "ERROR: FindPulseTDSNRZ::process: "
	     << "FindPulseTDSNRZ has not been initialised"
	     << endl;
	return returnFeatures;
    }

    for (unsigned i=0; i < m_blockSize; ++i) {
        // grab sample as complex I/Q pair
        
        float I = inputBuffers[0][i];
        float Q = inputBuffers[1][i];

        m_dc_offset_I = (2047.0 * m_dc_offset_I + I) / 2048.0;
        m_dc_offset_Q = (2047.0 * m_dc_offset_Q + Q) / 2048.0;

        std::complex < float > sample (I - m_dc_offset_I, Q - m_dc_offset_Q);

        // buffer it

        m_samp_buff.push_back(sample);

        float power = sample.real() * sample.real() + sample.imag() * sample.imag();
        // send it to the pulse detector

        if ((*m_fpd) (power)) {
            // found a pulse

            // how many samples back was the centre of the pulse?
            // this is only as precise as the fft step (fft_size - overlap)

            int centre = (int) i - (int) m_fpd->location() + m_plen_samples / 2.0;
            Vamp::RealTime ts = timestamp + Vamp::RealTime::frame2RealTime(centre, (size_t) m_inputSampleRate);
 
            // dump the feature
            Feature feature;
            feature.hasTimestamp = true;
            feature.hasDuration = false;
            
            // The pulse timestamp is taken to be the centre of the pulse window
            
            feature.timestamp = ts;
     
            float sig = dB (m_fpd->signal() - m_fpd->bkgd()); // + m_power_scale_dB;
            float noise = dB (m_fpd->bkgd()); // + m_power_scale_dB;

            // get an estimate of frequency offset
            auto a1 = m_samp_buff.array_one();
            auto a2 = m_samp_buff.array_two();

            int n1, n2;
            n1 = std::min(m_plen_samples, (int) a1.second);
            n2 = m_plen_samples - n1;

            float freq = m_fest->get(a1.first, n1, a2.first, n2) / (double) DEFAULT_FFT_PADDING;

            freq *= m_inputSampleRate / (1000.0 * m_plen_samples);

            if (freq < m_min_freq || freq > m_max_freq)
                continue;
               
            if (m_batch_host) {
                feature.values.push_back(freq);
                feature.values.push_back(sig);
                feature.values.push_back(noise);
            } else {
                std::stringstream ss;
                ss.precision(5);
                    
                ss << "freq: " << freq << " (kHz)"
                   << "; pwr: " << sig << " dB"
                   << "; bgkd: " << noise  << " dB"
                   << "; Z: " << m_fpd->Z();
            
                feature.label = ss.str();
            }
            returnFeatures[0].push_back(feature);
        }
    }
    return returnFeatures;
}

FindPulseTDSNRZ::FeatureSet
FindPulseTDSNRZ::getRemainingFeatures()
{
    return FeatureSet();
}

float FindPulseTDSNRZ::m_default_plen = 2.5; // milliseconds
float FindPulseTDSNRZ::m_default_bkgd = 12.5; // milliseconds

double FindPulseTDSNRZ::m_default_min_SNR_dB = 6; // minimum SNR 
double FindPulseTDSNRZ::m_default_min_Z = 150; // z-score
double FindPulseTDSNRZ::m_default_max_noise_for_Z_dB = -70; // z-score
float FindPulseTDSNRZ::m_default_min_freq = -5.0; // -4 kHz
float FindPulseTDSNRZ::m_default_max_freq =  5.0; // +4 kHz
