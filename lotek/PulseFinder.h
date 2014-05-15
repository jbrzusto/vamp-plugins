/* -*- mode: c++ -*- */
/*
    Copyright 2012 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _PULSE_FINDER_H
#define _PULSE_FINDER_H

// Pulsefinder: search an input sequence for pulses of width n, using this probe function:
//
//            1/n|              ---- (n) ----
//               |
//               |
//               |
//        1/(2m) | ---- (m) ----             ---- (m) ----
//               |
//
// A pulse is flagged every time the probe function convolved with the
// input sequence reaches a maximum over a window of 2k+1 consecutive
// slots.  The convolution with the central high segment gives the "signal"
// power, while that with the two low side segments gives the "noise".

#include <boost/circular_buffer.hpp>
#include "KahanSum.h"

template < typename DATATYPE >
class PulseFinder {
 public:
  PulseFinder (int recalc_interval=131072, size_t n = 1, size_t m = 1, size_t k = 1) :
    m_pulse_width(n),
    m_noise_width(m),
    m_pulse_sep(k),
    m_back(m+n),
    m_div_signal (n),
    m_div_noise (2 * m),
    m_mult_SNR (2 * m / (DATATYPE) n),
    m_sample_buf (n + 2*m),
    m_probe_signal_buf (2 * k + 1),
    m_probe_noise_buf (2 * k + 1),
    m_signal (0),
    m_noise (0),
    m_noise_floor(2.511886432E-10 * 2 * m),
    m_max_probe_index (-1),
    m_got_pulse(false),
    m_recalc_countdown(recalc_interval),
    m_recalc_interval(recalc_interval)
  {
  };

  bool got_pulse () {  // have we seen a pulse? True if the probe value was maximum at the centre of the k-sized window.
    return m_got_pulse;
  }

  DATATYPE pulse_signal () { // the signal value for the pulse; only valid if got_pulse() is true
    return m_pulse_signal / m_div_signal;
  }

  DATATYPE pulse_noise () { // the signal value for the pulse; only valid if got_pulse() is true
    return m_pulse_noise / m_div_noise;
  }

  DATATYPE pulse_SNR () { // the SNR for the pulse; only valid if got_pulse() is true
    return SNR(m_pulse_signal, m_pulse_noise);
  }

  inline DATATYPE SNR(DATATYPE signal, DATATYPE noise) const {
    // signal to noise ratio; note that noise is never allowed to reach 0
    // by logic in process()
    return (signal * m_div_noise / m_div_signal - noise) / noise;
  };

  void process(DATATYPE d) { // process a value from the data stream
    m_got_pulse = false;

    if (m_sample_buf.full())
      // the first sample is moving from the left noise zone out of the probe window
      m_noise -= m_sample_buf[0];

    int n = m_sample_buf.size();

    if (n >= m_noise_width) {
      // a sample is moving from the right noise zone to the central signal zone
      m_noise  -= m_sample_buf[n - m_noise_width];
      m_signal += m_sample_buf[n - m_noise_width];
    }
    if (n >= m_back) {
      // a sample is moving from the central signal zone to the left noise zone
      m_signal -= m_sample_buf[n - m_back];
      m_noise  += m_sample_buf[n - m_back];
    }

    // the new sample moves into the right noise zone
    m_noise += d;
    m_sample_buf.push_back(d);

    // we might need to do a full recalculation of running sums,
    // since pulse has just left the entire window.

      -- m_recalc_countdown;
      if (m_recalc_countdown <= 0) {
	m_recalc_countdown = m_recalc_interval;
	m_signal = m_noise = 0.0;
	for (int j = 0; j < m_noise_width; ++j)
	  m_noise += m_sample_buf[j] + m_sample_buf[j + m_back];
	for (int j = m_noise_width; j < m_back; ++j)
	  m_signal += m_sample_buf[j];
      }
    

    // Due to rounding errors, m_signal or m_noise might be negative,
    // even though incoming data are all non-negative.
    // fix this!

      m_signal = std::max((DATATYPE) 0.0, (DATATYPE) m_signal);
    m_noise = std::max(m_noise_floor,  (DATATYPE) m_noise);

    if (m_sample_buf.full()) {
      // we have a full probe value; push it into the probe buffer
      // and recalculate the buffer maximum

      // note whether the current max position must move
      bool max_will_move = m_probe_signal_buf.full();

      m_probe_signal_buf.push_back(m_signal);
      m_probe_noise_buf.push_back(m_noise);

      DATATYPE snr = SNR(m_signal, m_noise);

      if (m_max_probe_index < 0 || snr >= m_max_probe_SNR) {
	m_max_probe_SNR = snr;
	m_max_probe_index = m_probe_signal_buf.size() - 1;
      } else if (max_will_move) {
	// the location of the max will change
	if (m_max_probe_index == 0) {
	  // the max was shifted out so rescan
	  m_max_probe_index = m_probe_signal_buf.size() - 1;
	  m_max_probe_SNR = snr;
	  for (int i = m_probe_signal_buf.size() - 2; i >= 0; --i) {
	    snr = SNR(m_probe_signal_buf[i], m_probe_noise_buf[i]);
	    if (snr > m_max_probe_SNR) {
	      m_max_probe_index = i;
	      m_max_probe_SNR = snr;
	    }
	  }
	} else {
	  // the max value is unchanged, but the window has advanced one slot
	  -- m_max_probe_index;
	}
      }
      if (m_max_probe_index == m_pulse_sep) {
	// the central value in the probe buffer is the maximum,
	// so indicate we have a pulse.
	m_got_pulse = true;
	m_pulse_signal = m_probe_signal_buf[m_max_probe_index];
	m_pulse_noise = m_probe_noise_buf[m_max_probe_index];
      }
    }
  }
      
protected:

  int m_pulse_width;
  int m_noise_width;
  int m_pulse_sep;
  int m_back;
  DATATYPE m_div_signal;
  DATATYPE m_div_noise;
  DATATYPE m_mult_SNR;
  boost::circular_buffer < DATATYPE > m_sample_buf;
  boost::circular_buffer < DATATYPE > m_probe_signal_buf;
  boost::circular_buffer < DATATYPE > m_probe_noise_buf;
  KahanSum < DATATYPE > m_signal;
  KahanSum < DATATYPE > m_noise;
  DATATYPE m_noise_floor;
  int m_max_probe_index;
  DATATYPE m_max_probe_SNR;
  bool m_got_pulse;
  DATATYPE m_pulse_signal;
  DATATYPE m_pulse_noise;
  int m_recalc_countdown;
  int m_recalc_interval;
};

#endif //  _PULSE_FINDER_H
