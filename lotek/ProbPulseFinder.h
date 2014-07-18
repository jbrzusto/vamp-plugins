/* -*- mode: c++ -*- */
/*
    Copyright 2012 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _PROB_PULSE_FINDER_H
#define _PROB_PULSE_FINDER_H

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
// power, while that with the two low side segments gives the "noise" or background.

#include <boost/circular_buffer.hpp>
#include "KahanSum.h"

template < typename DATATYPE >
class ProbPulseFinder {
 public:
  ProbPulseFinder (int recalc_interval=131072, size_t n = 1, size_t m = 1, size_t k = 1) :
    m_pulse_width(n),
    m_bkgd_width(m),
    m_pulse_sep(k),
    m_back(m+n),
    m_sig_scale (1.0 / n),
    m_bkgd_scale (1.0 / (2 * m)),
    m_sample_buf (n + 2*m),
    m_probe_signal_buf (2 * k + 1),
    m_probe_bkgd_buf (2 * k + 1),
    m_probe_signal2_buf (2 * k + 1),
    m_probe_bkgd2_buf (2 * k + 1),
    m_signal (0),
    m_bkgd (0),
    m_signal2 (0),
    m_bkgd2 (0),
    m_bkgd_floor(2.511886432E-10 * 2 * m),
    m_max_probe_index (-1),
    m_recalc_countdown(recalc_interval),
    m_recalc_interval(recalc_interval)
  {
  };

  DATATYPE pulse_signal () { // the signal value for the pulse; 
    // only valid immediately after operator() returns true
    return m_pulse_signal * m_sig_scale;
	m_pulse_signal = m_probe_signal_buf[m_max_probe_index];
	m_pulse_bkgd = m_probe_bkgd_buf[m_max_probe_index];
	m_pulse_signal2 = m_probe_signal2_buf[m_max_probe_index];
	m_pulse_bkgd2 = m_probe_bkgd2_buf[m_max_probe_index];
  }

  DATATYPE pulse_bkgd () { // the signal value for the pulse; only valid if got_pulse() is true
    return m_pulse_bkgd * m_bkgd_scale;
  }

  DATATYPE pulse_Z () { // the Z score for the pulse; only valid if got_pulse() is true
    return Z(m_pulse_signal, m_pulse_bkgd, m_pulse_signal2, m_pulse_bkgd2);
  }

  inline DATATYPE Z(DATATYPE signal, DATATYPE bkgd, DATATYPE signal2, DATATYPE bkgd2) const {
    // z score for signal - bkgd
    // signal2 and bkgd2 are sums of squares of power values

    // Use pooled variance:
    // DATATYPE N = m_pulse_width + m_bkgd_width;
    // DATATYPE s2 = signal2 + bkgd2;
    // DATATYPE s1 = signal + bkgd;
    // DATATYPE sd = sqrt(N * s2 - s1 * s1) / N;

    // Use variance of bkgd window only:
    double sd = sqrt(2.0 * m_bkgd_width * bkgd2 - bkgd * bkgd) * m_bkgd_scale;
    return (signal * m_sig_scale - bkgd * m_bkgd_scale) / sd;
  };

  bool operator () (DATATYPE d) { // process a value from the data stream
    bool got_pulse = false;

    if (m_sample_buf.full()) {
      // the first sample is moving from the left bkgd zone out of the probe window
      m_bkgd -= m_sample_buf[0];
      m_bkgd2 -= m_sample_buf[0] * m_sample_buf[0];
    }

    int n = m_sample_buf.size();

    if (n >= m_bkgd_width) {
      // a sample is moving from the right bkgd zone to the central signal zone
      DATATYPE mov = m_sample_buf[n - m_bkgd_width];
      DATATYPE mov2 = mov * mov;
        
      m_bkgd  -= mov;
      m_bkgd2  -= mov2;
      m_signal += mov;
      m_signal2 += mov2;
    }
    if (n >= m_back) {
      // a sample is moving from the central signal zone to the left bkgd zone
      DATATYPE mov = m_sample_buf[n - m_back];
      DATATYPE mov2 = mov * mov;

      m_signal -= mov;
      m_signal2 -= mov2;
      m_bkgd  += mov;
      m_bkgd2  += mov2;
    }

    // the new sample moves into the right bkgd zone
    m_bkgd += d;
    m_bkgd2 += d * d;
    m_sample_buf.push_back(d);

    // we might need to do a full recalculation of running sums,
    // since pulse has just left the entire window.

      -- m_recalc_countdown;
      if (m_recalc_countdown <= 0) {
	m_recalc_countdown = m_recalc_interval;
	m_signal = m_bkgd = m_signal2 = m_bkgd2 = 0.0;
	for (int j = 0; j < m_bkgd_width; ++j) {
	  m_bkgd += m_sample_buf[j] + m_sample_buf[j + m_back];
          m_bkgd2 += m_sample_buf[j] * m_sample_buf[j] + m_sample_buf[j + m_back] * m_sample_buf[j + m_back];
        }
	for (int j = m_bkgd_width; j < m_back; ++j) {
	  m_signal += m_sample_buf[j];
          m_signal2 += m_sample_buf[j] * m_sample_buf[j];
        }
      }
    

    // Due to rounding errors, m_signal or m_bkgd might be negative,
    // even though incoming data are all non-negative.
    // fix this!

      m_signal = std::max((DATATYPE) 0.0, (DATATYPE) m_signal);
      m_bkgd = std::max(m_bkgd_floor,  (DATATYPE) m_bkgd);

    if (m_sample_buf.full()) {
      // we have a full probe value; push it into the probe buffer
      // and recalculate the buffer maximum

      // get the squares of signal and bkgd

      // note whether the current max position must move
      bool max_will_move = m_probe_signal_buf.full();

      m_probe_signal_buf.push_back(m_signal);
      m_probe_bkgd_buf.push_back(m_bkgd);
      m_probe_signal2_buf.push_back(m_signal2);
      m_probe_bkgd2_buf.push_back(m_bkgd2);

      DATATYPE z = Z(m_signal, m_bkgd, m_signal*m_signal, m_bkgd*m_bkgd);

      if (m_max_probe_index < 0 || z >= m_max_probe_Z) {
	m_max_probe_Z = z;
	m_max_probe_index = m_probe_signal_buf.size() - 1;
      } else if (max_will_move) {
	// the location of the max will change
	if (m_max_probe_index == 0) {
	  // the max was shifted out so rescan
	  m_max_probe_index = m_probe_signal_buf.size() - 1;
	  m_max_probe_Z = z;
	  for (int i = m_probe_signal_buf.size() - 2; i >= 0; --i) {
	    z = Z(m_probe_signal_buf[i], m_probe_bkgd_buf[i], 
                  m_probe_signal2_buf[i], m_probe_bkgd2_buf[i]);
	    if (z > m_max_probe_Z) {
	      m_max_probe_index = i;
	      m_max_probe_Z = z;
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
	got_pulse = true;
      }
    }
    return got_pulse;
  }
      
protected:

  int m_pulse_width;
  int m_bkgd_width;
  int m_pulse_sep;
  int m_back;
  DATATYPE m_sig_scale;
  DATATYPE m_bkgd_scale;
  boost::circular_buffer < DATATYPE > m_sample_buf;
  boost::circular_buffer < DATATYPE > m_probe_signal_buf;
  boost::circular_buffer < DATATYPE > m_probe_bkgd_buf;
  boost::circular_buffer < DATATYPE > m_probe_signal2_buf;
  boost::circular_buffer < DATATYPE > m_probe_bkgd2_buf;
  KahanSum < DATATYPE > m_signal;
  KahanSum < DATATYPE > m_bkgd;
  KahanSum < DATATYPE > m_signal2;
  KahanSum < DATATYPE > m_bkgd2;
  DATATYPE m_bkgd_floor;
  int m_max_probe_index;
  DATATYPE m_max_probe_Z;
  bool m_got_pulse;
  int m_recalc_countdown;
  int m_recalc_interval;
};

#endif //  _PROB_PULSE_FINDER_H
