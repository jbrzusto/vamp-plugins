/* -*- mode: c++ -*- */
/*
    Copyright 2011 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _PEAK_FINDER_H
#define _PEAK_FINDER_H

// PeakFinder: given a window of odd-size 2n+1, search an input sequence for 
//             peaks in its moving average.  A peak occurs at slot i when 
//             MA(i-n) < MA(i-n+1) <= ... <= MA(i) >= MA(i+1) >= ... >= MA(i+n)
//             where MA(m) = (x[m-n] + x[m-n+1] + ... + x[m+n]) / (2n+1)

#include "MovingAverager.h"

template < typename DATATYPE >
class PeakFinder : public MovingAverager < DATATYPE, DATATYPE > {
 public:
  PeakFinder (size_t n = 1) :
    MovingAverager < DATATYPE, DATATYPE > (n | 1),
    m_win_size(n | 1),
    m_num_rises(0),
    m_num_falls(0),
    m_rising(false),
    m_falling(false),
    m_have_prev(false),
    m_got_peak(false)
  {
  };

  bool got_peak () {  // have we seen a peak?  True if the moving average went up then down in the most recent three consecutive calls to process()
    return m_got_peak;
  }

  DATATYPE peak_val () { // the value of peak; only valid if got_peak() is true
    return m_peak_val;
  }

  void process(DATATYPE d) { // process a value from the data stream
    
    MovingAverager < DATATYPE, DATATYPE > :: process(d);
    m_got_peak = false;
    if (MovingAverager < DATATYPE, DATATYPE > :: have_average()) {
      DATATYPE cur_val =  MovingAverager < DATATYPE, DATATYPE > :: get_average();
      if (m_have_prev) {
	if (m_rising) {
	  if (cur_val >= m_prev) {
	    if (m_num_rises < m_win_size / 2)
	      ++m_num_rises;
	  } else {
	    m_peak_val = m_prev;
	    m_rising = false;
	    if (m_num_rises >= m_win_size / 2) {
	      m_falling = true;
	      m_num_falls = 1;
	    }
	  }
	} else if (m_falling) {
	  if (cur_val <= m_prev) {
	    if (++m_num_falls >= m_win_size / 2) {
	      // we got a peak
	      m_got_peak = true;
	      m_falling = false;
	      m_rising = true;
	    }
	  } else {
	    m_falling = false;
	    if (cur_val > m_prev) {
	      m_rising = true;
	      m_num_rises = 1;
	    }
	  }
	} else {
	  // both m_falling and m_rising are false
	  if (cur_val > m_prev) {
	    m_rising = true;
	    m_num_rises = 1;
	  }
	}
      }
      m_prev = cur_val;
      m_have_prev = true;
    }
  }

 protected:
  int m_win_size;
  int m_num_rises;
  int m_num_falls;
  bool m_rising; // true when we're rising
  bool m_falling; // true when we've risen sufficiently and hit a peak
  DATATYPE m_prev;
  bool m_have_prev;
  bool m_got_peak;
  DATATYPE m_peak_val;
};

#endif //  _PEAK_FINDER_H
