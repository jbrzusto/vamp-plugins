/* -*- mode: c++ -*- */
/*
    Copyright 2011-2014 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _KAHAN_MOVING_AVERAGER_H
#define _KAHAN_MOVING_AVERAGER_H

#include <boost/circular_buffer.hpp>
#include "KahanSum.h"

template < typename DATATYPE, typename AVGTYPE > 
class KahanMovingAverager {
 public:
  KahanMovingAverager (size_t n = 1) :
  buf(n),
  m_total(0)
    {
    };

  void clear() {
    m_total = 0;
    buf.clear();
  }

  operator AVGTYPE () { // get the moving average
    if (buf.empty())
      throw std::runtime_error("No data in buffer");
    return m_total / buf.size();
  };

  AVGTYPE sum() { // get the sum
    return m_total;
  };

  AVGTYPE rem() { // get the remainder being carried
    return m_total.rem();
  };

  bool operator() (const DATATYPE & d) { 
    // process a value from the data stream
    // return true if a (full window) moving average is available
    // (always true once we've seen n samples)
    if (buf.full())
      m_total -= buf.front();
    buf.push_back(d);
    m_total += d;
    return buf.full();
  };

  bool is_linearized() {
    return buf.is_linearized();
  };

 protected:
  boost::circular_buffer < DATATYPE > buf;
  KahanSum < AVGTYPE > m_total;
};

#endif //  _KAHAN_MOVING_AVERAGER_H
