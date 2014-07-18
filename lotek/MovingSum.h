/* -*- mode: c++ -*- */
/*
    Copyright 2014 John Brzustowski

    Licence: GPL version 2 or later

    This version does a full recalculation each time
    the buffer wraps.  Seems like a reasonable trade-off
    between efficiency and protection against catastrophic
    loss of precision.

*/

#ifndef _MOVING_SUM_H
#define _MOVING_SUM_H

#include <boost/circular_buffer.hpp>

template < typename DATATYPE> 
class MovingSum {
 public:
  MovingSum (size_t n = 1) :
    n(n),
    buf(n),
    m_total(0)
  {
  };
  
  void clear() {
    m_total = 0;
    buf.clear();
  };

  operator DATATYPE () { // get the moving sum
    return m_total;
  };

  DATATYPE mean() {
    return m_total / buf.size();
  };

  bool operator () (DATATYPE d) {
    // process a value from the data stream
    // return true if a full window has been summed

    if (buf.full()) {
      if (buf.array_one().second == 1) {
        // the new element will go into the last physical slot, creating
        // a linearized buffer.  Insert the new element and do a full recalc
        buf.push_back(d);
        DATATYPE total = 0.0;
        DATATYPE * p = & buf.front();
        for (int i = 0; i < n; ++i)
          total += p[i];
        m_total = total;
        return true;
        
      }
      // first remove the oldest element
      m_total -= buf.front();
    }
    buf.push_back(d);
    m_total += d;
    return buf.full();
  }

 protected:
  int n;
  boost::circular_buffer < DATATYPE > buf;
  DATATYPE m_total;
};

#endif //  _MOVING_SUM_H
