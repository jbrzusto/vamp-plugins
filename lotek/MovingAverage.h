/* -*- mode: c++ -*- */
/*
    Copyright 2011 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _MOVING_AVERAGE_H
#define _MOVING_AVERAGE_H

#include <boost/circular_buffer.hpp>

template < typename DATATYPE, typename AVGTYPE > 
class MovingAverage {
 public:
  MovingAverage (size_t n = 1) :
  buf(n),
  m_total(0)
    {
    };

  void clear() {
    m_total = 0;
    buf.clear();
  }

  bool have_average () {  // is a moving average available yet?
    return buf.full();
  }

  AVGTYPE get_average () { // get the moving average
    if (! buf.full())
      throw std::runtime_error("Not enough data for moving average");
    return m_total / buf.capacity();
  };

  AVGTYPE get_buffer_average () { // get the average of items seen so far
    if (buf.empty())
      throw std::runtime_error("No data in buffer");
    return m_total / (static_cast < AVGTYPE > (buf.size()));
  };

  void process(DATATYPE d) { // process a value from the data stream
    if (buf.full())
      m_total -= buf.front();
    buf.push_back(d);
    m_total += d;
  }


 protected:
  boost::circular_buffer < DATATYPE > buf;
  AVGTYPE m_total;
};

#endif //  _MOVING_AVERAGE_H
