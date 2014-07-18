/* -*- mode: c++ -*- */
/*
    Copyright 2014 John Brzustowski

    Licence: GPL version 2 or later

    Ring buffer where each element is a fixed-size buffer of given type;
    Fixed storage for the vectors is allocated at creation.
    Analogues to the circular_buffer operation take/return pointers to
    fixed-size buffers, either internal or external.
    FIXME: only push_back() and operator[] are implemented.

*/

#ifndef _VECTOR_RING_BUFFER_H
#define _VECTOR_RING_BUFFER_H

#include <boost/circular_buffer.hpp>
#include <memory.h>

template < typename DATATYPE > 

class VectorRingBuffer {
  
public:
  VectorRingBuffer (int numvec, int vecsize) :
    numvec (numvec),
    vecsize (vecsize),
    eltbuf (new DATATYPE [numvec * vecsize]),
    ptrbuf (numvec),
    count (0)
  {
    // initialize circular buffer of pointers to DATATYPE
    // each pointer is to the first element of the vector in
    // slot i.
    // These vectors are simply rotated in the circular buffer.
    for (int i = 0; i < numvec; ++i)
      ptrbuf.push_back(& eltbuf[vecsize * i]);
  };
  
  ~VectorRingBuffer () {
    delete [] eltbuf;
  };

  bool push_back (DATATYPE * in) {
    // copy in a new vector; return true
    // when buffer is full
    memcpy (ptrbuf[0], in, vecsize * sizeof(DATATYPE) );

    // rotate the pointer buffer
    ptrbuf.rotate(++ ptrbuf.begin());

    if (full())
      return true;
    ++count;
    return false;
  };

  DATATYPE * operator[] (int i) {
    return ptrbuf[i];
  };

  bool full() {
    return count == numvec;
  };

protected:
  int numvec;
  int vecsize;
  DATATYPE * eltbuf;
  boost::circular_buffer < DATATYPE * > ptrbuf;
  int count; // actual number of vectors in buffer
};
    
#endif // _VECTOR_RING_BUFFER_H
