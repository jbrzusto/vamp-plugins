/* -*- mode: c++ -*- */
/*
    Copyright 2014 John Brzustowski

    Licence: GPL version 2 or later

    Keep a running sum of floating point numbers, using the
    Kahan algorithm to reduce effects of catastrophic loss
    of precision.

    Reference: http://en.wikipedia.org/wiki/Kahan_summation

*/

#ifndef _KAHAN_SUM_H
#define _KAHAN_SUM_H

template < typename DATATYPE >
class KahanSum {
 public:
  KahanSum (DATATYPE initial = 0.0) {
    sum = initial;
    carry = 0.0;
  };

#pragma GCC push_options
#pragma GCC optimize ("O0")


    // NB: don't allow GCC O3 optimization, since that will
    // optimize out the Kahan summation

  KahanSum & operator += (const DATATYPE & x) {

    if (std::abs(sum) >= std::abs(x)) {
      carry += x;
      DATATYPE newsum = sum + carry;
      carry = carry - (newsum - sum);
      sum = newsum;
    } else {
      // new entry is larger than sum, swap roles.
      DATATYPE newsum = x + sum;
      carry += sum - (newsum - x);
      sum = newsum;
    }
    return * this;
  };

#pragma GCC pop_options

  KahanSum & operator -= (const DATATYPE & x) {
    return *this += (-x);
  };

  KahanSum & operator = (const DATATYPE & x) {
    sum = x;
    carry = 0.0;
    return *this;
  };


  operator DATATYPE () {
    return sum;
  };

  DATATYPE rem() {
    return carry;
  };

 private:
  DATATYPE sum;
  DATATYPE carry;
};

#endif // _KAHAN_SUM_H
