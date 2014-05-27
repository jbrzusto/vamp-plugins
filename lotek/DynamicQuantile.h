/* -*- mode: c++ -*- */
/*
    Copyright 2014 John Brzustowski

    Licence: GPL version 2 or later
*/

#ifndef _DYNAMIC_QUANTILE_H
#define _DYNAMIC_QUANTILE_H

#include <math.h>

#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/tail.hpp>
#include <boost/accumulators/statistics/tail_quantile.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>


using namespace boost::accumulators;
using namespace boost;

// DynamicQuantile: maintain an estimate of a quantile over the
// (approximately) n most recent items in a datastream.
//
// We implement this using the boost p_square_quantile accumulator.
// Because p_square_quantile does not allow 'unlearning' old data, we
// maintain a staggered, overlapping set of p_square_quantiles,
// switching to the 2nd oldest one and emptying the oldest one at
// regular intervals.  We provide a parameter 'overlap' which specifies the
// fraction of overlap in samples between the oldest and 2nd oldest
// p_square_quantile when this switchover occurs.  The closer p is to
// 1, the more 'smooth' the switchover, but the larger the storage and
// (especially) processing costs.
//
// Storage: O(n) * sizeof (TYPE)
// Complexity: process() is O(log(n)); other operations constant.

typedef accumulator_set < float, stats < tag::p_square_quantile > > quant_est;

template < typename TYPE >

class DynamicQuantile {

 public:

  DynamicQuantile (double prob, size_t n = 100000, double overlap = 0.8) :
    prob(prob),
    n(n),
    m(n * (1. - overlap)),
    rot_count(0),
    num_estimators(ceil(1. / (1. - overlap))),
    oldest(0)
  {
    quant_estimators.push_back (quant_est(quantile_probability = prob));
  };

  void operator() (TYPE d) { // process a value from the data stream
    if (rot_count == m) {
      // it's time to rotate estimators, or at least to start a fresh one
      rot_count = 0;

      // if we have a full set of active estimators, rotate 
      if (quant_estimators.size() == num_estimators) {
        // start a fresh estimator in the slot of the oldest
        quant_estimators[oldest] = quant_est(quantile_probability = prob);
        ++ oldest;
        if (oldest == num_estimators)
          oldest = 0;
      } else {
        // otherwise, add a fresh estimator to the back of the vector
        quant_estimators.push_back (quant_est(quantile_probability = prob));
      }
    } else {
      ++ rot_count;
    }
      
    // add the new element to all active estimators

    for (typename std::vector < quant_est > ::iterator qe = quant_estimators.begin(); qe != quant_estimators.end(); ++qe)
      (*qe)(d);
  }

  operator TYPE () {
    // return the quantile estimate from the oldest estimator
    return p_square_quantile(quant_estimators[oldest]);
  };

 protected:
  double prob;
  int n; // approximate number of (most recent) samples to maintain quantile for
  int m; // number of samples after which we rotate quantile estimators, if possible
  int rot_count; // count samples up to a rotation
  int num_estimators; // number of estimators needed to achieve desired overlap
  std::vector < quant_est > quant_estimators;
  int oldest; // index of oldest quantile in quant_estimators
};

#endif //  _DYNAMIC_QUANTILE_H
