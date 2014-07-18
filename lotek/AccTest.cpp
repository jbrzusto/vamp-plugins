/*
  use the boost::accumulators library to find pulses by accumulated stats

  idea: detect pulse 'on' and 'off' edges separately, constraining that a pulse must
  have an on followed by an off within a constrained range, with no intervening on or off.
  Edges are places where
     TotalPower(x[m+1]...x[m+n]) - TotalPower(x[m-(n-1)]...x[m]) has sufficiently
     low P value and is locally maximal (minimal) over a window of size 2n.
     (so we compare power in n samples to the power in the next n samples, look for
     a difference with low probability (so there's an edge-like feature) and seek
     nearby for a potentially better difference (so we match the edge as precisely
     as possible).

     The sequence of values of TotalPower(x[m+1]...x[m+n]) - TotalPower(x[m-(n-1)]...x[m])
     is kept in a distribution, and quantiles are maintained for it.

*/

///////////////////////////////////////////////////////////////////////////////
// main.hpp
//
//  Copyright 2005 Eric Niebler. Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <iostream>
#include <algorithm>
#include <boost/ref.hpp>
#include <boost/bind.hpp>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/p_square_quantile.hpp>
#include <boost/accumulators/statistics/tail_quantile.hpp>
#include <boost/accumulators/statistics/tail.hpp>

using namespace boost;
using namespace boost::accumulators;

void example4() {
  typedef accumulator_set<double, stats<tag::p_square_quantile> > accumulator_t;

 accumulator_t acc0(quantile_probability = 0.001);

 double sample;
 for (;;) {
   std::cin >> sample;
   if (std::cin.eof()) break;
   acc0(sample);
   std::cout << p_square_quantile(acc0) << std::endl;
 }
}

void example5() {
  std::size_t c =  10000; // cache size

  typedef accumulator_set<double, stats<tag::tail_quantile<left> > > accumulator_t_left;

  accumulator_t_left  acc0( tag::tail<left>::cache_size = c );

 double sample;
 unsigned i;
 for (unsigned i = 1;;++i) {
   std::cin >> sample;
   if (std::cin.eof()) break;
   acc0(sample);
   if (i % c == 0)
     std::cout << quantile(acc0, quantile_probability  = 0.001) << std::endl;
 }
}

void example6() {
  std::size_t c =  10000; // cache size

  typedef accumulator_set<double, stats<tag::tail<left> >, features< tag::max >  > accumulator_t_left;

  accumulator_t_left  acc0( tag::tail<left>::cache_size = c );

 double sample;
 unsigned i;
 for (unsigned i = 1;;++i) {
   std::cin >> sample;
   if (std::cin.eof()) break;
   acc0(sample);
   if (i % c == 0) {
     std::cout << extract_result< tag::max > ( acc0 ) << std::endl;
   }
 }
}


///////////////////////////////////////////////////////////////////////////////
// main
int main()
{
  example6();
}
