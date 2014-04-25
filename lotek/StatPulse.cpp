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
#include <boost/random.hpp>
#include <boost/test/unit_test.hpp>

using namespace boost;
using namespace boost::accumulators;

// Helper that uses BOOST_FOREACH to display a range of doubles
template<typename Range>
void output_range(Range const &rng)
{
    bool first = true;
    BOOST_FOREACH(double d, rng)
    {
        if(!first) std::cout << ", ";
        std::cout << d;
        first = false;
    }
    std::cout << '\n';
}

///////////////////////////////////////////////////////////////////////////////
// example1
//
//  Calculate some useful stats using accumulator_set<> and std::for_each()
//
void example1()
{
    accumulator_set<
        double
      , stats<tag::min, tag::mean(immediate), tag::sum, tag::moment<2> >
    > acc;

    boost::array<double, 4> data = {0., 1., -1., 3.14159};

    // std::for_each pushes each sample into the accumulator one at a
    // time, and returns a copy of the accumulator.
    acc = std::for_each(data.begin(), data.end(), acc);

    // The following would be equivalent, and could be more efficient
    // because it doesn't pass and return the entire accumulator set
    // by value.
    //std::for_each(data.begin(), data.end(), bind<void>(ref(acc), _1));

    std::cout << "  min""(acc)        = " << (min)(acc) << std::endl; // Extra quotes are to prevent complaints from Boost inspect tool
    std::cout << "  mean(acc)       = " << mean(acc) << std::endl;

    // since mean depends on count and sum, we can get their results, too.
    std::cout << "  count(acc)      = " << count(acc) << std::endl;
    std::cout << "  sum(acc)        = " << sum(acc) << std::endl;
    std::cout << "  moment<2>(acc)  = " << accumulators::moment<2>(acc) << std::endl;
}

///////////////////////////////////////////////////////////////////////////////
// example2
//
//  Calculate some tail statistics. This demonstrates how to specify
//  constructor and accumulator parameters. Note that the tail statistics
//  return multiple values, which are returned in an iterator_range.
//
//  It pushes data in and displays the intermediate results to demonstrate
//  how the tail statistics are updated.
void example2()
{
    // An accumulator which tracks the right tail (largest N items) and
    // some data that are covariate with them. N == 4.
    accumulator_set<
        double
      , stats<tag::tail_variate<double, tag::covariate1, right> >
    > acc(tag::tail<right>::cache_size = 4);

    acc(2.1, covariate1 = .21);
    acc(1.1, covariate1 = .11);
    acc(2.1, covariate1 = .21);
    acc(1.1, covariate1 = .11);

    std::cout << "  tail            = "; output_range(tail(acc));
    std::cout << "  tail_variate    = "; output_range(tail_variate(acc));
    std::cout << std::endl;

    acc(21.1, covariate1 = 2.11);
    acc(11.1, covariate1 = 1.11);
    acc(21.1, covariate1 = 2.11);
    acc(11.1, covariate1 = 1.11);

    std::cout << "  tail            = "; output_range(tail(acc));
    std::cout << "  tail_variate    = "; output_range(tail_variate(acc));
    std::cout << std::endl;

    acc(42.1, covariate1 = 4.21);
    acc(41.1, covariate1 = 4.11);
    acc(42.1, covariate1 = 4.21);
    acc(41.1, covariate1 = 4.11);

    std::cout << "  tail            = "; output_range(tail(acc));
    std::cout << "  tail_variate    = "; output_range(tail_variate(acc));
    std::cout << std::endl;

    acc(32.1, covariate1 = 3.21);
    acc(31.1, covariate1 = 3.11);
    acc(32.1, covariate1 = 3.21);
    acc(31.1, covariate1 = 3.11);

    std::cout << "  tail            = "; output_range(tail(acc));
    std::cout << "  tail_variate    = "; output_range(tail_variate(acc));
}

///////////////////////////////////////////////////////////////////////////////
// example3
//
//  Demonstrate how to calculate weighted statistics. This example demonstrates
//  both a simple weighted statistical calculation, and a more complicated
//  calculation where the weight statistics are calculated and stored in an
//  external weight accumulator.
void example3()
{
    // weight == double
    double w = 1.;

    // Simple weighted calculation
    {
        // stats that depend on the weight are made external
        accumulator_set<double, stats<tag::mean>, double> acc;

        acc(0., weight = w);
        acc(1., weight = w);
        acc(-1., weight = w);
        acc(3.14159, weight = w);

        std::cout << "  mean(acc)       = " << mean(acc) << std::endl;
    }

    // Weighted calculation with an external weight accumulator
    {
        // stats that depend on the weight are made external
        accumulator_set<double, stats<tag::mean>, external<double> > acc;

        // Here's an external weight accumulator
        accumulator_set<void, stats<tag::sum_of_weights>, double> weight_acc;

        weight_acc(weight = w); acc(0., weight = w);
        weight_acc(weight = w); acc(1., weight = w);
        weight_acc(weight = w); acc(-1., weight = w);
        weight_acc(weight = w); acc(3.14159, weight = w);

        std::cout << "  mean(acc)       = " << mean(acc, weights = weight_acc) << std::endl;
    }
}

void example4() {
typedef accumulator_set<double, stats<tag::p_square_quantile> > accumulator_t;

// tolerance in %
double epsilon = 1;

// a random number generator
boost::lagged_fibonacci607 rng;

accumulator_t acc0(quantile_probability = 0.001);
accumulator_t acc1(quantile_probability = 0.01 );
accumulator_t acc2(quantile_probability = 0.1  );
accumulator_t acc3(quantile_probability = 0.25 );
accumulator_t acc4(quantile_probability = 0.5  );
accumulator_t acc5(quantile_probability = 0.75 );
accumulator_t acc6(quantile_probability = 0.9  );
accumulator_t acc7(quantile_probability = 0.99 );
accumulator_t acc8(quantile_probability = 0.999);

for (int i=0; i<100000; ++i)
{
    double sample = rng();
    acc0(sample);
    acc1(sample);
    acc2(sample);
    acc3(sample);
    acc4(sample);
    acc5(sample);
    acc6(sample);
    acc7(sample);
    acc8(sample);
}

BOOST_CHECK_CLOSE( p_square_quantile(acc0), 0.001, 15*epsilon );
BOOST_CHECK_CLOSE( p_square_quantile(acc1), 0.01 , 5*epsilon );
BOOST_CHECK_CLOSE( p_square_quantile(acc2), 0.1  , epsilon );
BOOST_CHECK_CLOSE( p_square_quantile(acc3), 0.25 , epsilon );
BOOST_CHECK_CLOSE( p_square_quantile(acc4), 0.5  , epsilon );
BOOST_CHECK_CLOSE( p_square_quantile(acc5), 0.75 , epsilon );
BOOST_CHECK_CLOSE( p_square_quantile(acc6), 0.9  , epsilon );
BOOST_CHECK_CLOSE( p_square_quantile(acc7), 0.99 , epsilon );
BOOST_CHECK_CLOSE( p_square_quantile(acc8), 0.999, epsilon );
}
///////////////////////////////////////////////////////////////////////////////
// main
int main()
{
    std::cout << "Example 1:\n";
    example1();

    std::cout << "\nExample 2:\n";
    example2();

    std::cout << "\nExample 3:\n";
    example3();

    std::cout << "\nExample 4:\n";
    return 0;
}
