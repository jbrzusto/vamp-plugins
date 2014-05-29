#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "PulseDetector.h"

int
main (int argc, char *argv[])
{
  if (argc < 6) {
    std::cout << 
"Usage: testPulseDetector MINWIDTH MAXWIDTH MINDIFF MAXPROB NUMSAMPLES\n\
Find large magnitude or low probability edges in a data stream.\n\
MINWIDTH: minimum pulse width, in samples\n\
MAXWIDTH: maximum pulse width, in samples\n\
MINDIFF: minimum difference between left and right window averages to be an edge\n\
MAXPROB: maximum probability of difference between left and right window\n\
        averages, to count as an edge\n\
NUMSAMPLES: number of random samples in [0,1] to generate; 0 means read from stdin\n\n\
An edge is detected if either the MINDIFF or the MAXPROB criterion is satisfied.\n\
Pulses will be detected no closer than N samples apart, and an edge will\n\
only be detected when it is the best edge for N samples in either direction.\n";
    exit(1);
  }
    
  srand48(time(0));

  int min_width  = atoi(argv[1]);
  int max_width  = atoi(argv[2]);
  double mindiff = atof(argv[3]);
  double maxprob = atof(argv[4]);
  int m          = atoi(argv[5]);

  PulseDetector < float > pd (min_width, max_width, mindiff, maxprob);
  unsigned long long count = 0;

  int i;
  for (i = 0; m == 0 || i < m; ++i) {
    float val;
    if (m == 0) {
      std::cin >> val;
      if (std::cin.eof())
        break;
      ++count;
    } else {
      val = (float) drand48();
      std::cout << val << std::endl;
    }
    if (pd(val)) {
      std::cout << count - pd.location() << ',' << pd.width() << ',';
      auto a1 = pd.pulse_array1();
      auto a2 = pd.pulse_array2();
      if (a1.first) {
        std::cout << * a1.first;
      } else {
        std::cout << * a2.first;
      }
      std::cout << ',';
      if (a2.first) {
        std::cout << * (a2.first + (a2.second - 1));
      } else {
        std::cout << * (a1.first + (a1.second - 1));
      }
      std::cout << std::endl;
    }
  }
}
