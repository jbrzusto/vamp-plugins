#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "FixedPulseDetector.h"

int
main (int argc, char *argv[])
{
  if (argc != 5) {
    std::cout << 
"Usage: testPulseDetector WIDTH MINDIFF MAXPROB NUMSAMPLES\n\
Find large magnitude or low probability edges in a data stream.\n\
WIDTH: pulse width, in samples\n\
MINDIFF: minimum difference between left and right window averages to be an edge\n\
MAXPROB: maximum probability of difference between left and right window\n\
        averages, to count as an edge\n\
NUMSAMPLES: number of random samples in [0,1] to generate; 0 means read from stdin\n\n\
A pulse is detected if either the MINDIFF or the MAXPROB criterion is satisfied.\n\
Pulses will be detected no closer than WIDTH samples apart\n";
    exit(1);
  }
    
  srand48(time(0));

  int width      = atoi(argv[1]);
  double mindiff = atof(argv[2]);
  double maxprob = atof(argv[3]);
  int m          = atoi(argv[4]);

  FixedPulseDetector < float > pd (width, mindiff, maxprob);
  unsigned long long count = 0;

  int i;
  for (i = 0; m == 0 || i < m; ++i) {
    float val;
    if (m == 0) {
      std::cin >> val;
      if (std::cin.eof())
        break;
    } else {
      val = (float) drand48();
      std::cout << val << std::endl;
    }
    ++count;
    if (pd(val)) {
      std::cout << count - pd.location() << ',' << pd.signal() << ',' << pd.bkgd() << ',';
      if (pd.big()) 
        std::cout << " diff " << pd.signal() - pd.bkgd() << " > " << mindiff << " ";
      if (pd.unlikely())
        std::cout << " prob(diff) <= " << maxprob << " quantile = " << pd.quantile();

      std::cout << std::endl;
    }
  }
}
