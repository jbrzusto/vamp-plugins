#include "PeakFinder.h"
#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>

int
main (int argc, char *argv[])
{
  if (argc < 3) {
    std::cout << 
"Usage: testPeakFinder winsize numsamples\n\
look for local strict maxima (over a window of size 2*winsize + 1)\
in a stream of numsamples random floating point numbers in [0, 1].\n";
    exit(1);
  }
    
  srand48(time(0));

  int m = atoi(argv[1]);
  PeakFinder < float > pf(m);

  int n = atoi(argv[2]);

  int i;
  for (i = 0; i < n; ++i) {
    float val = (float) drand48();
    bool peak = pf.process(val);
    if (i < m)
      std::cout << "  " << pf[i] << std::endl;
    else if (pf.full())
      std::cout << (peak ? "^ " : "  ") << pf[m] << std::endl;
  }
  for (i = 1; i <= m; ++i)
    std::cout << "  " << pf[m + i] << std::endl;

}
