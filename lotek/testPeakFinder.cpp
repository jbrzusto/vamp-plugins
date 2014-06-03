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
"Usage: testPeakFinder WINSIZE NUMSAMPLES\n\
look for local strict maxima (over a window of size 2*WINSIZE + 1)\n\
in a stream of NUMSAMPLES random floating point numbers in [0, 1].\n\
Specify NUMSAMPLES as 0 to read samples from stdin.\n\
";
    exit(1);
  }
    
  srand48(time(0));

  int m = atoi(argv[1]);
  PeakFinder < float > pf(m, true, true);

  int n = atoi(argv[2]);

  int i;
  for (i = 0; n == 0 || i < n; ++i) {
    float val;
    if (n == 0) {
      std::cin >> val;
    } else {
      val = (float) drand48();
      std::cout << val << std::endl;
    }
    if (pf(val))
      std::cout << "Peak (back " << m << ") " <<  (pf ? "^ " : "v ") << " = " << (double) pf << std::endl;
  }
}
