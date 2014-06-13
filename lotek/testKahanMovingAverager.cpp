#include "KahanMovingAverager.h"
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
  if (argc != 3) {
    std::cout << 
"Usage: testKahanMovingAverager WINSIZE NUMSAMPLES\n\
calculate a running moving average using a window of size WINSIZE\n\
in a stream of NUMSAMPLES random floating point numbers in [0, 1].\n\
Specify NUMSAMPLES as 0 to read samples from stdin.\n\
";
    exit(1);
  }
    
  srand48(time(0));

  int m = atoi(argv[1]);
  KahanMovingAverager < float, double > ma(m);

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
    if (ma(val))
      std::cout << "av: " << ma << std::endl;
  }
}
