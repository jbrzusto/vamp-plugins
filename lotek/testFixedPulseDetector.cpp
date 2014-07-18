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
  if (argc != 8) {
    std::cout << 
"Usage: testPulseDetector WIDTH BGWIDTH AVOID MINSNRDB MINZ MAXNOISEZDB NUMSAMPLES\n\
Find large magnitude or low probability edges in a data stream.\n\
WIDTH: pulse width, in samples\n\
BGWIDTH: background width, in samples (on each side of pulse)\n\
AVOID: avoidance zone, in samples (on each side of pulse between pulse and BKGD)\n\
to remove leakage between pulse and correlated background samples\n\
MINSNRDB: minimum pulse signal to noise ration, in dB\n\
MINZ: minimum z-score for difference between signal and noise\n\
MAXNOISEZDB: maximum noise at which z-score is used to accept pulse\n\
NUMSAMPLES: number of random samples in [0,1] to generate; 0 means read from stdin\n\n\
A pulse is detected if either the MINSNR or the MINZ criterion is satisfied.\n\
Pulses will be detected no closer than WIDTH samples apart\n";
    exit(1);
  }
    
  srand48(time(0));

  int width      = atoi(argv[1]);
  int bgwidth    = atoi(argv[2]);
  int avoid      = atoi(argv[3]);
  double minsnrdB = atof(argv[4]);
  double minsnr  = undB(minsnrdB);
  double minz    = atof(argv[5]);
  double maxnoisez = undB(atof(argv[6]));
  int m          = atoi(argv[7]);

  FixedPulseDetector < float > pd (width, bgwidth, avoid, minsnr, minz, maxnoisez);
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
        std::cout << " SNR " << 10 * log10(pd.SNR()) << " > " << minsnrdB << " ";
      if (pd.unlikely())
        std::cout << " Z = " << pd.Z();
      std::cout << std::endl;
    }
  }
}
