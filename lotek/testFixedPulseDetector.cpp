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
  if (argc != 6) {
    std::cout << 
"Usage: testPulseDetector WIDTH BGWIDTH MINSNRDB MINZ NUMSAMPLES\n\
Find large magnitude or low probability edges in a data stream.\n\
WIDTH: pulse width, in samples\n\
BGWIDTH: background width, in samples (on each side of pulse)\n\
MINSNRDB: minimum pulse signal to noise ration, in dB\n\
MINZ: minimum z-score for difference between signal and noise\n\
NUMSAMPLES: number of random samples in [0,1] to generate; 0 means read from stdin\n\n\
A pulse is detected if either the MINSNR or the MINZ criterion is satisfied.\n\
Pulses will be detected no closer than WIDTH samples apart\n";
    exit(1);
  }
    
  srand48(time(0));

  int width      = atoi(argv[1]);
  int bgwidth    = atoi(argv[2]);
  double minsnrdB = atof(argv[3]);
  double minsnr  = exp10(minsnrdB / 10.);
  double minz    = atof(argv[4]);
  int m          = atoi(argv[5]);

  FixedPulseDetector < float > pd (width, bgwidth, minsnr, minz);
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
