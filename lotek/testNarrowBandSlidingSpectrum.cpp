#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "NarrowBandSlidingSpectrum.h"

int
main (int argc, char *argv[])
{
  if (argc != 5) {
    std::cout << 
"Usage: testNarrowBandSlidingSpectrum RATE SIZE FREQLOW FREQHIGH\n\
Generate sliding FFTs of an input stream.\n\
RATE: sampling rate (Hz) \n\
SIZE: window size for DFT\n\
MINFREQ: min freq. of interest, in kHz\n\
MAXFREQ: max freq. of interest, in kHz\n\
\n\
Input stream is interleaved I/Q floats.\n\
\n\
Output is lines with REAL,IMAGINARY,POWER\n\
for each bin for each frame.\n";

      exit(1);
  }
  
  int rate = atoi(argv[1]);
  int size = atoi(argv[2]);
  float freqlow = 1000.0 * atof(argv[3]);
  float freqhigh = 1000.0 * atof(argv[4]);
  float freqstep = rate / (float) size;
  int numbins = 1 + ceil((freqhigh - freqlow) / freqstep);
  
  NarrowBandSlidingSpectrum ss (size, freqlow, freqstep, numbins);

  for (;;) {
    float x_i, x_q;
    std::cin >> x_i >> x_q;
    if (std::cin.eof())
      break;
    std::complex < float > x (x_i, x_q);
    
    if (ss(x)) {
      for (int i = 0; i < numbins; ++i) {
        auto X = ss[i];
        std::cout << X.real() << ',' << X.imag() << ',' << std::abs(X) * std::abs(X) << std::endl;
      }
    }
  }
}
