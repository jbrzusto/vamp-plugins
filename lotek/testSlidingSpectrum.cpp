#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "SlidingSpectrum.h"

int
main (int argc, char *argv[])
{
  if (argc != 4 && argc != 6 && argc != 7) {
    std::cout << 
"Usage: testSlidingSpectrum SIZE PAD OVERLAP [MINBIN MAXBIN [WINTYPE]]\n\
Generate sliding FFTs of an input stream.\n\
SIZE: samples (not including zero-padding) per FFT\n\
PAD: number of zero samples to include in FFT\n\
OVERLAP: number of samples to overlap between consecutive FFTs\n\
MINBIN: min bin number to report (taking into account padding)\n\
MAXBIN: min bin number to report (taking into account padding)\n\
WINTYPE: defaults to 0 for Hamming; 1 means rectangular\n\
\n\
Input stream is interleaved I/Q floats.\n\
\n\
Output is lines with REAL,IMAGINARY,POWER\n\
for each bin for each frame.\n";

      exit(1);
  }
  
  int size = atoi(argv[1]);
  int pad = atoi(argv[2]);
  int overlap = atoi(argv[3]);
  int wintype = 0;
  int minbin = 0;
  int maxbin = pad + size - 1;
  
  if (argc >= 6) {
    minbin = atoi(argv[4]);
    maxbin = atoi(argv[5]);
    if (argc == 7)
      wintype = atoi(argv[6]);
  };
  
  SlidingSpectrum ss (size, pad, overlap, wintype == 0 ? SlidingSpectrum::WINDOW_HAMMING : SlidingSpectrum::WINDOW_RECTANGULAR);

  for (;;) {
    float x_i, x_q;
    std::cin >> x_i >> x_q;
    if (std::cin.eof())
      break;
    std::complex < float > x (x_i, x_q);
    
    if (ss(x)) {
      for (int i = minbin; i <= maxbin; ++i) {
        int ii = i;
        if (ii < 0)
          ii += pad + size;
        auto X = ss[ii];
        std::cout << X.real() << ',' << X.imag() << ',' << std::abs(X) * std::abs(X) << std::endl;
      }
    }
  }
}
