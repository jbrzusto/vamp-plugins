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
  if (argc < 4 || argc > 5) {
    std::cout << 
"Usage: testSlidingSpectrum SIZE PAD OVERLAP [WINTYPE]\n\
Generate sliding FFTs of an input stream.\n\
SIZE: samples (not including zero-padding) per FFT\n\
PAD: number of zero samples to include in FFT\n\
OVERLAP: number of samples to overlap between consecutive FFTs\n\
WINTYPE: defaults to 0 for Hamming; 1 means rectangular\n\
\n\
Input stream is interleaved I/Q floats.\n";
      exit(1);
  }
  
  int size = atoi(argv[1]);
  int pad = atoi(argv[2]);
  int overlap = atoi(argv[3]);
  int wintype = 0;
  if (argc > 4)
    wintype = atoi(argv[4]);

  SlidingSpectrum ss (size, pad, overlap, wintype == 0 ? SlidingSpectrum::WINDOW_HAMMING : SlidingSpectrum::WINDOW_RECTANGULAR);

  for (;;) {
    float x_i, x_q;
    std::cin >> x_i >> x_q;
    if (std::cin.eof())
      break;
    std::complex < float > x (x_i, x_q);
    
    if (ss(x)) {
      for (int i = 0; i < size + pad; ++i) {
        auto X = ss[i];
        std::cout << X.real() << ',' << X.imag() << ',' << std::abs(X) * std::abs(X) << std::endl;
      }
    }
  }
}
