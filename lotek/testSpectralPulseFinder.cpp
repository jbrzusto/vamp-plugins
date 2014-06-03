#include <iostream>
#include <fstream>
#include <cstdio>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "SpectralPulseFinder.h"

int
main (int argc, char *argv[])
{
  if (argc < 10 || argc > 11) {
    std::cout << 
"Usage: testSpectralPulseFinder WIDTH BKGDWIDTH WINSIZE PAD OVERLAP MINBIN MAXBIN MINSNRDB MINZ [FILE]\n\
Find fixed-width pulses in the spectrum of an input stream of complex numbers.\n\
WIDTH: width of pulses to find, in samples\n\
BKGWIDTH: width of background window on each side, in samples\n\
WINSIZE: number of data samples in each FFT \n\
PAD: number of zero (padding) samples for each FFT\n\
OVERLAP: number of samples of overlap between consecutive FFTs\n\
MINBIN: minimum bin number of interest in finding pulses (-WINSIZE / 2 ... WINSIZE / 2)\n\
MAXBIN: maximum bin number of interest in finding pulses (-WINSIZE / 2 ... WINSIZE / 2)\n\
MINSNRDB: min signal to noise ratio of a pulse, in dB\n\
MINZ: min z-score for a pulse's sig - noise value.\n\
A pulse is accepted if it satisfies either the MINSNRDB or the MINZ criterion.\n\
Data are read from FILE, if specified, or stdin as R1 I1 R2 I2 R3 I3 ...\n\
(i.e. interleaved real and imaginary parts of each sample).\n";
      exit(1);
  }
  
  int width = atoi(argv[1]);
  int bkgd = atoi(argv[2]);
  int winsize = atoi(argv[3]);
  int pad = atoi(argv[4]);
  int overlap = atoi(argv[5]);
  int minbin = atoi(argv[6]);
  int maxbin = atoi(argv[7]);
  double minsnrdb = atof(argv[8]);
  double minz = atof(argv[9]);

  std::istream * in = (argc > 10) ? new std::ifstream(argv[10]) : & std::cin;

  SpectralPulseFinder spf (width, bkgd, winsize, pad, overlap, minbin, maxbin, minsnrdb, minz);

  unsigned long long count = 0;

  for (;;) {
    float x_i, x_q;
    *in >> x_i >> x_q;
    if (in->eof() || in->fail())
      break;
    std::complex < float > x (x_i, x_q);

    if (spf(x)) {
      for (auto biniter = spf.beginbin(); biniter != spf.endbin(); ++biniter) {
        int i = *biniter;
        std::cout << count << ',' << i << ',' << spf.signal(i) << ',' << spf.noise(i) << ',' << spf.Z(i) << std::endl;
      }
    }
    ++count;
  }
  if (argc > 9)
    delete in;
}
