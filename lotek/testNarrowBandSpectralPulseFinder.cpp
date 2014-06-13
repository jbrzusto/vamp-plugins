#include <iostream>
#include <fstream>
#include <cstdio>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "NarrowBandSpectralPulseFinder.h"

int
main (int argc, char *argv[])
{
  if (argc < 8 || argc > 9) {
    std::cout << 
"Usage: testSpectralPulseFinder RATE WIDTH BINLOW BINHIGH MINSNRDB MINZ MAXNOISEZ [FILE]\n\
RATE: sample rate (samples per second)\n\
Find fixed-width pulses in the spectrum of an input stream of complex numbers.\n\
WIDTH: width of pulses to find, in samples\n\
FREQLWO: lowest frequency of interest, in kHz\n\
FREQHIGH: highest frequency of interest, in kHz\n\
MINSNRDB: min signal to noise ratio of a pulse, in dB\n\
MINZ: min z-score for a pulse's sig - noise value.\n\
MAXNOISEZ: max noise (dB) at which the MINZ criterion is used.\n\
A pulse is accepted if it satisfies either the MINSNRDB or the MINZ criterion.\n\
Data are read from FILE, if specified, or stdin as R1 I1 R2 I2 R3 I3 ...\n\
(i.e. interleaved real and imaginary parts of each sample).\n";
      exit(1);
  }
  
  int rate          = atoi(argv[1]);
  int width         = atoi(argv[2]);
  double minfreq    = atof(argv[3]);
  double maxfreq    = atof(argv[4]);
  double minsnrdb   = atof(argv[5]);
  double minz       = atof(argv[6]);
  double maxnoisez  = atof(argv[7]);

  double freqstep = rate / (float) width;
  int minbin        = floor(minfreq * 1000.0 / freqstep);
  int maxbin        = ceil(maxfreq * 1000.0 / freqstep);

  std::istream * in = (argc > 8) ? new std::ifstream(argv[8]) : & std::cin;

  NarrowBandSpectralPulseFinder nspf (width, minbin, maxbin, minsnrdb, minz, maxnoisez);

  unsigned long long count = 0;

  // scale for the fact samples are -32767..32767 instead of -1..1
  float offset = dB(1.0 / (32767 * 32767));
  for (;;) {
    float x_i, x_q;
    *in >> x_i >> x_q;
    if (in->eof() || in->fail())
      break;
    std::complex < float > x (x_i, x_q);

    if (nspf(x)) {
      std::cout << (count - nspf.location()) / (float) rate << ',' << nspf.freq() * rate / (float) width << ',' << dB(nspf.signal()) + offset << ',' << dB(nspf.noise()) + offset << ',' << dB(nspf.SNR()) << /* ',' << nspf.Z() << */ std::endl;
    }
    ++count;
  }
  if (argc > 8)
    delete in;
}
