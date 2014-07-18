#include <complex>

typedef std::complex < float > C;

C cumsum (C * in, int n, C sum) {
  for (int i = 0; i < n; ++i)
    sum += in[i];
  return sum;
}


float cumsum (float * in, int n) {
  float sum = 0;
  for (int i = 0; i < n; ++i)
    sum += in[i];
  return sum;
}


