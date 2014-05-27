#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "DynamicQuantile.h"

int
main (int argc, char *argv[])
{
  if (argc < 3) {
    std::cout << 
"Usage: testDynamicQuantile PROB N OVERLAP M\n\
Maintain an estimate of the PROB quantile of a distribution of M random numbers\n\
using staggered windows of size N, where consecutive windows have\n\
proportion OVERLAP samples in common.\n\
\n\
e.g. testDynamicQuantile 0.90 1000 0.8 100000\n";

    exit(1);
  }
    
  srand48(time(0));

  double prob    = atof (argv[1]);
  int n          = atoi(argv[2]);
  double overlap = atof (argv[3]);
  int m          = atoi(argv[4]);

  DynamicQuantile < float > dq (prob, n, overlap);

  int i;
  std::cout << "Estimator#,NewValue,Quantile\n";
  for (i = 0; i < m; ++i) {
    float val = (float) drand48();
    dq(val);
    std::cout << dq.current() << ',' << val << ',' << dq << std::endl;
  }
}
