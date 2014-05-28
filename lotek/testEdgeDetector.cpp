#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "EdgeDetector.h"

int
main (int argc, char *argv[])
{
  if (argc < 5) {
    std::cout << 
"Usage: testEdgeDetector N MINDIFF MAXPROB NUMSAMPLES\n\
Find large magnitude or low probability edges in a data stream.\n\
N: size of window.  An edge consists of a significant change in mean\n\
  value between N samples ('left window') and the subsequent N samples ('right window').\n\
MINDIFF: minimum difference between left and right window averages to be an edge\n\
MAXPROB: maximum probability of difference between left and right window\n\
        averages, to count as an edge\n\
NUMSAMPLES: number of random samples in [0,1] to generate; 0 means read from stdin\n\n\
An edge is detected if either the MINDIFF or the MAXPROB criterion is satisfied.\n\
Edges will be detected no closer than N samples apart, and an edge will\n\
only be detected when it is the best edge for N samples in either direction.\n";
    exit(1);
  }
    
  srand48(time(0));

  int n          = atoi(argv[1]);
  double mindiff = atof(argv[2]);
  double maxprob = atof(argv[3]);
  int m          = atoi(argv[4]);

  EdgeDetector < float > ej (n, mindiff, maxprob);

  int i;
  for (i = 0; m == 0 || i < m; ++i) {
    float val;
    if (m == 0) {
      std::cin >> val;
    } else {
      val = (float) drand48();
      std::cout << val << std::endl;
    }
    if (ej(val)) {
      std::cout << "Edge (back " << 2 * n << "): " << (ej.rising() ? "^ " : "v ") << ej.diff() 
                << ',' << ej.left() << ',' << ej.right() 
                << ',' << ej.quantile() << std::endl;
    }
  }
}
