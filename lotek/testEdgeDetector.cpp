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
NUMSAMPLES: number of random samples in [0,1] to generate\n\n\
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

  boost::circular_buffer < float > sampbuf (n);
  int i;
  for (i = 0; i < m; ++i) {
    float val = (float) drand48();
    if (sampbuf.full()) {
      std::cout << sampbuf[0] << std::endl;
    }
    sampbuf.push_back(val);
    if (ej(val)) {
      std::cout << "Edge: " << (ej.rising() ? "^ " : "v ") << ej.diff() 
                << ',' << ej.left() << ',' << ej.right() 
                << ',' << ej.quantile() << std::endl;
    }
  }
}
