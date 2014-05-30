#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include "Accordion.h"

Accordion < float > xx (20, 5);

int
main (int argc, char *argv[])
{
  if (argc != 4) {
    std::cout << 
"Usage: testAccordion WIDTH OVERLAP NUMSAMPLES\n\
Replicate an input sequence so that the output is given as windows of\n\
WIDTH samples at a time, with OVERLAP samples shared between\n\
consecutive windows.  The input sequence is read from stdin, \n\
if NUMSAMPLES == 0, or is 1 ... NUMSAMPLES otherwise.\n";

    exit(1);
  }
    
  size_t width   = atoi(argv[1]);
  size_t overlap = atof(argv[2]);
  int m          = atoi(argv[3]);

  Accordion < float > ac (width, overlap);

  int i;
  for (i = 1; m == 0 || i <= m; ++i) {
    float val;
    if (m == 0) {
      std::cin >> val;
      if (std::cin.eof())
        break;
    } else {
      val = i;
    }
    if (ac(val)) {
      auto a1 = ac.array_one();
      for (size_t j = 0; j < a1.second; ++j)
        std::cout << a1.first[j] << std::endl;

      auto a2 = ac.array_two();
      for (size_t j = 0; j < a2.second; ++j)
        std::cout << a2.first[j] << std::endl;
    }
  }
  // output any items which have not yet been seen

  auto a1 = ac.partial_array_one();
  for (size_t j = 0; j < a1.second; ++j)
    std::cout << a1.first[j] << std::endl;

  auto a2 = ac.partial_array_two();
  for (size_t j = 0; j < a2.second; ++j)
    std::cout << a2.first[j] << std::endl;
}
