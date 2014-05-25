#include "FreqEstimator.h"
#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <cmath>
#include <vector>
#include <stdlib.h>

#define DEFAULT_TEST_FREQ 19.12345
#define DEFAULT_NUM_SAMPLES 120
#define MULT1 5
#define MULT2 10
#define MULT3 20

int
main (int argc, char *argv[]) {

  double test_freq = DEFAULT_TEST_FREQ;
  int num_samples = DEFAULT_NUM_SAMPLES;
  if (argc > 1)
    test_freq = strtod (argv[1], 0);

  if (argc > 2)
    num_samples = strtol (argv[2], 0, 0);

  int16_t inbuff[2 * num_samples];

  float true_power = 0.0;
  for (int i = 0; i < num_samples; ++i) {
    inbuff[2*i]   = round(32767 * cos(2 * M_PI * test_freq * i / num_samples));
    inbuff[2*i+1] = round(32767 * sin(2 * M_PI * test_freq * i / num_samples));
    true_power += inbuff[2*i]*inbuff[2*i] + inbuff[2*i+1]*inbuff[2*i+1];
  }
  true_power = 10 * log10f(true_power / 2 / 32767 / 32767);

  FreqEstimator fe0(num_samples), fe1(num_samples * MULT1), fe2(num_samples * MULT2), fe3(num_samples * MULT3), fe4(4096), fe5(2048), fe6(1024);

  float pwr;

  std::cout  << "True Freq: " << test_freq << ", power: " << true_power << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe0.get(inbuff, num_samples, &pwr) << " without zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe1.get(inbuff, num_samples, &pwr) << " with " << MULT1 << " times zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe2.get(inbuff, num_samples, &pwr) << " with " << MULT2 << " times zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe3.get(inbuff, num_samples, &pwr) << " with " << MULT3 << " times zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe6.get(inbuff, num_samples, &pwr) << " with padding to 1024 samples";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe5.get(inbuff, num_samples, &pwr) << " with padding to 2048 samples";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe4.get(inbuff, num_samples, &pwr) << " with padding to 4096 samples";
  std::cout << ", power: " << pwr << "\n";

  fe0.generateWindowingCoefficients();
  fe1.generateWindowingCoefficients();
  fe2.generateWindowingCoefficients();
  fe3.generateWindowingCoefficients();
  fe4.generateWindowingCoefficients();
  fe5.generateWindowingCoefficients();
  fe6.generateWindowingCoefficients();

  std::cout << " After Windowing" << std::endl;
  std::cout  << "True Freq: " << test_freq << ", power: " << true_power << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe0.get(inbuff, num_samples, &pwr) << " without zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe1.get(inbuff, num_samples, &pwr) << " with " << MULT1 << " times zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe2.get(inbuff, num_samples, &pwr) << " with " << MULT2 << " times zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe3.get(inbuff, num_samples, &pwr) << " with " << MULT3 << " times zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe6.get(inbuff, num_samples, &pwr) << " with padding to 1024 samples";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe5.get(inbuff, num_samples, &pwr) << " with padding to 2048 samples";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << num_samples * fe4.get(inbuff, num_samples, &pwr) << " with padding to 4096 samples";
  std::cout << ", power: " << pwr << "\n";
}
