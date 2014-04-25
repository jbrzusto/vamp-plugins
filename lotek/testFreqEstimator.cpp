#include "FreqEstimator.h"
#include <iostream>
#include <cstdio>
#include <stdint.h>
#include <cmath>

#define NUM_SAMPLES 240
#define TEST_FREQ 119.765
#define MULT1 5
#define MULT2 10
#define MULT3 20

main (int argc, char *argv[]) {

  int16_t inbuff[2 * NUM_SAMPLES];
  
  float true_power = 0.0;
  for (int i = 0; i < NUM_SAMPLES; ++i) {
    inbuff[2*i]   = round(32767 * cos(2 * M_PI * TEST_FREQ * i / NUM_SAMPLES));
    inbuff[2*i+1] = round(32767 * sin(2 * M_PI * TEST_FREQ * i / NUM_SAMPLES));
    true_power += inbuff[2*i]*inbuff[2*i] + inbuff[2*i+1]*inbuff[2*i+1];
  }

  FreqEstimator fe0(NUM_SAMPLES), fe1(NUM_SAMPLES * MULT1), fe2(NUM_SAMPLES * MULT2), fe3(NUM_SAMPLES * MULT3), fe4(4096), fe5(2048), fe6(1024);

  float pwr;

  std::cout  << "True Freq: " << TEST_FREQ << ", power: " << true_power << "\n";
  std::cout  << "Est. Freq: " << NUM_SAMPLES * fe0.get(inbuff, NUM_SAMPLES, &pwr) << " without zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << NUM_SAMPLES * fe1.get(inbuff, NUM_SAMPLES, &pwr) << " with " << MULT1 << " times zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << NUM_SAMPLES * fe2.get(inbuff, NUM_SAMPLES, &pwr) << " with " << MULT2 << " times zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << NUM_SAMPLES * fe3.get(inbuff, NUM_SAMPLES, &pwr) << " with " << MULT3 << " times zero padding";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << NUM_SAMPLES * fe6.get(inbuff, NUM_SAMPLES, &pwr) << " with padding to 1024 samples";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << NUM_SAMPLES * fe5.get(inbuff, NUM_SAMPLES, &pwr) << " with padding to 2048 samples";
  std::cout << ", power: " << pwr << "\n";
  std::cout  << "Est. Freq: " << NUM_SAMPLES * fe4.get(inbuff, NUM_SAMPLES, &pwr) << " with padding to 4096 samples";
  std::cout << ", power: " << pwr << "\n";
}
