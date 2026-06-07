#ifndef WAVELET_HPP
#define WAVELET_HPP

#include "utils.h"

#include "../src/par.h"

class Wavelet 
{
private:
  const config_t &c;
public:
  Wavelet(const config_t& config) 
    : c(config) 
  {
    wavelet = allocf(c.ntw);
    wavelet_dt = allocf(c.ntw);
  }

  float* wavelet;
  float* wavelet_dt;

  void get();
  void second_derivative();
};

#endif
