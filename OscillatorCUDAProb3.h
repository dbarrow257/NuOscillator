#ifndef __OSCILLATOR_CUDAPROB3_H__
#define __OSCILLATOR_CUDAPROB3_H__

//C++ Includes
#include <iostream>

//MaCh3 Code Includes
#include "OscillatorBase.h"

//Oscillation Implementation Includes
#include "CUDAProb3/constants.hpp"
#include "CUDAProb3/propagator.hpp"
#include "CUDAProb3/physics.hpp"

class OscillatorCUDAProb3 : virtual OscillatorBase {
 public:
  OscillatorCUDAProb3();
  void WhoAmI();
};

#endif
