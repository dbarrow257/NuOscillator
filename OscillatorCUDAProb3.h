#ifndef __OSCILLATOR_CUDAPROB3_H__
#define __OSCILLATOR_CUDAPROB3_H__

#include "OscillatorBase.h"

#include <memory>

#include "CUDAProb3/constants.hpp"
#include "CUDAProb3/propagator.hpp"
#include "CUDAProb3/physics.hpp"

#include "CUDAProb3/cpupropagator.hpp"

using namespace cudaprob3;

class OscillatorCUDAProb3 : virtual OscillatorBase {
 public:
  OscillatorCUDAProb3();
  void WhoAmI();

 private:
  void InitialisePropagator();
  std::unique_ptr< Propagator< FLOAT_T > > propagator;
};

#endif
