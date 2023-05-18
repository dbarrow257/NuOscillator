#ifndef __OSCILLATOR_CUDAPROB3_H__
#define __OSCILLATOR_CUDAPROB3_H__

#include "OscillatorBase.h"

#include <memory>

template<typename T> class Propagator;

class OscillatorCUDAProb3 : public OscillatorBase {
 public:
  OscillatorCUDAProb3();
  void WhoAmI();

 private:
  void InitialisePropagator();
  std::unique_ptr< Propagator< FLOAT_T > > propagator;
};

#endif
