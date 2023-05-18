#include "OscillatorCUDAProb3.h"

#include "CUDAProb3/constants.hpp"
#include "CUDAProb3/propagator.hpp"
#include "CUDAProb3/physics.hpp"

#include "CUDAProb3/cpupropagator.hpp"

using namespace cudaprob3;

OscillatorCUDAProb3::OscillatorCUDAProb3() : OscillatorBase(){}

void OscillatorCUDAProb3::WhoAmI() {
  std::cout << "OscillatorCUDAProb3" << std::endl;
}

void OscillatorCUDAProb3::InitialisePropagator() {
  int nThreads = 1;
  propagator = std::unique_ptr< Propagator< FLOAT_T > > ( new CpuPropagator<FLOAT_T>(nCosZPoints, nEnergyPoints, nThreads)); // MultiThread CPU propagator
}
