#include "OscillatorBase.h"

#ifdef UseCUDAProb3
#include "OscProbCalcer_CUDAProb3.h"
#endif

#ifdef UseProb3ppLinear
#include "OscProbCalcer_Prob3ppLinear.h"
#endif

#ifdef UseProbGPULinear
#include "OscProbCalcer_ProbGPULinear.h"
#endif

#include <iostream>

OscillatorBase::OscillatorBase() {
  InitialiseOscProbCalcer();
}

void OscillatorBase::InitialiseOscProbCalcer() {

#ifdef UseCUDAProb3
  OscProbCalcerCUDAProb3* CUDAProb3 = new OscProbCalcerCUDAProb3();
  OPCalcer = (OscProbCalcerBase*)CUDAProb3;
#endif

#ifdef UseProb3ppLinear
  OscProbCalcerProb3ppLinear* Prob3ppLinear = new OscProbCalcerProb3ppLinear();
  OPCalcer = (OscProbCalcerBase*)Prob3ppLinear;
#endif

#ifdef UseProbGPULinear
  OscProbCalcerProbGPULinear* ProbGPULinear = new OscProbCalcerProbGPULinear();
  OPCalcer = (OscProbCalcerBase*)ProbGPULinear;
#endif
}

void OscillatorBase::ImplementationName() {
  std::cout << OPCalcer->ReturnImplementationName() << std::endl;
}
