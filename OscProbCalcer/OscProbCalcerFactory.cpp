#include "OscProbCalcerFactory.h"

#if UseCUDAProb3==1
#include "OscProbCalcer/OscProbCalcer_CUDAProb3.h"
#endif

#if UseCUDAProb3Linear==1
#include "OscProbCalcer/OscProbCalcer_CUDAProb3Linear.h"
#endif

#if UseProb3ppLinear==1
#include "OscProbCalcer/OscProbCalcer_Prob3ppLinear.h"
#endif

#if UseProbGPULinear==1
#include "OscProbCalcer/OscProbCalcer_ProbGPULinear.h"
#endif

#if UseNuFASTLinear==1
#include "OscProbCalcer/OscProbCalcer_NuFASTLinear.h"
#endif

#include <iostream>

OscProbCalcerFactory::OscProbCalcerFactory() {
}

OscProbCalcerFactory::~OscProbCalcerFactory() {
}

OscProbCalcerBase* OscProbCalcerFactory::CreateOscProbCalcer(std::string OscProbCalcerImplementationToCreate, std::string OscProbCalcerConfigName, int Instance, int Verbose) {
  OscProbCalcerBase* Calcer = NULL;
  
  if (OscProbCalcerImplementationToCreate == "CUDAProb3") {
#if UseCUDAProb3==1
    OscProbCalcerCUDAProb3* CUDAProb3 = new OscProbCalcerCUDAProb3(OscProbCalcerConfigName,Instance);
    Calcer = (OscProbCalcerBase*)CUDAProb3;
    if (Verbose >= INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }

  else if (OscProbCalcerImplementationToCreate == "CUDAProb3Linear") {
#if UseCUDAProb3Linear==1
    OscProbCalcerCUDAProb3Linear* CUDAProb3Linear = new OscProbCalcerCUDAProb3Linear(OscProbCalcerConfigName,Instance);
    Calcer = (OscProbCalcerBase*)CUDAProb3Linear;
    if (Verbose >= INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else if (OscProbCalcerImplementationToCreate == "Prob3ppLinear") {
#if UseProb3ppLinear==1
    OscProbCalcerProb3ppLinear* Prob3ppLinear = new OscProbCalcerProb3ppLinear(OscProbCalcerConfigName,Instance);
    Calcer = (OscProbCalcerBase*)Prob3ppLinear;
    if (Verbose >= INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else if (OscProbCalcerImplementationToCreate == "ProbGPULinear") {
#if UseProbGPULinear==1
    OscProbCalcerProbGPULinear* ProbGPULinear = new OscProbCalcerProbGPULinear(OscProbCalcerConfigName,Instance);
    Calcer = (OscProbCalcerBase*)ProbGPULinear;
    if (Verbose >= INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }

  else if (OscProbCalcerImplementationToCreate == "NuFASTLinear") {
#if UseNuFASTLinear==1
    OscProbCalcerNuFASTLinear* NuFASTLinear = new OscProbCalcerNuFASTLinear(OscProbCalcerConfigName,Instance);
    Calcer = (OscProbCalcerBase*)NuFASTLinear;
    if (Verbose >= INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else {
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but this is not implemented within " << __FILE__ << std::endl;
    std::cerr << "Please correct the mistake or implement the new OscProbCalcer" << std::endl;
    throw;
  }

  return Calcer;
}
