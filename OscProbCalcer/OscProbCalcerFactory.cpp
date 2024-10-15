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

#if UseOscProb==1
#include "OscProbCalcer/OscProbCalcer_OscProb.h"
#endif

#include <iostream>

OscProbCalcerFactory::OscProbCalcerFactory() {
}

OscProbCalcerFactory::~OscProbCalcerFactory() {
}

OscProbCalcerBase* OscProbCalcerFactory::CreateOscProbCalcer(std::string OscProbCalcerConfigName) {
  YAML::Node Config = YAML::LoadFile(OscProbCalcerConfigName);
  return CreateOscProbCalcer(Config);
}

OscProbCalcerBase* OscProbCalcerFactory::CreateOscProbCalcer(YAML::Node OscProbCalcerConfig) {
  OscProbCalcerBase* Calcer = NULL;

  std::string OscProbCalcerImplementationToCreate = OscProbCalcerConfig["OscProbCalcerSetup"]["ImplementationName"].as<std::string>();
  std::string Verbosity = OscProbCalcerConfig["General"]["Verbosity"].as<std::string>();
  int Verbose = Verbosity_StrToInt(Verbosity);

  if (OscProbCalcerImplementationToCreate == "CUDAProb3") {
#if UseCUDAProb3==1
    OscProbCalcerCUDAProb3* CUDAProb3 = new OscProbCalcerCUDAProb3(OscProbCalcerConfig);
    Calcer = (OscProbCalcerBase*)CUDAProb3;
    if (Verbose >= NuOscillator::INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscProbCalcerFactory object" << std::endl;}
#else
    std::cerr << "OscProbCalcerFactory was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }

  else if (OscProbCalcerImplementationToCreate == "CUDAProb3Linear") {
#if UseCUDAProb3Linear==1
    OscProbCalcerCUDAProb3Linear* CUDAProb3Linear = new OscProbCalcerCUDAProb3Linear(OscProbCalcerConfig);
    Calcer = (OscProbCalcerBase*)CUDAProb3Linear;
    if (Verbose >= NuOscillator::INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscProbCalcerFactory object" << std::endl;}
#else
    std::cerr << "OscProbCalcerFactory was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else if (OscProbCalcerImplementationToCreate == "Prob3ppLinear") {
#if UseProb3ppLinear==1
    OscProbCalcerProb3ppLinear* Prob3ppLinear = new OscProbCalcerProb3ppLinear(OscProbCalcerConfig);
    Calcer = (OscProbCalcerBase*)Prob3ppLinear;
    if (Verbose >= NuOscillator::INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscProbCalcerFactory object" << std::endl;}
#else
    std::cerr << "OscProbCalcerFactory was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else if (OscProbCalcerImplementationToCreate == "ProbGPULinear") {
#if UseProbGPULinear==1
    OscProbCalcerProbGPULinear* ProbGPULinear = new OscProbCalcerProbGPULinear(OscProbCalcerConfig);
    Calcer = (OscProbCalcerBase*)ProbGPULinear;
    if (Verbose >= NuOscillator::INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscProbCalcerFactory object" << std::endl;}
#else
    std::cerr << "OscProbCalcerFactory was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }

  else if (OscProbCalcerImplementationToCreate == "NuFASTLinear") {
#if UseNuFASTLinear==1
    OscProbCalcerNuFASTLinear* NuFASTLinear = new OscProbCalcerNuFASTLinear(OscProbCalcerConfig);
    Calcer = (OscProbCalcerBase*)NuFASTLinear;
    if (Verbose >= NuOscillator::INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscProbCalcerFactory object" << std::endl;}
#else
    std::cerr << "OscProbCalcerFactory was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }

  else if (OscProbCalcerImplementationToCreate == "OscProb") {
#if UseOscProb==1
    OscProbCalcerOscProb* OscProb = new OscProbCalcerOscProb(OscProbCalcerConfig);
    Calcer = (OscProbCalcerBase*)OscProb;
    if (Verbose >= NuOscillator::INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscProbCalcerFactory object" << std::endl;}
#else
    std::cerr << "OscProbCalcerFactory was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else {
    std::cerr << "OscProbCalcerFactory was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but this is not implemented within " << __FILE__ << std::endl;
    std::cerr << "Please correct the mistake or implement the new OscProbCalcer" << std::endl;
    throw;
  }

  return Calcer;
}
