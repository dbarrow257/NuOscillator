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

OscillatorBase::OscillatorBase(std::vector<std::string> OscProbCalcerImplementationToCreate_) {
  //DB Grab OscProbCalcerImplementationToCreate and fVerbose from config manager once implemented (And store yaml config such that it can be used to setup OscProbCalcer objects as well)
  fOscProbCalcerImplementationToCreate = OscProbCalcerImplementationToCreate_;
  fVerbose = INFO;
  fCosineZIgnored = false;

  fOscProbCalcers = std::vector<OscProbCalcerBase*>();
  fOscProbCalcerSet = false;
  
  //DB Move to seperate function
  fNCalcers = fOscProbCalcerImplementationToCreate.size();
  if (fNCalcers <= 0) {
    std::cerr << "Number of OscProbCalcerBase objects to be initialised is unreasonable" << std::endl;
    std::cout << "fNCalcers:" << fNCalcers << std::endl;
    throw;
  }

  for (size_t iOsc=0;iOsc<fOscProbCalcerImplementationToCreate.size();iOsc++) {
    OscProbCalcerBase* Calcer = InitialiseOscProbCalcer(fOscProbCalcerImplementationToCreate[iOsc]);
    if (Calcer == NULL) {
      std::cerr << "OscProbCalcer was not correctly initialised. Please check the setup" << std::endl;
      std::cerr << "OscProbCalcerImplementationToCreate:" << fOscProbCalcerImplementationToCreate[iOsc] << std::endl;
      std::cerr << "iOsc:" << iOsc << std::endl;
      throw;
    }
    fOscProbCalcers.push_back(Calcer);
  }
  fOscProbCalcerSet = true;
}


OscProbCalcerBase* OscillatorBase::InitialiseOscProbCalcer(std::string OscProbCalcerImplementationToCreate) {
  OscProbCalcerBase* Calcer;

  if (OscProbCalcerImplementationToCreate == "CUDAProb3") {
#ifdef UseCUDAProb3
    OscProbCalcerCUDAProb3* CUDAProb3 = new OscProbCalcerCUDAProb3("",fVerbose);
    Calcer = (OscProbCalcerBase*)CUDAProb3;
    if (fVerbose >= INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else if (OscProbCalcerImplementationToCreate == "Prob3ppLinear") {
#ifdef UseProb3ppLinear
    OscProbCalcerProb3ppLinear* Prob3ppLinear = new OscProbCalcerProb3ppLinear(fVerbose);
    Calcer = (OscProbCalcerBase*)Prob3ppLinear;
    if (fVerbose >= INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else if (OscProbCalcerImplementationToCreate == "ProbGPULinear") {
#ifdef UseProbGPULinear
    OscProbCalcerProbGPULinear* ProbGPULinear = new OscProbCalcerProbGPULinear(fVerbose);
    Calcer = (OscProbCalcerBase*)ProbGPULinear;
    if (fVerbose >= INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
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

void OscillatorBase::SetEnergyArrayInCalcer(std::vector<FLOAT_T> Array, int CalcerIndex) {
  if (CalcerIndex < 0 || CalcerIndex >= fNCalcers) {
    std::cerr << "Requested to set Energy array at invalid index within fOscProbCalcers array" << std::endl;
    std::cerr << "CalcerIndex:" << CalcerIndex << std::endl;
    std::cerr << "fNCalcers:" << fNCalcers << std::endl;
    throw;
  }
  if (fVerbose >= INFO) {std::cout << "Setting Energy array in OscProbCalcer Implementation:" << fOscProbCalcers[CalcerIndex]->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
  fOscProbCalcers[CalcerIndex]->SetEnergyArray(Array);
}

void OscillatorBase::SetCosineZArrayInCalcer(std::vector<FLOAT_T> Array, int CalcerIndex) {
  if (CalcerIndex < 0 || CalcerIndex >= fNCalcers) {
    std::cerr << "Requested to set CosineZ array at invalid index within fOscProbCalcers array" << std::endl;
    std::cerr << "CalcerIndex:"<< CalcerIndex << std::endl;
    std::cerr << "fNCalcers:" << fNCalcers << std::endl;
    throw;
  }
  if (fVerbose >= INFO) {std::cout << "Setting CosineZ array in OscProbCalcer Implementation:" << fOscProbCalcers[CalcerIndex]->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
  fOscProbCalcers[CalcerIndex]->SetCosineZArray(Array);
}

void OscillatorBase::CalculateProbabilities(std::vector<FLOAT_T> OscParams) {
  for (int CalcerIndex=0;CalcerIndex<fNCalcers;CalcerIndex++) {
    if (fVerbose >= INFO) {std::cout << "Calculating oscillation probabilities using OscProbCalcer Implementation:" << fOscProbCalcers[CalcerIndex]->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
    fOscProbCalcers[CalcerIndex]->Reweight(OscParams);
  }
}

void OscillatorBase::Setup() {
  for (int CalcerIndex=0;CalcerIndex<fNCalcers;CalcerIndex++) {
    if (fVerbose >= INFO) {std::cout << "Setting up OscProbCalcer Implementation:" << fOscProbCalcers[CalcerIndex]->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
    fOscProbCalcers[CalcerIndex]->Setup();
  }

  SanityCheck();
}

int OscillatorBase::ReturnNOscParams(int CalcerIndex) {
  if (CalcerIndex < 0 || CalcerIndex >= fNCalcers) {
    std::cerr << "Requested to Return NOscParams at invalid index within fOscProbCalcers array" << std::endl;
    std::cerr << "CalcerIndex:"<< CalcerIndex << std::endl;
    std::cerr << "fNCalcers:" << fNCalcers << std::endl;
    throw;
  }

  return fOscProbCalcers[CalcerIndex]->ReturnNOscParams();
}

void OscillatorBase::SanityCheck() {
  bool IsSane = fOscProbCalcerSet;

  for (int iCalcer=0;iCalcer<fNCalcers;iCalcer++) {
    IsSane = IsSane && fOscProbCalcers[iCalcer]->SanityCheck();
    IsSane = IsSane && (fCosineZIgnored == fOscProbCalcers[iCalcer]->ReturnCosineZIgnored());
  }

  if (!IsSane) {
    std::cerr << "OscillatorBase object has been found to not be 'sane' - The following booleans were expected to be true" << std::endl;
    std::cerr << "fOscProbCalcerSet:" << fOscProbCalcerSet << std::endl;
    std::cerr << "fCosineZIgnored:" << fCosineZIgnored << std::endl;
    for (int iCalcer=0;iCalcer<fNCalcers;iCalcer++) {
      std::cerr << "fOscProbCalcers[iCalcer]->SanityCheck():" << fOscProbCalcers[iCalcer]->SanityCheck() << std::endl;
      std::cerr << "fOscProbCalcers[iCalcer]->ReturnCosineZIgnored():" << fOscProbCalcers[iCalcer]->ReturnCosineZIgnored() << std::endl;
      std::cerr << "iCalcer:" << iCalcer << std::endl;
    }
    throw;
  } else {
    if (fVerbose >= INFO) {std::cout << "OscillatorBase instance passed SanityCheck" << std::endl;}
  }
}

void OscillatorBase::PrintWeights(int CalcerIndex) {
  if (CalcerIndex < 0 || CalcerIndex >= fNCalcers) {
    std::cerr << "Requested to PrintWeights at invalid index within fOscProbCalcers array" << std::endl;
    std::cerr << "CalcerIndex:"<< CalcerIndex << std::endl;
    std::cerr << "fNCalcers:" << fNCalcers << std::endl;
    throw;
  }

  fOscProbCalcers[CalcerIndex]->PrintWeights();
}
