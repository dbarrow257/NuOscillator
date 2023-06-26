#include "OscillatorBase.h"

#if UseCUDAProb3==1
#include "OscProbCalcer_CUDAProb3.h"
#endif

#if UseCUDAProb3Linear==1
#include "OscProbCalcer_CUDAProb3Linear.h"
#endif

#if UseProb3ppLinear==1
#include "OscProbCalcer_Prob3ppLinear.h"
#endif

#if UseProbGPULinear==1
#include "OscProbCalcer_ProbGPULinear.h"
#endif

#include <iostream>

OscillatorBase::OscillatorBase(std::string ConfigName_) {
  fOscProbCalcerImplementationToCreate = std::vector<std::string>();
  fVerbose = NONE;
  fCosineZIgnored = false;

  fEvalPointsSetInConstructor = false;
  fCalculationTypeName = "";

  fOscProbCalcers = std::vector<OscProbCalcerBase*>();
  fOscProbCalcerSet = false;

  // Create config manager
  std::cout << "Reading config in OscillatorBase: " << ConfigName_ << std::endl;
  Config = YAML::LoadFile(ConfigName_);

  std::string Verbosity = Config["General"]["Verbosity"].as<std::string>();
  fVerbose = Verbosity_StrToInt(Verbosity);

  for (auto const &OscProbCalcer : Config["OscProbCalcer"]) {
    fOscProbCalcerImplementationToCreate.push_back(OscProbCalcer["Config"].as<std::string>());
  }
  if (fVerbose >= INFO) {
    std::cout << "Size of fOscProbCalcerImplementationToCreate:" << fOscProbCalcerImplementationToCreate.size() << std::endl;
    for (size_t i=0;i<fOscProbCalcerImplementationToCreate.size();i++) {
      std::cout << i << " " << fOscProbCalcerImplementationToCreate[i] << std::endl;
    }
  }

  fCosineZIgnored = Config["General"]["CosineZIgnored"].as<bool>();

  if (fVerbose >= INFO) {std::cout << "Read:" << ConfigName_ << "\n" << Config << std::endl;}
  InitialiseOscProbCalcers();
}

void OscillatorBase::InitialiseOscProbCalcers() {
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

OscProbCalcerBase* OscillatorBase::InitialiseOscProbCalcer(std::string OscProbCalcerImplementationToCreateString) {
  OscProbCalcerBase* Calcer;

  std::string OscProbCalcerImplementationToCreate = "";
  std::string OscProbCalcerConfigname = "";
  std::string Instance_Str = "";
  int Instance = -1;

  size_t Delimiter = OscProbCalcerImplementationToCreateString.find(":");
  if (Delimiter != std::string::npos) {
    OscProbCalcerImplementationToCreate = OscProbCalcerImplementationToCreateString.substr(0,Delimiter);

    std::string SubStr = OscProbCalcerImplementationToCreateString.substr(Delimiter+1,OscProbCalcerImplementationToCreateString.size());
    Delimiter = SubStr.find(":");
    OscProbCalcerConfigname = SubStr.substr(0,Delimiter);

    if (Delimiter == std::string::npos) {
      Instance_Str = "0";
    } else {
      Instance_Str = SubStr.substr(Delimiter+1,SubStr.size());
    }
    Instance = std::stoi(Instance_Str);
    
  } else {
    std::cerr << "Expected a string formatted as: 'OscProbCalcerImplementation:PathToYAMLConfig:InstanceNumber'" << std::endl;
    std::cerr << "Recieved:" << OscProbCalcerImplementationToCreateString << std::endl;
    throw;
  }

  if (fVerbose >= INFO) {
    std::cout << "Parsed " << OscProbCalcerImplementationToCreateString << " and found:" << std::endl;
    std::cout << "\tOscProbCalcerImplementationToCreate:" << OscProbCalcerImplementationToCreate << std::endl;
    std::cout << "\tOscProbCalcerConfigname:" << OscProbCalcerConfigname << std::endl;
    std::cout << "\tInstance:" << Instance << std::endl;
  }

  if (OscProbCalcerImplementationToCreate == "CUDAProb3") {
#if UseCUDAProb3==1
    OscProbCalcerCUDAProb3* CUDAProb3 = new OscProbCalcerCUDAProb3(OscProbCalcerConfigname,Instance);
    Calcer = (OscProbCalcerBase*)CUDAProb3;
    if (fVerbose >= INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }

  else if (OscProbCalcerImplementationToCreate == "CUDAProb3Linear") {
#if UseCUDAProb3Linear==1
    OscProbCalcerCUDAProb3Linear* CUDAProb3Linear = new OscProbCalcerCUDAProb3Linear(OscProbCalcerConfigname,Instance);
    Calcer = (OscProbCalcerBase*)CUDAProb3Linear;
    if (fVerbose >= INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else if (OscProbCalcerImplementationToCreate == "Prob3ppLinear") {
#if UseProb3ppLinear==1
    OscProbCalcerProb3ppLinear* Prob3ppLinear = new OscProbCalcerProb3ppLinear(OscProbCalcerConfigname,Instance);
    Calcer = (OscProbCalcerBase*)Prob3ppLinear;
    if (fVerbose >= INFO) {std::cout << "Initalised OscProbCalcer Implementation:" << Calcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else if (OscProbCalcerImplementationToCreate == "ProbGPULinear") {
#if UseProbGPULinear==1
    OscProbCalcerProbGPULinear* ProbGPULinear = new OscProbCalcerProbGPULinear(OscProbCalcerConfigname,Instance);
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
  std::cout << "CalcerIndex:" << CalcerIndex << std::endl;
  if (CalcerIndex < 0 || CalcerIndex >= fNCalcers) {
    std::cerr << "Requested to set Energy array at invalid index within fOscProbCalcers array" << std::endl;
    std::cerr << "CalcerIndex:" << CalcerIndex << std::endl;
    std::cerr << "fNCalcers:" << fNCalcers << std::endl;
    throw;
  }
  if (fOscProbCalcers[CalcerIndex]->ReturnHasSetEnergyArray()) {
    std::cerr << "Have already set the Energy array in the requested OscProbCalcer" << std::endl;
    std::cerr << "This seems like a fault in the setup" << std::endl;
    throw;
  }
  if (fVerbose >= INFO) {
    std::cout << "Setting Energy array in OscProbCalcer Implementation:" << fOscProbCalcers[CalcerIndex]->ReturnImplementationName() << " in OscillatorBase object" << std::endl;
    std::cout << "Using CalcerIndex:" << CalcerIndex << std::endl;
  }
  fOscProbCalcers[CalcerIndex]->SetEnergyArray(Array);
}

void OscillatorBase::SetCosineZArrayInCalcer(std::vector<FLOAT_T> Array, int CalcerIndex) {
  if (CalcerIndex < 0 || CalcerIndex >= fNCalcers) {
    std::cerr << "Requested to set CosineZ array at invalid index within fOscProbCalcers array" << std::endl;
    std::cerr << "CalcerIndex:"<< CalcerIndex << std::endl;
    std::cerr << "fNCalcers:" << fNCalcers << std::endl;
    throw;
  }
  if (fOscProbCalcers[CalcerIndex]->ReturnHasSetCosineZArray()) {
    std::cerr << "Have already set the CosineZ array in the requested OscProbCalcer" << std::endl;
    std::cerr << "This seems like a fault in the setup"<< std::endl;
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

int OscillatorBase::ReturnNEnergyPoints(int CalcerIndex) {
  if (CalcerIndex < 0 || CalcerIndex >= fNCalcers) {
    std::cerr << "Requested to Return NEnergyPoints at invalid index within fOscProbCalcers array" << std::endl;
    std::cerr << "CalcerIndex:"<< CalcerIndex << std::endl;
    std::cerr << "fNCalcers:" << fNCalcers << std::endl;
    throw;
  }

  return fOscProbCalcers[CalcerIndex]->ReturnNEnergyPoints();
}

const FLOAT_T* OscillatorBase::ReturnPointerToWeightinCalcer(int CalcerIndex, int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal) {
  if (CalcerIndex < 0 || CalcerIndex >= fNCalcers) {
    std::cerr << "Requested to ReturnPointerToWeightinCalcer at invalid index within fOscProbCalcers array" << std::endl;
    std::cerr << "CalcerIndex:"<< CalcerIndex << std::endl;
    std::cerr << "fNCalcers:" << fNCalcers << std::endl;
    throw;
  }

  return fOscProbCalcers[CalcerIndex]->ReturnPointerToWeight(InitNuFlav,FinalNuFlav,EnergyVal,CosineZVal);
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

std::string OscillatorBase::ReturnImplementationName() {
  std::string ReturnString = "";
  
  ReturnString += fCalculationTypeName;
  for (int iCalcer=0;iCalcer<fNCalcers;iCalcer++) {
    ReturnString += "_"+fOscProbCalcers[iCalcer]->ReturnImplementationName();
  }
  
  return ReturnString;
}

bool OscillatorBase::HasOscProbCalcerGotOscillationChannel(int GeneratedFlavour, int DetectedFlavour, int CalcerIndex) {
  if (CalcerIndex < 0 || CalcerIndex >= fNCalcers) {
    std::cerr << "Requested to PrintWeights at invalid index within fOscProbCalcers array" << std::endl;
    std::cerr << "CalcerIndex:"<< CalcerIndex << std::endl;
    std::cerr << "fNCalcers:" << fNCalcers << std::endl;
    throw;
  }

  return fOscProbCalcers[CalcerIndex]->HasOscillationChannel(GeneratedFlavour,DetectedFlavour);
}
