#include "Oscillator/OscillatorBase.h"

#include "OscProbCalcer/OscProbCalcerFactory.h"

#if UseOscProb==1
#include "OscProbCalcer/OscProbCalcer_OscProb.h"
#endif

#include <iostream>

OscillatorBase::OscillatorBase(std::string ConfigName_) {
  fVerbose = NONE;
  fCosineZIgnored = false;

  fEvalPointsSetInConstructor = false;
  fCalculationTypeName = "";

  fOscProbCalcerSet = false;

  // Create config manager
  std::cout << "Reading config in OscillatorBase: " << ConfigName_ << std::endl;
  Config = YAML::LoadFile(ConfigName_);

  std::string Verbosity = Config["General"]["Verbosity"].as<std::string>();
  fVerbose = Verbosity_StrToInt(Verbosity);

  fCosineZIgnored = Config["General"]["CosineZIgnored"].as<bool>();

  InitialiseOscProbCalcer();
  
  if (fVerbose >= INFO) {std::cout << "Read:" << ConfigName_ << "\n" << Config << std::endl;}
}

OscillatorBase::~OscillatorBase() {
  delete fOscProbCalcer;
}

void OscillatorBase::InitialiseOscProbCalcer() {
  OscProbCalcerFactory* OscProbCalcFactory = new OscProbCalcerFactory();
  fOscProbCalcer = OscProbCalcFactory->CreateOscProbCalcer(Config);
  fOscProbCalcerSet = true;
}

void OscillatorBase::SetEnergyArrayInCalcer(std::vector<FLOAT_T> Array) {
  if (fOscProbCalcer->ReturnHasSetEnergyArray()) {
    std::cerr << "Have already set the Energy array in the requested OscProbCalcer" << std::endl;
    std::cerr << "This seems like a fault in the setup" << std::endl;
    throw;
  }
  if (fVerbose >= INFO) {
    std::cout << "Setting Energy array in OscProbCalcer Implementation:" << fOscProbCalcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;
  }
  fOscProbCalcer->SetEnergyArray(Array);
}

void OscillatorBase::SetCosineZArrayInCalcer(std::vector<FLOAT_T> Array) {
  if (fOscProbCalcer->ReturnHasSetCosineZArray()) {
    std::cerr << "Have already set the CosineZ array in the requested OscProbCalcer" << std::endl;
    std::cerr << "This seems like a fault in the setup"<< std::endl;
    throw;
  }
  if (fVerbose >= INFO) {std::cout << "Setting CosineZ array in OscProbCalcer Implementation:" << fOscProbCalcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
  fOscProbCalcer->SetCosineZArray(Array);
}

// It is assumed that all OscProbCalcers within a single instance of an OscillatorBase object take the same oscillation probabilities
void OscillatorBase::CalculateProbabilities(std::vector<FLOAT_T> OscParams) {
  if (fVerbose >= INFO) {std::cout << "Calculating oscillation probabilities using OscProbCalcer Implementation:" << fOscProbCalcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
  fOscProbCalcer->Reweight(OscParams);
}

void OscillatorBase::Setup() {
  if (fVerbose >= INFO) {std::cout << "Setting up OscProbCalcer Implementation:" << fOscProbCalcer->ReturnImplementationName() << " in OscillatorBase object" << std::endl;}
  fOscProbCalcer->Setup();

  SanityCheck();
}

int OscillatorBase::ReturnNOscParams() {
  return fOscProbCalcer->ReturnNOscParams();
}

int OscillatorBase::ReturnNEnergyPoints() {
  return fOscProbCalcer->ReturnNEnergyPoints();
}

const FLOAT_T* OscillatorBase::ReturnPointerToWeightinCalcer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal) {
  const FLOAT_T* Pointer = fOscProbCalcer->ReturnPointerToWeight(InitNuFlav,FinalNuFlav,EnergyVal,CosineZVal);
  return Pointer;
}

void OscillatorBase::SanityCheck() {
  bool IsSane = fOscProbCalcerSet;

  IsSane = IsSane && fOscProbCalcer->SanityCheck();
  IsSane = IsSane && (fCosineZIgnored == fOscProbCalcer->ReturnCosineZIgnored());

  if (!IsSane) {
    std::cerr << "OscillatorBase object has been found to not be 'sane' - The following booleans were expected to be true" << std::endl;
    std::cerr << "fOscProbCalcerSet:" << fOscProbCalcerSet << std::endl;
    std::cerr << "fCosineZIgnored:" << fCosineZIgnored << std::endl;
    std::cerr << "fOscProbCalcer->SanityCheck():" << fOscProbCalcer->SanityCheck() << std::endl;
    std::cerr << "fOscProbCalcer->ReturnCosineZIgnored():" << fOscProbCalcer->ReturnCosineZIgnored() << std::endl;
    throw;
  } else {
    if (fVerbose >= INFO) {std::cout << "OscillatorBase instance passed SanityCheck" << std::endl;}
  }
}

void OscillatorBase::PrintWeights() {
  fOscProbCalcer->PrintWeights();
}

std::string OscillatorBase::ReturnImplementationName() {
  std::string ReturnString = fCalculationTypeName+"_"+fOscProbCalcer->ReturnImplementationName();
  return ReturnString;
}

bool OscillatorBase::HasOscProbCalcerGotOscillationChannel(int GeneratedFlavour, int DetectedFlavour) {
  return fOscProbCalcer->HasOscillationChannel(GeneratedFlavour,DetectedFlavour);
}
