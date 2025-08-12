#include "OscProbCalcer/OscProbCalcerBase.h"

#include <cmath>
#include <math.h>

#include <iostream>
#include <iomanip>

OscProbCalcerBase::OscProbCalcerBase(YAML::Node InputConfig_) {
  // Set default values of all variables within this base object
  fVerbose = NuOscillator::NONE;

  fNNeutrinoTypes = DUMMYVAL;
  fNeutrinoTypes = std::vector<int>();

  fNOscillationChannels = DUMMYVAL;
  fOscillationChannels = std::vector<NuOscillator::OscillationChannel>();

  fNEnergyPoints = DUMMYVAL;
  fNCosineZPoints = DUMMYVAL;
  fEnergyArray = std::vector<FLOAT_T>();
  fCosineZArray = std::vector<FLOAT_T>();

  fNWeights = DUMMYVAL;
  fWeightArray = std::vector<FLOAT_T>();

  fNOscParams = DUMMYVAL;
  fOscParamsCurr = std::vector<FLOAT_T>();
  fExpectedOscillationParameterNames = std::vector<std::string>();
  fOscParams = std::vector<FLOAT_T*>();

  fCosineZIgnored = false;

  fEnergyArraySet = false;
  fCosineZArraySet = false;
  fPropagatorSet = false;
  fWeightArrayInit = false;
  fNuMappingSet = false;

  PrecisionLimit = 1.0e-6;
  
  Config = InputConfig_;

  std::string Verbosity = Config["General"]["Verbosity"].as<std::string>();
  fVerbose = Verbosity_StrToInt(Verbosity);

  if (!Config["OscProbCalcerSetup"]) {
    std::cerr << "Did not find the 'OscProbCalcerSetup' Node within the config" << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  
  if (Config["OscProbCalcerSetup"]["PrecisionLimit"]) {
    //Update the default valye of 1.0e-6
    PrecisionLimit = Config["OscProbCalcerSetup"]["PrecisionLimit"].as<FLOAT_T>();
  }

  if (!Config["OscProbCalcerSetup"]["ImplementationName"]) {
    std::cerr << "Did not find the 'OscProbCalcerSetup''ImplementationName' Node within the config" << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  fImplementationName = Config["OscProbCalcerSetup"]["ImplementationName"].as<std::string>();
  if (fVerbose >= NuOscillator::INFO) {
    std::cout << "From config, found implementation:" << fImplementationName << std::endl;
  }
  
  if (!Config["OscProbCalcerSetup"]["OscChannelMapping"]) {
    std::cerr << "Expected to find a 'OscChannelMapping' Node within the 'OscProbCalcerSetup'Node" << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  InitialiseOscillationChannelMapping();

  fUseLegacyMode = false;
  fUseLegacyMode_OscParsSet = false;
  if (Config["OscProbCalcerSetup"]["UseLegacyMode"]) {
    fUseLegacyMode = Config["OscProbCalcerSetup"]["UseLegacyMode"].as<bool>();
    if (fUseLegacyMode) {
      std::cerr << "WARNING - Using NuOscillator in Legacy mode! Expecting user to use OscProbCalcerBase::Reweight(const std::vector<FLOAT>& OscParams)" << std::endl;
    }
  }
  
}

OscProbCalcerBase::~OscProbCalcerBase() {
}

void OscProbCalcerBase::SetEnergyArray(std::vector<FLOAT_T> EnergyArray) {
  if (fEnergyArraySet) {
    // Already defined the Energy array, or the implementation is designed such not to care about it
    if (fVerbose >= NuOscillator::INFO) {std::cout << "EnergyArray was already found to be set in implementation:" << fImplementationName << std::endl;}
    return;
  }

  if (EnergyArray.size()==0) {
    std::cerr << "Invalid array passed to OscProbCalcerBase::SetEnergyArray as it has size 0" << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  fEnergyArray = EnergyArray;
  for (size_t iEnergy=0;iEnergy<fEnergyArray.size();iEnergy++) {
    if (fEnergyArray[iEnergy] <= 0) {
      std::cerr << "Found a negative neutrino energy. This indicates a problem in the array passed to: void OscProbCalcerBase::SetEnergyArray(std::vector<FLOAT_T> EnergyArray)" << std::endl;
      std::cerr << "iEnergy:" << iEnergy << std::endl;
      std::cerr << "fEnergyArray[iEnergy]:" << fEnergyArray[iEnergy] << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }

  for (size_t iEnergy=1;iEnergy<fEnergyArray.size();iEnergy++) {
    if (fEnergyArray[iEnergy-1] > fEnergyArray[iEnergy]) {
      std::cerr << "Found that the Energy array given to OscProbCalcerBase object is not ordered. This will become a problem when performing the binary search for indices" << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }

  fNEnergyPoints = EnergyArray.size();
  fEnergyArraySet = true;
  if (fVerbose >= NuOscillator::INFO) {
    std::cout << "Set EnergyArray in implementation:" << fImplementationName << " -" << std::endl;
    for (size_t iEnergy=1;iEnergy<fEnergyArray.size();iEnergy++) {
      std::cout << fEnergyArray[iEnergy] << ", ";
    }
    std::cout << std::endl;
  }
}

void OscProbCalcerBase::SetCosineZArray(std::vector<FLOAT_T> CosineZArray) {
  if (fCosineZIgnored) {
    // The implementation is designed such not to care about it
    if (fVerbose >= NuOscillator::INFO) {std::cout << "CosineZArray is ignored in implementation:" << fImplementationName << std::endl;}
    return;
  }

  if (fCosineZArraySet) {
    // Already defined the CosineZ array, or the implementation is designed such not to care about it
    if (fVerbose >= NuOscillator::INFO) {std::cout << "CosineZArray was already found to be set in implementation:" << fImplementationName << std::endl;}
    return;
  }

  if (CosineZArray.size()==0) {
    std::cerr << "Invalid array passed to OscProbCalcerBase::SetCosineZArray as it has size 0" << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  fCosineZArray = CosineZArray;
  for (size_t iCosineZ=0;iCosineZ<fCosineZArray.size();iCosineZ++) {
    if (fCosineZArray[iCosineZ] < -1.0 || fCosineZArray[iCosineZ] > 1.0) {
      std::cerr<< "Found a CosineZ outside of [-1.0,1.0]. This indicates a problem in the array passed to: void OscProbCalcerBase::SetCosineZArray(std::vector<FLOAT_T> CosineZArray)" << std::endl;
      std::cerr << "iCosineZ:" << iCosineZ << std::endl;
      std::cerr << "fCosineZArray[iCosineZ]:" << fCosineZArray[iCosineZ] << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }

  for (size_t iCosineZ=1;iCosineZ<fCosineZArray.size();iCosineZ++) {
    if (fCosineZArray[iCosineZ-1] > fCosineZArray[iCosineZ]) {
      std::cerr << "Found that the CosineZ array given to OscProbCalcerBase object is not ordered. This will become a problem when performing the binary search for indices" << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }

  fNCosineZPoints = CosineZArray.size();
  fCosineZArraySet = true;
  if (fVerbose >= NuOscillator::INFO) {std::cout << "Set CosineZArray in implementation:" << fImplementationName << std::endl;}
}

void OscProbCalcerBase::IgnoreCosineZBinning(bool Ignore) {
  if (Ignore) {
    fCosineZArraySet = true;
    fCosineZIgnored = true;
    if (fVerbose >= NuOscillator::INFO) {std::cout << "Ignoring CosineZArray in implementation:" << fImplementationName << std::endl;}
  }
}

void OscProbCalcerBase::Setup() {
  if (fVerbose>=NuOscillator::INFO) {std::cout << "About to setup OscProbCalcerBase implementation:" << fImplementationName << std::endl;}

  if (!fEnergyArraySet) {
    std::cerr << "Must call OscProbCalcerBase::SetEnergyArray(std::vector<FLOAT_T> EnergyArray) before trying to initialise propagator" << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  if (!fCosineZArraySet) {
    std::cerr << "Must call OscProbCalcerBase::SetCosineZArray(std::vector<FLOAT_T> CosineZArray) before trying to initialise propagator" << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  ResetCurrOscParams();
  if (fVerbose>=NuOscillator::INFO) {std::cout << "Reset Saved OscParams in OscProbCalcerBase implementation:" << fImplementationName << std::endl;}
  IntialiseWeightArray();
  if (fVerbose>=NuOscillator::INFO) {std::cout << "Initialised fWeightArray in OscProbCalcerBase implementation:" << fImplementationName << std::endl;}

  SetupPropagator();
  fPropagatorSet = true;

  CheckOscillationParametersDefined();
  CheckNuFlavourMapping();

  SanityCheck();
}

// Neutrinos and antineutrinos are separated based on the sign of the flavour (Thus need to check whether the sign of both flavours is consistent)
// No other requirements are made based on the flavours
const FLOAT_T* OscProbCalcerBase::ReturnPointerToWeight(int InitNuFlav, int FinalNuFlav, FLOAT_T Energy, FLOAT_T CosineZ) {
  int Product = InitNuFlav*FinalNuFlav;
  if (Product < 0) {
    std::cerr << "Initial neutrino flavour and final neutrino flavour are different Neutrino types (one is positive integer and the other is negative)" << std::endl;
    std::cerr << "InitNuFlav:" << InitNuFlav << std::endl;
    std::cerr << "FinalNuFlav:" << FinalNuFlav << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  int NuTypeIndex = ReturnNuTypeFromFlavour(InitNuFlav);
  int OscChanIndex = ReturnOscChannelIndexFromFlavours(std::abs(InitNuFlav),std::abs(FinalNuFlav));
  int CosineZIndex = -1;
  if (!ReturnCosineZIgnored()) {
    CosineZIndex = ReturnCosineZIndexFromValue(CosineZ);
  }
  int EnergyIndex = ReturnEnergyIndexFromValue(Energy);

  int WeightArrayIndex;
  if (!ReturnCosineZIgnored()) {
    WeightArrayIndex = ReturnWeightArrayIndex(NuTypeIndex,OscChanIndex,EnergyIndex,CosineZIndex);
  } else {
    WeightArrayIndex = ReturnWeightArrayIndex(NuTypeIndex,OscChanIndex,EnergyIndex);
  }
  
  if (WeightArrayIndex < 0 || WeightArrayIndex >= (int)fWeightArray.size()) {
    std::cerr << "Array index in fWeightArray is outside of the array size. This indicates that the implementation of ReturnWeightArrayIndex is incorrect." << std::endl;
    std::cerr << "NuTypeIndex:" << NuTypeIndex << std::endl;
    std::cerr << "OscChanIndex:" << OscChanIndex << std::endl;
    std::cerr << "CosineZIndex:" << CosineZIndex << std::endl;
    std::cerr << "EnergyIndex:" << EnergyIndex << std::endl;
    std::cerr << "WeightArrayIndex:" << WeightArrayIndex << std::endl;
    std::cerr << "fWeightArray.size():" << fWeightArray.size() << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  if (fVerbose >= NuOscillator::VERBOSE) {std::cout << "Implementation:" << fImplementationName << " returned pointer to index " << WeightArrayIndex << std::endl;}
  return &(fWeightArray[WeightArrayIndex]);
}

std::vector<NuOscillator::OscillationProbability> OscProbCalcerBase::ReturnProbabilities() {
  std::vector<NuOscillator::OscillationProbability> ReturnVec(fNWeights);
  std::vector<bool> CheckVec(fNWeights,false);

  int Index;
  int NuType;
  NuOscillator::OscillationChannel OscChan;
  FLOAT_T Energy;
  FLOAT_T CosineZ;
  FLOAT_T Weight;
 
  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
    NuType = fNeutrinoTypes[iNuType];

    for (int iOscChan=0;iOscChan<fNOscillationChannels;iOscChan++) {
      OscChan = fOscillationChannels[iOscChan];

      if (fCosineZIgnored) {
	for (int iEnergy=0;iEnergy<fNEnergyPoints;iEnergy++) {
	  Energy = fEnergyArray[iEnergy];
	  CosineZ = DUMMYVAL;

	  Index = ReturnWeightArrayIndex(iNuType, iOscChan, iEnergy);
	  Weight = fWeightArray[Index];
	  
	  NuOscillator::OscillationProbability OscProb = {NuType,OscChan,Energy,CosineZ,Weight};
	  ReturnVec[Index] = OscProb;
	  CheckVec[Index] = true;
	}	
      } else {
	for (int iCosZ=0;iCosZ<fNCosineZPoints;iCosZ++) {
          CosineZ = fCosineZArray[iCosZ];

	  for (int iEnergy=0;iEnergy<fNEnergyPoints;iEnergy++) {
	    Energy = fEnergyArray[iEnergy];	 

	    Index = ReturnWeightArrayIndex(iNuType, iOscChan, iEnergy, iCosZ);
	    Weight = fWeightArray[Index];
	    
	    NuOscillator::OscillationProbability OscProb = {NuType,OscChan,Energy,CosineZ,Weight};
	    ReturnVec[Index] = OscProb;
	    CheckVec[Index] = true;
	  }
	}
      }
      
    }
  }
  for (int iPoint=0;iPoint<fNWeights;iPoint++) {
    if (CheckVec[iPoint] == false) {
      std::cerr << "Index:" << iPoint << " has not been filled within the returning vector of std::vector<NuOscillator::OscillationProbability> OscProbCalcerBase::ReturnProbabilities()" << std::endl;
      std::cerr << "Indicates a problem!" << std::endl;
      throw std::runtime_error("Invalid setup");
    } 
  }

  return ReturnVec;
}

/* 
 * fUseLegacyMode_OscParsSet is only true if OscProbCalcerBase::Reweight() called from OscProbCalcerBase::Reweight(const std::vector<FLOAT_T>& OscParams_)
 *
 * fUseLegacyMode = true,  fUseLegacyMode_OscParsSet = true   --- Valid; Use Reweight(const std::vector<FLOAT_T>& OscParams_), not Reweight()
 * fUseLegacyMode = true,  fUseLegacyMode_OscParsSet = false  --- Invalid; Calling Reweight() without calling DefineParameter(std::string ParName, FLOAT_T* ParVal)
 * fUseLegacyMode = false, fUseLegacyMode_OscParsSet = true   --- Invalid; Calling Reweight(const std::vector<FLOAT_T>& OscParams_) after calling DefineParameter(std::string ParName, FLOAT_T* ParVal)
 * fUseLegacyMode = false, fUseLegacyMode_OscParsSet = false  --- Valid; Use Reweight(), not Reweight(const std::vector<FLOAT_T>& OscParams_)
 *
 * Consequently, fUseLegacyMode == fUseLegacyMode_OscParsSet to be Valid
 */ 

void OscProbCalcerBase::Reweight() {
  if (fVerbose >= NuOscillator::INFO) {std::cout << "Implementation:" << fImplementationName << " starting reweight" << std::endl;}

  if (!(fUseLegacyMode == fUseLegacyMode_OscParsSet)) {
    std::cerr << "Invalid running mode -" << std::endl;
    std::cerr << "fUseLegacyMode:" << fUseLegacyMode << std::endl;
    std::cerr << "fUseLegacyMode_OscParsSet:" << fUseLegacyMode_OscParsSet << std::endl;

    //fUseLegacyMode = true, fUseLegacyMode_OscParsSet = false
    if (fUseLegacyMode) {
      std::cerr << "Running in Legacy Mode but trying to call Reweight() - This will not work. Call OscProbCalcerBase::Reweight(const std::vector<FLOAT_T>& OscParams) instead!" << std::endl;
      throw std::runtime_error("Invalid used of Legacy Mode and OscProbCalcerBase::Reweight() call");
    }

    if (fUseLegacyMode_OscParsSet) {
      std::cerr << "Not running in Legacy Mode but trying to call Reweight(const std::vector<FLOAT_T>& OscParams_) - This will not work. Call OscProbCalcerBase::Reweight() instead!" << std::endl;
      throw std::runtime_error("Invalid used of Legacy Mode and OscProbCalcerBase::Reweight() call");
    }
  }
  
  if (!AreOscParamsChanged()) {
    return;
  }
  SetCurrOscParams();

  CalculateProbabilities();
  SanitiseProbabilities();
  if (fVerbose >= NuOscillator::INFO) {std::cout << "Implementation:" << fImplementationName << " completed reweight and was found to have sensible oscillation weights" << std::endl;}
}

void OscProbCalcerBase::Reweight(const std::vector<FLOAT_T>& OscParams_) {
  if ((int)OscParams_.size() != fNOscParams) {
    std::cerr << "Number of oscillation parameters passed to calculater does not match that expected by the implementation" << std::endl;
    std::cerr << "OscParams.size():" << OscParams_.size() << std::endl;
    std::cerr << "fNOscParams:" << fNOscParams << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  
  for (int iOscPar=0;iOscPar<fNOscParams;iOscPar++) {
    *(fOscParams[iOscPar]) = OscParams_[iOscPar];
  }

  fUseLegacyMode_OscParsSet = true;
  Reweight();
  fUseLegacyMode_OscParsSet = false;
}

void OscProbCalcerBase::CheckOscillationParametersDefined() {
  // Using Legacy mode where oscillation parameters passed through function argument
  if (fUseLegacyMode) return;

  for (size_t iOscPar=0;iOscPar<fNOscParams;iOscPar++) {
    if (!fOscillationParametersSetCheck[iOscPar]) {
      std::cerr << "Oscillation parameter missing from definition!" << std::endl;
      std::cerr << "Missing definition of fExpectedOscillationParameterNames[iOscPar]:" << fExpectedOscillationParameterNames[iOscPar] << std::endl;
      throw std::runtime_error("Missing oscillation parameter definition");
    }
  }
}

void OscProbCalcerBase::DefineParameter(std::string ParName_, FLOAT_T* ParValue_) {
  for (size_t iOscPar=0;iOscPar<fNOscParams;iOscPar++) {
    if (fExpectedOscillationParameterNames[iOscPar] == ParName_) {

      if (!fOscillationParametersSetCheck[iOscPar]) {
	fOscParams[iOscPar] = ParValue_;
	fOscillationParametersSetCheck[iOscPar] = true;
	return;
      } else {
	// If we hit here, that means the parameter was set twice which we will treat as an error
	std::cerr << "Parameter: " << ParName_ << " was found to be set twice which is treated as an error!" << std::endl;
	throw std::runtime_error("Parameter was set twice:"+ParName_);
      }
      
    }
  }

  // If you hit this point - then the provided oscillation parameter name is not defined in the engine/model configuration and thus is not valid
  throw std::runtime_error("Invalid oscillation parameter: "+ParName_);
}

FLOAT_T OscProbCalcerBase::GetOscillationParameter(int Index) {
  // Check if the requested index is appropriate value
  if (!(Index >= 0 && Index < fNOscParams)) {
    std::cerr << "Oscillation parameter requested - Invalid index:" << Index << std::endl;
    std::cerr	<< "fNOscParams:" << fNOscParams << std::endl;
  }

  // Because NuOscillator doesn't own the memory, check it's atleast not a nullptr
  if (fOscParams[Index] != nullptr) {
    return *fOscParams[Index];
  }
  
  // If we hit here, the pointer of the requested oscillation parameter is null
  throw std::runtime_error("Requested oscillation parameter pointer is nullptr");
}

void OscProbCalcerBase::SanitiseProbabilities() {

  // Precompute these here
  const double lower_limit = -1.0*PrecisionLimit;
  const double upper_limit = 1.0 + PrecisionLimit;

  for (int iWeight=0;iWeight<fNWeights;++iWeight) {
    if (std::isnan(fWeightArray[iWeight])) {
      std::cerr << "Found nan probability in fWeightArray" << std::endl;
      std::cerr << "iWeight:" << iWeight << std::endl;
      throw std::runtime_error("Invalid probability");
    }

    if ((fWeightArray[iWeight] >= 0.0) && (fWeightArray[iWeight] <= 1.0)) {
      //Check if it's between 0.0 and 1.0, if so continue to next event
      continue;
    }
    if (fWeightArray[iWeight] < lower_limit) {
      //Check if it's below 0.0-PrecisionLimit
      std::cerr << "Found probability which is below the allowable precision of: 0.0-" << PrecisionLimit << std::endl;
      std::cerr << "iWeight:" << iWeight << std::endl;
      std::cerr << "fWeightArray[iWeight]:" << fWeightArray[iWeight] << std::endl;
      throw std::runtime_error("Probability below zero");
    }
    if ((fWeightArray[iWeight] > lower_limit) && (fWeightArray[iWeight] < 0)) {
      //Check if it's just below 0 (within some precision) and set to 0 if it is
      fWeightArray[iWeight] = 0.;
      continue;
    }
    if ((fWeightArray[iWeight] > 1.0) && (fWeightArray[iWeight] < (upper_limit))) {
      //Check if it's just above 1 (within some precision) and set to 1 if it is
      fWeightArray[iWeight] = 1.;
      continue;
    }
    if (fWeightArray[iWeight] > (upper_limit)) {
      //Check if it's above 1.0+PrecisionLimit
      std::cerr << "Found probability which is above the allowable precision of: 1.0+" << PrecisionLimit << std::endl;
      std::cerr << "iWeight:" << iWeight << std::endl;
      std::cerr << "fWeightArray[iWeight]:" << fWeightArray[iWeight] << std::endl;
      throw std::runtime_error("Probability below zero");
    }
  }
  
}

bool OscProbCalcerBase::AreOscParamsChanged() {
  std::cout << "fOscParams.size():" << fOscParams.size() << std::endl;
  
  for (int iParam=0;iParam<fNOscParams;++iParam) {
    std::cout << "iParam:" << iParam << std::endl;
    std::cout << "fOscParams[iParam]:" << fOscParams[iParam] << std::endl;
    std::cout << "*fOscParams[iParam]:" << *fOscParams[iParam] << std::endl;
    std::cout << "fOscParamsCurr[iParam]:" << fOscParamsCurr[iParam] << std::endl;
    
    if (*fOscParams[iParam] != fOscParamsCurr[iParam]) {
      if (fVerbose >= NuOscillator::INFO) {std::cout << "Implementation:" << fImplementationName << " was found to have different oscillation parameters than the previous calculation" << std::endl;}
      return true;
    }
  }
  if (fVerbose >= NuOscillator::INFO) {std::cout << "Implementation:" << fImplementationName << " was found to have the same oscillation parameters than the previous calculation" << std::endl;}
  return false;
}

void OscProbCalcerBase::ResetCurrOscParams() {
  fOscParamsCurr = std::vector<FLOAT_T>(fNOscParams,DUMMYVAL);
}

void OscProbCalcerBase::SetCurrOscParams() {
  for (int iParam=0;iParam<fNOscParams;iParam++) {
    fOscParamsCurr[iParam] = *fOscParams[iParam];
  }
  if (fVerbose >= NuOscillator::INFO) {std::cout << "Saved oscillation parameters in Implementation:" << fImplementationName << std::endl;}
}

void OscProbCalcerBase::PrintWeights() {
  if (fVerbose >= NuOscillator::INFO) {std::cout << "Printing weights in Implementation:" << fImplementationName << std::endl;}
  /*
  for (size_t i=0;i<fWeightArray.size();i++) {
    std::cout << std::setw(10) << i << " | " << fWeightArray[i] << std::endl;
    if (fWeightArray[i] == DUMMYVAL) {
      std::cerr << "Found oscillation probability which has not been correctly calculated!" << std::endl;
      std::cerr << "This indicates that the mapping between the propagator and fWeightArray is incorrect" << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }
  */

  std::vector<NuOscillator::OscillationProbability> Probabilities = ReturnProbabilities();
  for (size_t iProb=0;iProb<Probabilities.size();iProb++) {
    std::cout << "Index:" << std::setw(10) << iProb << " | NuType:" << std::setw(3) << Probabilities[iProb].NuType << " | OscChan:" << std::setw(3) << Probabilities[iProb].OscChan.GeneratedFlavour << " -> " << std::setw(3) << Probabilities[iProb].OscChan.DetectedFlavour << " | Energy:" << std::setw(10) << Probabilities[iProb].Energy << " | CosZ:" << std::setw(10) << Probabilities[iProb].CosineZ << " | Prob:" << std::setw(10) << Probabilities[iProb].Probability << std::endl;
    if (Probabilities[iProb].Probability == DUMMYVAL) {
      std::cerr << "Found oscillation probability which has not been correctly calculated!" << std::endl;
      std::cerr << "This indicates that the mapping between the propagator and fWeightArray is incorrect" << std::endl;
      throw std::runtime_error("Invalid setup");
    } 
  }
  
}

int OscProbCalcerBase::ReturnEnergyIndexFromValue(FLOAT_T EnergyVal) {
  if (!fEnergyArraySet) {
    std::cerr << "Can not find Energy index as Energy array has not been set" << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  int EnergyIndex = -1;
  auto it = std::lower_bound(fEnergyArray.begin(), fEnergyArray.end(), EnergyVal);
  if (it == fEnergyArray.end() || *it != EnergyVal) {
    EnergyIndex = -1;
  } else {
    EnergyIndex = std::distance(fEnergyArray.begin(), it);
  } 
  
  if (EnergyIndex == -1) {
    std::cerr << "Did not find Energy in the array used in calculating oscillation probabilities" << std::endl;
    std::cerr << "Requested Energy:" << EnergyVal << std::endl;
    std::cerr << "Binning - " << std::endl;
    for (size_t iEnergy=0;iEnergy<fEnergyArray.size();iEnergy++) {
      std::cerr << fEnergyArray[iEnergy] << ", ";
    }
    std::cerr << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  if (fVerbose >= NuOscillator::VERBOSE) {std::cout << "Returning Energy index:" << EnergyIndex << " for Energy value:" << EnergyVal << " in Implementation:" << fImplementationName << std::endl;}
  return EnergyIndex;
}

int OscProbCalcerBase::ReturnCosineZIndexFromValue(FLOAT_T CosineZVal) {
  if (!fCosineZArraySet) {
    std::cerr << "Can not find CosineZ index as CosineZ array has not been set" << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  int CosineZIndex = -1;
  auto it = std::lower_bound(fCosineZArray.begin(), fCosineZArray.end(), CosineZVal);
  if (it == fCosineZArray.end() || *it != CosineZVal) {
    CosineZIndex = -1;
  } else {
    CosineZIndex = std::distance(fCosineZArray.begin(), it);
  } 
  
  if (CosineZIndex == -1) {
    std::cerr << "Did not find CosineZ in the array used in calculating oscillation probabilities" << std::endl;
    std::cerr << "Requested CosineZ:" << CosineZVal << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  if (fVerbose >= NuOscillator::VERBOSE) {std::cout << "Returning CosineZ index:" << CosineZIndex << " for CosineZ value:" << CosineZVal << " in Implementation:" << fImplementationName << std::endl;}
  return CosineZIndex;
}

int OscProbCalcerBase::ReturnOscChannelIndexFromFlavours(int InitFlav, int FinalFlav) {
  for (int iOscChan=0;iOscChan<fNOscillationChannels;iOscChan++) {
    if (InitFlav == fOscillationChannels[iOscChan].GeneratedFlavour && FinalFlav == fOscillationChannels[iOscChan].DetectedFlavour) {
      return iOscChan;
    }
  }

  std::cerr << "Did not find reasonable oscillation channel index for the requested generated flavour: " << InitFlav << " and detected flavour: " << FinalFlav << std::endl;
  PrintKnownOscillationChannels();
  throw std::runtime_error("Invalid setup");
}

int OscProbCalcerBase::ReturnNuTypeFromFlavour(int NuFlav) {
  int NuType = (NuFlav > 0) - (NuFlav < 0); // Calculates the sign of NuFlav
  
  for (int iType=0;iType<fNNeutrinoTypes;iType++) {
    if (NuType == fNeutrinoTypes[iType]) {
      if (fVerbose >= NuOscillator::VERBOSE) {std::cout << "Returning type:" << iType << " for NuFlav:" << NuFlav << " in Implementation:" << fImplementationName << std::endl;}
      return iType;
    }
  }

  std::cerr << "Requested Neutrino type is not defined within the NeutrinoType map!" << std::endl;
  std::cerr << "NuFlav:" << NuFlav << std::endl;
  std::cerr << "Associated NuType:" << NuType << std::endl;
  throw std::runtime_error("Invalid setup");
}

void OscProbCalcerBase::IntialiseWeightArray() {
  if (fVerbose >= NuOscillator::INFO) {std::cout << "Asking OscProbCalcerBase implementation:" << fImplementationName << " for the size" << std::endl;}
  fNWeights = DefineWeightArraySize();
  if (fNWeights <= 0) {
    std::cerr << "Number of weights which will be stored is less than 0. This indicates a fault in the calculation specific code: DefineWeightArraySize() is incorrect or has overflow-ed the 'long' type used as the return type" << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  if (fVerbose >= NuOscillator::INFO) {std::cout << "Asked OscProbCalcerBase implementation:" << fImplementationName << " for the size and got " << fNWeights << std::endl;}
  fWeightArray = std::vector<FLOAT_T>(fNWeights,DUMMYVAL);  
  fWeightArrayInit = true;
  if (fVerbose >= NuOscillator::INFO) {std::cout << "Initialising fWeightArray to be of size:" << fNWeights << " in Implementation:" << fImplementationName << std::endl;}
}

void OscProbCalcerBase::InitialiseNeutrinoTypesArray(int Size) {
  if (Size <= 0) {
    std::cerr << "Attempting to initialise fNeutrinoTypes array with size:" << Size << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  fNeutrinoTypes = std::vector<int>(Size,DUMMYVAL);
  if (fVerbose >= NuOscillator::INFO) {std::cout << "Initialising fNeutrinoTypes to be of size:" << Size << " in Implementation:" << fImplementationName << std::endl;}
}

void OscProbCalcerBase::InitialiseOscillationChannelMapping() {
  for (auto const &OscChannel : Config["OscProbCalcerSetup"]["OscChannelMapping"]) {
    NuOscillator::OscillationChannel myOscChan = ReturnOscillationChannel(OscChannel["Entry"].as<std::string>());
    fOscillationChannels.push_back(myOscChan);
  }
  fNOscillationChannels = fOscillationChannels.size();

  if (fVerbose >= NuOscillator::INFO) {PrintKnownOscillationChannels();}
}

void OscProbCalcerBase::CheckNuFlavourMapping() {
  if (fNNeutrinoTypes == DUMMYVAL || fNOscillationChannels == DUMMYVAL || fNNeutrinoTypes <= 0 || fNOscillationChannels <= 0) {
    std::cerr << "Number of neutrino types or flavours have not been correctly defined:" << std::endl;
    std::cerr << "fNNeutrinoTypes:" << fNNeutrinoTypes << std::endl;
    std::cerr << "fNOscillationChannels:" << fNOscillationChannels << std::endl;
    std::cerr << "DUMMYVAL:" << DUMMYVAL << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  if (fNNeutrinoTypes != (int)fNeutrinoTypes.size()) {
    std::cerr << "fNeutrinoTypes array not equal in size to fNNeutrinoTypes" << std::endl;
    std::cerr << "fNNeutrinoTypes:" << fNNeutrinoTypes << std::endl;
    std::cerr << "fNeutrinoTypes.size():" << fNeutrinoTypes.size() << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
    if (fNeutrinoTypes[iNuType]==DUMMYVAL) {
      std::cerr << "Found DUMMYVAL in fNeutrinoTypes" << std::endl;
      std::cerr << "iNuType:" << iNuType << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }

  if (fNOscillationChannels != (int)fOscillationChannels.size()) {
    std::cerr << "fOscillationChannels array not equal in size to fNOscillationChannels" << std::endl;
    std::cerr << "fNOscillationChannels:" << fNOscillationChannels << std::endl;
    std::cerr << "fOscillationChannels.size():" << fOscillationChannels.size() << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  for (int iNuFlav=0;iNuFlav<fNOscillationChannels;iNuFlav++) {
    if (fOscillationChannels[iNuFlav].DetectedFlavour == 0 || fOscillationChannels[iNuFlav].GeneratedFlavour == 0) {
      std::cerr << "Found DUMMYVAL in fOscillationChannels" << std::endl;
      std::cerr << "iNuFlav:" << iNuFlav << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }

  if (fVerbose >= NuOscillator::INFO) {std::cout << "NeutrinoType and NeutrinoFlavour mapping was found to be OK in Implementation:" << fImplementationName << std::endl;}
  fNuMappingSet = true;
}

bool OscProbCalcerBase::SanityCheck() {
  bool IsSane = fEnergyArraySet && fCosineZArraySet && fPropagatorSet && fWeightArrayInit && fNuMappingSet;
  
  if (!IsSane) {
    std::cerr << "OscProbCalcerBase object has been found to not be 'sane' - The following booleans were expected to be true" << std::endl;
    std::cerr << "fEnergyArraySet:" << fEnergyArraySet << std::endl;
    std::cerr << "fCosineZArraySet:" << fCosineZArraySet << std::endl;
    std::cerr << "fPropagatorSet:" << fPropagatorSet << std::endl;
    std::cerr << "fWeightArrayInit:" << fWeightArrayInit << std::endl;
    std::cerr << "fNuMappingSet:" << fNuMappingSet << std::endl;
    throw std::runtime_error("Invalid setup"); 
  } else {
    if (fVerbose >= NuOscillator::INFO) {std::cout << "Implementation:" << fImplementationName << " passed SanityCheck" << std::endl;}
  }

  return IsSane;
}

void OscProbCalcerBase::PrintKnownOscillationChannels() {
  std::cout << "Number of requested oscillation channels:" << fOscillationChannels.size() << std::endl;
  for (int i=0;i<fNOscillationChannels;i++) {
    std::cout << "\t" << i << " : " << fOscillationChannels[i].GeneratedFlavour << " (" << std::setw(10) << NeutrinoFlavour_IntToStr(fOscillationChannels[i].GeneratedFlavour) << ") -> " << fOscillationChannels[i].DetectedFlavour << " (" << std::setw(10) << NeutrinoFlavour_IntToStr(fOscillationChannels[i].DetectedFlavour) << ")" << std::endl;
  }
}

bool OscProbCalcerBase::HasOscillationChannel(int GeneratedFlavour, int DetectedFlavour) {
  for (int iOscChan=0;iOscChan<fNOscillationChannels;iOscChan++) {
    if (fOscillationChannels[iOscChan].GeneratedFlavour == GeneratedFlavour && fOscillationChannels[iOscChan].DetectedFlavour == DetectedFlavour) {
      return true;
    }
  }
  return false;
}
