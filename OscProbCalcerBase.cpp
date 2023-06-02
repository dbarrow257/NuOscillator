#include "OscProbCalcerBase.h"

#include <cmath>
#include <math.h>

#include <iostream>
#include <iomanip>

OscProbCalcerBase::OscProbCalcerBase() {
  // Set deafult values of all variables within this base object
  fVerbose = NONE;
  fImplementationName = std::string();

  fNNeutrinoTypes = DUMMYVAL;
  fNeutrinoTypes = std::vector<int>();
  fNInitialFlavours = DUMMYVAL;
  fInitialFlavours = std::vector<int>();
  fNFinalFlavours = DUMMYVAL;
  fFinalFlavours = std::vector<int>();

  fNEnergyPoints = DUMMYVAL;
  fNCosineZPoints = DUMMYVAL;
  fEnergyArray = std::vector<FLOAT_T>();
  fCosineZArray = std::vector<FLOAT_T>();

  fNWeights = DUMMYVAL;
  fWeightArray = std::vector<FLOAT_T>();

  fNOscParams = DUMMYVAL;
  fOscParamsCurr = std::vector<FLOAT_T>();

  fCosineZIgnored = false;

  fEnergyArraySet = false;
  fCosineZArraySet = false;
  fPropagatorSet = false;
  fWeightArrayInit = false;
  fNuMappingSet = false;
}

void OscProbCalcerBase::SetEnergyArray(std::vector<FLOAT_T> EnergyArray) {
  if (fEnergyArraySet) {
    // Already defined the Energy array, or the implementation is designed such not to care about it
    if (fVerbose >= INFO) {std::cout << "EnergyArray was already found to be set in implementation:" << fImplementationName << std::endl;}
    return;
  }

  fEnergyArray = EnergyArray;
  for (size_t iEnergy=0;iEnergy<fEnergyArray.size();iEnergy++) {
    if (fEnergyArray[iEnergy] <= 0) {
      std::cerr << "Found a negative neutrino energy. This indicates a problem in the array passed to: void OscProbCalcerBase::SetEnergyArray(std::vector<FLOAT_T> EnergyArray)" << std::endl;
      std::cerr << "iEnergy:" << iEnergy << std::endl;
      std::cerr << "fEnergyArray[iEnergy]:" << fEnergyArray[iEnergy] << std::endl;
      throw;
    }
  }

  fNEnergyPoints = EnergyArray.size();
  fEnergyArraySet = true;
  if (fVerbose >= INFO) {std::cout << "Set EnergyArray in implementation:" << fImplementationName << std::endl;}
}

void OscProbCalcerBase::SetCosineZArray(std::vector<FLOAT_T> CosineZArray) {
  if (fCosineZIgnored) {
    // The implementation is designed such not to care about it
    if (fVerbose >= INFO) {std::cout << "CosineZArray is ignored in implementation:" << fImplementationName << std::endl;}
    return;
  }

  if (fCosineZArraySet) {
    // Already defined the CosineZ array, or the implementation is designed such not to care about it
    if (fVerbose >= INFO) {std::cout << "CosineZArray was already found to be set in implementation:" << fImplementationName << std::endl;}
    return;
  }

  fCosineZArray = CosineZArray;
  for (size_t iCosineZ=0;iCosineZ<fCosineZArray.size();iCosineZ++) {
    if (fCosineZArray[iCosineZ] < -1.0 || fCosineZArray[iCosineZ] > 1.0) {
      std::cerr<< "Found a CosineZ outside of [-1.0,1.0]. This indicates a problem in the array passed to: void OscProbCalcerBase::SetCosineZArray(std::vector<FLOAT_T> CosineZArray)" << std::endl;
      std::cerr << "iCosineZ:" << iCosineZ << std::endl;
      std::cerr << "fCosineZArray[iCosineZ]:" << fCosineZArray[iCosineZ] << std::endl;
      throw;
    }
  }

  fNCosineZPoints = CosineZArray.size();
  fCosineZArraySet = true;
  if (fVerbose >= INFO) {std::cout << "Set CosineZArray in implementation:" << fImplementationName << std::endl;}
}

void OscProbCalcerBase::IgnoreCosineZBinning(bool Ignore) {
  //DB Add an actual ignore flag which can then be interrogated from a getter

  if (Ignore) {
    fCosineZArraySet = true;
    fCosineZIgnored = true;
    if (fVerbose >= INFO) {std::cout << "Ignoring CosineZArray in implementation:" << fImplementationName << std::endl;}
  }
}

void OscProbCalcerBase::Setup() {
  if (!fEnergyArraySet) {
    std::cerr << "Must call OscProbCalcerBase::SetEnergyArray(std::vector<FLOAT_T> EnergyArray) before trying to initialise propagator" << std::endl;
    throw;
  }

  if (!fCosineZArraySet) {
    std::cerr << "Must call OscProbCalcerBase::SetCosineZArray(std::vector<FLOAT_T> CosineZArray) before trying to initialise propagator" << std::endl;
    throw;
  }

  ResetCurrOscParams();
  IntialiseWeightArray();

  SetupPropagator();
  fPropagatorSet = true;

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
    throw;
  }

  int NuTypeIndex = ReturnNuTypeFromFlavour(InitNuFlav);
  int InitNuIndex = ReturnInitialIndexFromFlavour(InitNuFlav);
  int FinalNuIndex = ReturnFinalIndexFromFlavour(FinalNuFlav);
  int CosineZIndex = ReturnCosineZIndexFromValue(CosineZ);
  int EnergyIndex = ReturnEnergyIndexFromValue(Energy);

  int WeightArrayIndex = ReturnWeightArrayIndex(NuTypeIndex,InitNuIndex,FinalNuIndex,EnergyIndex,CosineZIndex);
  if (WeightArrayIndex < 0 || WeightArrayIndex >= (int)fWeightArray.size()) {
    std::cerr << "Array index in fWeightArray is outside of the array size. This indicates that the implementation of ReturnWeightArrayIndex is incorrect." << std::endl;
    std::cerr << "WeightArrayIndex:" << WeightArrayIndex << std::endl;
    std::cerr << "fWeightArray.size():" << fWeightArray.size() << std::endl;
    throw;
  }

  if (fVerbose >= INFO) {std::cout << "Implementation:" << fImplementationName << " returned pointer to index " << WeightArrayIndex << std::endl;}
  return &(fWeightArray[WeightArrayIndex]);
}

void OscProbCalcerBase::Reweight(std::vector<FLOAT_T> OscParams) {
  if (fVerbose >= INFO) {std::cout << "Implementation:" << fImplementationName << " starting reweight" << std::endl;}

  if ((int)OscParams.size() != fNOscParams) {
    std::cerr << "Number of oscillation parameters passed to calculater does not match that expected by the implementation" << std::endl;
    std::cerr << "OscParams.size():" << OscParams.size() << std::endl;
    std::cerr << "fNOscParams:" << fNOscParams << std::endl;
    throw;
  }

  if (!AreOscParamsChanged(OscParams)) {
    return;
  }
  SetCurrOscParams(OscParams);

  CalculateProbabilities(OscParams);
  SanitiseProbabilities();
  if (fVerbose >= INFO) {std::cout << "Implementation:" << fImplementationName << " completed reweight and was found to have sensible oscillation weights" << std::endl;}
}

void OscProbCalcerBase::SanitiseProbabilities() {
  for (int iWeight=0;iWeight<fNWeights;iWeight++) {
    if (std::isnan(fWeightArray[iWeight]) || fWeightArray[iWeight] < 0.0 || fWeightArray[iWeight] > 1.0) {
      std::cerr << "Found unreasonable weight in fWeightArray" << std::endl;
      std::cerr << "iWeight:" << iWeight << std::endl;
      std::cerr << "fWeightArray[iWeight]:" << fWeightArray[iWeight] << std::endl;
      throw;
    }
  }
}

bool OscProbCalcerBase::AreOscParamsChanged(std::vector<FLOAT_T> OscParamsToCheck) {
  for (int iParam=0;iParam<fNOscParams;iParam++) {
    if (OscParamsToCheck[iParam] != fOscParamsCurr[iParam]) {
      if (fVerbose >= INFO) {std::cout << "Implementation:" << fImplementationName << " was found to have different oscillation parameters than the previous calculation" << std::endl;}
      return true;
    }
  }
  if (fVerbose >= INFO) {std::cout << "Implementation:" << fImplementationName << " was found to have the same oscillation parameters than the previous calculation" << std::endl;}
  return false;
}

void OscProbCalcerBase::ResetCurrOscParams() {
  fOscParamsCurr = std::vector<FLOAT_T>(fNOscParams,DUMMYVAL);
}

void OscProbCalcerBase::SetCurrOscParams(std::vector<FLOAT_T> OscParamsToSave) {
  for (int iParam=0;iParam<fNOscParams;iParam++) {
    fOscParamsCurr[iParam] = OscParamsToSave[iParam];
  }
  if (fVerbose >= INFO) {std::cout << "Saved oscillation parameters in Implementation:" << fImplementationName << std::endl;}
}

void OscProbCalcerBase::PrintWeights() {
  if (fVerbose >= INFO) {std::cout << "Printing weights in Implementation:" << fImplementationName << std::endl;}
  for (size_t i=0;i<fWeightArray.size();i++) {
    std::cout << std::setw(10) << i << " | " << fWeightArray[i] << std::endl;
    if (fWeightArray[i] == DUMMYVAL) {
      std::cerr << "Found oscillation probability which has not been correctly calculated!" << std::endl;
      std::cerr << "This indicates that the mapping between the propagator and fWeightArray is incorrect" << std::endl;
      throw;
    }
  }
}

int OscProbCalcerBase::ReturnEnergyIndexFromValue(FLOAT_T EnergyVal) {
  if (!fEnergyArraySet) {
    std::cerr << "Can not find Energy index as Energy array has not been set" << std::endl;
    throw;
  }

  int EnergyIndex = -1;
  for (size_t iEnergy=0;iEnergy<fEnergyArray.size();iEnergy++) {
    if (EnergyVal == fEnergyArray[iEnergy]) {
      EnergyIndex = iEnergy;
      break;
    }
  }
  if (EnergyIndex == -1) {
    std::cerr << "Did not find Energy in the array used in calculating oscillation probabilities" << std::endl;
    std::cerr << "Requested Energy:" << EnergyVal << std::endl;
    throw;
  }

  if (fVerbose >= INFO) {std::cout << "Returning Energy index:" << EnergyIndex << " for Energy value:" << EnergyVal << " in Implementation:" << fImplementationName << std::endl;}
  return EnergyIndex;
}

int OscProbCalcerBase::ReturnCosineZIndexFromValue(FLOAT_T CosineZVal) {
  if (!fCosineZArraySet) {
    std::cerr << "Can not find CosineZ index as CosineZ array has not been set" << std::endl;
    throw;
  }

  int CosineZIndex = -1;
  for (size_t iCosineZ=0;iCosineZ<fCosineZArray.size();iCosineZ++) {
    if (CosineZVal == fCosineZArray[iCosineZ]) {
      CosineZIndex = iCosineZ;
      break;
    }
  }
  if (CosineZIndex == -1) {
    std::cerr << "Did not find CosineZ in the array used in calculating oscillation probabilities" << std::endl;
    std::cerr << "Requested CosineZ:" << CosineZVal << std::endl;
    throw;
  }

  if (fVerbose >= INFO) {std::cout << "Returning CosineZ index:" << CosineZIndex << " for CosineZ value:" << CosineZVal << " in Implementation:" << fImplementationName << std::endl;}
  return CosineZIndex;
}

int OscProbCalcerBase::ReturnInitialIndexFromFlavour(int InitFlav) {
  for (size_t iFlav=0;iFlav<fInitialFlavours.size();iFlav++) {
    if (fabs(InitFlav) == fInitialFlavours[iFlav]) {
      if (fVerbose >= INFO) {std::cout << "Returning index:" << iFlav << " for InitFlav:" << InitFlav << " in Implementation:" << fImplementationName << std::endl;}
      return iFlav;
    }
  }

  std::cerr << "Requested Initial Neutrino flavour is not defined within the InitialFlavour map!" << std::endl;
  std::cerr << "InitNuFlav:" << InitFlav << std::endl;
  std::cerr << "fNInitialFlavours:" << fNInitialFlavours << std::endl;
  throw;
}

int OscProbCalcerBase::ReturnFinalIndexFromFlavour(int FinalFlav) {
  for (size_t iFlav=0;iFlav<fFinalFlavours.size();iFlav++) {
    if (fabs(FinalFlav) == fFinalFlavours[iFlav]) {
      if (fVerbose >= INFO) {std::cout << "Returning index:" << iFlav << " for FinalFlav:" << FinalFlav << " in Implementation:" << fImplementationName << std::endl;}
      return iFlav;
    }
  }

  std::cerr << "Requested Final Neutrino flavour is not defined within the FinalFlavour map!" << std::endl;
  std::cerr << "InitNuFlav:" << FinalFlav << std::endl;
  std::cerr << "fNFinalFlavours:" << fNFinalFlavours << std::endl;
  throw;
}

int OscProbCalcerBase::ReturnNuTypeFromFlavour(int NuFlav) {
  int NuType = (NuFlav > 0) - (NuFlav < 0); // Calculates the sign of NuFlav
  
  for (int iType=0;iType<fNNeutrinoTypes;iType++) {
    if (NuType == fNeutrinoTypes[iType]) {
      if (fVerbose >= INFO) {std::cout << "Returning type:" << iType << " for NuFlav:" << NuFlav << " in Implementation:" << fImplementationName << std::endl;}
      return iType;
    }
  }

  std::cerr << "Requested Neutrino type is not defined within the NeutrinoType map!" << std::endl;
  std::cerr << "NuFlav:" << NuFlav << std::endl;
  std::cerr << "Associated NuType:" << NuType << std::endl;
  throw;
}

void OscProbCalcerBase::IntialiseWeightArray() {
  fNWeights = DefineWeightArraySize();
  fWeightArray = std::vector<FLOAT_T>(fNWeights,DUMMYVAL);  
  fWeightArrayInit = true;
  if (fVerbose >= INFO) {std::cout << "Initialising fWeightArray to be of size:" << fNWeights << " in Implementation:" << fImplementationName << std::endl;}
}

void OscProbCalcerBase::InitialiseNeutrinoTypesArray(int Size) {
  if (Size <= 0) {
    std::cerr << "Attempting to initialise fNeutrinoTypes array with size:" << Size << std::endl;
    throw;
  }
  fNeutrinoTypes = std::vector<int>(Size,DUMMYVAL);
  if (fVerbose >= INFO) {std::cout << "Initialising fNeutrinoTypes to be of size:" << Size << " in Implementation:" << fImplementationName << std::endl;}
}

void OscProbCalcerBase::InitialiseInitialFlavoursArray(int Size) {
  if (Size <= 0) {
    std::cerr << "Attempting to initialise fInitialFlavours array with size:" << Size << std::endl;
    throw;
  }
  if (fVerbose >= INFO) {std::cout << "Initialising fInitialFlavours to be of size:" << Size << " in Implementation:" << fImplementationName << std::endl;}
  fInitialFlavours = std::vector<int>(Size,DUMMYVAL);
}

void OscProbCalcerBase::InitialiseFinalFlavoursArray(int Size) {
  if (Size <= 0) {
    std::cerr << "Attempting to initialise fFinalFlavours array with size:" << Size << std::endl;
    throw;
  }
  if (fVerbose >= INFO) {std::cout << "Initialising fFinalFlavours to be of size:" << Size << " in Implementation:" << fImplementationName << std::endl;}
  fFinalFlavours = std::vector<int>(Size,DUMMYVAL);
}

void OscProbCalcerBase::CheckNuFlavourMapping() {
  if (fNNeutrinoTypes == DUMMYVAL || fNInitialFlavours == DUMMYVAL || fNFinalFlavours == DUMMYVAL) {
    std::cerr << "Number of neutrino types or flavours have not been correctly defined:" << std::endl;
    std::cerr << "fNNeutrinoTypes:" << fNNeutrinoTypes << std::endl;
    std::cerr << "fNInitialFlavours:" << fNInitialFlavours << std::endl;
    std::cerr << "fNFinalFlavours:" << fNFinalFlavours << std::endl;
    std::cerr << "DUMMYVAL:" << DUMMYVAL << std::endl;
    throw;
  }

  if (fNNeutrinoTypes != (int)fNeutrinoTypes.size()) {
    std::cerr << "fNeutrinoTypes array not equal in size to fNNeutrinoTypes" << std::endl;
    std::cerr << "fNNeutrinoTypes:" << fNNeutrinoTypes << std::endl;
    std::cerr << "fNeutrinoTypes.size():" << fNeutrinoTypes.size() << std::endl;
    throw;
  }
  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
    if (fNeutrinoTypes[iNuType]==DUMMYVAL) {
      std::cerr << "Found DUMMYVAL in fNeutrinoTypes" << std::endl;
      std::cerr << "iNuType:" << iNuType << std::endl;
      throw;
    }
  }

  if (fNInitialFlavours != (int)fInitialFlavours.size()) {
    std::cerr << "fInitialFlavours array not equal in size to fNInitialFlavours" << std::endl;
    std::cerr << "fNInitialFlavours:" << fNInitialFlavours << std::endl;
    std::cerr << "fInitialFlavours.size():" << fInitialFlavours.size() << std::endl;
    throw;
  }
  for (int iNuFlav=0;iNuFlav<fNInitialFlavours;iNuFlav++) {
    if (fInitialFlavours[iNuFlav]==DUMMYVAL) {
      std::cerr << "Found DUMMYVAL in fInitialFlavours" << std::endl;
      std::cerr << "iNuFlav:" << iNuFlav << std::endl;
      throw;
    }
  }

  if (fNFinalFlavours != (int)fFinalFlavours.size()) {
    std::cerr << "fFinalFlavours array not equal in size to fNFinalFlavours" << std::endl;
    std::cerr << "fNFinalFlavours:" << fNFinalFlavours << std::endl;
    std::cerr << "fFinalFlavours.size():" << fFinalFlavours.size() << std::endl;
    throw;
  }
  for (int iNuFlav=0;iNuFlav<fNFinalFlavours;iNuFlav++) {
    if (fFinalFlavours[iNuFlav]==DUMMYVAL) {
      std::cerr << "Found DUMMYVAL in fFinalFlavours" << std::endl;
      std::cerr << "iNuFlav:" << iNuFlav << std::endl;
      throw;
    }
  }

  if (fVerbose >= INFO) {std::cout << "NeutrinoType and NeutrinoFlavour mapping was found to be OK in Implementation:" << fImplementationName << std::endl;}
  fNuMappingSet = true;
}

void OscProbCalcerBase::SanityCheck() {
  bool IsSane = fEnergyArraySet && fCosineZArraySet && fPropagatorSet && fWeightArrayInit && fNuMappingSet;
  
  if (!IsSane) {
    std::cerr << "OscProbCalcerBase object has been found to not be 'sane' - The following booleans were expected to be true" << std::endl;
    std::cerr << "fEnergyArraySet:" << fEnergyArraySet << std::endl;
    std::cerr << "fCosineZArraySet:" << fCosineZArraySet << std::endl;
    std::cerr << "fPropagatorSet:" << fPropagatorSet << std::endl;
    std::cerr << "fWeightArrayInit:" << fWeightArrayInit << std::endl;
    std::cerr << "fNuMappingSet:" << fNuMappingSet << std::endl;
    throw; 
  } else {
    if (fVerbose >= INFO) {std::cout << "Implementation:" << fImplementationName << " passed SanityCheck" << std::endl;}
  }
}
