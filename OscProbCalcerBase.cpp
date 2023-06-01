#include "OscProbCalcerBase.h"

#include <math.h>

#include <iostream>
#include <iomanip>

OscProbCalcerBase::OscProbCalcerBase() {
  fEnergyArraySet = false;
  fCosineZArraySet = false;
  fPropagatorSet = false;
  fWeightArrayInit = false;

  

  fNEnergyPoints = -1;
  fNCosineZPoints = -1;
  fEnergyArray = std::vector<FLOAT_T>();
  fCosineZArray = std::vector<FLOAT_T>();

  nWeights = -1;
  fWeightArray = std::vector<FLOAT_T>();

  fNOscParams = -1;
  fOscParamsCurr = std::vector<FLOAT_T>();
}

void OscProbCalcerBase::SetEnergyArray(std::vector<FLOAT_T> EnergyArray) {
  if (fEnergyArraySet) {
    // Already defined the Energy array, or the implementation is designed such not to care about it
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
}

void OscProbCalcerBase::SetCosineZArray(std::vector<FLOAT_T> CosineZArray) {
  if (fCosineZArraySet) {
    // Already defined the CosineZ array, or the implementation is designed such not to care about it
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
}

void OscProbCalcerBase::IgnoreCosineZBinning(bool Ignore) {
  if (Ignore) {
    fCosineZArraySet = true;
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

  return &(fWeightArray[WeightArrayIndex]);
}

void OscProbCalcerBase::Reweight(std::vector<FLOAT_T> OscParams) {
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
}

bool OscProbCalcerBase::AreOscParamsChanged(std::vector<FLOAT_T> OscParamsToCheck) {
  for (int iParam=0;iParam<fNOscParams;iParam++) {
    if (OscParamsToCheck[iParam] != fOscParamsCurr[iParam]) {
      return true;
    }
  }
  return false;
}

void OscProbCalcerBase::ResetCurrOscParams() {
  fOscParamsCurr = std::vector<FLOAT_T>(fNOscParams,DUMMYVAL);
}

void OscProbCalcerBase::SetCurrOscParams(std::vector<FLOAT_T> OscParamsToSave) {
  for (int iParam=0;iParam<fNOscParams;iParam++) {
    fOscParamsCurr[iParam] = OscParamsToSave[iParam];
  }
}

void OscProbCalcerBase::PrintWeights() {
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

  return CosineZIndex;
}

int OscProbCalcerBase::ReturnInitialIndexFromFlavour(int InitFlav) {
  for (size_t iFlav=0;iFlav<InitialFlavours.size();iFlav++) {
    if (fabs(InitFlav) == InitialFlavours[iFlav]) {
      return iFlav;
    }
  }

  std::cerr << "Requested Initial Neutrino flavour is not defined within the InitialFlavour map!" << std::endl;
  std::cerr << "InitNuFlav:" << InitFlav << std::endl;
  std::cerr << "nInitialFlavours:" << nInitialFlavours << std::endl;
  throw;
}

int OscProbCalcerBase::ReturnFinalIndexFromFlavour(int FinalFlav) {
  for (size_t iFlav=0;iFlav<FinalFlavours.size();iFlav++) {
    if (fabs(FinalFlav) == FinalFlavours[iFlav]) {
      return iFlav;
    }
  }

  std::cerr << "Requested Final Neutrino flavour is not defined within the FinalFlavour map!" << std::endl;
  std::cerr << "InitNuFlav:" << FinalFlav << std::endl;
  std::cerr << "nFinalFlavours:" << nFinalFlavours << std::endl;
  throw;
}

int OscProbCalcerBase::ReturnNuTypeFromFlavour(int NuFlav) {
  int NuType = (NuFlav > 0) - (NuFlav < 0); // Calculates the sign of NuFlav
  
  for (int iType=0;iType<nNeutrinoTypes;iType++) {
    if (NuType == NeutrinoTypes[iType]) {
      return iType;
    }
  }

  std::cerr << "Requested Neutrino type is not defined within the NeutrinoType map!" << std::endl;
  std::cerr << "NuFlav:" << NuFlav << std::endl;
  std::cerr << "Associated NuType:" << NuType << std::endl;
  throw;
}

void OscProbCalcerBase::IntialiseWeightArray() {
  nWeights = DefineWeightArraySize();
  std::cout << "Creating weight array with " << nWeights << " entries" << std::endl;
  fWeightArray = std::vector<FLOAT_T>(nWeights,DUMMYVAL);  
  fWeightArrayInit = true;
}

void OscProbCalcerBase::SanityCheck() {
  //DB Check that required variables are sensible

  bool IsSane = fEnergyArraySet && fCosineZArraySet && fPropagatorSet && fWeightArrayInit;
  
  if (!IsSane) {
    std::cerr << "OscProbCalcerBase object has been found to not be 'sane' - The following booleans were expected to be true" << std::endl;
    std::cerr << "fEnergyArraySet:" << fEnergyArraySet << std::endl;
    std::cerr << "fCosineZArraySet:" << fCosineZArraySet << std::endl;
    std::cerr << "fPropagatorSet:" << fPropagatorSet << std::endl;
    std::cerr << "fWeightArrayInit:" << fWeightArrayInit << std::endl;
    throw; 
  } else {
    std::cout << "OscProbCalcerBase object has been found to be 'sane'" << std::endl;
  }
}
