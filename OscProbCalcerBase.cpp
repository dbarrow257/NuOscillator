#include "OscProbCalcerBase.h"

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

  fNOscParams = -1;
  fOscParamsCurr = std::vector<FLOAT_T>();
}

void OscProbCalcerBase::SetEnergyArray(std::vector<FLOAT_T> EnergyArray) {
  if (fEnergyArraySet) {
    // Already defined the Energy array, or the implementation is designed such not to care about it
    return;
  }

  fEnergyArray = EnergyArray;
  fNEnergyPoints = EnergyArray.size();

  //DB Could probably do some manipulation of the vector to see if they are meaningful values (i.e. positive)

  fEnergyArraySet = true;
}

void OscProbCalcerBase::SetCosineZArray(std::vector<FLOAT_T> CosineZArray) {
  if (fCosineZArraySet) {
    // Already defined the CosineZ array, or the implementation is designed such not to care about it
    return;
  }

  fCosineZArray = CosineZArray;
  fNCosineZPoints = CosineZArray.size();

  //DB Could probably do some manipulation of the vector to see if they are meaningful values (i.e. between -1. and 1.)

  fCosineZArraySet = true;
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
  IntiailiseWeightArray();
  SetupPropagator();
  SanityCheck();
}

const FLOAT_T* OscProbCalcerBase::ReturnPointerToWeight(FLOAT_T Energy, FLOAT_T CosineZ) {
  return ReturnPointer(Energy,CosineZ);
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

void OscProbCalcerBase::SanityCheck() {
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
