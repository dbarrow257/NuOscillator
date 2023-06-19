#include "OscillatorUnbinned.h"

#include <iostream>

OscillatorUnbinned::OscillatorUnbinned(std::string ConfigName_) : OscillatorBase(ConfigName_) {
  fCalculationTypeName = "Unbinned";

  //=======
  // Grab the following from config manager - Currently brought through via constructor
  //DB
  //fOscProbCalcerImplementationToCreate = OscProbCalcerImplementationToCreate_;
  //=======

  InitialiseOscProbCalcers();
}

const FLOAT_T* OscillatorUnbinned::ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal) {
  int CalcerIndex = 0;
  return ReturnPointerToWeightinCalcer(CalcerIndex,InitNuFlav,FinalNuFlav,EnergyVal,CosineZVal);
}

void OscillatorUnbinned::SetEnergyArray(std::vector<FLOAT_T> Array) {
  SetEnergyArrayInCalcer(Array);
}

void OscillatorUnbinned::SetCosineZArray(std::vector<FLOAT_T> Array) {
  SetCosineZArrayInCalcer(Array);
}
