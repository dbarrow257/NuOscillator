#include "OscillatorUnbinned.h"

#include <iostream>

OscillatorUnbinned::OscillatorUnbinned(std::vector<std::string> OscProbCalcerImplementationToCreate_, int Verbosity_, bool CosineZIgnored_) : OscillatorBase() {
  fCalculationTypeName = "Unbinned";

  //=======
  //DB Grab the following from config manager - Currently brought through via constructor
  fOscProbCalcerImplementationToCreate = OscProbCalcerImplementationToCreate_;
  fVerbose = Verbosity_;
  fCosineZIgnored = CosineZIgnored_;
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
