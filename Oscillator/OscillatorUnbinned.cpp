#include "OscillatorUnbinned.h"

#include <iostream>

OscillatorUnbinned::OscillatorUnbinned(std::string ConfigName_) : OscillatorBase(ConfigName_) {
  fCalculationTypeName = "Unbinned";
}

const FLOAT_T* OscillatorUnbinned::ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal) {
  int CalcerIndex = 0;
  const FLOAT_T* Pointer = ReturnPointerToWeightinCalcer(CalcerIndex,InitNuFlav,FinalNuFlav,EnergyVal,CosineZVal);
  return Pointer;
}
