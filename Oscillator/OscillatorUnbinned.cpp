#include "Oscillator/OscillatorUnbinned.h"

#include <iostream>

OscillatorUnbinned::OscillatorUnbinned(std::string ConfigName_) : OscillatorBase(ConfigName_) {
  fCalculationTypeName = "Unbinned";
}

OscillatorUnbinned::~OscillatorUnbinned() {

}

const FLOAT_T* OscillatorUnbinned::ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal) {
  const FLOAT_T* Pointer = ReturnPointerToWeightinCalcer(InitNuFlav,FinalNuFlav,EnergyVal,CosineZVal);
  return Pointer;
}
