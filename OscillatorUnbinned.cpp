#include "OscillatorUnbinned.h"

#include <iostream>

OscillatorUnbinned::OscillatorUnbinned(std::vector<std::string> OscProbCalcerImplementationToCreate, int Verbosity_, bool CosineZIgnored_) : OscillatorBase(OscProbCalcerImplementationToCreate) {
  //=======
  //DB Grab the following from config manager - Currently brought through via constructor
  fVerbose = Verbosity_;
  fCosineZIgnored = CosineZIgnored_;
  //=======
}

const FLOAT_T* OscillatorUnbinned::ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal) {
  int CalcerIndex = 0;
  if (CalcerIndex < 0 || CalcerIndex >= fNCalcers) {
    std::cerr << "Requested to Return NOscParams at invalid index within OPCalcers array" << std::endl;
    std::cerr << "CalcerIndex:"<< CalcerIndex << std::endl;
    std::cerr << "fNCalcers:" << fNCalcers << std::endl;
    throw;
  }

  return OPCalcers[CalcerIndex]->ReturnPointerToWeight(InitNuFlav,FinalNuFlav,EnergyVal,CosineZVal);
}

void OscillatorUnbinned::SetEnergyArray(std::vector<FLOAT_T> Array) {
  SetEnergyArrayInCalcer(Array);
}

void OscillatorUnbinned::SetCosineZArray(std::vector<FLOAT_T> Array) {
  SetCosineZArrayInCalcer(Array);
}
