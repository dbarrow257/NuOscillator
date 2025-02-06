#include "Oscillator/OscillatorUnbinned.h"

#include <iostream>

OscillatorUnbinned::OscillatorUnbinned(std::string ConfigName_) : OscillatorBase(ConfigName_) {
  Initialise();
}

OscillatorUnbinned::OscillatorUnbinned(YAML::Node Config_) : OscillatorBase(Config_) {
  Initialise();
}

void OscillatorUnbinned::Initialise() {
  fCalculationTypeName = "Unbinned";
}

OscillatorUnbinned::~OscillatorUnbinned() {

}

const FLOAT_T* OscillatorUnbinned::ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal) {
  const FLOAT_T* Pointer = ReturnPointerToWeightinCalcer(InitNuFlav,FinalNuFlav,EnergyVal,CosineZVal);
  return Pointer;
}

std::vector<FLOAT_T> OscillatorUnbinned::ReturnBinEdgesForPlotting(bool ReturnEnergy) {
  std::vector<FLOAT_T> BinEdges;

  std::vector<FLOAT_T> EvalPoints;
  if (ReturnEnergy) {
    EvalPoints = fOscProbCalcer->ReturnEnergyArray();
  } else {
    EvalPoints = fOscProbCalcer->ReturnCosineZArray();
  }

  BinEdges.resize(EvalPoints.size()+1);
  for (int i=0;i<(EvalPoints.size()-1);i++) {
    BinEdges[i] = EvalPoints[i]-(EvalPoints[i+1]-EvalPoints[i])/2.0;
  }
  BinEdges[EvalPoints.size()-1] = EvalPoints[EvalPoints.size()-1]-(EvalPoints[EvalPoints.size()-1]-EvalPoints[EvalPoints.size()-2])/2.0;
  BinEdges[EvalPoints.size()] = EvalPoints[EvalPoints.size()-1]+(EvalPoints[EvalPoints.size()-1]-EvalPoints[EvalPoints.size()-2])/2.0;
  
  return BinEdges;
}
