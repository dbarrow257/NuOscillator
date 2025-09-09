#include "OscProbCalcer_OscLib.h"

#include "OscLib/OscCalcPMNS.h"

OscProbCalcerOscLib::OscProbCalcerOscLib(YAML::Node Config_) : OscProbCalcerBase(Config_) {
  fNOscParams = 6;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;
  
  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);
  
  osc::_OscCalcPMNS<FLOAT_T>* OscLib = new osc::_OscCalcPMNS<FLOAT_T>();
}

OscProbCalcerOscLib::~OscProbCalcerOscLib() {
}

void OscProbCalcerOscLib::SetupPropagator() {
}

void OscProbCalcerOscLib::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
}

int OscProbCalcerOscLib::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = 1;
  return IndexToReturn;
}

long OscProbCalcerOscLib::DefineWeightArraySize() {
  int nCalculationPoints = 100;
  return nCalculationPoints;
}
