#include "OscProbCalcer_NuSQUIDSLinear.h"

#include <iostream>

#include "nuSQuIDS/nuSQuIDS.h"
#include "Decoherence/nuSQUIDSDecoh.h"

OscProbCalcerNuSQUIDSLinear::OscProbCalcerNuSQUIDSLinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config

  //=======

  fNOscParams = kNOscParams;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);

  squids::Const units;
  unsigned int numneu = 3;
  nusquids::nuSQUIDSDecoh nus(nusquids::logspace(1.e2*units.GeV,1.e6*units.GeV,200),numneu,nusquids::neutrino,true);
}

OscProbCalcerNuSQUIDSLinear::~OscProbCalcerNuSQUIDSLinear() {
}

void OscProbCalcerNuSQUIDSLinear::SetupPropagator() {
}

void OscProbCalcerNuSQUIDSLinear::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
}

int OscProbCalcerNuSQUIDSLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerNuSQUIDSLinear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
