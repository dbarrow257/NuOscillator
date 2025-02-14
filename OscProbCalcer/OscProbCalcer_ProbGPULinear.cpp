#include "OscProbCalcer_ProbGPULinear.h"

extern "C" void setMNS(double x12, double x13, double x23, double m21, double m23, double Delta, bool kSquared);
extern "C" void GetProb(int Alpha, int Beta, double Path, double Density, double *Energy, int n, double *oscw);

#include <iostream>

OscProbCalcerProbGPULinear::OscProbCalcerProbGPULinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
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

  // Implementation specific variables
  doubled_angle = true;
}

OscProbCalcerProbGPULinear::~OscProbCalcerProbGPULinear() {

}

void OscProbCalcerProbGPULinear::SetupPropagator() {
  // This implementation doesn't really need to do anything in the setup due to probGPU's horrific implementation
}

void OscProbCalcerProbGPULinear::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
  setMNS(OscParams[kTH12], OscParams[kTH13], OscParams[kTH23], OscParams[kDM12], OscParams[kDM23], OscParams[kDCP], doubled_angle);

  // ProbGPULinear calculates oscillation probabilities for each NeutrinoType, so need to copy them from the calculator into fWeightArray
  int CopyArrSize = fNEnergyPoints;
  FLOAT_T* CopyArr = new FLOAT_T[CopyArrSize];

  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
    for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {
      GetProb(fNeutrinoTypes[iNuType]*fOscillationChannels[iOscChannel].GeneratedFlavour, fNeutrinoTypes[iNuType]*fOscillationChannels[iOscChannel].DetectedFlavour, OscParams[kPATHL], OscParams[kDENS], fEnergyArray.data(), fNEnergyPoints, CopyArr);
      
      // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
      int IndexToFill = iNuType*fNOscillationChannels*CopyArrSize + iOscChannel*CopyArrSize;
      for (int iOscProb=0;iOscProb<CopyArrSize;iOscProb++) {
        fWeightArray[IndexToFill+iOscProb] = CopyArr[iOscProb];
      }
    }
  }
  delete[] CopyArr;
}

int OscProbCalcerProbGPULinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerProbGPULinear::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
