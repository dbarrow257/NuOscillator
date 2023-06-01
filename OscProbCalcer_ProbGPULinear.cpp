#include "OscProbCalcer_ProbGPULinear.h"

extern "C" void setMNS(double x12, double x13, double x23, double m21, double m23, double Delta, bool kSquared);
extern "C" void GetProb(int Alpha, int Beta, double Path, double Density, double *Energy, int n, double *oscw);

#include <iostream>

OscProbCalcerProbGPULinear::OscProbCalcerProbGPULinear() : OscProbCalcerBase()
{
  // Required variables
  fNOscParams = kNOscParams;

  nNeutrinoTypes = 2;
  NeutrinoTypes.resize(nNeutrinoTypes);
  NeutrinoTypes[0] = Nu;
  NeutrinoTypes[1] = Nubar;

  nInitialFlavours = 2;
  InitialFlavours.resize(nInitialFlavours);
  InitialFlavours[0] = Electron;
  InitialFlavours[1] = Muon;

  nFinalFlavours = 2;
  FinalFlavours.resize(nFinalFlavours);
  FinalFlavours[0] = Electron;
  FinalFlavours[1] = Muon;

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);

  // Implementation specific variables
  doubled_angle = true;
}

void OscProbCalcerProbGPULinear::SetupPropagator() {
  // This implementation doesn't really need to do anything in the setup due to probGPU's horrific implementation
}

void OscProbCalcerProbGPULinear::CalculateProbabilities(std::vector<FLOAT_T> OscParams) {
  setMNS(OscParams[kTH12], OscParams[kTH23], OscParams[kTH13], OscParams[kDM12], OscParams[kDM23], OscParams[kDCP], doubled_angle);

  // ProbGPULinear calculates oscillation probabilites for each NeutrinoType, so need to copy them from the calculator into fWeightArray
  int CopyArrSize = fNEnergyPoints;
  FLOAT_T* CopyArr = new FLOAT_T[CopyArrSize];

  for (int iNuType=0;iNuType<nNeutrinoTypes;iNuType++) {
    for (int iInitFlav=0;iInitFlav<nInitialFlavours;iInitFlav++) {
      for (int iFinalFlav=0;iFinalFlav<nFinalFlavours;iFinalFlav++) {
	GetProb(NeutrinoTypes[iNuType]*iInitFlav, NeutrinoTypes[iNuType]*iFinalFlav, OscParams[kPATHL], OscParams[kDENS], fEnergyArray.data(), fNEnergyPoints, CopyArr);

        // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
        int IndexToFill = iNuType*nInitialFlavours*nFinalFlavours*CopyArrSize + iInitFlav*nFinalFlavours*CopyArrSize + iFinalFlav*CopyArrSize;
        for (int iOscProb=0;iOscProb<CopyArrSize;iOscProb++) {
          fWeightArray[IndexToFill+iOscProb] = CopyArr[iOscProb];
        }
      }
    }
  }
}

int OscProbCalcerProbGPULinear::ReturnWeightArrayIndex(int NuTypeIndex, int InitNuIndex, int FinalNuIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*nInitialFlavours*nFinalFlavours*fNEnergyPoints + InitNuIndex*nFinalFlavours*fNEnergyPoints + FinalNuIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerProbGPULinear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * nInitialFlavours * nFinalFlavours * nNeutrinoTypes;
  return nCalculationPoints;
}
