#include "OscProbCalcer_ProbGPULinear.h"

extern "C" void setMNS(double x12, double x13, double x23, double m21, double m23, double Delta, bool kSquared);
extern "C" void GetProb(int Alpha, int Beta, double Path, double Density, double *Energy, int n, double *oscw);

#include <iostream>

OscProbCalcerProbGPULinear::OscProbCalcerProbGPULinear(std::string ConfigName_) : OscProbCalcerBase(ConfigName_)
{
  // Required variables
  fImplementationName = "ProbGPU";

  fNOscParams = kNOscParams;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  fNInitialFlavours = 3;
  InitialiseInitialFlavoursArray(fNInitialFlavours);
  fInitialFlavours[0] = Electron;
  fInitialFlavours[1] = Muon;
  fInitialFlavours[2] = Tau;

  fNFinalFlavours = 3;
  InitialiseFinalFlavoursArray(fNFinalFlavours);
  fFinalFlavours[0] = Electron;
  fFinalFlavours[1] = Muon;
  fFinalFlavours[2] = Tau;

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);

  // Implementation specific variables
  doubled_angle = true;
}

void OscProbCalcerProbGPULinear::SetupPropagator() {
  // This implementation doesn't really need to do anything in the setup due to probGPU's horrific implementation
}

void OscProbCalcerProbGPULinear::CalculateProbabilities(std::vector<FLOAT_T> OscParams) {
  setMNS(OscParams[kTH12], OscParams[kTH13], OscParams[kTH23], OscParams[kDM12], OscParams[kDM23], OscParams[kDCP], doubled_angle);

  // ProbGPULinear calculates oscillation probabilites for each NeutrinoType, so need to copy them from the calculator into fWeightArray
  int CopyArrSize = fNEnergyPoints;
  FLOAT_T* CopyArr = new FLOAT_T[CopyArrSize];

  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
    for (int iInitFlav=0;iInitFlav<fNInitialFlavours;iInitFlav++) {
      for (int iFinalFlav=0;iFinalFlav<fNFinalFlavours;iFinalFlav++) {
	GetProb(fNeutrinoTypes[iNuType]*fInitialFlavours[iInitFlav], fNeutrinoTypes[iNuType]*fFinalFlavours[iFinalFlav], OscParams[kPATHL], OscParams[kDENS], fEnergyArray.data(), fNEnergyPoints, CopyArr);

        // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
        int IndexToFill = iNuType*fNInitialFlavours*fNFinalFlavours*CopyArrSize + iInitFlav*fNFinalFlavours*CopyArrSize + iFinalFlav*CopyArrSize;
        for (int iOscProb=0;iOscProb<CopyArrSize;iOscProb++) {
          fWeightArray[IndexToFill+iOscProb] = CopyArr[iOscProb];
        }
      }
    }
  }
}

int OscProbCalcerProbGPULinear::ReturnWeightArrayIndex(int NuTypeIndex, int InitNuIndex, int FinalNuIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNInitialFlavours*fNFinalFlavours*fNEnergyPoints + InitNuIndex*fNFinalFlavours*fNEnergyPoints + FinalNuIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerProbGPULinear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNInitialFlavours * fNFinalFlavours * fNNeutrinoTypes;
  return nCalculationPoints;
}
