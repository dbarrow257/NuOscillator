#include "OscProbCalcer_Prob3ppLinear.h"

#include <iostream>

OscProbCalcerProb3ppLinear::OscProbCalcerProb3ppLinear() : OscProbCalcerBase()
{
  // Required variables
  fVerbose = INFO;
  fImplementationName = "Prob3pp";

  fNOscParams = kNOscParams;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  fNInitialFlavours = 2;
  InitialiseInitialFlavoursArray(fNInitialFlavours);
  fInitialFlavours[0] = Electron;
  fInitialFlavours[1] = Muon;

  fNFinalFlavours = 2;
  InitialiseFinalFlavoursArray(fNFinalFlavours);
  fFinalFlavours[0] = Electron;
  fFinalFlavours[1] = Muon;

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);

  // Implementation specific variables
  doubled_angle = true;
}

void OscProbCalcerProb3ppLinear::SetupPropagator() {
   bNu = new BargerPropagator();
   bNu->UseMassEigenstates(false);
   bNu->SetOneMassScaleMode(false);
   bNu->SetWarningSuppression(true);
}

void OscProbCalcerProb3ppLinear::CalculateProbabilities(std::vector<FLOAT_T> OscParams) {
  // Prob3++ calculates oscillation probabilites for each NeutrinoType and each energy, so need to copy them from the calculator into fWeightArray
  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
    for (int iInitFlav=0;iInitFlav<fNInitialFlavours;iInitFlav++) {
      for (int iFinalFlav=0;iFinalFlav<fNFinalFlavours;iFinalFlav++) {

        // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
        int IndexToFill = iNuType*fNInitialFlavours*fNFinalFlavours*fNEnergyPoints + iInitFlav*fNFinalFlavours*fNEnergyPoints + iFinalFlav*fNEnergyPoints;

        for (int iOscProb=0;iOscProb<fNEnergyPoints;iOscProb++) {
	  bNu->SetMNS(OscParams[kTH12], OscParams[kTH13], OscParams[kTH23], OscParams[kDM12], OscParams[kDM23], OscParams[kDCP], fEnergyArray[iOscProb], doubled_angle, fNeutrinoTypes[iNuType]);
	  bNu->propagateLinear(fNeutrinoTypes[iNuType]*fInitialFlavours[iInitFlav], OscParams[kPATHL], OscParams[kDENS]);
          fWeightArray[IndexToFill+iOscProb] = bNu->GetProb(fNeutrinoTypes[iNuType]*fInitialFlavours[iInitFlav], fNeutrinoTypes[iNuType]*fFinalFlavours[iFinalFlav]);
        }
      }
    }
  }

}

int OscProbCalcerProb3ppLinear::ReturnWeightArrayIndex(int NuTypeIndex, int InitNuIndex, int FinalNuIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNInitialFlavours*fNFinalFlavours*fNEnergyPoints + InitNuIndex*fNFinalFlavours*fNEnergyPoints + FinalNuIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerProb3ppLinear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNInitialFlavours * fNFinalFlavours * fNNeutrinoTypes;
  return nCalculationPoints;
}
