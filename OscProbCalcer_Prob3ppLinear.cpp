#include "OscProbCalcer_Prob3ppLinear.h"

#include <iostream>

OscProbCalcerProb3ppLinear::OscProbCalcerProb3ppLinear() : OscProbCalcerBase()
{
  //Base variables
  fNOscParams = kNOscParams;

  //This implementation only considers linear propagation, thus no requirement to set cosineZ array
  fCosineZArraySet = true;

  //Implementation specific variables
  doubled_angle = true;

  nNeutrinoSigns = kNNeutrinoTypes;
  nInitialFlavours = 2; // =2 if excluding taus, =3 if including taus
  nFinalFlavours = 3; // =2 if excluding taus, =3 if including taus

  NeutrinoTypes.resize(nNeutrinoSigns);
  NeutrinoTypes[0] = Neutrino;
  NeutrinoTypes[1] = AntiNeutrino;
}

void OscProbCalcerProb3ppLinear::SetupPropagator() {
   bNu = new BargerPropagator();
   bNu->UseMassEigenstates(false);
   bNu->SetOneMassScaleMode(false);
   bNu->SetWarningSuppression(true); 

   fPropagatorSet = true;
}

void OscProbCalcerProb3ppLinear::CalculateProbabilities(std::vector<FLOAT_T> OscParams) {
  // Prob3++ calculates oscillation probabilites for each NeutrinoType and each energy, so need to copy them from the calculator into fWeightArray

  for (int iNuType=0;iNuType<nNeutrinoSigns;iNuType++) {
    for (int iInitFlav=0;iInitFlav<nInitialFlavours;iInitFlav++) {
      for (int iFinalFlav=0;iFinalFlav<nFinalFlavours;iFinalFlav++) {

        // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
        int IndexToFill = iNuType*nInitialFlavours*nFinalFlavours*fNEnergyPoints + iInitFlav*nFinalFlavours*fNEnergyPoints + iFinalFlav*fNEnergyPoints;

        for (int iOscProb=0;iOscProb<fNEnergyPoints;iOscProb++) {
	  bNu->SetMNS(OscParams[kTH12], OscParams[kTH23], OscParams[kTH13], OscParams[kDM12], OscParams[kDM23], OscParams[kDCP], fEnergyArray[iOscProb], doubled_angle);
	  bNu->propagateLinear(NeutrinoTypes[iNuType]*iInitFlav, OscParams[kPATHL], OscParams[kDENS]);
          fWeightArray[IndexToFill+iOscProb] = bNu->GetProb(NeutrinoTypes[iNuType]*iInitFlav, NeutrinoTypes[iNuType]*iFinalFlav);
        }
      }
    }
  }

}

const FLOAT_T* OscProbCalcerProb3ppLinear::ReturnPointer(FLOAT_T Energy, FLOAT_T CosineZ) {
  return NULL;
}

void OscProbCalcerProb3ppLinear::IntiailiseWeightArray() {
  int nCalculationPoints = fNEnergyPoints * nInitialFlavours * nFinalFlavours * nNeutrinoSigns;
  std::cout << "Creating weight array with " << nCalculationPoints << " entries" << std::endl;

  fWeightArray = std::vector<FLOAT_T>(nCalculationPoints,DUMMYVAL);
  fWeightArrayInit = true;
}
