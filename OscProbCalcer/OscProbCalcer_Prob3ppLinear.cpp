#include "OscProbCalcer_Prob3ppLinear.h"

#include <iostream>

OscProbCalcerProb3ppLinear::OscProbCalcerProb3ppLinear(std::string ConfigName_, int Instance_) : OscProbCalcerBase(ConfigName_,"Prob3ppLinear",Instance_)
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

OscProbCalcerProb3ppLinear::~OscProbCalcerProb3ppLinear() {

  if(bNu != nullptr) delete bNu;
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
    for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {
    
      // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
      int IndexToFill = iNuType*fNOscillationChannels*fNEnergyPoints + iOscChannel*fNEnergyPoints;
      
      for (int iOscProb=0;iOscProb<fNEnergyPoints;iOscProb++) {
        bNu->SetMNS(OscParams[kTH12], OscParams[kTH13], OscParams[kTH23], OscParams[kDM12], OscParams[kDM23], OscParams[kDCP], fEnergyArray[iOscProb], doubled_angle, fNeutrinoTypes[iNuType]);
        bNu->propagateLinear(fNeutrinoTypes[iNuType]*fOscillationChannels[iOscChannel].GeneratedFlavour, OscParams[kPATHL], OscParams[kDENS]);
        fWeightArray[IndexToFill+iOscProb] = bNu->GetProb(fNeutrinoTypes[iNuType]*fOscillationChannels[iOscChannel].GeneratedFlavour, fNeutrinoTypes[iNuType]*fOscillationChannels[iOscChannel].DetectedFlavour);
      }
    }
  }
}

int OscProbCalcerProb3ppLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerProb3ppLinear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
