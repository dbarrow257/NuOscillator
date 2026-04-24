#include "OscProbCalcer_Prob3ppLinear.h"

#include <iostream>

OscProbCalcerProb3ppLinear::OscProbCalcerProb3ppLinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config

  //=======
  std::vector<std::string> OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp","path_length","matter_density"};
  SetExpectedParameterNames(OscParNames);
  
  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);

  // Implementation specific variables
  doubled_angle = true;

  bNu = nullptr;
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

void OscProbCalcerProb3ppLinear::CalculateProbabilities() {
  // Prob3++ calculates oscillation probabilities for each NeutrinoType and each energy, so need to copy them from the calculator into fWeightArray
  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
    for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {
    
      // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
      int IndexToFill = iNuType*fNOscillationChannels*fNEnergyPoints + iOscChannel*fNEnergyPoints;
      
      for (int iOscProb=0;iOscProb<fNEnergyPoints;iOscProb++) {
        bNu->SetMNS(GetOscillationParameter(kTH12), GetOscillationParameter(kTH13), GetOscillationParameter(kTH23), GetOscillationParameter(kDM12), GetOscillationParameter(kDM23), GetOscillationParameter(kDCP), fEnergyArray[iOscProb], doubled_angle, fNeutrinoTypes[iNuType]);
        bNu->propagateLinear(fNeutrinoTypes[iNuType]*fOscillationChannels[iOscChannel].GeneratedFlavour, GetOscillationParameter(kPATHL), GetOscillationParameter(kDENS));
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
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
