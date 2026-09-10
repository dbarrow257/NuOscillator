#include "OscProbCalcer_Prob3ppLinear.h"

#include <iostream>

OscProbCalcerProb3ppLinear::OscProbCalcerProb3ppLinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config
  IgnoreCosineZBinning(Config_["General"]["CosineZIgnored"].as<bool>());

  if (!fCosineZIgnored) {
    if (!Config_["OscProbCalcerSetup"]["PREMFile"]) {
      std::cerr << "Expected to find a 'PREMFile' Node within the 'OscProbCalcerSetup' Node" << std::endl;
      throw std::runtime_error("YAML node not found");
    }

    fPremFile = Config_["OscProbCalcerSetup"]["PREMFile"].as<std::string>();
  }

  if (fCosineZIgnored) fImplementationName += "Linear";

  //=======
  std::vector<std::string> OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};
  if (fCosineZIgnored){
    OscParNames.push_back("path_length");
    OscParNames.push_back("matter_density");
  } else {
    OscParNames.push_back("production_height");
  }
  SetExpectedParameterNames(OscParNames);
  
  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  // Implementation specific variables
  doubled_angle = true;

  bNu = nullptr;
}

OscProbCalcerProb3ppLinear::~OscProbCalcerProb3ppLinear() {
  if(bNu != nullptr) delete bNu;
}

void OscProbCalcerProb3ppLinear::SetupPropagator() {
  if (fCosineZIgnored){
    bNu = new BargerPropagator();
  } else {
    bNu = new BargerPropagator(fPremFile.c_str());
  }
  bNu->UseMassEigenstates(false);
  bNu->SetOneMassScaleMode(false);
  bNu->SetWarningSuppression(true);
}

void OscProbCalcerProb3ppLinear::CalculateProbabilities() {
  if (fCosineZIgnored) {
    CalculateProbabilitiesBeam();
  } else {
    CalculateProbabilitiesAtm();
  }
}

void OscProbCalcerProb3ppLinear::CalculateProbabilitiesBeam() {
  const double Baseline = GetOscillationParameter(ReturnNOscParams() - 2); // km
  const double rho = GetOscillationParameter(ReturnNOscParams() - 1); // g/cc

  // Prob3++ calculates oscillation probabilities for each NeutrinoType and each energy, so need to copy them from the calculator into fWeightArray
  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
    for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {

      // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
      int IndexToFill = iNuType*fNOscillationChannels*fNEnergyPoints + iOscChannel*fNEnergyPoints;

      for (int iOscProb=0;iOscProb<fNEnergyPoints;iOscProb++) {
        bNu->SetMNS(GetOscillationParameter(kTH12), GetOscillationParameter(kTH13), GetOscillationParameter(kTH23), GetOscillationParameter(kDM12), GetOscillationParameter(kDM23), GetOscillationParameter(kDCP), fEnergyArray[iOscProb], doubled_angle, fNeutrinoTypes[iNuType]);
        bNu->propagateLinear(fNeutrinoTypes[iNuType]*fOscillationChannels[iOscChannel].GeneratedFlavour, Baseline, rho);
        fWeightArray[IndexToFill+iOscProb] = bNu->GetProb(fNeutrinoTypes[iNuType]*fOscillationChannels[iOscChannel].GeneratedFlavour, fNeutrinoTypes[iNuType]*fOscillationChannels[iOscChannel].DetectedFlavour);
      }
    }
  }
}

void OscProbCalcerProb3ppLinear::CalculateProbabilitiesAtm() {
  const double productionHeight = GetOscillationParameter(ReturnNOscParams() - 1);
  for (int iNuType = 0; iNuType < fNNeutrinoTypes; ++iNuType)
  {
    for (int iCosZ = 0; iCosZ < fNCosineZPoints; ++iCosZ)
    {
      const double cosZ = fCosineZArray[iCosZ];

      bNu->DefinePath(cosZ, productionHeight);
      for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; ++iOscChannel)
      {
        for (int iEnergy = 0; iEnergy < fNEnergyPoints; ++iEnergy)
        {
          const int IndexToFill = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosZ);
          const double energy = fEnergyArray[iEnergy];
          bNu->SetMNS(GetOscillationParameter(kTH12), GetOscillationParameter(kTH13),
                      GetOscillationParameter(kTH23), GetOscillationParameter(kDM12),
                      GetOscillationParameter(kDM23), GetOscillationParameter(kDCP),
                      energy, doubled_angle, fNeutrinoTypes[iNuType]);
          const int generatedFlavour = fNeutrinoTypes[iNuType] * fOscillationChannels[iOscChannel].GeneratedFlavour;

          const int detectedFlavour = fNeutrinoTypes[iNuType] * fOscillationChannels[iOscChannel].DetectedFlavour;
          bNu->propagate(generatedFlavour);
          fWeightArray[IndexToFill] = bNu->GetProb(generatedFlavour, detectedFlavour);
        }
      }
    }
  }
}

int OscProbCalcerProb3ppLinear::GetNCosineZ() {
  return fCosineZIgnored ? 1 : fNCosineZPoints;
}

int OscProbCalcerProb3ppLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = ((NuTypeIndex *  fNOscillationChannels + OscChanIndex) * GetNCosineZ() + std::max(CosineZIndex,0)) * fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerProb3ppLinear::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * GetNCosineZ() * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
