#include "OscProbCalcer_nuOscillationLinear.h"
#include <iostream>

#include "header/vOscillation.hh"

#if UseMultithreading == 1
#include "omp.h"
#endif

OscProbCalcernuOscillationLinear::OscProbCalcernuOscillationLinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  fOscillator = std::make_unique<vOscillation>();
  //=======
  std::vector<std::string> OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp","path_length"};
  SetExpectedParameterNames(OscParNames);
  
  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  #if UseMultithreading == 1
  fImplementationName += "-CPU-"+std::to_string(omp_get_max_threads());
  #else
  fImplementationName += "-CPU-"+std::to_string(1);
  #endif

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);
}

OscProbCalcernuOscillationLinear::~OscProbCalcernuOscillationLinear() {
}

void OscProbCalcernuOscillationLinear::SetupPropagator() {
}

TString FlavourToVOscillationString(int Flavour) {
  switch (Flavour) {
    case NuOscillator::kElectron:
      return "e";
    case NuOscillator::kMuon:
      return "mu";
    case NuOscillator::kTau:
      return "tau";
    default:
      throw std::invalid_argument("vOscillation supports only Electron, Muon and Tau flavours. "
                                  "Received flavour index: " + std::to_string(Flavour));
  }
}

void OscProbCalcernuOscillationLinear::CalculateProbabilities() {
  // Mixing parameters
  const double s12sq = GetOscillationParameter(kTH12);
  const double s23sq = GetOscillationParameter(kTH23);
  const double s13sq = GetOscillationParameter(kTH13);

  // Mass-squared differences
  const double dm21 = GetOscillationParameter(kDM12);
  const double dm32 = GetOscillationParameter(kDM23);
  const double dm31 = dm21 + dm32;

  // CP phase
  const double delta = GetOscillationParameter(kDCP);
  const double delta_degrees = delta * 180.0 / TMath::Pi();
  // Configure vOscillation
  fOscillator->Setsin(std::sqrt(s12sq), std::sqrt(s13sq), std::sqrt(s23sq));
  fOscillator->SetDm2(dm21, dm31);
  fOscillator->Setdelta_CP(delta_degrees);

  // Experimental baseline, need to covnert KM to M
  const double L_m = 1000*GetOscillationParameter(kPATHL);
  #if UseMultithreading == 1
  #pragma omp parallel for collapse(2)
  #endif
  for (int iNuType = 0; iNuType < fNNeutrinoTypes; ++iNuType) {
    for (int iEnergy = 0; iEnergy < fNEnergyPoints; ++iEnergy) {
      const bool anti = (fNeutrinoTypes[iNuType] == Nubar);
      // vOscillation uses MeV
      const double E_MeV = fEnergyArray[iEnergy] * 1000.0;

      for (int iChannel = 0; iChannel < fNOscillationChannels; ++iChannel) {
        const int InitialFlavour = fOscillationChannels[iChannel].GeneratedFlavour;
        const int FinalFlavour = fOscillationChannels[iChannel].DetectedFlavour;

        // Convert flavour integers to: 1 = e, 2 = mu, 3 = tau
        // Passing strings isn't super smart...
        const TString Initial = FlavourToVOscillationString(InitialFlavour);
        const TString Final = FlavourToVOscillationString(FinalFlavour);

        const double Probability = fOscillator->GetProbability(L_m, E_MeV, Initial, Final, anti);
        const int Index = ReturnWeightArrayIndex(iNuType, iChannel, iEnergy, -1);
        fWeightArray[Index] = Probability;
      }
    }
  }
}

int OscProbCalcernuOscillationLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcernuOscillationLinear::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
