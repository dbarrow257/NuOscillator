#include "OscProbCalcer_CHIC.h"

#include <iostream>
#include "CHIC.h"
#include "CHIC_EARTH.h"

OscProbCalcerCHIC::OscProbCalcerCHIC(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config
  //=======
  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  IgnoreCosineZBinning(Config_["General"]["CosineZIgnored"].as<bool>());
  if (!fCosineZIgnored) {
    if (!Config_["OscProbCalcerSetup"]["PremName"]) {
      std::cerr << "Expected to find a 'PremName' Node within the 'OscProbCalcerSetup' Node" << std::endl;
      throw std::runtime_error("YAML node not found");
    }
    fPremName = Config_["OscProbCalcerSetup"]["PremName"].as<std::string>();
    if (!Config_["OscProbCalcerSetup"]["DetDepth"]) {
      std::cerr << "Expected to find a 'DetDepth' Node within the 'OscProbCalcerSetup' Node" << std::endl;
      throw std::runtime_error("YAML node not found");
    }

    fDetDepth = Config_["OscProbCalcerSetup"]["DetDepth"].as<double>();
  }

  std::vector<std::string> OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};
  if (fCosineZIgnored){
    OscParNames.push_back("path_length");
    OscParNames.push_back("matter_density");
  }

  SetExpectedParameterNames(OscParNames);
}

OscProbCalcerCHIC::~OscProbCalcerCHIC() {

}

void OscProbCalcerCHIC::SetupPropagator() {
 if (fCosineZIgnored) {
  // Initialise engine for nu and nubar
  chic_nu    = std::make_unique<CHIC>("neutrino");
  chic_nubar = std::make_unique<CHIC>("antineutrino");
 } else {
   chicearth_nu = std::make_unique<CHICEARTH>(
     "neutrino",
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // irrelevant values, we will override later
     fPremName, fDetDepth);

   chicearth_nubar = std::make_unique<CHICEARTH>(
     "antineutrino",
     0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // irrelevant values, we will override later
     fPremName, fDetDepth );
  }
}


template <typename Propagator>
void OscProbCalcerCHIC::SetPMNSParameters(Propagator* p) {
  // CHIC expects angles, not sin^2
  p->update_th12(std::asin(std::sqrt(GetOscillationParameter(kTH12))));
  p->update_th13(std::asin(std::sqrt(GetOscillationParameter(kTH13))));
  p->update_th23(std::asin(std::sqrt(GetOscillationParameter(kTH23))));

  p->update_dcp(GetOscillationParameter(kDCP));
  p->update_dm221(GetOscillationParameter(kDM12));
  p->update_dm231(GetOscillationParameter(kDM23) + GetOscillationParameter(kDM12));
}


void OscProbCalcerCHIC::CalculateProbabilitiesBeam() {
  for (int iNuType = 0; iNuType < fNNeutrinoTypes; ++iNuType) {
    // KS: CHIC sets Nu and NuBar based on constructor those we switch between both
    CHIC* chic_propagator = (fNeutrinoTypes[iNuType] == Nu) ? chic_nu.get() : chic_nubar.get();
    // ------------------------------- //
    // Set the experimental parameters //
    // ------------------------------- //
    const double Baseline = GetOscillationParameter(ReturnNOscParams() - 2); // km
    const double rho = GetOscillationParameter(ReturnNOscParams() - 1); // g/cc
    chic_propagator->update_density(rho);

    SetPMNSParameters(chic_propagator);

    // KS: Do not multithread, compute_oscillations is mutable as it stores neutrino energy
    // multithreading would break physics
    for (int iOscProb = 0; iOscProb < fNEnergyPoints; ++iOscProb) {
      const double Energy = fEnergyArray[iOscProb];

      Eigen::Matrix3d prob = chic_propagator->compute_oscillations(Energy, Baseline);
      for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; ++iOscChannel) {
        const int from = fOscillationChannels[iOscChannel].GeneratedFlavour - 1;
        const int to   = fOscillationChannels[iOscChannel].DetectedFlavour - 1;

        const double Weight = prob(to, from);

        const int IndexToFill = iNuType * fNOscillationChannels * fNEnergyPoints +
        iOscChannel * fNEnergyPoints + iOscProb;
        fWeightArray[IndexToFill] = Weight;
      }
    }
  }
}

void OscProbCalcerCHIC::CalculateProbabilitiesAtm() {
  for (int iNuType = 0; iNuType < fNNeutrinoTypes; ++iNuType) {
    // KS: CHIC sets Nu and NuBar based on constructor those we switch between both
    CHICEARTH* chic_propagator = (fNeutrinoTypes[iNuType] == Nu) ? chicearth_nu.get() : chicearth_nubar.get();
    SetPMNSParameters(chic_propagator);

    for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {
      for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {
        const double Energy = fEnergyArray[iEnergy];
        const double CosineZ = fCosineZArray[iCosineZ];
        Eigen::Matrix3d prob = chic_propagator->compute_oscillations(Energy, CosineZ);
        for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; ++iOscChannel) {
          const int from = fOscillationChannels[iOscChannel].GeneratedFlavour - 1;
          const int to   = fOscillationChannels[iOscChannel].DetectedFlavour - 1;
          const double Weight = prob(to, from);

          const int IndexToFill = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosineZ);
          fWeightArray[IndexToFill] = Weight;
        }
      }
    }
  }
}

void OscProbCalcerCHIC::CalculateProbabilities() {
  // Oscpars, as given from MaCh3, expresses the mixing angles in sin^2(theta). This propagator expects them in theta
  for (int iOscPar = 0;iOscPar <= kTH13; iOscPar++) {
    if (GetOscillationParameter(iOscPar) < 0) {
      std::cerr << "Invalid oscillation parameter (Can not sqrt this value)!:" << GetOscillationParameter(iOscPar) << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }

  if (fCosineZIgnored) {
    CalculateProbabilitiesBeam();
  } else {
    CalculateProbabilitiesAtm();
  }
}

int OscProbCalcerCHIC::GetNCosineZ() {
  return fCosineZIgnored ? 1 : fNCosineZPoints;
}

int OscProbCalcerCHIC::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = ((NuTypeIndex *  fNOscillationChannels + OscChanIndex) * GetNCosineZ() + std::max(CosineZIndex,0)) * fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerCHIC::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * GetNCosineZ() * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
