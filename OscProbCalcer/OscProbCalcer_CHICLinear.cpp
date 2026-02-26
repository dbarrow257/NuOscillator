#include "OscProbCalcer_CHICLinear.h"

#include <iostream>
#include "CHIC.h"

OscProbCalcerCHICLinear::OscProbCalcerCHICLinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
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
}

OscProbCalcerCHICLinear::~OscProbCalcerCHICLinear() {

}

void OscProbCalcerCHICLinear::SetupPropagator() {
  // Initialise engine for nu and nubar
  chic_nu    = std::make_unique<CHIC>("neutrino");
  chic_nubar = std::make_unique<CHIC>("antineutrino");
}

void OscProbCalcerCHICLinear::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
  // Oscpars, as given from MaCh3, expresses the mixing angles in sin^2(theta). This propagator expects them in theta
  for (int iOscPar = 0;iOscPar <= kTH13; iOscPar++) {
    if (OscParams[iOscPar] < 0) {
      std::cerr << "Invalid oscillation parameter (Can not sqrt this value)!:" << OscParams[iOscPar] << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }

  for (int iNuType = 0; iNuType < fNNeutrinoTypes; ++iNuType) {
    // KS: CHIC sets Nu and NuBar based on constructor those we switch between both
    CHIC* chic_propagator = (fNeutrinoTypes[iNuType] == Nu) ? chic_nu.get() : chic_nubar.get();
    // ------------------------------- //
    // Set the experimental parameters //
    // ------------------------------- //
    const double Baseline = OscParams[kPATHL]; // km
    const double rho = OscParams[kDENS]; // g/cc

    // CHIC expects angles, not sin^2
    chic_propagator->update_th12(std::asin(std::sqrt(OscParams[kTH12])));
    chic_propagator->update_th13(std::asin(std::sqrt(OscParams[kTH13])));
    chic_propagator->update_th23(std::asin(std::sqrt(OscParams[kTH23])));

    chic_propagator->update_dcp(OscParams[kDCP]);
    chic_propagator->update_dm221(OscParams[kDM12]);
    chic_propagator->update_dm231(OscParams[kDM23] + OscParams[kDM12]);
    chic_propagator->update_density(rho);

    // KS: Do not multithread, compute_oscillations is mutable as it stores neutrino energy
    // multihreading would break physics
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

int OscProbCalcerCHICLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerCHICLinear::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
