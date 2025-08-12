#include "OscProbCalcer_NuFASTLinear.h"

#include <iostream>

#include "c++/NuFast.cpp"

#if UseMultithreading == 1
#include "omp.h"
#endif

OscProbCalcerNuFASTLinear::OscProbCalcerNuFASTLinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config

  // --------------------------------------------------------------------- //
  // Set the number of Newton-Raphson iterations which sets the precision. //
  // 0 is close to the single precision limit and is better than DUNE/HK   //
  // in the high statistics regime. Increasig N_Newton to 1,2,... rapidly  //
  // improves the precision at a modest computational cost                 //
  // --------------------------------------------------------------------- //
  N_Newton = 3;
  if (Config_["OscProbCalcerSetup"]["nNewtonIter"]) {
    N_Newton = Config_["OscProbCalcerSetup"]["nNewtonIters"].as<int>();
  }
  
  //=======
  std::vector<std::string> OscParNames = {"sin2_th12","sin2_th13","sin2_th23","dm2_12","dm2_23","delta_cp","path_length","matter_density","electron_density"};
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

OscProbCalcerNuFASTLinear::~OscProbCalcerNuFASTLinear() {

}

void OscProbCalcerNuFASTLinear::SetupPropagator() {
}

void OscProbCalcerNuFASTLinear::CalculateProbabilities() {
  // ------------------------------- //
  // Set the experimental parameters //
  // ------------------------------- //
  const double L = *fOscParams[kPATHL]; // km
  const double rho = *fOscParams[kDENS]; // g/cc
  const double Ye = *fOscParams[kELECDENS];
  
  // ------------------------------------- //
  // Set the vacuum oscillation parameters //
  // ------------------------------------- //
  const double s12sq = *fOscParams[kTH12];
  const double s13sq = *fOscParams[kTH13];
  const double s23sq = *fOscParams[kTH23];
  const double delta = *fOscParams[kDCP];
  const double Dmsq21 = *fOscParams[kDM12];

  //Need to convert fOscParams[kDM23] to kDM31
  const double Dmsq31 = *fOscParams[kDM23]+*fOscParams[kDM12]; // eV^2
  
  double probs_returned[3][3];
  // ------------------------------------------ //
  // Calculate all 9 oscillations probabilities //
  // ------------------------------------------ //
  #if UseMultithreading == 1
  #pragma omp parallel for collapse(2) private(probs_returned)
  #endif
  for (int iOscProb=0;iOscProb<fNEnergyPoints;iOscProb++) {
    for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
      //+ve energy for neutrinos, -ve energy for antineutrinos
      const double E = fEnergyArray[iOscProb] * fNeutrinoTypes[iNuType];

      Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, &probs_returned);

      #if UseMultithreading == 1
      #pragma omp simd
      #endif
      for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {
        // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
        const int IndexToFill = iNuType*fNOscillationChannels*fNEnergyPoints + iOscChannel*fNEnergyPoints;

        const double Weight = probs_returned[fOscillationChannels[iOscChannel].GeneratedFlavour-1][fOscillationChannels[iOscChannel].DetectedFlavour-1];
        fWeightArray[IndexToFill+iOscProb] = Weight;
      }
    }
  }

}

int OscProbCalcerNuFASTLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerNuFASTLinear::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
