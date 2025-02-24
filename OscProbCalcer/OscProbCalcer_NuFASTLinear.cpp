#include "OscProbCalcer_NuFASTLinear.h"

#include <iostream>

#include "c++/NuFast.cpp"

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

  fNOscParams = kNOscParams;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);
}

OscProbCalcerNuFASTLinear::~OscProbCalcerNuFASTLinear() {

}

void OscProbCalcerNuFASTLinear::SetupPropagator() {
}

void OscProbCalcerNuFASTLinear::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
  double L, E, rho, Ye, probs_returned[3][3];
  double s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31;

  // ------------------------------- //
  // Set the experimental parameters //
  // ------------------------------- //
  L = OscParams[kPATHL]; // km
  rho = OscParams[kDENS]; // g/cc
  Ye = OscParams[kELECDENS];
  
  // ------------------------------------- //
  // Set the vacuum oscillation parameters //
  // ------------------------------------- //
  s12sq = OscParams[kTH12];
  s13sq = OscParams[kTH13];
  s23sq = OscParams[kTH23];
  delta = OscParams[kDCP];
  Dmsq21 = OscParams[kDM12];

  //Need to convert OscParams[kDM23] to kDM31
  Dmsq31 = OscParams[kDM23]+OscParams[kDM12]; // eV^2
  
  // ------------------------------------------ //
  // Calculate all 9 oscillationa probabilities //
  // ------------------------------------------ //

  for (int iOscProb=0;iOscProb<fNEnergyPoints;iOscProb++) {
    for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
      
      //+ve energy for neutrinos, -ve energy for antineutrinos
      E = fEnergyArray[iOscProb] * fNeutrinoTypes[iNuType];

      Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, &probs_returned);

      for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {
	// Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
	int IndexToFill = iNuType*fNOscillationChannels*fNEnergyPoints + iOscChannel*fNEnergyPoints;

	double Weight = probs_returned[fOscillationChannels[iOscChannel].GeneratedFlavour-1][fOscillationChannels[iOscChannel].DetectedFlavour-1];
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
