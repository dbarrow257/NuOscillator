#include "OscProbCalcer_NuFASTLinear.h"

#include <iostream>

#include "c++/NuFast.cpp"

OscProbCalcerNuFASTLinear::OscProbCalcerNuFASTLinear(std::string ConfigName_, int Instance_) : OscProbCalcerBase(ConfigName_,"NuFASTLinear",Instance_)
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

void OscProbCalcerNuFASTLinear::SetupPropagator() {
}

void OscProbCalcerNuFASTLinear::CalculateProbabilities(std::vector<FLOAT_T> OscParams) {
  double L, E, rho, Ye, probs_returned[3][3];
  double s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31;
  int N_Newton;

  // ------------------------------- //
  // Set the experimental parameters //
  // ------------------------------- //
  L = 1300; // km
  E = 2.5; // GeV
  rho = 3; // g/cc
  Ye = 0.5;
  
  // --------------------------------------------------------------------- //
  // Set the number of Newton-Raphson iterations which sets the precision. //
  // 0 is close to the single precision limit and is better than DUNE/HK   //
  // in the high statistics regime. Increasig N_Newton to 1,2,... rapidly  //
  // improves the precision at a modest computational cost                 //
  // --------------------------------------------------------------------- //
  N_Newton = 0;

  // ------------------------------------- //
  // Set the vacuum oscillation parameters //
  // ------------------------------------- //
  s12sq = 0.31;
  s13sq = 0.02;
  s23sq = 0.55;
  delta = 0.7 * M_PI;
  Dmsq21 = 7.5e-5; // eV^2
  Dmsq31 = 2.5e-3; // eV^2
  
  // ------------------------------------------ //
  // Calculate all 9 oscillationa probabilities //
  // ------------------------------------------ //
  Probability_Matter_LBL(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, L, E, rho, Ye, N_Newton, &probs_returned);

  // --------------------------- //
  // Print out the probabilities //
  // --------------------------- //
  printf("L = %g E = %g rho = %g\n", L, E, rho);
  printf("Probabilities:\n");
  printf("alpha beta P(nu_alpha -> nu_beta)\n");
  for (int alpha = 0; alpha < 3; alpha++)
    {
      for (int beta = 0; beta < 3; beta++)
	{
	  printf("%d %d %g\n", alpha, beta, probs_returned[alpha][beta]);
	} // beta, 3
    } // alpha, 3
}

int OscProbCalcerNuFASTLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerNuFASTLinear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
