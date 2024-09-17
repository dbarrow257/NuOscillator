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
  L = OscParams[kPATHL]; // km
  rho = OscParams[kDENS]; // g/cc

  //Ye = OscParams[kELECDENS];
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
  s12sq = OscParams[kTH12];
  s13sq = OscParams[kTH13];
  s23sq = OscParams[kTH23];
  delta = OscParams[kDCP];
  Dmsq21 = OscParams[kDM12];

  //Need to convert OscParams[kDM23] to kDM31
  Dmsq31 = 2.5e-3; // eV^2
  
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

	double Weight = probs_returned[fOscillationChannels[iOscChannel].GeneratedFlavour][fOscillationChannels[iOscChannel].DetectedFlavour];

	//Cancel floating point precision
	if (Weight<0. && Weight>-1e-6) {Weight = 0.;}

	if (Weight<0. || Weight > 1.) {
	  std::cout << "s12sq:" << s12sq << " s13sq:" << s13sq << " s23sq:" << s23sq << " delta:" << delta << " Dmsq21:" << Dmsq21 << " Dmsq31:" << Dmsq31 << " L:" << L << " E:" << E << " rho:" << rho << " Ye:" << Ye << " N_Newton:" << N_Newton << std::endl;
	  std::cout << "iOscProb:" << iOscProb << " iNuType:" << iNuType << " iOscChannel:" << iOscChannel << " IndexToFill:" << IndexToFill << " fWeightArray[IndexToFill+iOscProb]:" << Weight << std::endl;
	  throw;
	}
	
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
  long nCalculationPoints = fNEnergyPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
