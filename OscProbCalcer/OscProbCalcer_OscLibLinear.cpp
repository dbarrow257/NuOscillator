#include "OscProbCalcer_OscLibLinear.h"

OscProbCalcerOscLibLinear::OscProbCalcerOscLibLinear(YAML::Node Config_) : OscProbCalcerBase(Config_) {
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
}

OscProbCalcerOscLibLinear::~OscProbCalcerOscLibLinear() {
  if (OscLib) {delete OscLib;}
}

void OscProbCalcerOscLibLinear::SetupPropagator() {
  OscLib = new osc::_OscCalcPMNS<FLOAT_T>();
}

void OscProbCalcerOscLibLinear::CalculateProbabilities() {
  for (int iOscPar=0;iOscPar<=kTH13;iOscPar++) {
    if (GetOscillationParameter(iOscPar) < 0) {
      std::cerr << "Invalid oscillation parameter (Can not sqrt this value)!:" << GetOscillationParameter(iOscPar) << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }

  const FLOAT_T theta12 = asin(sqrt(GetOscillationParameter(kTH12)));
  const FLOAT_T theta23 = asin(sqrt(GetOscillationParameter(kTH23)));
  const FLOAT_T theta13 = asin(sqrt(GetOscillationParameter(kTH13)));
  const FLOAT_T Dmsq21 = GetOscillationParameter(kDM12);
  const FLOAT_T Dmsq32 = GetOscillationParameter(kDM23);
  const FLOAT_T delta = GetOscillationParameter(kDCP);
  const FLOAT_T L = GetOscillationParameter(kPATHL);
  const FLOAT_T rho = GetOscillationParameter(kDENS);
  
  OscLib->SetL(L);
  OscLib->SetRho(rho);
  OscLib->SetDmsq21(Dmsq21);
  OscLib->SetDmsq32(Dmsq32);
  OscLib->SetTh12(theta12);
  OscLib->SetTh13(theta13);
  OscLib->SetTh23(theta23);
  OscLib->SetdCP(delta);
  
  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
    for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {
      
      int IndexToFill = iNuType*fNOscillationChannels*fNEnergyPoints + iOscChannel*fNEnergyPoints;

      for (int iOscProb=0;iOscProb<fNEnergyPoints;iOscProb++) {
	FLOAT_T Energy = fEnergyArray[iOscProb];

	int GenFlav = -1;
	switch(fOscillationChannels[iOscChannel].GeneratedFlavour) {
	case 1:
	  GenFlav = 12;
	  break;
	case 2:
	  GenFlav = 14;
	  break;
	case 3:
	  GenFlav = 16;
	  break;
	default:
	  std::cerr << "Invalid flav" << std::endl;
	  throw;
	}
	GenFlav *= fNeutrinoTypes[iNuType];

	int DetFlav = -1;
	switch(fOscillationChannels[iOscChannel].DetectedFlavour) {
	case 1:
	  DetFlav = 12;
	  break;
	case 2:
	  DetFlav = 14;
	  break;
	case 3:
	  DetFlav = 16;
	  break;
	default:
	  std::cerr << "Invalid flav" << std::endl;
	  throw;
	}
	DetFlav *= fNeutrinoTypes[iNuType];

	fWeightArray[IndexToFill+iOscProb] = OscLib->P(GenFlav,DetFlav,Energy);
      }
    }
  }
}

int OscProbCalcerOscLibLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerOscLibLinear::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
