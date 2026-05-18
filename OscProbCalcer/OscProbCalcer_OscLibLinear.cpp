#include "OscProbCalcer_OscLibLinear.h"

OscProbCalcerOscLibLinear::OscProbCalcerOscLibLinear(YAML::Node Config_) : OscProbCalcerBase(Config_) {
  //=======
  if (!Config_["OscProbCalcerSetup"]["PMNSType"]) {
    std::cerr << "Expected to find a 'PMNSType' Node within the 'OscProbCalcerSetup' Node" << std::endl;
    throw std::runtime_error("YAML node not found");
  }
  std::string PMNSType = Config_["OscProbCalcerSetup"]["PMNSType"].as<std::string>();
  fImplementationName += "-" + PMNSType;
  fOscType = PMNS_StrToInt(PMNSType);
  SetOscParams();
  
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
  if(fOscType == kPMNS) OscLib = new osc::_OscCalcPMNS<FLOAT_T>();
  if(fOscType == kNSI) OscLib = new osc::OscCalcPMNS_NSI();
}

int OscProbCalcerOscLibLinear::PMNS_StrToInt(const std::string& PMNSType) {
  if (PMNSType == "PMNS" || PMNSType == "pmns")          return kPMNS;
  if (PMNSType == "NSI" || PMNSType == "nsi")            return kNSI;

  std::cerr << "Invalid PMNS matrix type provided:" << PMNSType << std::endl;
  throw std::runtime_error("Invalid setup");

  return -1;
}

int GetGeneratedFlavour(int generatedFlavour, int neutrinoType) {
  int GenFlav = -1;
  switch (generatedFlavour) {
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
      throw std::invalid_argument("Invalid generated flavour");
  }
  return GenFlav * neutrinoType;
}

void OscProbCalcerOscLibLinear::CalculateProbabilities() {
  for (int iOscPar=0;iOscPar<=kTH13;iOscPar++) {
    if (GetOscillationParameter(iOscPar) < 0) {
      std::cerr << "Invalid oscillation parameter (Can not sqrt this value)!:" << GetOscillationParameter(iOscPar) << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }

  const FLOAT_T theta12 = std::asin(std::sqrt(GetOscillationParameter(kTH12)));
  const FLOAT_T theta23 = std::asin(std::sqrt(GetOscillationParameter(kTH23)));
  const FLOAT_T theta13 = std::asin(std::sqrt(GetOscillationParameter(kTH13)));
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
  
  if(fOscType == kNSI) {
    auto OscLib_NSI = static_cast<osc::OscCalcPMNS_NSI*>(OscLib);
    OscLib_NSI->SetEps_ee(GetOscillationParameter(kEps_ee));
    OscLib_NSI->SetEps_emu(GetOscillationParameter(kEps_emu));
    OscLib_NSI->SetEps_etau(GetOscillationParameter(kEps_etau));
    OscLib_NSI->SetEps_mumu(GetOscillationParameter(kEps_mumu));
    OscLib_NSI->SetEps_mutau(GetOscillationParameter(kEps_mutau));
    OscLib_NSI->SetEps_tautau(GetOscillationParameter(kEps_tautau));
    OscLib_NSI->SetDelta_emu(GetOscillationParameter(kDelta_emu));
    OscLib_NSI->SetDelta_etau(GetOscillationParameter(kDelta_etau));
    OscLib_NSI->SetDelta_mutau(GetOscillationParameter(kDelta_mutau));
  }
  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
    for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {
      
      int IndexToFill = iNuType*fNOscillationChannels*fNEnergyPoints + iOscChannel*fNEnergyPoints;

      for (int iOscProb=0;iOscProb<fNEnergyPoints;iOscProb++) {
        FLOAT_T Energy = fEnergyArray[iOscProb];
        int GenFlav = GetGeneratedFlavour(fOscillationChannels[iOscChannel].GeneratedFlavour, fNeutrinoTypes[iNuType]);
        int DetFlav = GetGeneratedFlavour(fOscillationChannels[iOscChannel].DetectedFlavour, fNeutrinoTypes[iNuType]);

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


void OscProbCalcerOscLibLinear::SetOscParams() {
  std::vector<std::string> OscParNames;

  switch (fOscType) {
    case kPMNS:
      OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp","path_length","matter_density"};
      break;
    case kNSI:
      OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp","path_length","matter_density",
                     "eps_ee", "eps_emu", "eps_etau", "eps_mumu", "eps_mutau", "eps_tautau", "delta_emu", "delta_etau", "delta_mutau"};
      break;
    default:
      std::cerr << "Invalid Oscillation Model type provided:" << fOscType << std::endl;
      throw std::runtime_error("Invalid Oscillation Model");
  }

  SetExpectedParameterNames(OscParNames);
}
