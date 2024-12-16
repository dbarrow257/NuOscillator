#include "OscProbCalcer_OscProb.h"

#include "TMath.h"

OscProbCalcerOscProb::OscProbCalcerOscProb(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config
  if (!Config_["OscProbCalcerSetup"]["PMNSType"]) {
    std::cerr << "Expected to find a 'PMNSType' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  std::string OscMatrix = Config_["OscProbCalcerSetup"]["PMNSType"].as<std::string>();

  if (!Config_["OscProbCalcerSetup"]["PREMFile"]) {
    std::cerr << "Expected to find a 'PREMFile' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  premfile = Config_["OscProbCalcerSetup"]["PREMFile"].as<std::string>();
  //=======

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  fOscType = PMNS_StrToInt(OscMatrix);
  fNOscParams = GetNOscParams(fOscType);
  
  std::cout << "PMNS Type : " << fOscType << std::endl;
  std::cout << "Number of parameters : " << fNOscParams << std::endl;
  std::cout << "PREM Model : " << premfile << std::endl;
}

OscProbCalcerOscProb::~OscProbCalcerOscProb() {
}

void OscProbCalcerOscProb::SetupPropagator() {
  
  PremModel = OscProb::PremModel(premfile);

}

void OscProbCalcerOscProb::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {

  double det_radius = 6371. - OscParams[kDetDepth];

  PremModel.SetDetPos(det_radius);

  switch(fOscType) {
    case kFast: 
      CalcProbPMNS_Fast(OscParams);
      break;

    case kPMNSSterile1:
      CalcProbPMNS_Sterile(OscParams, 4);
      break;

    case kPMNSSterile2:
      CalcProbPMNS_Sterile(OscParams, 5);
      break;

    case kPMNSSterile3:
      CalcProbPMNS_Sterile(OscParams, 6);
      break;

    case kDecay:
      CalcProbPMNS_Decay(OscParams);
      break;
    
    case kDeco:
      CalcProbPMNS_Deco(OscParams);
      break;

    case kNSI:
      CalcProbPMNS_NSI(OscParams);
      break;

    case kSNSI:
      CalcProbPMNS_SNSI(OscParams);
      break;

    case kIter:
      CalcProbPMNS_Iter(OscParams);
      break;

    case kNUNM:
      CalcProbPMNS_NUNM(OscParams);
      break;

    default:
      break;
  }
}

void OscProbCalcerOscProb::CalcProbPMNS_Fast(const std::vector<FLOAT_T>& OscParams) {
  double energy, cosZ;
  double weight;
  int index;

  OscProb::PMNS_Fast myPMNS;

  SetPMNSParams_Fast(&myPMNS, OscParams);
  

  for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

    cosZ = fCosineZArray[iCosineZ];

    PremModel.FillPath(cosZ);

    myPMNS.SetPath(PremModel.GetNuPath());
    
    for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {

      energy = fEnergyArray[iEnergy];

      for (int iNuType = 0; iNuType < fNNeutrinoTypes; iNuType++) {

        myPMNS.SetIsNuBar(fNeutrinoTypes[iNuType]==Nubar);

        for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
          
          weight = myPMNS.Prob(fOscillationChannels[iOscChannel].GeneratedFlavour-1, fOscillationChannels[iOscChannel].DetectedFlavour-1, energy);
          index = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosineZ);

          fWeightArray[index] = weight;

        }

      }
    }
  }

}

void OscProbCalcerOscProb::CalcProbPMNS_Sterile(const std::vector<FLOAT_T>& OscParams, int neutrino_number) {
  double energy, cosZ;
  double weight;
  int index;

  OscProb::PMNS_Sterile myPMNS(neutrino_number);

  SetPMNSParams_Sterile(&myPMNS, OscParams);
  
  std::cout << "Proba with PMNS Sterile and " << neutrino_number << " neutrinos" << std::endl;

  for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

    cosZ = fCosineZArray[iCosineZ];

    PremModel.FillPath(cosZ);

    myPMNS.SetPath(PremModel.GetNuPath());
    
    for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {

      energy = fEnergyArray[iEnergy];

      for (int iNuType = 0; iNuType < fNNeutrinoTypes; iNuType++) {

        myPMNS.SetIsNuBar(fNeutrinoTypes[iNuType]==Nubar);

        for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
          
          weight = myPMNS.Prob(fOscillationChannels[iOscChannel].GeneratedFlavour-1, fOscillationChannels[iOscChannel].DetectedFlavour-1, energy);
          index = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosineZ);

          fWeightArray[index] = weight;

        }

      }
    }
  }
}

void OscProbCalcerOscProb::CalcProbPMNS_Decay(const std::vector<FLOAT_T>& OscParams) {
  double energy, cosZ;
  double weight;
  int index;

  OscProb::PMNS_Decay myPMNS;

  SetPMNSParams_Decay(&myPMNS, OscParams);

  for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

    cosZ = fCosineZArray[iCosineZ];

    PremModel.FillPath(cosZ);

    myPMNS.SetPath(PremModel.GetNuPath());
    
    for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {

      energy = fEnergyArray[iEnergy];

      for (int iNuType = 0; iNuType < fNNeutrinoTypes; iNuType++) {

        myPMNS.SetIsNuBar(fNeutrinoTypes[iNuType]==Nubar);

        for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
          
          weight = myPMNS.Prob(fOscillationChannels[iOscChannel].GeneratedFlavour-1, fOscillationChannels[iOscChannel].DetectedFlavour-1, energy);
          index = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosineZ);

          fWeightArray[index] = weight;

        }

      }
    }
  }
}

void OscProbCalcerOscProb::CalcProbPMNS_Deco(const std::vector<FLOAT_T>& OscParams) {
  double energy, cosZ;
  double weight;
  int index;

  OscProb::PMNS_Deco myPMNS;

  SetPMNSParams_Deco(&myPMNS, OscParams);

  for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

    cosZ = fCosineZArray[iCosineZ];

    PremModel.FillPath(cosZ);

    myPMNS.SetPath(PremModel.GetNuPath());
    
    for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {

      energy = fEnergyArray[iEnergy];

      for (int iNuType = 0; iNuType < fNNeutrinoTypes; iNuType++) {

        myPMNS.SetIsNuBar(fNeutrinoTypes[iNuType]==Nubar);

        for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
          
          weight = myPMNS.Prob(fOscillationChannels[iOscChannel].GeneratedFlavour-1, fOscillationChannels[iOscChannel].DetectedFlavour-1, energy);
          index = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosineZ);

          fWeightArray[index] = weight;

        }

      }
    }
  }
}

void OscProbCalcerOscProb::CalcProbPMNS_NSI(const std::vector<FLOAT_T>& OscParams) {
  double energy, cosZ;
  double weight;
  int index;

  OscProb::PMNS_NSI myPMNS;

  SetPMNSParams_NSI(&myPMNS, OscParams);

  for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

    cosZ = fCosineZArray[iCosineZ];

    PremModel.FillPath(cosZ);

    myPMNS.SetPath(PremModel.GetNuPath());
    
    for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {

      energy = fEnergyArray[iEnergy];

      for (int iNuType = 0; iNuType < fNNeutrinoTypes; iNuType++) {

        myPMNS.SetIsNuBar(fNeutrinoTypes[iNuType]==Nubar);

        for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
          
          weight = myPMNS.Prob(fOscillationChannels[iOscChannel].GeneratedFlavour-1, fOscillationChannels[iOscChannel].DetectedFlavour-1, energy);
          index = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosineZ);

          fWeightArray[index] = weight;

        }

      }
    }
  }
}

void OscProbCalcerOscProb::CalcProbPMNS_SNSI(const std::vector<FLOAT_T>& OscParams) {
  double energy, cosZ;
  double weight;
  int index;

  OscProb::PMNS_SNSI myPMNS;

  SetPMNSParams_SNSI(&myPMNS, OscParams);

  for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

    cosZ = fCosineZArray[iCosineZ];

    PremModel.FillPath(cosZ);

    myPMNS.SetPath(PremModel.GetNuPath());
    
    for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {

      energy = fEnergyArray[iEnergy];

      for (int iNuType = 0; iNuType < fNNeutrinoTypes; iNuType++) {

        myPMNS.SetIsNuBar(fNeutrinoTypes[iNuType]==Nubar);

        for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
          
          weight = myPMNS.Prob(fOscillationChannels[iOscChannel].GeneratedFlavour-1, fOscillationChannels[iOscChannel].DetectedFlavour-1, energy);
          index = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosineZ);

          fWeightArray[index] = weight;

        }

      }
    }
  }
}

void OscProbCalcerOscProb::CalcProbPMNS_Iter(const std::vector<FLOAT_T>& OscParams) {
  double energy, cosZ;
  double weight;
  int index;

  OscProb::PMNS_Iter myPMNS;

  SetPMNSParams_Iter(&myPMNS, OscParams);

  for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

    cosZ = fCosineZArray[iCosineZ];

    PremModel.FillPath(cosZ);

    myPMNS.SetPath(PremModel.GetNuPath());
    
    for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {

      energy = fEnergyArray[iEnergy];

      for (int iNuType = 0; iNuType < fNNeutrinoTypes; iNuType++) {

        myPMNS.SetIsNuBar(fNeutrinoTypes[iNuType]==Nubar);

        for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
          
          weight = myPMNS.Prob(fOscillationChannels[iOscChannel].GeneratedFlavour-1, fOscillationChannels[iOscChannel].DetectedFlavour-1, energy);
          index = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosineZ);

          fWeightArray[index] = weight;

        }

      }
    }
  }
}

void OscProbCalcerOscProb::CalcProbPMNS_NUNM(const std::vector<FLOAT_T>& OscParams) {
  double energy, cosZ;
  double weight;
  int index;

  OscProb::PMNS_NUNM myPMNS;

  SetPMNSParams_NUNM(&myPMNS, OscParams);

  for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

    cosZ = fCosineZArray[iCosineZ];

    PremModel.FillPath(cosZ);

    myPMNS.SetPath(PremModel.GetNuPath());
    
    for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {

      energy = fEnergyArray[iEnergy];

      for (int iNuType = 0; iNuType < fNNeutrinoTypes; iNuType++) {

        myPMNS.SetIsNuBar(fNeutrinoTypes[iNuType]==Nubar);

        for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
          
          weight = myPMNS.Prob(fOscillationChannels[iOscChannel].GeneratedFlavour-1, fOscillationChannels[iOscChannel].DetectedFlavour-1, energy);
          index = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosineZ);

          fWeightArray[index] = weight;

        }

      }
    }
  }
}

int OscProbCalcerOscProb::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNCosineZPoints*fNEnergyPoints + OscChanIndex*fNCosineZPoints*fNEnergyPoints + CosineZIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerOscProb::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNCosineZPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}

void OscProbCalcerOscProb::SetPMNSParams_Fast(OscProb::PMNS_Fast *Fast, const std::vector<FLOAT_T>& OscParams) {
  double th12, th13, th23;
  double dm21, dm31;
  double dcp;

  th12 = asin(sqrt(OscParams[kTH12]));
  th13 = asin(sqrt(OscParams[kTH13]));
  th23 = asin(sqrt(OscParams[kTH23]));
  dcp = OscParams[kDCP];
  dm21 = OscParams[kDM12];

  //Need to convert OscParams[kDM23] to kDM31
  dm31 = OscParams[kDM23]+OscParams[kDM12]; // eV^2

  // Set PMNS parameters
  Fast->SetDm(2, dm21);
  Fast->SetDm(3, dm31);
  Fast->SetAngle(1,2, th12);
  Fast->SetAngle(1,3, th13);
  Fast->SetAngle(2,3, th23);
  Fast->SetDelta(1,3, dcp);
}

void OscProbCalcerOscProb::SetPMNSParams_Sterile(OscProb::PMNS_Sterile *Sterile, const std::vector<FLOAT_T>& OscParams) {
  double th12, th13, th23;
  double dm21, dm31;
  double dcp;

  th12 = asin(sqrt(OscParams[kTH12]));
  th13 = asin(sqrt(OscParams[kTH13]));
  th23 = asin(sqrt(OscParams[kTH23]));
  dcp = OscParams[kDCP];
  dm21 = OscParams[kDM12];

  //Need to convert OscParams[kDM23] to kDM31
  dm31 = OscParams[kDM23]+OscParams[kDM12]; // eV^2

  // Set PMNS parameters for standard oscillation
  Sterile->SetDm(2, dm21);
  Sterile->SetDm(3, dm31);
  Sterile->SetAngle(1,2, th12);
  Sterile->SetAngle(1,3, th13);
  Sterile->SetAngle(2,3, th23);
  Sterile->SetDelta(1,3, dcp);

  //Set PMNS parameters for first sterile state
  Sterile->SetDm(4, OscParams[kDM14]);
  Sterile->SetAngle(1,4, asin(sqrt(OscParams[kTH14])));
  Sterile->SetAngle(2,4, asin(sqrt(OscParams[kTH24])));
  Sterile->SetAngle(3,4, asin(sqrt(OscParams[kTH34])));
  Sterile->SetDelta(1,4, OscParams[kDelta14]);
  Sterile->SetDelta(2,4, OscParams[kDelta24]);

  //Set PMNS parameters for second sterile state if there is one
  if(fOscType == kPMNSSterile2 || fOscType == kPMNSSterile3) {
    Sterile->SetDm(5, OscParams[kDM15]);
    Sterile->SetAngle(1,5, asin(sqrt(OscParams[kTH15])));
    Sterile->SetAngle(2,5, asin(sqrt(OscParams[kTH25])));
    Sterile->SetAngle(3,5, asin(sqrt(OscParams[kTH35])));
    Sterile->SetAngle(4,5, asin(sqrt(OscParams[kTH45])));
    Sterile->SetDelta(1,5, OscParams[kDelta15]);
    Sterile->SetDelta(2,5, OscParams[kDelta25]);
    Sterile->SetDelta(3,5, OscParams[kDelta35]);
  }

  //Set PMNS parameters for third sterile state if there is one
  if(fOscType == kPMNSSterile3) {
    Sterile->SetDm(6, OscParams[kDM16]);
    Sterile->SetAngle(1,6, asin(sqrt(OscParams[kTH16])));
    Sterile->SetAngle(2,6, asin(sqrt(OscParams[kTH26])));
    Sterile->SetAngle(3,6, asin(sqrt(OscParams[kTH36])));
    Sterile->SetAngle(4,6, asin(sqrt(OscParams[kTH46])));
    Sterile->SetAngle(5,6, asin(sqrt(OscParams[kTH56])));
    Sterile->SetDelta(1,6, OscParams[kDelta16]);
    Sterile->SetDelta(2,6, OscParams[kDelta26]);
    Sterile->SetDelta(3,6, OscParams[kDelta36]);
    Sterile->SetDelta(4,6, OscParams[kDelta46]);
  }

}

void OscProbCalcerOscProb::SetPMNSParams_Decay(OscProb::PMNS_Decay *Decay, const std::vector<FLOAT_T>& OscParams) {
  double th12, th13, th23;
  double dm21, dm31;
  double dcp;
  double alpha2, alpha3;

  th12 = asin(sqrt(OscParams[kTH12]));
  th13 = asin(sqrt(OscParams[kTH13]));
  th23 = asin(sqrt(OscParams[kTH23]));
  dcp = OscParams[kDCP];
  dm21 = OscParams[kDM12];
  alpha2 = OscParams[kAlpha2];
  alpha3 = OscParams[kAlpha3];

  //Need to convert OscParams[kDM23] to kDM31
  dm31 = OscParams[kDM23]+OscParams[kDM12]; // eV^2

  // Set PMNS parameters
  Decay->SetDm(2, dm21);
  Decay->SetDm(3, dm31);
  Decay->SetAngle(1,2, th12);
  Decay->SetAngle(1,3, th13);
  Decay->SetAngle(2,3, th23);
  Decay->SetDelta(1,3, dcp);
  Decay->SetAlpha2(alpha2);
  Decay->SetAlpha3(alpha3);
}

void OscProbCalcerOscProb::SetPMNSParams_Deco(OscProb::PMNS_Deco *Deco, const std::vector<FLOAT_T>& OscParams) {
  double th12, th13, th23;
  double dm21, dm31;
  double dcp;
  double gamma21, gamma31;
  double deco_angle, power;

  th12 = asin(sqrt(OscParams[kTH12]));
  th13 = asin(sqrt(OscParams[kTH13]));
  th23 = asin(sqrt(OscParams[kTH23]));
  dcp = OscParams[kDCP];
  dm21 = OscParams[kDM12];
  gamma21 = OscParams[kGamma21];
  gamma31 = OscParams[kGamma31];
  deco_angle = OscParams[kDecoAngle];
  power = OscParams[kPower];

  //Need to convert OscParams[kDM23] to kDM31
  dm31 = OscParams[kDM23]+OscParams[kDM12]; // eV^2

  // Set PMNS parameters
  Deco->SetDm(2, dm21);
  Deco->SetDm(3, dm31);
  Deco->SetAngle(1,2, th12);
  Deco->SetAngle(1,3, th13);
  Deco->SetAngle(2,3, th23);
  Deco->SetDelta(1,3, dcp);
  Deco->SetGamma(2, gamma21);
  Deco->SetGamma(3, gamma31);
  Deco->SetDecoAngle(deco_angle);
  Deco->SetPower(power);
}

void OscProbCalcerOscProb::SetPMNSParams_NSI(OscProb::PMNS_NSI *NSI, const std::vector<FLOAT_T>& OscParams) {
  double th12, th13, th23;
  double dm21, dm31;
  double dcp;

  th12 = asin(sqrt(OscParams[kTH12]));
  th13 = asin(sqrt(OscParams[kTH13]));
  th23 = asin(sqrt(OscParams[kTH23]));
  dcp = OscParams[kDCP];
  dm21 = OscParams[kDM12];

  //Need to convert OscParams[kDM23] to kDM31
  dm31 = OscParams[kDM23]+OscParams[kDM12]; // eV^2

  // Set PMNS parameters
  NSI->SetDm(2, dm21);
  NSI->SetDm(3, dm31);
  NSI->SetAngle(1,2, th12);
  NSI->SetAngle(1,3, th13);
  NSI->SetAngle(2,3, th23);
  NSI->SetDelta(1,3, dcp);
  NSI->SetNSI(OscParams[kEps_ee], OscParams[kEps_emu], OscParams[kEps_etau], OscParams[kEps_mumu], OscParams[kEps_mutau], OscParams[kEps_tautau],
              OscParams[kDelta_emu], OscParams[kDelta_etau], OscParams[kDelta_mutau]);
  NSI->SetFermCoup(OscParams[kElecCoup], OscParams[kUpCoup], OscParams[kDownCoup]);
}

void OscProbCalcerOscProb::SetPMNSParams_SNSI(OscProb::PMNS_SNSI *SNSI, const std::vector<FLOAT_T>& OscParams) {
  double th12, th13, th23;
  double dm21, dm31;
  double dcp;

  th12 = asin(sqrt(OscParams[kTH12]));
  th13 = asin(sqrt(OscParams[kTH13]));
  th23 = asin(sqrt(OscParams[kTH23]));
  dcp = OscParams[kDCP];
  dm21 = OscParams[kDM12];

  //Need to convert OscParams[kDM23] to kDM31
  dm31 = OscParams[kDM23]+OscParams[kDM12]; // eV^2

  // Set PMNS parameters
  SNSI->SetDm(2, dm21);
  SNSI->SetDm(3, dm31);
  SNSI->SetAngle(1,2, th12);
  SNSI->SetAngle(1,3, th13);
  SNSI->SetAngle(2,3, th23);
  SNSI->SetDelta(1,3, dcp);
  SNSI->SetNSI(OscParams[kEps_ee], OscParams[kEps_emu], OscParams[kEps_etau], OscParams[kEps_mumu], OscParams[kEps_mutau], OscParams[kEps_tautau],
              OscParams[kDelta_emu], OscParams[kDelta_etau], OscParams[kDelta_mutau]);
  SNSI->SetFermCoup(OscParams[kElecCoup], OscParams[kUpCoup], OscParams[kDownCoup]);
  SNSI->SetLowestMass(OscParams[kLightMass]);
}

void OscProbCalcerOscProb::SetPMNSParams_Iter(OscProb::PMNS_Iter *Iter, const std::vector<FLOAT_T>& OscParams) {
  double th12, th13, th23;
  double dm21, dm31;
  double dcp;
  double prec;

  prec = OscParams[kPrec];

  th12 = asin(sqrt(OscParams[kTH12]));
  th13 = asin(sqrt(OscParams[kTH13]));
  th23 = asin(sqrt(OscParams[kTH23]));
  dcp = OscParams[kDCP];
  dm21 = OscParams[kDM12];

  //Need to convert OscParams[kDM23] to kDM31
  dm31 = OscParams[kDM23]+OscParams[kDM12]; // eV^2

  // Set PMNS parameters
  Iter->SetDm(2, dm21);
  Iter->SetDm(3, dm31);
  Iter->SetAngle(1,2, th12);
  Iter->SetAngle(1,3, th13);
  Iter->SetAngle(2,3, th23);
  Iter->SetDelta(1,3, dcp);
  Iter->SetPrec(prec);
}

void OscProbCalcerOscProb::SetPMNSParams_NUNM(OscProb::PMNS_NUNM *NUNM, const std::vector<FLOAT_T>& OscParams) {
  double th12, th13, th23;
  double dm21, dm31;
  double dcp;
 
  th12 = asin(sqrt(OscParams[kTH12]));
  th13 = asin(sqrt(OscParams[kTH13]));
  th23 = asin(sqrt(OscParams[kTH23]));
  dcp = OscParams[kDCP];
  dm21 = OscParams[kDM12];

  //Need to convert OscParams[kDM23] to kDM31
  dm31 = OscParams[kDM23]+OscParams[kDM12]; // eV^2

  // Set PMNS parameters
  NUNM->SetDm(2, dm21);
  NUNM->SetDm(3, dm31);
  NUNM->SetAngle(1,2, th12);
  NUNM->SetAngle(1,3, th13);
  NUNM->SetAngle(2,3, th23);
  NUNM->SetDelta(1,3, dcp);
  NUNM->SetAlpha_11(OscParams[kAlpha11]);
  NUNM->SetAlpha_22(OscParams[kAlpha22]);
  NUNM->SetAlpha_33(OscParams[kAlpha33]);
  NUNM->SetAlpha_21(OscParams[kAlpha21], OscParams[kPhi21]);
  NUNM->SetAlpha_31(OscParams[kAlpha31], OscParams[kPhi31]);
  NUNM->SetAlpha_32(OscParams[kAlpha32], OscParams[kPhi32]);
  NUNM->SetFracVnc(OscParams[kFracVnc]);
}

int OscProbCalcerOscProb::PMNS_StrToInt(std::string PMNSType) {
  if (PMNSType == "Fast" || PMNSType == "fast") {
    return kFast;
  } 
  else if (PMNSType == "Sterile+1" || PMNSType == "Sterile1") {
    return kPMNSSterile1;
  } 
  else if (PMNSType == "Sterile+2" || PMNSType == "Sterile2") {
    return kPMNSSterile2;
  }
  else if (PMNSType == "Sterile+3" || PMNSType == "Sterile3") {
    return kPMNSSterile3;
  }
  else if (PMNSType == "Decay" || PMNSType == "decay") {
    return kDecay;
  }
  else if (PMNSType == "Deco" || PMNSType == "deco") {
    return kDeco;
  }
  else if (PMNSType == "NSI" || PMNSType == "nsi") {
    return kNSI;
  }
  else if (PMNSType == "LIV" || PMNSType == "liv") {
    return kLIV;
  }
  else if (PMNSType == "NUNM" || PMNSType == "nunm") {
    return kNUNM;
  }
  else if (PMNSType == "SNSI" || PMNSType == "snsi") {
    return kSNSI;
  }
  else if (PMNSType == "Iter" || PMNSType == "iter") {
    return kIter;
  }
  else {
    std::cerr << "Invalid PMNS matrix type provided:" << PMNSType << std::endl;
    throw;
  }
  
  return -1;
}


int OscProbCalcerOscProb::GetNOscParams(int OscType) {
  if (OscType == kFast) {
    return kNOscParams;
  }
  else if (OscType == kPMNSSterile1) {
    return kNOscParams+6;
  }
  else if (OscType == kPMNSSterile2) {
    return kNOscParams+14;
  }
  else if (OscType == kPMNSSterile3) {
    return kNOscParams+24;
  }
  else if (OscType == kDecay) {
    return kNOscParams+2;
  }
  else if (OscType == kDeco) {
    return kNOscParams+4;
  }
  else if (OscType == kNSI) {
    return kNOscParams+12;
  }
  else if (OscType == kIter) {
    return kNOscParams+1;
  }
  else if (OscType == kNUNM) {
    return kNOscParams+10;
  }
  else if (OscType == kLIV) {
    std::cerr << "LIV PMNS matrix not implemented yet" << std::endl;
    throw;
  }
  else if (OscType == kSNSI) {
    return kNOscParams+13;
  }
  else {
    std::cerr << "Invalid PMNS matrix type provided:" << OscType << std::endl;
    throw;
  }

  return -1;
}
