#include "OscProbCalcer_OscProb.h"
#include <fstream>

OscProbCalcerOscProb::OscProbCalcerOscProb(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config
  if (!Config_["OscProbCalcerSetup"]["PMNSType"]) {
    std::cerr << "Expected to find a 'PMNSType' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw std::runtime_error("YAML node not found");
  }

  std::string OscMatrix = Config_["OscProbCalcerSetup"]["PMNSType"].as<std::string>();

  if (!Config_["OscProbCalcerSetup"]["PREMFile"]) {
    std::cerr << "Expected to find a 'PREMFile' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw std::runtime_error("YAML node not found");
  }

  premfile = Config_["OscProbCalcerSetup"]["PREMFile"].as<std::string>();

  if (!Config_["OscProbCalcerSetup"]["DetDepth"]) {
    std::cerr << "Expected to find a 'DetDepth' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw std::runtime_error("YAML node not found");
  }

  fDetDepth = Config_["OscProbCalcerSetup"]["DetDepth"].as<double>();
  //=======

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  fOscType = PMNS_StrToInt(OscMatrix);
  fNOscParams = GetNOscParams(fOscType);

  fMaxGenFlavour = 1;
  fMaxDetFlavour = 1;
  for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
    int gflv = fOscillationChannels[iOscChannel].GeneratedFlavour;
    int dflv = fOscillationChannels[iOscChannel].DetectedFlavour;
    if(gflv>fMaxGenFlavour) fMaxGenFlavour = gflv;
    if(dflv>fMaxDetFlavour) fMaxDetFlavour = dflv;
  }

  std::cout << "PMNS Type : " << fOscType << std::endl;
  std::cout << "Number of parameters : " << fNOscParams << std::endl;
  std::cout << "PREM Model : " << premfile << std::endl;
  std::cout << "Detector depth : " << fDetDepth << "km" << std::endl;
}

OscProbCalcerOscProb::~OscProbCalcerOscProb() {
}

void OscProbCalcerOscProb::SetupPropagator() {

  std::ifstream file(premfile);
  if(!file) throw std::runtime_error("could not open PREM file " + premfile);

  PremModel = OscProb::PremModel(premfile);

}

void OscProbCalcerOscProb::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {

  double det_radius = 6371. - fDetDepth;

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

    case kLIV:
      CalcProbPMNS_LIV(OscParams);
      break;

    default:
      break;
  }
}

void OscProbCalcerOscProb::CalcProbPMNS(OscProb::PMNS_Base* myPMNS) {

  for (int iNuType = 0; iNuType < fNNeutrinoTypes; iNuType++) {

    myPMNS->SetIsNuBar(fNeutrinoTypes[iNuType]==Nubar);

    for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

      PremModel.FillPath(fCosineZArray[iCosineZ]);

      myPMNS->SetPath(PremModel.GetNuPath());

      for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {

        OscProb::matrixD probMatrix = myPMNS->ProbMatrix(fMaxGenFlavour,
                                                         fMaxDetFlavour,
                                                         fEnergyArray[iEnergy]);

        for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {

          int gflv = fOscillationChannels[iOscChannel].GeneratedFlavour-1;
          int dflv = fOscillationChannels[iOscChannel].DetectedFlavour-1;
          double weight = probMatrix[gflv][dflv];
          int index = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosineZ);

          fWeightArray[index] = weight;

        }

      }
    }
  }

}

void OscProbCalcerOscProb::CalcProbPMNS_Fast(const std::vector<FLOAT_T>& OscParams) {

  OscProb::PMNS_Fast myPMNS;
  SetPMNSParams(&myPMNS, OscParams);
  return CalcProbPMNS(&myPMNS);

}

void OscProbCalcerOscProb::CalcProbPMNS_Sterile(const std::vector<FLOAT_T>& OscParams, int neutrino_number) {

  OscProb::PMNS_Sterile myPMNS(neutrino_number);
  SetPMNSParams_Sterile(&myPMNS, OscParams);

  std::cout << "Proba with PMNS Sterile and " << neutrino_number << " neutrinos" << std::endl;

  return CalcProbPMNS(&myPMNS);

}

void OscProbCalcerOscProb::CalcProbPMNS_Decay(const std::vector<FLOAT_T>& OscParams) {

  OscProb::PMNS_Decay myPMNS;
  SetPMNSParams_Decay(&myPMNS, OscParams);
  return CalcProbPMNS(&myPMNS);

}

void OscProbCalcerOscProb::CalcProbPMNS_Deco(const std::vector<FLOAT_T>& OscParams) {

  OscProb::PMNS_Deco myPMNS;
  SetPMNSParams_Deco(&myPMNS, OscParams);
  return CalcProbPMNS(&myPMNS);

}

void OscProbCalcerOscProb::CalcProbPMNS_NSI(const std::vector<FLOAT_T>& OscParams) {

  OscProb::PMNS_NSI myPMNS;
  SetPMNSParams_NSI(&myPMNS, OscParams);
  return CalcProbPMNS(&myPMNS);

}

void OscProbCalcerOscProb::CalcProbPMNS_SNSI(const std::vector<FLOAT_T>& OscParams) {

  OscProb::PMNS_SNSI myPMNS;
  SetPMNSParams_SNSI(&myPMNS, OscParams);
  return CalcProbPMNS(&myPMNS);

}

void OscProbCalcerOscProb::CalcProbPMNS_Iter(const std::vector<FLOAT_T>& OscParams) {

  OscProb::PMNS_Iter myPMNS;
  SetPMNSParams_Iter(&myPMNS, OscParams);
  return CalcProbPMNS(&myPMNS);

}

void OscProbCalcerOscProb::CalcProbPMNS_NUNM(const std::vector<FLOAT_T>& OscParams) {

  OscProb::PMNS_NUNM myPMNS;
  SetPMNSParams_NUNM(&myPMNS, OscParams);
  return CalcProbPMNS(&myPMNS);

}

void OscProbCalcerOscProb::CalcProbPMNS_LIV(const std::vector<FLOAT_T>& OscParams) {

  OscProb::PMNS_LIV myPMNS;
  SetPMNSParams_LIV(&myPMNS, OscParams);
  return CalcProbPMNS(&myPMNS);

}

int OscProbCalcerOscProb::ReturnWeightArrayIndex(int NuTypeIndex,
                                                 int OscChanIndex,
                                                 int EnergyIndex,
                                                 int CosineZIndex) {
  int IndexToReturn = ((NuTypeIndex *  fNOscillationChannels +
                        OscChanIndex) * fNCosineZPoints +
                        CosineZIndex) * fNEnergyPoints +
                        EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerOscProb::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) *
                            fNCosineZPoints *
                            fNOscillationChannels *
                            fNNeutrinoTypes;
  return nCalculationPoints;
}

void OscProbCalcerOscProb::SetPMNSParams(OscProb::PMNS_Base* pmns,
                                         const std::vector<FLOAT_T>& OscParams) {

  double th12 = asin(sqrt(OscParams[kTH12]));
  double th13 = asin(sqrt(OscParams[kTH13]));
  double th23 = asin(sqrt(OscParams[kTH23]));
  double dcp  = OscParams[kDCP];
  double dm21 = OscParams[kDM12];

  //Need to convert OscParams[kDM23] to kDM31
  double dm31 = OscParams[kDM23] + OscParams[kDM12]; // eV^2

  // Set PMNS parameters
  pmns->SetDm(2, dm21);
  pmns->SetDm(3, dm31);
  pmns->SetAngle(1,2, th12);
  pmns->SetAngle(1,3, th13);
  pmns->SetAngle(2,3, th23);
  pmns->SetDelta(1,3, dcp);

}

void OscProbCalcerOscProb::SetPMNSParams_Sterile(OscProb::PMNS_Sterile* Sterile,
                                                 const std::vector<FLOAT_T>& OscParams) {
  SetPMNSParams(Sterile, OscParams);

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

void OscProbCalcerOscProb::SetPMNSParams_Decay(OscProb::PMNS_Decay* Decay,
                                               const std::vector<FLOAT_T>& OscParams) {

  SetPMNSParams(Decay, OscParams);

  // Set Decay parameters
  Decay->SetAlpha2(OscParams[kAlpha2]);
  Decay->SetAlpha3(OscParams[kAlpha3]);

}

void OscProbCalcerOscProb::SetPMNSParams_Deco(OscProb::PMNS_Deco* Deco,
                                              const std::vector<FLOAT_T>& OscParams) {

  SetPMNSParams(Deco, OscParams);

  // Set Deco parameters
  Deco->SetGamma(2, OscParams[kGamma21]);
  Deco->SetGamma(3, OscParams[kGamma31]);
  Deco->SetDecoAngle(OscParams[kDecoAngle]);
  Deco->SetPower(OscParams[kPower]);

}

void OscProbCalcerOscProb::SetPMNSParams_NSI(OscProb::PMNS_NSI* NSI,
                                             const std::vector<FLOAT_T>& OscParams) {
  SetPMNSParams(NSI, OscParams);

  // Set NSI parameters
  NSI->SetNSI(OscParams[kEps_ee],
              OscParams[kEps_emu],
              OscParams[kEps_etau],
              OscParams[kEps_mumu],
              OscParams[kEps_mutau],
              OscParams[kEps_tautau],
              OscParams[kDelta_emu],
              OscParams[kDelta_etau],
              OscParams[kDelta_mutau]);

  NSI->SetFermCoup(OscParams[kElecCoup],
                   OscParams[kUpCoup],
                   OscParams[kDownCoup]);

}

void OscProbCalcerOscProb::SetPMNSParams_SNSI(OscProb::PMNS_SNSI* SNSI,
                                              const std::vector<FLOAT_T>& OscParams) {

  SetPMNSParams_NSI(SNSI, OscParams);

  // Set SNSI parameters
  SNSI->SetLowestMass(OscParams[kLightMass]);

}

void OscProbCalcerOscProb::SetPMNSParams_Iter(OscProb::PMNS_Iter* Iter,
                                              const std::vector<FLOAT_T>& OscParams) {

  SetPMNSParams(Iter, OscParams);

  // Set Iter parameters
  Iter->SetPrec(OscParams[kPrec]);

}

void OscProbCalcerOscProb::SetPMNSParams_NUNM(OscProb::PMNS_NUNM* NUNM,
                                              const std::vector<FLOAT_T>& OscParams) {

  SetPMNSParams(NUNM, OscParams);

  // Set PMNS parameters
  NUNM->SetAlpha_11(OscParams[kAlpha11]);
  NUNM->SetAlpha_22(OscParams[kAlpha22]);
  NUNM->SetAlpha_33(OscParams[kAlpha33]);
  NUNM->SetAlpha_21(OscParams[kAlpha21], OscParams[kPhi21]);
  NUNM->SetAlpha_31(OscParams[kAlpha31], OscParams[kPhi31]);
  NUNM->SetAlpha_32(OscParams[kAlpha32], OscParams[kPhi32]);
  NUNM->SetFracVnc(OscParams[kFracVnc]);

}

void OscProbCalcerOscProb::SetPMNSParams_LIV(OscProb::PMNS_LIV *LIV, const std::vector<FLOAT_T>& OscParams) {

  SetPMNSParams(LIV, OscParams);

  // Set LIV parameters
  LIV->SetaT(0, 0, 3, OscParams[kaT_ee_3], 0.);
  LIV->SetaT(0, 1, 3, OscParams[kaT_emu_3], OscParams[kDelta_emu_3]);
  LIV->SetaT(0, 2, 3, OscParams[kaT_etau_3], OscParams[kDelta_etau_3]);
  LIV->SetaT(1, 1, 3, OscParams[kaT_mumu_3], 0.);
  LIV->SetaT(1, 2, 3, OscParams[kaT_mutau_3], OscParams[kDelta_mutau_3]);
  LIV->SetaT(2, 2, 3, OscParams[kaT_tautau_3], 0.);
  LIV->SetcT(0, 0, 4, OscParams[kcT_ee_4], 0.);
  LIV->SetcT(0, 1, 4, OscParams[kcT_emu_4], OscParams[kDelta_emu_4]);
  LIV->SetcT(0, 2, 4, OscParams[kcT_etau_4], OscParams[kDelta_etau_4]);
  LIV->SetcT(1, 1, 4, OscParams[kcT_mumu_4], 0.);
  LIV->SetcT(1, 2, 4, OscParams[kcT_mutau_4], OscParams[kDelta_mutau_4]);
  LIV->SetcT(2, 2, 4, OscParams[kcT_tautau_4], 0.);
  LIV->SetaT(0, 0, 5, OscParams[kaT_ee_5], 0.);
  LIV->SetaT(0, 1, 5, OscParams[kaT_emu_5], OscParams[kDelta_emu_5]);
  LIV->SetaT(0, 2, 5, OscParams[kaT_etau_5], OscParams[kDelta_etau_5]);
  LIV->SetaT(1, 1, 5, OscParams[kaT_mumu_5], 0.);
  LIV->SetaT(1, 2, 5, OscParams[kaT_mutau_5], OscParams[kDelta_mutau_5]);
  LIV->SetaT(2, 2, 5, OscParams[kaT_tautau_5], 0.);
  LIV->SetcT(0, 0, 6, OscParams[kcT_ee_6], 0.);
  LIV->SetcT(0, 1, 6, OscParams[kcT_emu_6], OscParams[kDelta_emu_6]);
  LIV->SetcT(0, 2, 6, OscParams[kcT_etau_6], OscParams[kDelta_etau_6]);
  LIV->SetcT(1, 1, 6, OscParams[kcT_mumu_6], 0.);
  LIV->SetcT(1, 2, 6, OscParams[kcT_mutau_6], OscParams[kDelta_mutau_6]);
  LIV->SetcT(2, 2, 6, OscParams[kcT_tautau_6], 0.);
  LIV->SetaT(0, 0, 7, OscParams[kaT_ee_7], 0.);
  LIV->SetaT(0, 1, 7, OscParams[kaT_emu_7], OscParams[kDelta_emu_7]);
  LIV->SetaT(0, 2, 7, OscParams[kaT_etau_7], OscParams[kDelta_etau_7]);
  LIV->SetaT(1, 1, 7, OscParams[kaT_mumu_7], 0.);
  LIV->SetaT(1, 2, 7, OscParams[kaT_mutau_7], OscParams[kDelta_mutau_7]);
  LIV->SetaT(2, 2, 7, OscParams[kaT_tautau_7], 0.);
  LIV->SetcT(0, 0, 8, OscParams[kcT_ee_8], 0.);
  LIV->SetcT(0, 1, 8, OscParams[kcT_emu_8], OscParams[kDelta_emu_8]);
  LIV->SetcT(0, 2, 8, OscParams[kcT_etau_8], OscParams[kDelta_etau_8]);
  LIV->SetcT(1, 1, 8, OscParams[kcT_mumu_8], 0.);
  LIV->SetcT(1, 2, 8, OscParams[kcT_mutau_8], OscParams[kDelta_mutau_8]);
  LIV->SetcT(2, 2, 8, OscParams[kcT_tautau_8], 0.);
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
    throw std::runtime_error("Invalid setup");
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
    return kNOscParams+54;
  }
  else if (OscType == kSNSI) {
    return kNOscParams+13;
  }
  else {
    std::cerr << "Invalid PMNS matrix type provided:" << OscType << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  return -1;
}
