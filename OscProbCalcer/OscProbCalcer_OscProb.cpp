#include "OscProbCalcer_OscProb.h"
#include <fstream>
#include "inc/PMNS_Fast.h"
#include "inc/PMNS_Sterile.h"
#include "inc/PMNS_Decay.h"
#include "inc/PMNS_Deco.h"
#include "inc/PMNS_NSI.h"
#include "inc/PMNS_Iter.h"
#include "inc/PMNS_NUNM.h"
#include "inc/PMNS_SNSI.h"
#include "inc/PMNS_LIV.h"

OscProbCalcerOscProb::OscProbCalcerOscProb(YAML::Node Config_) :
  OscProbCalcerBase(Config_),
  fPMNSObj(nullptr)
{
  //=======
  //Grab information from the config
  if (!Config_["OscProbCalcerSetup"]["PMNSType"]) {
    std::cerr << "Expected to find a 'PMNSType' Node within the 'OscProbCalcerSetup' Node" << std::endl;
    throw std::runtime_error("YAML node not found");
  }

  std::string PMNSType = Config_["OscProbCalcerSetup"]["PMNSType"].as<std::string>();

  if (!Config_["General"]["CosineZIgnored"]) {
    std::cerr << "Expected to find a 'CosineZIgnored' Node within the 'General' Node" << std::endl;
    throw std::runtime_error("YAML node not found");
  }

  IgnoreCosineZBinning(Config_["General"]["CosineZIgnored"].as<bool>());

  if (!fCosineZIgnored) {
    if (!Config_["OscProbCalcerSetup"]["PREMFile"]) {
      std::cerr << "Expected to find a 'PREMFile' Node within the 'OscProbCalcerSetup' Node" << std::endl;
      throw std::runtime_error("YAML node not found");
    }

    fPremFile = Config_["OscProbCalcerSetup"]["PREMFile"].as<std::string>();

    if (!Config_["OscProbCalcerSetup"]["DetDepth"]) {
      std::cerr << "Expected to find a 'DetDepth' Node within the 'OscProbCalcerSetup' Node" << std::endl;
      throw std::runtime_error("YAML node not found");
    }

    fDetDepth = Config_["OscProbCalcerSetup"]["DetDepth"].as<double>();
  }
  //=======

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  fOscType = PMNS_StrToInt(PMNSType);
  fNOscParams = GetNOscParams();
  if(fCosineZIgnored) fNOscParams += 3;

  fMaxGenFlavour = 1;
  fMaxDetFlavour = 1;
  for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
    int gflv = fOscillationChannels[iOscChannel].GeneratedFlavour;
    int dflv = fOscillationChannels[iOscChannel].DetectedFlavour;
    if(gflv>fMaxGenFlavour) fMaxGenFlavour = gflv;
    if(dflv>fMaxDetFlavour) fMaxDetFlavour = dflv;
  }

  std::cout << "PMNS Type : " << PMNSType << std::endl;
  std::cout << "Number of parameters : " << fNOscParams << std::endl;
  std::cout << "Propagation Method : ";
  if (fCosineZIgnored) std::cout << "Linear" << std::endl;
  else {
    std::cout << "Atmospheric" << std::endl;
    std::cout << "PREM Model : " << fPremFile << std::endl;
    std::cout << "Detector depth : " << fDetDepth << "km" << std::endl;
  }
}

OscProbCalcerOscProb::~OscProbCalcerOscProb() {
  delete fPMNSObj;
}

int OscProbCalcerOscProb::GetNCosineZ() {
  return fCosineZIgnored ? 1 : fNCosineZPoints;
}

void OscProbCalcerOscProb::SetupPropagator() {

  if(!fCosineZIgnored) {
    std::ifstream file(fPremFile);
    if(!file) throw std::runtime_error("could not open PREM file " + fPremFile);

    double det_radius = 6371. - fDetDepth;

    fPremModel.SetDetPos(det_radius);
    fPremModel.LoadModel(fPremFile);
  }

  if(fPMNSObj) delete fPMNSObj;
  fPMNSObj = GetPMNSObj();

}

OscProb::PMNS_Base* OscProbCalcerOscProb::GetPMNSObj() {

 if(fOscType==kPMNSSterile1) return new OscProb::PMNS_Sterile(4);
 if(fOscType==kPMNSSterile2) return new OscProb::PMNS_Sterile(5);
 if(fOscType==kPMNSSterile3) return new OscProb::PMNS_Sterile(6);
 if(fOscType==kDecay)        return new OscProb::PMNS_Decay();
 if(fOscType==kDeco)         return new OscProb::PMNS_Deco();
 if(fOscType==kNSI)          return new OscProb::PMNS_NSI();
 if(fOscType==kSNSI)         return new OscProb::PMNS_SNSI();
 if(fOscType==kIter)         return new OscProb::PMNS_Iter();
 if(fOscType==kNUNM)         return new OscProb::PMNS_NUNM();
 if(fOscType==kLIV)          return new OscProb::PMNS_LIV();

 return new OscProb::PMNS_Fast();

}

void OscProbCalcerOscProb::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {

  SetPMNSParams(OscParams);

  CalcProbPMNS(OscParams);

}

void OscProbCalcerOscProb::SetPath(const std::vector<FLOAT_T>& OscParams,
                                   int iCosineZ) {

  if(fCosineZIgnored) {
    fPMNSObj->SetLength ( OscParams[fNOscParams - 3] );
    fPMNSObj->SetDensity( OscParams[fNOscParams - 2] );
    fPMNSObj->SetZoA    ( OscParams[fNOscParams - 1] );
  }
  else {
    fPremModel.FillPath(fCosineZArray[iCosineZ]);
    fPMNSObj->SetPath(fPremModel.GetNuPath());
  }

}

void OscProbCalcerOscProb::CalcProbPMNS(const std::vector<FLOAT_T>& OscParams) {

  for (int iNuType = 0; iNuType < fNNeutrinoTypes; iNuType++) {

    fPMNSObj->SetIsNuBar(fNeutrinoTypes[iNuType]==Nubar);

    for (int iCosineZ = 0; iCosineZ < GetNCosineZ(); iCosineZ++) {

      SetPath(OscParams, iCosineZ);

      for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {

        OscProb::matrixD probMatrix = fPMNSObj->ProbMatrix(fMaxGenFlavour,
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

int OscProbCalcerOscProb::ReturnWeightArrayIndex(int NuTypeIndex,
                                                 int OscChanIndex,
                                                 int EnergyIndex,
                                                 int CosineZIndex) {
  int IndexToReturn = ((NuTypeIndex *  fNOscillationChannels +
                        OscChanIndex) * GetNCosineZ() +
                        std::max(CosineZIndex,0)) * fNEnergyPoints +
                        EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerOscProb::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) *
                            GetNCosineZ() *
                            fNOscillationChannels *
                            fNNeutrinoTypes;
  return nCalculationPoints;
}

void OscProbCalcerOscProb::SetPMNSParams(const std::vector<FLOAT_T>& OscParams) {

  // Set PMNS parameters
  fPMNSObj->SetDm(2, OscParams[kDM12]);
  fPMNSObj->SetDm(3, OscParams[kDM23] + OscParams[kDM12]);
  fPMNSObj->SetAngle(1,2, asin(sqrt(OscParams[kTH12])));
  fPMNSObj->SetAngle(1,3, asin(sqrt(OscParams[kTH13])));
  fPMNSObj->SetAngle(2,3, asin(sqrt(OscParams[kTH23])));
  fPMNSObj->SetDelta(1,3, OscParams[kDCP]);

  //Set PMNS parameters for first sterile state
  if(OscProb::PMNS_Sterile* Sterile = dynamic_cast<OscProb::PMNS_Sterile*>(fPMNSObj)) {
    Sterile->SetDm(4, OscParams[kDM14]);
    Sterile->SetAngle(1,4, asin(sqrt(OscParams[kTH14])));
    Sterile->SetAngle(2,4, asin(sqrt(OscParams[kTH24])));
    Sterile->SetAngle(3,4, asin(sqrt(OscParams[kTH34])));
    Sterile->SetDelta(1,4, OscParams[kDelta14]);
    Sterile->SetDelta(2,4, OscParams[kDelta24]);

    //Set PMNS parameters for second sterile state
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

    //Set PMNS parameters for third sterile state
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

  // Set Decay parameters
  if(OscProb::PMNS_Decay* Decay = dynamic_cast<OscProb::PMNS_Decay*>(fPMNSObj)) {
    Decay->SetAlpha2(OscParams[kAlpha2]);
    Decay->SetAlpha3(OscParams[kAlpha3]);
  }

  // Set Deco parameters
  if(OscProb::PMNS_Deco* Deco = dynamic_cast<OscProb::PMNS_Deco*>(fPMNSObj)) {
    Deco->SetGamma(2, OscParams[kGamma21]);
    Deco->SetGamma(3, OscParams[kGamma31]);
    Deco->SetDecoAngle(OscParams[kDecoAngle]);
    Deco->SetPower(OscParams[kPower]);
  }

  // Set NSI parameters
  if(OscProb::PMNS_NSI* NSI = dynamic_cast<OscProb::PMNS_NSI*>(fPMNSObj)) {
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

  // Set SNSI parameters
  if(OscProb::PMNS_SNSI* SNSI = dynamic_cast<OscProb::PMNS_SNSI*>(fPMNSObj)) {
    SNSI->SetLowestMass(OscParams[kLightMass]);
  }

  // Set Iter parameters
  if(OscProb::PMNS_Iter* Iter = dynamic_cast<OscProb::PMNS_Iter*>(fPMNSObj)) {
    Iter->SetPrec(OscParams[kPrec]);
  }

  // Set NUNM parameters
  if(OscProb::PMNS_NUNM* NUNM = dynamic_cast<OscProb::PMNS_NUNM*>(fPMNSObj)) {
    NUNM->SetAlpha_11(OscParams[kAlpha11]);
    NUNM->SetAlpha_22(OscParams[kAlpha22]);
    NUNM->SetAlpha_33(OscParams[kAlpha33]);
    NUNM->SetAlpha_21(OscParams[kAlpha21], OscParams[kPhi21]);
    NUNM->SetAlpha_31(OscParams[kAlpha31], OscParams[kPhi31]);
    NUNM->SetAlpha_32(OscParams[kAlpha32], OscParams[kPhi32]);
    NUNM->SetFracVnc(OscParams[kFracVnc]);
  }

  // Set LIV parameters
  if(OscProb::PMNS_LIV* LIV = dynamic_cast<OscProb::PMNS_LIV*>(fPMNSObj)) {
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

}

int OscProbCalcerOscProb::PMNS_StrToInt(std::string PMNSType) {

  if (PMNSType == "Fast" || PMNSType == "fast")          return kFast;
  if (PMNSType == "Sterile+1" || PMNSType == "Sterile1") return kPMNSSterile1;
  if (PMNSType == "Sterile+2" || PMNSType == "Sterile2") return kPMNSSterile2;
  if (PMNSType == "Sterile+3" || PMNSType == "Sterile3") return kPMNSSterile3;
  if (PMNSType == "Decay" || PMNSType == "decay")        return kDecay;
  if (PMNSType == "Deco" || PMNSType == "deco")          return kDeco;
  if (PMNSType == "NSI" || PMNSType == "nsi")            return kNSI;
  if (PMNSType == "LIV" || PMNSType == "liv")            return kLIV;
  if (PMNSType == "NUNM" || PMNSType == "nunm")          return kNUNM;
  if (PMNSType == "SNSI" || PMNSType == "snsi")          return kSNSI;
  if (PMNSType == "Iter" || PMNSType == "iter")          return kIter;

  std::cerr << "Invalid PMNS matrix type provided:" << PMNSType << std::endl;
  throw std::runtime_error("Invalid setup");

  return -1;

}


int OscProbCalcerOscProb::GetNOscParams() {

  if (fOscType == kPMNSSterile1) return kNOscParams_Sterile1;
  if (fOscType == kPMNSSterile2) return kNOscParams_Sterile2;
  if (fOscType == kPMNSSterile3) return kNOscParams_Sterile3;
  if (fOscType == kDecay)        return kNOscParams_Decay;
  if (fOscType == kDeco)         return kNOscParams_Deco;
  if (fOscType == kNSI)          return kNOscParams_NSI;
  if (fOscType == kIter)         return kNOscParams_Iter;
  if (fOscType == kNUNM)         return kNOscParams_NUNM;
  if (fOscType == kLIV)          return kNOscParams_LIV;
  if (fOscType == kSNSI)         return kNOscParams_SNSI;

  return kNOscParams;

}
