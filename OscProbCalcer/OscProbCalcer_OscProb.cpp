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

OscProbCalcerOscProb::OscProbCalcerOscProb(YAML::Node Config_) : OscProbCalcerBase(Config_), fPMNSObj(nullptr) {
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
  SetOscParams();

  fMaxGenFlavour = 1;
  fMaxDetFlavour = 1;
  for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
    int gflv = fOscillationChannels[iOscChannel].GeneratedFlavour;
    int dflv = fOscillationChannels[iOscChannel].DetectedFlavour;
    if(gflv>fMaxGenFlavour) fMaxGenFlavour = gflv;
    if(dflv>fMaxDetFlavour) fMaxDetFlavour = dflv;
  }

  if (fCosineZIgnored) fImplementationName += "Linear";
  fImplementationName += "-" + PMNSType;

  std::cout << "PMNS Type : " << PMNSType << std::endl;
  std::cout << "Number of parameters : " << ReturnNOscParams() << std::endl;
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

void OscProbCalcerOscProb::CalculateProbabilities() {
  SetPMNSParams();
  CalcProbPMNS();
}

void OscProbCalcerOscProb::SetPath(int iCosineZ) {
  if(fCosineZIgnored) {
    fPMNSObj->SetLength (GetOscillationParameter(ReturnNOscParams() - 3));
    fPMNSObj->SetDensity(GetOscillationParameter(ReturnNOscParams() - 2));
    fPMNSObj->SetZoA    (GetOscillationParameter(ReturnNOscParams() - 1));
  }
  else {
    fPremModel.SetTopLayerSize(GetOscillationParameter(ReturnNOscParams() - 1));
    fPremModel.FillPath(fCosineZArray[iCosineZ]);
    fPMNSObj->SetPath(fPremModel.GetNuPath());
  }
}

void OscProbCalcerOscProb::CalcProbPMNS() {
  for (int iNuType = 0; iNuType < fNNeutrinoTypes; iNuType++) {
    fPMNSObj->SetIsNuBar(fNeutrinoTypes[iNuType]==Nubar);

    for (int iCosineZ = 0; iCosineZ < GetNCosineZ(); iCosineZ++) {
      SetPath(iCosineZ);

      for (int iEnergy = 0; iEnergy < fNEnergyPoints; iEnergy++) {
        OscProb::matrixD probMatrix = fPMNSObj->ProbMatrix(fMaxGenFlavour, fMaxDetFlavour, fEnergyArray[iEnergy]);

        #if UseMultithreading == 1
        #pragma omp simd
        #endif
        for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {

          const int gflv = fOscillationChannels[iOscChannel].GeneratedFlavour-1;
          const int dflv = fOscillationChannels[iOscChannel].DetectedFlavour-1;
          const double weight = probMatrix[gflv][dflv];
          const int index = ReturnWeightArrayIndex(iNuType, iOscChannel, iEnergy, iCosineZ);

          fWeightArray[index] = weight;
        }

      }
    }
  }
}

int OscProbCalcerOscProb::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = ((NuTypeIndex *  fNOscillationChannels + OscChanIndex) * GetNCosineZ() + std::max(CosineZIndex,0)) * fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerOscProb::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * GetNCosineZ() * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}

void OscProbCalcerOscProb::SetPMNSParams() {
  // Set PMNS parameters
  fPMNSObj->SetDm(2, GetOscillationParameter(kDM12));
  fPMNSObj->SetDm(3, GetOscillationParameter(kDM23) + GetOscillationParameter(kDM12));
  fPMNSObj->SetAngle(1,2, asin(sqrt(GetOscillationParameter(kTH12))));
  fPMNSObj->SetAngle(1,3, asin(sqrt(GetOscillationParameter(kTH13))));
  fPMNSObj->SetAngle(2,3, asin(sqrt(GetOscillationParameter(kTH23))));
  fPMNSObj->SetDelta(1,3, GetOscillationParameter(kDCP));

  //Set PMNS parameters for first sterile state
  if(OscProb::PMNS_Sterile* Sterile = dynamic_cast<OscProb::PMNS_Sterile*>(fPMNSObj)) {
    Sterile->SetDm(4, GetOscillationParameter(kDM14));
    Sterile->SetAngle(1,4, asin(sqrt(GetOscillationParameter(kTH14))));
    Sterile->SetAngle(2,4, asin(sqrt(GetOscillationParameter(kTH24))));
    Sterile->SetAngle(3,4, asin(sqrt(GetOscillationParameter(kTH34))));
    Sterile->SetDelta(1,4, GetOscillationParameter(kDelta14));
    Sterile->SetDelta(2,4, GetOscillationParameter(kDelta24));

    //Set PMNS parameters for second sterile state
    if(fOscType == kPMNSSterile2 || fOscType == kPMNSSterile3) {
      Sterile->SetDm(5, GetOscillationParameter(kDM15));
      Sterile->SetAngle(1,5, asin(sqrt(GetOscillationParameter(kTH15))));
      Sterile->SetAngle(2,5, asin(sqrt(GetOscillationParameter(kTH25))));
      Sterile->SetAngle(3,5, asin(sqrt(GetOscillationParameter(kTH35))));
      Sterile->SetAngle(4,5, asin(sqrt(GetOscillationParameter(kTH45))));
      Sterile->SetDelta(1,5, GetOscillationParameter(kDelta15));
      Sterile->SetDelta(2,5, GetOscillationParameter(kDelta25));
      Sterile->SetDelta(3,5, GetOscillationParameter(kDelta35));
    }

    //Set PMNS parameters for third sterile state
    if(fOscType == kPMNSSterile3) {
      Sterile->SetDm(6, GetOscillationParameter(kDM16));
      Sterile->SetAngle(1,6, asin(sqrt(GetOscillationParameter(kTH16))));
      Sterile->SetAngle(2,6, asin(sqrt(GetOscillationParameter(kTH26))));
      Sterile->SetAngle(3,6, asin(sqrt(GetOscillationParameter(kTH36))));
      Sterile->SetAngle(4,6, asin(sqrt(GetOscillationParameter(kTH46))));
      Sterile->SetAngle(5,6, asin(sqrt(GetOscillationParameter(kTH56))));
      Sterile->SetDelta(1,6, GetOscillationParameter(kDelta16));
      Sterile->SetDelta(2,6, GetOscillationParameter(kDelta26));
      Sterile->SetDelta(3,6, GetOscillationParameter(kDelta36));
      Sterile->SetDelta(4,6, GetOscillationParameter(kDelta46));
    }
  }

  // Set Decay parameters
  if(OscProb::PMNS_Decay* Decay = dynamic_cast<OscProb::PMNS_Decay*>(fPMNSObj)) {
    Decay->SetAlpha2(GetOscillationParameter(kAlpha2));
    Decay->SetAlpha3(GetOscillationParameter(kAlpha3));
  }

  // Set Deco parameters
  if(OscProb::PMNS_Deco* Deco = dynamic_cast<OscProb::PMNS_Deco*>(fPMNSObj)) {
    Deco->SetGamma(2, GetOscillationParameter(kGamma21));
    Deco->SetGamma(3, GetOscillationParameter(kGamma31));
    Deco->SetDecoAngle(GetOscillationParameter(kDecoAngle));
    Deco->SetPower(GetOscillationParameter(kPower));
  }

  // Set NSI parameters
  if(OscProb::PMNS_NSI* NSI = dynamic_cast<OscProb::PMNS_NSI*>(fPMNSObj)) {
    NSI->SetNSI(GetOscillationParameter(kEps_ee),
                GetOscillationParameter(kEps_emu),
                GetOscillationParameter(kEps_etau),
                GetOscillationParameter(kEps_mumu),
                GetOscillationParameter(kEps_mutau),
                GetOscillationParameter(kEps_tautau),
                GetOscillationParameter(kDelta_emu),
                GetOscillationParameter(kDelta_etau),
                GetOscillationParameter(kDelta_mutau));

    NSI->SetFermCoup(GetOscillationParameter(kElecCoup),
                     GetOscillationParameter(kUpCoup),
                     GetOscillationParameter(kDownCoup));
  }

  // Set SNSI parameters
  if(OscProb::PMNS_SNSI* SNSI = dynamic_cast<OscProb::PMNS_SNSI*>(fPMNSObj)) {
    SNSI->SetLowestMass(GetOscillationParameter(kLightMass));
  }

  // Set Iter parameters
  if(OscProb::PMNS_Iter* Iter = dynamic_cast<OscProb::PMNS_Iter*>(fPMNSObj)) {
    Iter->SetPrec(GetOscillationParameter(kPrec));
  }

  // Set NUNM parameters
  if(OscProb::PMNS_NUNM* NUNM = dynamic_cast<OscProb::PMNS_NUNM*>(fPMNSObj)) {
    NUNM->SetAlpha_11(GetOscillationParameter(kAlpha11));
    NUNM->SetAlpha_22(GetOscillationParameter(kAlpha22));
    NUNM->SetAlpha_33(GetOscillationParameter(kAlpha33));
    NUNM->SetAlpha_21(GetOscillationParameter(kAlpha21), GetOscillationParameter(kPhi21));
    NUNM->SetAlpha_31(GetOscillationParameter(kAlpha31), GetOscillationParameter(kPhi31));
    NUNM->SetAlpha_32(GetOscillationParameter(kAlpha32), GetOscillationParameter(kPhi32));
    NUNM->SetFracVnc(GetOscillationParameter(kFracVnc));
  }

  // Set LIV parameters
  if(OscProb::PMNS_LIV* LIV = dynamic_cast<OscProb::PMNS_LIV*>(fPMNSObj)) {
    LIV->SetaT(0, 0, 3, GetOscillationParameter(kaT_ee_3), 0.);
    LIV->SetaT(0, 1, 3, GetOscillationParameter(kaT_emu_3), GetOscillationParameter(kDelta_emu_3));
    LIV->SetaT(0, 2, 3, GetOscillationParameter(kaT_etau_3), GetOscillationParameter(kDelta_etau_3));
    LIV->SetaT(1, 1, 3, GetOscillationParameter(kaT_mumu_3), 0.);
    LIV->SetaT(1, 2, 3, GetOscillationParameter(kaT_mutau_3), GetOscillationParameter(kDelta_mutau_3));
    LIV->SetaT(2, 2, 3, GetOscillationParameter(kaT_tautau_3), 0.);
    LIV->SetcT(0, 0, 4, GetOscillationParameter(kcT_ee_4), 0.);
    LIV->SetcT(0, 1, 4, GetOscillationParameter(kcT_emu_4), GetOscillationParameter(kDelta_emu_4));
    LIV->SetcT(0, 2, 4, GetOscillationParameter(kcT_etau_4), GetOscillationParameter(kDelta_etau_4));
    LIV->SetcT(1, 1, 4, GetOscillationParameter(kcT_mumu_4), 0.);
    LIV->SetcT(1, 2, 4, GetOscillationParameter(kcT_mutau_4), GetOscillationParameter(kDelta_mutau_4));
    LIV->SetcT(2, 2, 4, GetOscillationParameter(kcT_tautau_4), 0.);
    LIV->SetaT(0, 0, 5, GetOscillationParameter(kaT_ee_5), 0.);
    LIV->SetaT(0, 1, 5, GetOscillationParameter(kaT_emu_5), GetOscillationParameter(kDelta_emu_5));
    LIV->SetaT(0, 2, 5, GetOscillationParameter(kaT_etau_5), GetOscillationParameter(kDelta_etau_5));
    LIV->SetaT(1, 1, 5, GetOscillationParameter(kaT_mumu_5), 0.);
    LIV->SetaT(1, 2, 5, GetOscillationParameter(kaT_mutau_5), GetOscillationParameter(kDelta_mutau_5));
    LIV->SetaT(2, 2, 5, GetOscillationParameter(kaT_tautau_5), 0.);
    LIV->SetcT(0, 0, 6, GetOscillationParameter(kcT_ee_6), 0.);
    LIV->SetcT(0, 1, 6, GetOscillationParameter(kcT_emu_6), GetOscillationParameter(kDelta_emu_6));
    LIV->SetcT(0, 2, 6, GetOscillationParameter(kcT_etau_6), GetOscillationParameter(kDelta_etau_6));
    LIV->SetcT(1, 1, 6, GetOscillationParameter(kcT_mumu_6), 0.);
    LIV->SetcT(1, 2, 6, GetOscillationParameter(kcT_mutau_6), GetOscillationParameter(kDelta_mutau_6));
    LIV->SetcT(2, 2, 6, GetOscillationParameter(kcT_tautau_6), 0.);
    LIV->SetaT(0, 0, 7, GetOscillationParameter(kaT_ee_7), 0.);
    LIV->SetaT(0, 1, 7, GetOscillationParameter(kaT_emu_7), GetOscillationParameter(kDelta_emu_7));
    LIV->SetaT(0, 2, 7, GetOscillationParameter(kaT_etau_7), GetOscillationParameter(kDelta_etau_7));
    LIV->SetaT(1, 1, 7, GetOscillationParameter(kaT_mumu_7), 0.);
    LIV->SetaT(1, 2, 7, GetOscillationParameter(kaT_mutau_7), GetOscillationParameter(kDelta_mutau_7));
    LIV->SetaT(2, 2, 7, GetOscillationParameter(kaT_tautau_7), 0.);
    LIV->SetcT(0, 0, 8, GetOscillationParameter(kcT_ee_8), 0.);
    LIV->SetcT(0, 1, 8, GetOscillationParameter(kcT_emu_8), GetOscillationParameter(kDelta_emu_8));
    LIV->SetcT(0, 2, 8, GetOscillationParameter(kcT_etau_8), GetOscillationParameter(kDelta_etau_8));
    LIV->SetcT(1, 1, 8, GetOscillationParameter(kcT_mumu_8), 0.);
    LIV->SetcT(1, 2, 8, GetOscillationParameter(kcT_mutau_8), GetOscillationParameter(kDelta_mutau_8));
    LIV->SetcT(2, 2, 8, GetOscillationParameter(kcT_tautau_8), 0.);
  }

}

int OscProbCalcerOscProb::PMNS_StrToInt(const std::string& PMNSType) {
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


void OscProbCalcerOscProb::SetOscParams() {
  std::vector<std::string> OscParNames;  

  //DB TODO CORRECTLY IMPLEMENT THE BELOW
  
  switch (fOscType) {
  case kFast:
    OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};
    break;
  case kPMNSSterile1:
    OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};
    break;
  case kPMNSSterile2:
    OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};
    break;
  case kPMNSSterile3:
    OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};
    break;
  case kDecay:
    OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};    
    break;
  case kDeco:
    OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};    
    break;
  case kNSI:
    OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp","eps_ee", "eps_emu", "eps_etau", "eps_mumu", "eps_mutau", "eps_tautau", "delta_emu", "delta_etau", "delta_mutau", "elec_coup", "up_coup", "down_coup"};    
    break;
  case kIter:
    OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};    
    break;
  case kNUNM:
    OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};    
    break;
  case kLIV:
    OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};    
    break;
  case kSNSI:
    OscParNames = {"sin2_th12","sin2_th23","sin2_th13","dm2_12","dm2_23","delta_cp"};    
    break;
  default:
    std::cerr << "Invalid Oscillation Model type provided:" << fOscType << std::endl;
    throw std::runtime_error("Invalid Oscillation Model");
  }
  
  // This is needed to treat neutrino path parameters due to how these have
  // been used so far with classes with less flexibility
  if (fCosineZIgnored){
    OscParNames.push_back("path_length");
    OscParNames.push_back("matter_density");
    OscParNames.push_back("electron_density");
  } else {
    OscParNames.push_back("production_height");
  }

  SetExpectedParameterNames(OscParNames);
}
