#include "OscProbCalcer_OscProb.h"

#include "inc/PremModel.h"

#include "TMath.h"

OscProbCalcerOscProb::OscProbCalcerOscProb(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config

  //=======

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  if (!Config_["OscProbCalcerSetup"]["PMNSType"]) {
    std::cerr << "Expected to find a 'PMNSType' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  std::string OscMatrix = Config_["OscProbCalcerSetup"]["PMNSType"].as<std::string>();
  fOscType = PMNS_StrToInt(OscMatrix);

  std::cout << "PMNS Type : " << fOscType << std::endl;

  fNOscParams = GetNOscParams(fOscType);

  std::cout << "Number of parameters : " << fNOscParams << std::endl;

}

OscProbCalcerOscProb::~OscProbCalcerOscProb() {
}

void OscProbCalcerOscProb::SetupPropagator() {
}

std::string OscProbCalcerOscProb::SetupPREMModel(int model) {

  if (model > 3 || model < 0) std::cout << "Invalid prem model" << std::endl;

    // Get the model table paths
    std::string filename;

    switch (model) {
        case 0:
            filename = "./build/_deps/oscprob-src/PremTables/prem_default.txt";
            break;
        case 1:
            filename = "./build/_deps/oscprob-src/PremTables/prem_15layers.txt";
            break;
        case 2:
            filename = "./build/_deps/oscprob-src/PremTables/prem_44layers.txt";
            break;
        case 3:
            filename = "./build/_deps/oscprob-src/PremTables/prem_425layers.txt";
            break;
    }

    return filename;
}

void OscProbCalcerOscProb::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {

  fWeightArray.resize(DefineWeightArraySize());
  std::cout << fNEnergyPoints << " " << fNCosineZPoints << " " << fWeightArray.size() << std::endl;


  switch(fOscType) {
    case kFast: 
      CalcProbPMNS_Fast(OscParams);
      break;

    case kDecay:
      CalcProbPMNS_Decay(OscParams);

    case kIter:
      CalcProbPMNS_Iter(OscParams);
      break;

    default:
      break;
  }
  //std::cout << "OscProb Probability: Energy = " << Energy << " , Baseline:" << Baseline << " , Prob (mu->mu):" << myPMNS.Prob(1, 1, Energy, Baseline) << std::endl;

}

void OscProbCalcerOscProb::CalcProbPMNS_Fast(const std::vector<FLOAT_T>& OscParams) {

  int prem_model;

  double energy, cosZ;
  double weight;
  int index;

  prem_model = OscParams[kPREM];
  std::string premfile = SetupPREMModel(prem_model);
  OscProb::PremModel prem(premfile);

  std::cout << "Fast with PREM model : " << prem_model << std::endl;

  OscProb::PMNS_Fast myPMNS;

  SetPMNSParams(&myPMNS, OscParams);
  

  for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

    cosZ = fCosineZArray[iCosineZ];

    prem.FillPath(cosZ);

    myPMNS.SetPath(prem.GetNuPath());
    
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
  //std::cout << "OscProb Probability: Energy = " << Energy << " , Baseline:" << Baseline << " , Prob (mu->mu):" << myPMNS.Prob(1, 1, Energy, Baseline) << std::endl;

}

void OscProbCalcerOscProb::CalcProbPMNS_Decay(const std::vector<FLOAT_T>& OscParams) {

  int prem_model;

  double energy, cosZ;
  double weight;
  int index;

  prem_model = OscParams[kPREM];
  std::string premfile = SetupPREMModel(prem_model);
  OscProb::PremModel prem(premfile);

  std::cout << "Decay with PREM model : " << prem_model << std::endl;

  OscProb::PMNS_Decay myPMNS;

  SetPMNSParams(&myPMNS, OscParams);


  for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

    cosZ = fCosineZArray[iCosineZ];

    prem.FillPath(cosZ);

    myPMNS.SetPath(prem.GetNuPath());
    
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
  //std::cout << "OscProb Probability: Energy = " << Energy << " , Baseline:" << Baseline << " , Prob (mu->mu):" << myPMNS.Prob(1, 1, Energy, Baseline) << std::endl;

}

void OscProbCalcerOscProb::CalcProbPMNS_Iter(const std::vector<FLOAT_T>& OscParams) {

  int prem_model;

  double energy, cosZ;
  double weight;
  int index;

  prem_model = OscParams[kPREM];
  std::string premfile = SetupPREMModel(prem_model);
  OscProb::PremModel prem(premfile);

  std::cout << "Iter with PREM model : " << prem_model << std::endl;

  OscProb::PMNS_Iter myPMNS;

  SetPMNSParams(&myPMNS, OscParams);


  for (int iCosineZ = 0; iCosineZ < fNCosineZPoints; iCosineZ++) {

    cosZ = fCosineZArray[iCosineZ];

    prem.FillPath(cosZ);

    myPMNS.SetPath(prem.GetNuPath());
    
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
  //std::cout << "OscProb Probability: Energy = " << Energy << " , Baseline:" << Baseline << " , Prob (mu->mu):" << myPMNS.Prob(1, 1, Energy, Baseline) << std::endl;

}

int OscProbCalcerOscProb::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNCosineZPoints*fNEnergyPoints + OscChanIndex*fNCosineZPoints*fNEnergyPoints + CosineZIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerOscProb::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNCosineZPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}

void OscProbCalcerOscProb::SetPMNSParams(OscProb::PMNS_Fast *Fast, const std::vector<FLOAT_T>& OscParams) {
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

void OscProbCalcerOscProb::SetPMNSParams(OscProb::PMNS_Decay *Decay, const std::vector<FLOAT_T>& OscParams) {
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

void OscProbCalcerOscProb::SetPMNSParams(OscProb::PMNS_Iter *Iter, const std::vector<FLOAT_T>& OscParams) {
  double th12, th13, th23;
  double dm21, dm31;
  double dcp;
  double prec;

  prec = OscParams[kPrec];
  std::cout << "Precision for Iter PMNS class : " << prec << std::endl;
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
  if(OscType == kDecay) {
    return kNOscParams+2;
  }
  else if (OscType == kIter) {
    return kNOscParams+1;
  }
  else {
    std::cerr << "Invalid PMNS matrix type provided:" << OscType << std::endl;
    throw;
  }

  return -1;
}
