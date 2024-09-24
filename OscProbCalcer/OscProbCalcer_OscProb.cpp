#include "OscProbCalcer_OscProb.h"

#include "inc/PMNS_Fast.h"
#include "TMath.h"

OscProbCalcerOscProb::OscProbCalcerOscProb(std::string ConfigName_, int Instance_) : OscProbCalcerBase(ConfigName_,"OscProb",Instance_)
{
  //=======
  //Grab information from the config

  //=======

  fNOscParams = kNOscParams;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  OscProb::PMNS_Fast myPMNS;
  
  int mh = 1;
  double dm21 = 7.5e-5;
  double dm31 = mh>0 ? 2.457e-3 : -2.449e-3 + dm21;
  double th12 = asin(sqrt(0.304));
  double th13 = asin(sqrt(mh>0 ? 0.0218 : 0.0219));
  double th23 = asin(sqrt(mh>0 ? 0.452 : 0.579));
  double dcp  = (mh>0 ? 306 : 254)*TMath::Pi()/180;
  
  // Set PMNS parameters
  myPMNS.SetDm(2, dm21);
  myPMNS.SetDm(3, dm31);
  myPMNS.SetAngle(1,2, th12);
  myPMNS.SetAngle(1,3, th13);
  myPMNS.SetAngle(2,3, th23);
  myPMNS.SetDelta(1,3, dcp);

  double Energy = 1.5; //GeV
  double Baseline = 25.0; //km
  
  std::cout << "OscProb Probability: Energy = " << Energy << " , Baseline:" << Baseline << " , Prob (mu->mu):" << myPMNS.Prob(1, 1, Energy, Baseline) << std::endl;
}

OscProbCalcerOscProb::~OscProbCalcerOscProb() {
}

void OscProbCalcerOscProb::SetupPropagator() {
}

void OscProbCalcerOscProb::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
}

int OscProbCalcerOscProb::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNCosineZPoints*fNEnergyPoints + OscChanIndex*fNCosineZPoints*fNEnergyPoints + CosineZIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerOscProb::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNCosineZPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
