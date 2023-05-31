#include "OscProbCalcer_ProbGPULinear.h"

extern "C" void setMNS(double x12, double x13, double x23, double m21, double m23, double Delta, bool kSquared);
extern "C" void GetProb(int Alpha, int Beta, double Path, double Density, double *Energy, int n, double *oscw);

#include <iostream>

OscProbCalcerProbGPULinear::OscProbCalcerProbGPULinear() : OscProbCalcerBase()
{
  //Base variables
  fNOscParams = kNOscParams;

  //This implementation only considers linear propagation, thus no requirement to set cosineZ array
  fCosineZArraySet = true;

  //Implementation specific variables
  doubled_angle = true;

  nNeutrinoSigns = kNNeutrinoTypes;
  nInitialFlavours = 2; // =2 if excluding taus, =3 if including taus
  nFinalFlavours = 3; // =2 if excluding taus, =3 if including taus

  NeutrinoTypes.resize(nNeutrinoSigns);
  NeutrinoTypes[0] = Neutrino;
  NeutrinoTypes[1] = AntiNeutrino;
}

void OscProbCalcerProbGPULinear::SetupPropagator() {
  //This implementation doesn't really need to do anything in the setup due to probGPU's horrific implementation
  fPropagatorSet = true;  
}

void OscProbCalcerProbGPULinear::CalculateProbabilities(std::vector<FLOAT_T> OscParams) {
  setMNS(OscParams[kTH12], OscParams[kTH23], OscParams[kTH13], OscParams[kDM12], OscParams[kDM23], OscParams[kDCP], doubled_angle);

  // ProbGPULinear calculates oscillation probabilites for each NeutrinoType, so need to copy them from the calculator into fWeightArray
  int CopyArrSize = fNEnergyPoints;
  FLOAT_T* CopyArr = new FLOAT_T[CopyArrSize];

  for (int iNuType=0;iNuType<nNeutrinoSigns;iNuType++) {
    for (int iInitFlav=0;iInitFlav<nInitialFlavours;iInitFlav++) {
      for (int iFinalFlav=0;iFinalFlav<nFinalFlavours;iFinalFlav++) {
	GetProb(NeutrinoTypes[iNuType]*iInitFlav, NeutrinoTypes[iNuType]*iFinalFlav, OscParams[kPATHL], OscParams[kDENS], fEnergyArray.data(), fNEnergyPoints, CopyArr);

        // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
        int IndexToFill = iNuType*nInitialFlavours*nFinalFlavours*CopyArrSize + iInitFlav*nFinalFlavours*CopyArrSize + iFinalFlav*CopyArrSize;
        for (int iOscProb=0;iOscProb<CopyArrSize;iOscProb++) {
          fWeightArray[IndexToFill+iOscProb] = CopyArr[iOscProb];
        }
      }
    }
  }

}

const FLOAT_T* OscProbCalcerProbGPULinear::ReturnPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T Energy, FLOAT_T CosineZ) {
  return NULL;
}

void OscProbCalcerProbGPULinear::IntiailiseWeightArray() {
  int nCalculationPoints = fNEnergyPoints * nInitialFlavours * nFinalFlavours * nNeutrinoSigns;
  std::cout << "Creating weight array with " << nCalculationPoints << " entries" << std::endl;

  fWeightArray = std::vector<FLOAT_T>(nCalculationPoints,DUMMYVAL);
  fWeightArrayInit = true;
}
