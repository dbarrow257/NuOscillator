#include "OscProbCalcer_GLoBESLinear.h"

#include <iostream>

#include "globes/globes.h"
#include "globes/glb-modules.h"
// KS: Include the SNU header with C linkage to prevent C++ name mangling
// This ensures that the functions in snu.h can be linked correctly when compiled in C++.
extern "C" {
  #include "snu.h"
}

OscProbCalcerGLoBESLinear::OscProbCalcerGLoBESLinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  fNOscParams = kNOscParams;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  fImplementationName += "-CPU-"+std::to_string(1);
  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);
}

OscProbCalcerGLoBESLinear::~OscProbCalcerGLoBESLinear() {

}

void OscProbCalcerGLoBESLinear::SetupPropagator() {
  char name[] = "dummy";
  glbInit(name);
}

void OscProbCalcerGLoBESLinear::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
  // Set the experimental parameters
  const double L = OscParams[kPATHL]; // km
  const double rho = OscParams[kDENS];
  // Set the vacuum oscillation parameters
  const double theta12 = std::asin(std::sqrt(OscParams[kTH12]));
  const double theta13 = std::asin(std::sqrt(OscParams[kTH13]));
  const double theta23 = std::asin(std::sqrt(OscParams[kTH23]));
  const double delta = OscParams[kDCP];
  const double Dmsq21 = OscParams[kDM12];
  const double Dmsq31 = OscParams[kDM23] + OscParams[kDM12]; // eV^2

  // Set GLoBES oscillation parameters
  glb_params true_values = glbAllocParams();
  glbSetOscParams(true_values, theta12, GLB_THETA_12);
  glbSetOscParams(true_values, theta13, GLB_THETA_13);
  glbSetOscParams(true_values, theta23, GLB_THETA_23);
  glbSetOscParams(true_values, delta, GLB_DELTA_CP);
  glbSetOscParams(true_values, Dmsq21, GLB_DM_21);
  glbSetOscParams(true_values, Dmsq31, GLB_DM_31);
  glbSetOscillationParameters(true_values);

  // Loop over energy and oscillation channels
  // KS: WARNING according to manual it isn't thread safe...
  for (int iEnergy = 0; iEnergy < fNEnergyPoints; ++iEnergy) {
    const double E = fEnergyArray[iEnergy];
    for (int iChan = 0; iChan < fNOscillationChannels; ++iChan) {
      for (int iNuType = 0; iNuType < fNNeutrinoTypes; ++iNuType) {
        int alpha = fOscillationChannels[iChan].GeneratedFlavour;
        int beta  = fOscillationChannels[iChan].DetectedFlavour;
        int cp_sign = (iNuType == 0) ? +1 : -1; // neutrino vs antineutrino

        double P = glbConstantDensityProbability(alpha, beta, cp_sign, E, L, rho);

        int index = ReturnWeightArrayIndex(iNuType, iChan, iEnergy, 0);
        fWeightArray[index] = P;
      }
    }
  }
  glbFreeParams(true_values);
}

int OscProbCalcerGLoBESLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerGLoBESLinear::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
