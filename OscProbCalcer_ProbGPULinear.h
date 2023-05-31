#ifndef __OSCILLATOR_PROBGPULINEAR_H__
#define __OSCILLATOR_PROBGPULINEAR_H__

#include "OscProbCalcerBase.h"

class OscProbCalcerProbGPULinear : public OscProbCalcerBase {
 public:
  OscProbCalcerProbGPULinear();

  //========================================================================================================================================================================
  // Functions which need implementation specific code
  void SetupPropagator();
  void CalculateProbabilities(std::vector<FLOAT_T> OscParams);
  const FLOAT_T* ReturnPointer(FLOAT_T Energy, FLOAT_T CosineZ);
  void IntiailiseWeightArray();

  //========================================================================================================================================================================
  //Functions which help setup implementation specific code

  //========================================================================================================================================================================
  // Variables which are needed for implementation specific code
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPATHL, kDENS, kNOscParams};

  enum NeutrinoType{AntiNeutrino=-1, Neutrino=1, kNNeutrinoTypes=2};
  int nNeutrinoSigns;
  int nInitialFlavours;
  int nFinalFlavours;
  std::vector<int> NeutrinoTypes;

  bool doubled_angle;
};

#endif
