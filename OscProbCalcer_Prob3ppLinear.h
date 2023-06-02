#ifndef __OSCILLATOR_PROB3PPLINEAR_H__
#define __OSCILLATOR_PROB3PPLINEAR_H__

#include "OscProbCalcerBase.h"

#include "BargerPropagator.h"

class OscProbCalcerProb3ppLinear : public OscProbCalcerBase {
 public:
  OscProbCalcerProb3ppLinear(int Verbosity_=NONE);

  // ========================================================================================================================================================================
  // Functions which need implementation specific code
  void SetupPropagator();
  void CalculateProbabilities(std::vector<FLOAT_T> OscParams);
  int ReturnWeightArrayIndex(int NuTypeIndex, int InitNuIndex, int FinalNuIndex, int EnergyIndex, int CosineZIndex=-1);
  long DefineWeightArraySize();

  // ========================================================================================================================================================================
  // Functions which help setup implementation specific code

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPATHL, kDENS, kNOscParams};
  enum NuType{Nubar=-1, Nu=1};
  enum NuFlav{Electron=1, Muon=2, Tau=3};

  bool doubled_angle;

  BargerPropagator *bNu;
};

#endif
