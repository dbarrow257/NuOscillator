#ifndef __OSCILLATOR_CUDAPROB3_H__
#define __OSCILLATOR_CUDAPROB3_H__

#include "OscProbCalcerBase.h"

// ==================================================================
// Includes specific to CUDAProb3 implementation
#include <memory>
namespace cudaprob3 { template<typename T> class Propagator;}
// ==================================================================

class OscProbCalcerCUDAProb3 : public OscProbCalcerBase {
 public:
  OscProbCalcerCUDAProb3(std::string ConfigName_="", int Verbosity_=NONE);

 private:
  // ========================================================================================================================================================================
  // Functions which need implementation specific code
  void SetupPropagator();
  void CalculateProbabilities(std::vector<FLOAT_T> OscParams);
  int ReturnWeightArrayIndex(int NuTypeIndex, int InitNuIndex, int FinalNuIndex, int EnergyIndex, int CosineZIndex=-1);
  long DefineWeightArraySize();
  
  // ========================================================================================================================================================================
  //Functions which help setup implementation specific code

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPRODH, kNOscParams};
  enum NuType{Nubar=-1, Nu=1};
  enum NuFlav{Electron=1, Muon=2, Tau=3};

  std::string ConfigName;

  std::vector< std::vector<int> > OscChannels;
  int nThreads;
  std::unique_ptr< cudaprob3::Propagator< FLOAT_T > > propagator;
  std::string EarthDensityFile;
};

#endif
