#ifndef __OSCILLATOR_CUDAPROB3_H__
#define __OSCILLATOR_CUDAPROB3_H__

#include "OscProbCalcerBase.h"

//==================================================================
// Includes specific to CUDAProb3 implementation
#include <memory>
namespace cudaprob3 { template<typename T> class Propagator;}
//==================================================================

class OscProbCalcerCUDAProb3 : public OscProbCalcerBase {
 public:
  OscProbCalcerCUDAProb3();

 private:
  //========================================================================================================================================================================
  // Functions which need implementation specific code
  void SetupPropagator();
  void CalculateProbabilities(std::vector<FLOAT_T> OscParams);
  const FLOAT_T* ReturnPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T Energy, FLOAT_T CosineZ);
  void IntiailiseWeightArray();
  
  //========================================================================================================================================================================
  //Functions which help setup implementation specific code

  //========================================================================================================================================================================
  // Variables which are needed for implementation specific code
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPRODH, kNOscParams};

  int nNeutrinoSigns;
  int nInitialFlavours;
  int nFinalFlavours;

  std::vector< std::vector<int> > OscChannels;
  std::vector<int> NeutrinoTypes;

  int nThreads;

  std::unique_ptr< cudaprob3::Propagator< FLOAT_T > > propagator;
  std::string EarthDensityFile;
};

#endif
