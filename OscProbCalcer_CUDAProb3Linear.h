#ifndef __OSCILLATOR_CUDAPROB3LINEAR_H__
#define __OSCILLATOR_CUDAPROB3LINEAR_H__

#include "OscProbCalcerBase.h"
#include "OscillatorConstants.h"

#include <memory>

// ==================================================================
/**
 * @brief Includes specific to CUDAProb3 implementation
 */
namespace cudaprob3 { template<typename T> class Propagator;}
// ==================================================================

/**
 * @file OscProbCalcer_CUDAProb3Linear.h
 *
 * @class OscProbCalcerCUDAProb3Linear
 *
 * @brief Oscillation calculation engine for linear and atmospheric propagation in CUDAProb3.
 */
class OscProbCalcerCUDAProb3Linear : public OscProbCalcerBase {
 public:
    /**
   * @brief Default constructor
   *
   * @param ConfigName_ Name of config used to setup the OscProbCalcerCUDAProb3Linear() instance
   * @param Verbosity_ Verbosity of console output
   */
  OscProbCalcerCUDAProb3Linear(std::string ConfigName_="", int Verbosity_=NONE);

 private:
  // ========================================================================================================================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup CUDAProb3Linear specific variables
   *
   * Setup the cudaprob3::Propagator instance and set all of the variables that it needs like EarthModel etc.
   */
  void SetupPropagator();

  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities in CUDAProb3Linear. This links to Propagator->getProbability in CUDAProb3Linear. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray.
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalculateProbabilities(std::vector<FLOAT_T> OscParams);

  /**
   * @brief Return implementation specific index in the weight array for a specific combination of neutrino oscillation channel, energy and cosine zenith
   * 
   * @param NuTypeIndex The index in #fNeutrinoTypes (neutrino/antinuetrino) to return the pointer for 
   * @param InitNuIndex The index in #fInitialFlavours (electron/muon/tau) to return the pointer for 
   * @param FinalNuIndex The index in #fFinalFlavours (electron/muon/tau) to return the pointer for 
   * @param EnergyIndex The index in #fEnergyArray to return the pointer for 
   * @param CosineZIndex The index in #fCosineZArray to return the pointer for 
   *
   * @return Index in #fWeightArray which corresponds to the given inputs
   */
  int ReturnWeightArrayIndex(int NuTypeIndex, int InitNuIndex, int FinalNuIndex, int EnergyIndex, int CosineZIndex=-1);

  /**
   * @brief Define the size of fWeightArray
   *
   * This is implementation specific because because CUDAProb3Linear is setup to calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  long DefineWeightArraySize();
  
  // ========================================================================================================================================================================
  //Functions which help setup implementation specific code

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
   * @brief Definition of oscillation parameters which are expected in this ProbGPU implementation
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPATHL, kDENS, kNOscParams};

  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nubar=-1, Nu=1};

  /**
   * @brief Define the neutrino flavours expected by this implementation
   */
  enum NuFlav{Electron=1, Muon=2, Tau=3};

  /**
   * @brief The name of the config used to setup the particular instance of OscProbCalcerCUDAProb3Linear()
   */
  std::string ConfigName;

  /**
   * @brief The mapping of the oscillation channels defined in #fInitialFlavours and #fFinalFlavours to the CUDAProb3Linear constants
   */
  std::vector< std::vector<int> > OscChannels;

  /**
   * @brief The number of threads being used to perform the calculation
   */
  int nThreads;

  /**
   * @brief The instance of the CUDAProb3Linear Propagator being used in a particular instance of OscProbCalcerCUDAProb3Linear()
   */
  std::unique_ptr< cudaprob3::Propagator< FLOAT_T > > propagator;
};

#endif