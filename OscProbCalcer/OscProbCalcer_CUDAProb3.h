#ifndef __OSCILLATOR_CUDAPROB3_H__
#define __OSCILLATOR_CUDAPROB3_H__

#include "OscProbCalcerBase.h"

#include <memory>

// ==================================================================
/**
 * @brief Includes specific to CUDAProb3 implementation
 */
namespace cudaprob3 { template<typename T> class Propagator;}
// ==================================================================

/**
 * @file OscProbCalcer_CUDAProb3.h
 *
 * @class OscProbCalcerCUDAProb3
 *
 * @brief Oscillation calculation engine for linear and atmospheric propagation in CUDAProb3.
 */
class OscProbCalcerCUDAProb3 : public OscProbCalcerBase {
 public:
    /**
   * @brief Default constructor
   *
   * @param ConfigName_ Name of config used to setup the OscProbCalcerCUDAProb3() instance
   */
  //DB
  OscProbCalcerCUDAProb3(std::string ConfigName_="", int Instance_=0);

 private:
  // ========================================================================================================================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup CUDAProb3 specific variables
   *
   * Setup the cudaprob3::Propagator instance and set all of the variables that it needs like EarthModel etc.
   */
  void SetupPropagator();

  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities in CUDAProb3. This links to Propagator->getProbability in CUDAProb3. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray.
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalculateProbabilities(std::vector<FLOAT_T> OscParams);

  /**
   * @brief Return implementation specific index in the weight array for a specific combination of neutrino oscillation channel, energy and cosine zenith
   * 
   * @param NuTypeIndex The index in #fNeutrinoTypes (neutrino/antinuetrino) to return the pointer for 
   * @param OscChanIndex
   * @param EnergyIndex The index in #fEnergyArray to return the pointer for 
   * @param CosineZIndex The index in #fCosineZArray to return the pointer for 
   *
   * @return Index in #fWeightArray which corresponds to the given inputs
   */
  //DB
  int ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex=-1);

  /**
   * @brief Define the size of fWeightArray
   *
   * This is implementation specific because because CUDAProb3 is setup to calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  long DefineWeightArraySize();
  
  // ========================================================================================================================================================================
  //Functions which help setup implementation specific code

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
   * @brief Definition of oscillation parameters which are expected in this CUDAProb3 implementation
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPRODH, kNOscParams};

  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nubar=-1, Nu=1};

  /**
   * @brief The name of the config used to setup the particular instance of OscProbCalcerCUDAProb3()
   */
  std::string ConfigName;

  /**
   * @brief The mapping of the oscillation channels defined in #fInitialFlavours and #fFinalFlavours to the CUDAProb3 constants
   */
  //DB Check for all instances of fInitialFlavours and fFinalFlavours
  std::vector<int> OscChannels;

  /**
   * @brief The number of threads being used to perform the calculation
   */
  int nThreads;

  /**
   * @brief The instance of the CUDAProb3 Propagator being used in a particular instance of OscProbCalcerCUDAProb3()
   */
  std::unique_ptr< cudaprob3::Propagator< FLOAT_T > > propagator;

  /**
   * @brief The name of the Earth Density file being used in a particular instance of OscProbCalcerCUDAProb3()
   */
  std::string EarthDensityFile;
};

#endif
