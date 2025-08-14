#ifndef __OSCILLATOR_CUDAPROB3LINEAR_H__
#define __OSCILLATOR_CUDAPROB3LINEAR_H__

#include "OscProbCalcerBase.h"
#include "Constants/OscillatorConstants.h"

#include <memory>

// ==================================================================
/**
 * @brief Includes specific to CUDAProb3 implementation
 */
namespace cudaprob3linear { template<typename T> class Propagator;}
// ==================================================================

/**
 * @file OscProbCalcer_CUDAProb3Linear.h
 *
 * @class OscProbCalcerCUDAProb3Linear
 *
 * @brief Oscillation calculation engine for linear propagation in CUDAProb3.
 */
class OscProbCalcerCUDAProb3Linear : public OscProbCalcerBase {
 public:
  /**
   * @brief Default constructor
   *
   * @param Config_ YAML::Node to setup the OscProbCalcerCUDAProb3Linear() instance
   */
  OscProbCalcerCUDAProb3Linear(YAML::Node Config_);

  /**
   * @brief Constructor which takes a file path, creates a YAML::Node and calls the default constructor
   *
   * @param ConfigName_ File path to config
   */
  OscProbCalcerCUDAProb3Linear(std::string ConfigName_) : OscProbCalcerCUDAProb3Linear(YAML::LoadFile(ConfigName_)) {}

  /**
   * @brief Destructor
   */
  virtual ~OscProbCalcerCUDAProb3Linear();

 private:
  // ========================================================================================================================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup CUDAProb3Linear specific variables
   *
   * Setup the cudaprob3::Propagator instance and set all of the variables that it needs like EarthModel etc.
   */
  void SetupPropagator() override;

  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities in CUDAProb3Linear. This links to Propagator->getProbability in CUDAProb3Linear. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray.
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) override;

  /**
   * @brief Return implementation specific index in the weight array for a specific combination of neutrino oscillation channel, energy and cosine zenith
   * 
   * @param NuTypeIndex The index in #fNeutrinoTypes (neutrino/antinuetrino) to return the pointer for 
   * @param OscChanIndex The index in #fOscillationChannels to return the pointer for 
   * @param EnergyIndex The index in #fEnergyArray to return the pointer for 
   * @param CosineZIndex The index in #fCosineZArray to return the pointer for 
   *
   * @return Index in #fWeightArray which corresponds to the given inputs
   */
  int ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex=-1) override;

  /**
   * @brief Define the size of fWeightArray
   *
   * This is implementation specific because because CUDAProb3Linear is setup to calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  long DefineWeightArraySize() override;
  
  // ========================================================================================================================================================================
  //Functions which help setup implementation specific code

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
   * @brief Definition of oscillation parameters which are expected in this CUDAProb3Linear implementation
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPATHL, kDENS, kNOscParams};

  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nubar=-1, Nu=1};

  /**
   * @brief The mapping of the oscillation channels defined in #fOscillationChannels to the CUDAProb3Linear constants
   */
  std::vector<int> OscChannels;

  /**
   * @brief The number of threads being used to perform the calculation
   */
  int nThreads;

  /**
   * @brief The instance of the CUDAProb3Linear Propagator being used in a particular instance of OscProbCalcerCUDAProb3Linear()
   */
  std::unique_ptr< cudaprob3linear::Propagator< FLOAT_T > > propagator;
};

#endif
