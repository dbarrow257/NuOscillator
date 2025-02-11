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
   * @param Config_ YAML::Node to setup the OscProbCalcerNuCUDAProb3() instance
   */
  OscProbCalcerCUDAProb3(YAML::Node Config_);

  /**
   * @brief Constructor which takes a file path, creates a YAML::Node and calls the default constructor
   *
   * @param ConfigName_ File path to config
   */  
  OscProbCalcerCUDAProb3(std::string ConfigName_) : OscProbCalcerCUDAProb3(YAML::LoadFile(ConfigName_)) {}

  /**
   * @brief Destructor
   */
  virtual ~OscProbCalcerCUDAProb3();

 private:
  // ========================================================================================================================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup CUDAProb3 specific variables
   *
   * Setup the cudaprob3::Propagator instance and set all of the variables that it needs like EarthModel etc.
   */
  void SetupPropagator() override;

  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities in CUDAProb3. This links to Propagator->getProbability in CUDAProb3. This function both calculates and stores
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
   * This is implementation specific because because CUDAProb3 is setup to calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  long DefineWeightArraySize() override;
  
  // ========================================================================================================================================================================
  //Functions which help setup implementation specific code

  /**
 * @brief Return the PMNS Matrix type corresponding to a particular string 
 * 
 * @param PMNSType String to convert to enum value
 *
 * @return Enum value describing the PMNS Matrix to use
 */
  int PMNS_StrToInt(std::string PMNSType);

  /**
 * @brief Return number of parameters needed for a particular type of PMNS matrix
 * 
 * @param OscType int value corresponding to type of PMNS matrix
 *
 * @return number of parameters needed for the type of PMNS matrix
 */
  int GetNOscParams(int OscType);

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
  * @brief ...
  */
  enum PMNSMatrix{kStandard=0, k4layers=1};

  /**
   * @brief ...
   */
  int fOscType;

  /**
   * @brief Definition of oscillation parameters which are expected in this CUDAProb3 implementation
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPRODH, kNOscParams};

  /**
   * @brief ...
   */
  enum OscParams_PREM4layers{kBoundLayers12=kPRODH+1, kBoundLayers23=kPRODH+2, kBoundLayers34=kPRODH+3, kBoundLayers45=kPRODH+4, kWeightLayer1=kPRODH+5, kWeightLayer2=kPRODH+6, kWeightLayer3=kPRODH+7, kWeightLayer4=kPRODH+8};

  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nubar=-1, Nu=1};

  /**
   * @brief The name of the config used to setup the particular instance of OscProbCalcerCUDAProb3()
   */
  std::string ConfigName;

  /**
   * @brief The mapping of the oscillation channels defined in #fOscillationChannels to the CUDAProb3 constants
   */
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

  /**
   * @brief  Option to ... OscProbCalcerCUDAProb3()
   */
  bool UseEarthModelSystematics;
};

#endif
