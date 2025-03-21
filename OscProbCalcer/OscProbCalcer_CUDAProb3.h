#ifndef __OSCILLATOR_CUDAPROB3_H__
#define __OSCILLATOR_CUDAPROB3_H__

#include "OscProbCalcerBase.h"
#include "TH2D.h"
#include "TH3D.h"
#include "TFile.h"

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
  * @brief Set the variables needed to use the production heights averaging.
  */
  void SetProductionHeightsAveraging();
  
  /**
   * @brief Apply a new set of parameters to set the density model of the Earth in CUDAProb3
   *
   * @param OscParams The full parameter set to calculate oscillation probabilities at 
   */
  void ApplyEarthModelSystematics(const std::vector<FLOAT_T>& OscParams);
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
   * @brief  Option to use Production Heights averaging in a particular instance of OscProbCalcerCUDAProb3()
   */
  bool UseProductionHeightsAve;
  
  /**
   * @brief The name of the Production Heights file being used in a particular instance of OscProbCalcerCUDAProb3()
   *
   * The file should contain a collection of TH3Ds named "ProductionHeight_" + flavour suffixes, to be defined in [OscProbCalcerSetup][ProductionHeightsHistFlavourSuffixes].
   * The binning in x (Energy) and y (CosineZ) must be the same specified in the config file under [Binned][FileName].
   * The heights used for the average correspond to the lower edges of the z-axis bins.
   */
  std::string ProductionHeightsFile;

   /**
   * @brief The flavour suffixes strings of the TH3Ds in ProductionHeightsFile, read from the configuration file under [OscProbCalcerSetup][ProductionHeightsHistFlavourSuffixes] .
   */
  std::vector<std::string> ProductionHeightsHistFlavourSuffixes;

  /**
   * @brief  Option to apply density model systematics in OscProbCalcerCUDAProb3()
   */
  bool UseEarthModelSystematics;

  /**
   * @brief  Number of Earth layers
   */
  int nLayers;

  /**
   * @brief Size of the array to be copied.
   */
  int CopyArrSize;

  /**
   * @brief Pointer to the array used for copying oscillation probabilities.
   */
  FLOAT_T* CopyArr;

};

#endif
