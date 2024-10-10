#ifndef __OSCILLATOR_OSCPROB_H__
#define __OSCILLATOR_OSCPROB_H__

#include "OscProbCalcerBase.h"

#include "inc/PMNS_Fast.h"
#include "inc/PMNS_Decay.h"
#include "inc/PMNS_Iter.h"
#include "inc/PremModel.h"

/**
 * @file OscProbCalcer_OscProb.h
 *
 * @class OscProbCalcerOscProb
 *
 * @brief Oscillation calculation engine for linear propagation in NuFAST.
 */
class OscProbCalcerOscProb : public OscProbCalcerBase {
 public:

  /**
   * @brief Default constructor
   *
   * @param Config_ YAML::Node to setup the OscProbCalcerOscProb() instance
   */
  OscProbCalcerOscProb(YAML::Node Config_);

  /**
   * @brief Constructor which takes a file path, creates a YAML::Node and calls the default constructor
   *
   * @param ConfigName_ File path to config
   */
  OscProbCalcerOscProb(std::string ConfigName_) : OscProbCalcerOscProb(YAML::LoadFile(ConfigName_)) {}
  
  /**
   * @brief Destructor
   */
  virtual ~OscProbCalcerOscProb();

  // ========================================================================================================================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup NuFAST specific variables
   */  
  void SetupPropagator() override;

  /**
   * @brief Setup PREM Model with certain number of layers 
   *
   * Files with PREM Models already included in OscProb, used to compute neutrino path and densities
   *
   * @param model parameter to define which model to use. 0 -> default file, 1 -> 15 layers, 2 -> 44 layers, 3 -> 425 layers 
   */
  std::string SetupPREMModel(int model = 2);
  
  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calls particular version of CalcProbPMNS depending on PMNS type
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) override;

  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities with PMNS_Fast object. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray.
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalcProbPMNS_Fast(const std::vector<FLOAT_T>& OscParams);

   /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities with PMNS_Decay object. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray. 
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalcProbPMNS_Decay(const std::vector<FLOAT_T>& OscParams);

    /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities with PMNS_Iter object. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray. 
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalcProbPMNS_Iter(const std::vector<FLOAT_T>& OscParams);

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
  int ReturnWeightArrayIndex(int NuTypeIndex, int OscNuIndex, int EnergyIndex, int CosineZIndex=-1) override;
  
  /**
   * @brief Define the size of fWeightArray
   *
   * This is implementation specific because because NuFAST is setup to calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  long DefineWeightArraySize() override;

  // ========================================================================================================================================================================
  // Functions which help setup implementation specific code

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

  /**
 * @brief Set parameters for PMNS_Fast 
 * 
 * @param Fast object of class PMNS_Fast 
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams(OscProb::PMNS_Fast *Fast, const std::vector<FLOAT_T>& OscParams);

   /**
 * @brief Set parameters for PMNS_Decay
 * 
 * @param Iter object of class PMNS_Decay
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams(OscProb::PMNS_Decay *Decay, const std::vector<FLOAT_T>& OscParams);

  /**
 * @brief Set parameters for PMNS_Iter 
 * 
 * @param Iter object of class PMNS_Iter 
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams(OscProb::PMNS_Iter *Iter, const std::vector<FLOAT_T>& OscParams);

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
   * @brief Definition of oscillation parameters which are expected in this ProbGPU implementation. 
   * Base parameters for standard 3x3 oscillation matrix, works as is for PMNS_Fast class
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kNOscParams};

  /**
   * @brief Definition of extra oscillation parameters for PMNS_Decay class
   * kAlpha2 = m2/tau2, mass and lifetime of the 2nd state in the restframe
   * kAlpha3 = m3/tau3, mass and lifetime of the 3rd state in the restframe
   */
  enum OscParams_Decay{kAlpha2 = kDCP+1, kAlpha3 = kDCP+2};

   /**
   * @brief Definition of extra oscillation parameters for PMNS_Iter class
   * kPrec defines precision of the iterative method
   */
  enum OscParams_Iter{kPrec = kDCP+1};
  
  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nu=1,Nubar=-1};

  /**
  * @brief Different types of PMNS matrices currently supported within the analysis
  */
  enum PMNSMatrix{kFast=0, kPMNSSterile1=1, kPMNSSterile2=2, kPMNSSterile3=3, kDecay=4, kDeco=5, kNSI=6, kIter=7, kLIV=8, kNUNM=9, kSNSI=10};

  /**
   * @brief Define the type for the PMNS matrix
   */
  int fOscType;

private:
  
  OscProb::PremModel PremModel;
  
};

#endif
