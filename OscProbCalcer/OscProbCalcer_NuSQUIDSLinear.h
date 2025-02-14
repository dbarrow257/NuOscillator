#ifndef __OSCILLATOR_NUSQUIDSLINEAR_H__
#define __OSCILLATOR_NUSQUIDSLINEAR_H__

#include "OscProbCalcerBase.h"

#include "nuSQuIDS/nuSQuIDS.h"
#include "examples/Decoherence/nuSQUIDSDecoh.h"
#include "examples/NSI/NSI.h"
#include "examples/LV/LV.h"

/**
 * @file OscProbCalcer_NuSQUIDSLinear.h
 *
 * @class OscProbCalcerNuSQUIDSLinear
 *
 * @brief Oscillation calculation engine for linear propagation in NuSQUIDS.
 */
class OscProbCalcerNuSQUIDSLinear : public OscProbCalcerBase {
 public:

  /**
   * @brief Default constructor
   *
   * @param Config_ YAML::Node to setup the OscProbCalcerNuSQUIDSLinear() instance
   */
  OscProbCalcerNuSQUIDSLinear(YAML::Node Config_);

  /**
   * @brief Constructor which takes a file path, creates a YAML::Node and calls the default constructor
   *
   * @param ConfigName_ File path to config
   */
  OscProbCalcerNuSQUIDSLinear(std::string ConfigName_) : OscProbCalcerNuSQUIDSLinear(YAML::LoadFile(ConfigName_)) {}
  
  /**
   * @brief Destructor
   */
  virtual ~OscProbCalcerNuSQUIDSLinear();

  // ========================================================================================================================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup NuSQUIDS specific variables
   */  
  void SetupPropagator() override;
  
  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities in NuSQUIDS. This function both calculates and stores
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
  int ReturnWeightArrayIndex(int NuTypeIndex, int OscNuIndex, int EnergyIndex, int CosineZIndex=-1) override;
  
  /**
   * @brief Define the size of fWeightArray
   *
   * This is implementation specific because because NuSQUIDS is setup to calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  long DefineWeightArraySize() override;

  // ========================================================================================================================================================================
  // Functions which help setup implementation specific code

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
   * @brief Definition of oscillation parameters which are expected in this ProbGPU implementation
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPATHL, kDENS, kELECDENS, kEnergyStrength, kEnergyDep, kEnergyScale, kNOscParams};
 
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
  * @brief Different types of PMNS matrices currently supported within the analysis
  * LIV and SNSI still to be implemented at some point
  */
  enum PMNSMatrix{kFast=0, kPMNSSterile1=1, kPMNSSterile2=2, kPMNSSterile3=3, kDecay=4, kDeco=5, kNSI=6, kIter=7, kNUNM=8, kLIV=9, kSNSI=10, kPMNSSM=11};

  /**
   * @brief Define the type for the PMNS matrix
   */
  int fOscType;

  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nu=1,Nubar=-1};

  nusquids::nuSQUIDS* nus_base;
  nusquids::nuSQUIDS* nubars_base;
  
  nusquids::nuSQUIDS* nus_pmns;
  nusquids::nuSQUIDS* nubars_pmns;
  
  nusquids::nuSQUIDSDecoh* nus_decoh;
  nusquids::nuSQUIDSDecoh* nubars_decoh;

  nuSQUIDSNSI* nus_NSI;
  nuSQUIDSNSI* nubars_NSI;
  
  nusquids::nuSQUIDSLV* nus_LV;
  nusquids::nuSQUIDSLV* nubarss_LV;

  squids::Const units;

  double zenith_angle;
  double integration_step;
  double rel_error;
  double abs_error;
  std::string osc_model;

  std::string decoherence_model;
  nusquids::nuSQUIDSDecoh::DecoherenceModel nusquids_decoherence_model;
};

#endif
