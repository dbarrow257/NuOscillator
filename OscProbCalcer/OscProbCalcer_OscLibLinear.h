#ifndef __OSCILLATOR_OSCLIB_H__
#define __OSCILLATOR_OSCLIB_H__

#include "OscProbCalcerBase.h"

#include "OscLib/OscCalcPMNS.h"
#include "OscLib/OscCalcPMNS_NSI.h"

/**
 * @file OscProbCalcer_OscLibLinear.h
 *
 * @class OscProbCalcerOscLibLinear
 *
 * @brief Oscillation calculation engine for linear propagation in OscLibLinear.
 */
class OscProbCalcerOscLibLinear : public OscProbCalcerBase {
 public:

  /**
   * @brief Default constructor
   *
   * @param Config_ YAML::Node to setup the OscProbCalcerOscLibLinear() instance
   */
  OscProbCalcerOscLibLinear(YAML::Node Config_);

  /**
   * @brief Constructor which takes a file path, creates a YAML::Node and calls the default constructor
   *
   * @param ConfigName_ File path to config
   */
  OscProbCalcerOscLibLinear(std::string ConfigName_) : OscProbCalcerOscLibLinear(YAML::LoadFile(ConfigName_)) {}

  /**
   * @brief Destructor
   */
  virtual ~OscProbCalcerOscLibLinear();

  // ==========================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup Earth model
   */
  void SetupPropagator() final;

  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calls particular version of CalcProbPMNS depending on PMNS type
   */
  void CalculateProbabilities() final;

  /**
   * @brief Return implementation specific index in the weight array for a specific combination of neutrino oscillation channel, energy and cosine zenith
   *
   * @param NuTypeIndex The index in #fNeutrinoTypes (neutrino/antineutrino) to return the pointer for
   * @param OscChanIndex The index in #fOscillationChannels to return the pointer for
   * @param EnergyIndex The index in #fEnergyArray to return the pointer for
   * @param CosineZIndex The index in #fCosineZArray to return the pointer for
   *
   * @return Index in #fWeightArray which corresponds to the given inputs
   */
  int ReturnWeightArrayIndex(int NuTypeIndex, int OscNuIndex, int EnergyIndex, int CosineZIndex=-1) final;

  /**
   * @brief Define the size of fWeightArray
   *
   * This is implementation specific because because OscLibLinear is setup to
   * calculate all oscillation channels together, whilst others calculate
   * only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  long DefineWeightArraySize() final;

  /**
   * @brief Return the PMNS Matrix type corresponding to a particular string
   *
   * @param PMNSType String to convert to enum value
   *
   * @return Enum value describing the PMNS Matrix to use
   */
  int PMNS_StrToInt(const std::string& PMNSType);

  /**
   * @brief Set number of parameters needed for a particular type of PMNS matrix
   */
  void SetOscParams();

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code
  
  /**
   * @brief Definition of oscillation parameters which are expected in this OscLibLinear implementation
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPATHL, kDENS, kNOscParams};

  /**
   * @brief Definition of extra oscillation parameters for PMNS_NSI class
   * kEps quantity representing the intensity of NSI wrt standard matter effects between the different flavours
   * kDelta complex phases between the different flavours (non diagonal elements
   */
  enum OscParams_NSI{kEps_ee = kNOscParams, kEps_emu, kEps_etau, kEps_mumu,
    kEps_mutau, kEps_tautau, kDelta_emu, kDelta_etau, kDelta_mutau, kNOscParams_NSI};

  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nubar=-1, Nu=1};


  /**
   * @brief Different types of PMNS matrices currently supported within the analysis
   */
  enum PMNSMatrix{kPMNS, kNSI};

  /**
   * @brief OscLib object used to calculate the oscillation probabilities
   */
  osc::_IOscCalcAdjustable<FLOAT_T>* OscLib;

  /**
   * @brief Define the type for the PMNS matrix
   */
  int fOscType;
};

#endif
