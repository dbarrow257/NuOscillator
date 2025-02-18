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
   * @brief Definition of Standard Model oscillation parameters (same as BSM Non-Standard Interactions parameters) which are expected in this ProbGPU implementation
   */
  enum OscParams_PMNS{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPATHL, kDENS, kELECDENS, kNOscParams_PMNS};

  /**
   * @brief Definition of Standard Model oscillation parameters and additional parameters characterising the Beyond-the-Standard-Model physics (Decoherence)
   */
  enum OscParams_Decoh{kEnergyStrength=kNOscParams_PMNS, kEnergyDep, kEnergyScale, kNOscParams_Decoh};
  
  /**
   * @brief Definition of Standard Model oscillation parameters and additional parameters characterising the Beyond-the-Standard-Model physics (Lorentz-Violation Invariance)
   */
  enum OscParams_LIV{kEMuReal=kNOscParams_PMNS, kEMuImg, kMuTauReal, kMuTauImg, kEnergyPower, kNOscParams_LIV};

  /**
   * @brief Return the PMNS Matrix type corresponding to a particular string 
   * 
   * @param PMNSType String to convert to enum value
   *
   * @return Enum value describing the PMNS Matrix to use
   */
  int PMNS_StrToInt(std::string OscModel);
  
  /**
   * @brief Return number of parameters needed for a particular type of PMNS matrix
   * 
   * @param OscType int value corresponding to type of PMNS matrix
   *
   * @return number of parameters needed for the type of PMNS matrix
   */
  int GetNOscParams(int OscType);
 
  /**
   * @brief Set enums corresponding to BSM model 
   * 
   * @param BSMModel String to convert to enum value
   *
   * @return Enum value describing the BSMModel to use
   */
  enum OscModels{kSM=0, kDecoherence=1, kLIV=2, kNSI=3};
  
  /**
   * @brief Define the type for the PMNS matrix
   */
  int fOscModel;

  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nu=1,Nubar=-1};

  /**
   * @brief Define a Const object for handling the units
   */
  squids::Const units;
  
  /**
   * @brief Declaration of the NuSQUIDS base class object for neutrinos
   */
  nusquids::nuSQUIDS* nus_base;

  /**
   * @brief Declaration of the NuSQUIDS base class object for anti-neutrinos
   */
  nusquids::nuSQUIDS* nubars_base;

  /**
   * @brief Declaration of the NuSQUIDS SM (PMNS) derived class object for neutrinos
   */
  nusquids::nuSQUIDS* nus_PMNS;
  /**
   * @brief Declaration of the NuSQUIDS SM (PMNS) derived class object for anti-neutrinos
   */
  nusquids::nuSQUIDS* nubars_PMNS;
  
  /**
   * @brief Declaration of the NuSQUIDS BSM (Decoherence) derived class object for neutrinos
   */
  nusquids::nuSQUIDSDecoh* nus_decoh;

  /**
   * @brief Declaration of the NuSQUIDS BSM (Decoherence) derived class object for anti-neutrinos
   */
  nusquids::nuSQUIDSDecoh* nubars_decoh;

  /**
   * @brief Declaration of the string to choose the decoherence model
   */
  std::string decoherence_model;

  /**
   * @brief Declaration of the NuSQUIDS BSM (Decoherence) model object
   */
  nusquids::nuSQUIDSDecoh::DecoherenceModel nusquids_decoherence_model;

  /**
   * @brief Declaration of the NuSQUIDS BSM (Lorentz-Violation Invariance) derived class object for neutrinos
   */
  nusquids::nuSQUIDSLV* nus_LIV;
  /**
   * @brief Declaration of the NuSQUIDS BSM (Lorentz-Violation Invariance) derived class object for anti-neutrinos
   */
  nusquids::nuSQUIDSLV* nubars_LIV;

  /**
   * @brief Declaration of the NuSQUIDS BSM (Non-Standard Interactions) derived class object for neutrinos
   */
  nuSQUIDSNSI* nus_NSI;

  /**
   * @brief Declaration of the NuSQUIDS BSM (Non-Standard Interactions) derived class object for anti-neutrinos
   */
  nuSQUIDSNSI* nubars_NSI;

  /**
   * @brief Declaration of the muon-neutrino to tau-neutrino coupling for the NSI computation
   */
  FLOAT_T nsi_mutau_coupling;

  /**
   * @brief Declaration of the integration step
   */
  double integration_step;

  /**
   * @brief Declaration of the relative error
   */
  double rel_error;

  /**
   * @brief Declaration of the absolute error
   */
  double abs_error;

  /**
   * @brief Number of neutrino flavours considered in the analysis
   */
  int nNeutrinoFlavours;
};

#endif
