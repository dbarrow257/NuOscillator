#ifndef __OSCILLATOR_CHIC_H__
#define __OSCILLATOR_CHIC_H__

#include "OscProbCalcerBase.h"

class CHIC;
class CHICEARTH;
/**
 * @file OscProbCalcer_CHICLinear.h
 *
 * @class OscProbCalcerCHIC
 *
 * @brief Oscillation calculation engine for linear propagation in CHIC.
 */
class OscProbCalcerCHIC : public OscProbCalcerBase {
 public:

  /**
   * @brief Default constructor
   *
   * @param Config_ YAML::Node to setup the OscProbCalcerCHIC() instance
   */
  OscProbCalcerCHIC(YAML::Node Config_);

  /**
   * @brief Constructor which takes a file path, creates a YAML::Node and calls the default constructor
   *
   * @param ConfigName_ File path to config
   */
  OscProbCalcerCHIC(std::string ConfigName_) : OscProbCalcerCHIC(YAML::LoadFile(ConfigName_)) {}
  
  /**
   * @brief Destructor
   */
  virtual ~OscProbCalcerCHIC();

  // ========================================================================================================================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup CHIC specific variables
   */  
  void SetupPropagator() final;
  
  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities in CHIC. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray.
   */
  void CalculateProbabilities() final;

  /**
   * @brief Calculate oscillation probabilities using Beam
   * */
  void CalculateProbabilitiesBeam();

  /**
   * @brief Calculate oscillation probabilities using ATM
   * */
  void CalculateProbabilitiesAtm();

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
  int ReturnWeightArrayIndex(int NuTypeIndex, int OscNuIndex, int EnergyIndex, int CosineZIndex=-1) final;
  
  /**
   * @brief Define the size of fWeightArray
   *
   * This is implementation specific because because CHIC is setup to calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  long DefineWeightArraySize() final;

  /**
   * @brief Auxiliary function to handle ignored cosineZ cases
   */
  int GetNCosineZ();

  /**
   * @brief CHIC doesn't have cool inheritance so we need to use template to avoid copy pasting
   */
  template <typename Propagator>
  void SetPMNSParameters(Propagator* p);

  // ========================================================================================================================================================================
  // Functions which help setup implementation specific code

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
   * @brief Definition of oscillation parameters which are expected in this CHIC implementation
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kNOscParams};
  
  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nu=1,Nubar=-1};

  /**
   * @brief CHIC propagator for neutrinos
   */
  std::unique_ptr<CHIC> chic_nu;

  /**
   * @brief CHIC propagator for antineutrinos
   */
  std::unique_ptr<CHIC> chic_nubar;

  /**
   * @brief CHICEARTH propagator for neutrinos (atmospheric)
   */
  std::unique_ptr<CHICEARTH> chicearth_nu;

  /**
   * @brief CHICEARTH propagator for antineutrinos (atmospheric)
   */
  std::unique_ptr<CHICEARTH> chicearth_nubar;


  /**
   * @brief Name of earth model implemented in CHIC
   */
  std::string fPremName;

  /**
   * @brief Double storing the detector depth in km
   */
  double fDetDepth;
};

#endif
