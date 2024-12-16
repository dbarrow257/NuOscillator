#ifndef __OSCILLATOR_OSCPROB_H__
#define __OSCILLATOR_OSCPROB_H__

#include "OscProbCalcerBase.h"

#include "inc/PMNS_Fast.h"
#include "inc/PMNS_Sterile.h"
#include "inc/PMNS_Decay.h"
#include "inc/PMNS_Deco.h"
#include "inc/PMNS_NSI.h"
#include "inc/PMNS_Iter.h"
#include "inc/PMNS_NUNM.h"
#include "inc/PMNS_SNSI.h"
#include "inc/PMNS_LIV.h"
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
   * @brief Setup Earth model
   */  
  void SetupPropagator() override;
  
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
   * Calculator oscillation probabilities with PMNS_Sterile object for up to 3 additonal neutrino states. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray.
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   * @param neutrino_number Total number of neutrino states to consider
   */
  void CalcProbPMNS_Sterile(const std::vector<FLOAT_T>& OscParams, int neutrino_number=4);

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
   * Calculator oscillation probabilities with PMNS_Deco object. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray. 
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalcProbPMNS_Deco(const std::vector<FLOAT_T>& OscParams);

  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities with PMNS_NSI object. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray. 
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalcProbPMNS_NSI(const std::vector<FLOAT_T>& OscParams);

  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities with PMNS_SNSI object. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray. 
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalcProbPMNS_SNSI(const std::vector<FLOAT_T>& OscParams);

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
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities with PMNS_NUNM object. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray. 
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalcProbPMNS_NUNM(const std::vector<FLOAT_T>& OscParams);

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
  void SetPMNSParams_Fast(OscProb::PMNS_Fast *Fast, const std::vector<FLOAT_T>& OscParams);

  /**
 * @brief Set parameters for PMNS_Sterile with up to 3 additional neutrino states
 * 
 * @param Sterile object of class PMNS_Sterile
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_Sterile(OscProb::PMNS_Sterile *Sterile, const std::vector<FLOAT_T>& OscParams);

   /**
 * @brief Set parameters for PMNS_Decay
 * 
 * @param Decay object of class PMNS_Decay
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_Decay(OscProb::PMNS_Decay *Decay, const std::vector<FLOAT_T>& OscParams);

  /**
 * @brief Set parameters for PMNS_Deco
 * 
 * @param Deco object of class PMNS_Deco
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_Deco(OscProb::PMNS_Deco *Deco, const std::vector<FLOAT_T>& OscParams);

  /**
 * @brief Set parameters for PMNS_NSI
 * 
 * @param NSI object of class PMNS_NSI
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_NSI(OscProb::PMNS_NSI *NSI, const std::vector<FLOAT_T>& OscParams);

   /**
 * @brief Set parameters for PMNS_SNSI
 * 
 * @param SNSI object of class PMNS_SNSI
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_SNSI(OscProb::PMNS_SNSI *SNSI, const std::vector<FLOAT_T>& OscParams);

  /**
 * @brief Set parameters for PMNS_Iter 
 * 
 * @param Iter object of class PMNS_Iter 
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_Iter(OscProb::PMNS_Iter *Iter, const std::vector<FLOAT_T>& OscParams);

  /**
 * @brief Set parameters for PMNS_NUNM
 * 
 * @param NUNM object of class PMNS_NUNM
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_NUNM(OscProb::PMNS_NUNM *NUNM, const std::vector<FLOAT_T>& OscParams);

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
   * @brief Definition of oscillation parameters which are expected in this ProbGPU implementation. 
   * Base parameters for standard 3x3 oscillation matrix, works as is for PMNS_Fast class
   * DetDepth is the detector depth and is expected in km
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kDetDepth, kNOscParams};

  /**
   * @brief Definition of extra oscillation parameters for PMNS_Sterile class with up to 3 additional neutrinos
   * kTHij -> mixing angle between states i and j
   * kDeltaij -> CP violating phase between states i and j
   * kDM1j -> mass squared difference between states 1 and j in eV^2
   */
  enum OscParams_Sterile{kTH14=kDetDepth+1, kTH24=kDetDepth+2, kTH34=kDetDepth+3, kDM14=kDetDepth+4, kDelta14=kDetDepth+5, kDelta24=kDetDepth+6,
                         kTH15=kDetDepth+7, kTH25=kDetDepth+8, kTH35=kDetDepth+9, kTH45=kDetDepth+10, kDM15=kDetDepth+11, kDelta15=kDetDepth+12, kDelta25=kDetDepth+13, kDelta35=kDetDepth+14,
                         kTH16=kDetDepth+15, kTH26=kDetDepth+16, kTH36=kDetDepth+17, kTH46=kDetDepth+18, kTH56=kDetDepth+19, kDM16=kDetDepth+20, 
                         kDelta16=kDetDepth+21, kDelta26=kDetDepth+22, kDelta36=kDetDepth+23, kDelta46=kDetDepth+24};

  /**
   * @brief Definition of extra oscillation parameters for PMNS_Decay class
   * kAlpha2 = m2/tau2, mass and lifetime of the 2nd state in the restframe, must be positive and units are eV^2
   * kAlpha3 = m3/tau3, mass and lifetime of the 3rd state in the restframe, must be positive and units are eV^2
   */
  enum OscParams_Decay{kAlpha2 = kDetDepth+1, kAlpha3 = kDetDepth+2};

    /**
   * @brief Definition of extra oscillation parameters for PMNS_Deco class
   * kGamma21 decoherence parameter between mass states 1 and 2 
   * kGamma31 decoherence parameter between mass states 3 and 1
   * kDecoAngle decoherence angle
   * kPower power index of decoherence energy dependence
   */
  enum OscParams_Deco{kGamma21 = kDetDepth+1, kGamma31 = kDetDepth+2, kDecoAngle = kDetDepth+3, kPower = kDetDepth+4};

   /**
   * @brief Definition of extra oscillation parameters for PMNS_NSI class
   * kEps quantity reprensenting the intensity of NSI wrt standard matter effects between the different flavours 
   * kDelta complex phases between the different flavours (non diag elements)
   * kElecCoup electron coupling
   * kUpCoup u-quark coupling
   * kDownCoup d-quark coupling
   */
  enum OscParams_NSI{kEps_ee = kDetDepth+1, kEps_emu = kDetDepth+2, kEps_etau = kDetDepth+3, kEps_mumu = kDetDepth+4, kEps_mutau = kDetDepth+5, kEps_tautau= kDetDepth+6,
                     kDelta_emu = kDetDepth+7, kDelta_etau = kDetDepth+8, kDelta_mutau = kDetDepth+9, kElecCoup = kDetDepth+10, kUpCoup = kDetDepth+11, kDownCoup = kDetDepth+12};

  /**
   * @brief Definition of extra oscillation parameters for PMNS_SNSI class (scalar non standard interactions)
   * Takes class PMNS_NSI as base with all parameters associated to it
   * kLightMass mass of lightest neutrino
   * WARNING : kEps expected in units of 1/mass^2 in MeV^-2 and not dimensionless like for PMNS_NSI
   */
  enum OscParams_SNSI{kLightMass = kDownCoup+1};

   /**
   * @brief Definition of extra oscillation parameters for PMNS_Iter class
   * kPrec defines precision of the iterative method
   */
  enum OscParams_Iter{kPrec = kDetDepth+1};

  /**
   * @brief Definition of extra oscillation parameters for PMNS_NUNM class (non unitary neutrino mixing)
   * kAlphaij quantify deviation from unitary wrt standard mixing between states i and j
   * kPhiij complex phases for non diagonal elements between states i and j 
   * kFracVnc fraction of matter potential affecting NC
   */
  enum OscParams_NUNM{kAlpha11 = kDetDepth+1, kAlpha21 = kDetDepth+2, kAlpha31 = kDetDepth+3, kAlpha22 = kDetDepth+4, kAlpha32 = kDetDepth+5, kAlpha33 = kDetDepth+6,
                      kPhi21 = kDetDepth+7, kPhi31 = kDetDepth+8, kPhi32 = kDetDepth+9, kFracVnc = kDetDepth+10};
  
  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nu=1,Nubar=-1};

  /**
  * @brief Different types of PMNS matrices currently supported within the analysis
  * LIV and SNSI still to be implemented at some point
  */
  enum PMNSMatrix{kFast=0, kPMNSSterile1=1, kPMNSSterile2=2, kPMNSSterile3=3, kDecay=4, kDeco=5, kNSI=6, kIter=7, kNUNM=8, kLIV=9, kSNSI=10};

  /**
   * @brief Define the type for the PMNS matrix
   */
  int fOscType;

private:
  
  /**
   * @brief Object handling the Earth model and the estimation of neutrino paths through the Earth depending on their direction
   * Paths used by PMNS objects for neutrino propagation
   * Earth model implemented as succession of spherical shells with uniform density 
   */
  OscProb::PremModel PremModel;

  /**
   * @brief String storing the path of the density table file used to setup the Earth model
   */
  std::string premfile;
  
};

#endif
