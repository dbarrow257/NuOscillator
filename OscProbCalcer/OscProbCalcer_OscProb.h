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
   * Calculator oscillation probabilities with any PMNS object. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray.
   *
   * @param myPMNS The PMNS object to compute
   */
  void CalcProbPMNS(OscProb::PMNS_Base* myPMNS);

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
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities with PMNS_LIV object. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray.
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalcProbPMNS_LIV(const std::vector<FLOAT_T>& OscParams);

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
 * @brief Set parameters for PMNS_Base
 *
 * @param PMNS object of any class
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams(OscProb::PMNS_Base* pmns, const std::vector<FLOAT_T>& OscParams);

  /**
 * @brief Set parameters for PMNS_Sterile with up to 3 additional neutrino states
 *
 * @param Sterile object of class PMNS_Sterile
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_Sterile(OscProb::PMNS_Sterile* Sterile, const std::vector<FLOAT_T>& OscParams);

   /**
 * @brief Set parameters for PMNS_Decay
 *
 * @param Decay object of class PMNS_Decay
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_Decay(OscProb::PMNS_Decay* Decay, const std::vector<FLOAT_T>& OscParams);

  /**
 * @brief Set parameters for PMNS_Deco
 *
 * @param Deco object of class PMNS_Deco
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_Deco(OscProb::PMNS_Deco* Deco, const std::vector<FLOAT_T>& OscParams);

  /**
 * @brief Set parameters for PMNS_NSI
 *
 * @param NSI object of class PMNS_NSI
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_NSI(OscProb::PMNS_NSI* NSI, const std::vector<FLOAT_T>& OscParams);

   /**
 * @brief Set parameters for PMNS_SNSI
 *
 * @param SNSI object of class PMNS_SNSI
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_SNSI(OscProb::PMNS_SNSI* SNSI, const std::vector<FLOAT_T>& OscParams);

  /**
 * @brief Set parameters for PMNS_Iter
 *
 * @param Iter object of class PMNS_Iter
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_Iter(OscProb::PMNS_Iter* Iter, const std::vector<FLOAT_T>& OscParams);

  /**
 * @brief Set parameters for PMNS_NUNM
 *
 * @param NUNM object of class PMNS_NUNM
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_NUNM(OscProb::PMNS_NUNM* NUNM, const std::vector<FLOAT_T>& OscParams);

    /**
 * @brief Set parameters for PMNS_LIV
 *
 * @param LIV object of class PMNS_LIV
 * @param OscParams  The parameter set to calculate oscillation probabilities at
 *
 * @return Sets relevant parameters for PMNS Matrix
 */
  void SetPMNSParams_LIV(OscProb::PMNS_LIV* LIV, const std::vector<FLOAT_T>& OscParams);

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
   * @brief Definition of oscillation parameters which are expected in this ProbGPU implementation.
   * Base parameters for standard 3x3 oscillation matrix, works as is for PMNS_Fast class
   * DetDepth is the detector depth and is expected in km
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kNOscParams};

  /**
   * @brief Definition of extra oscillation parameters for PMNS_Sterile class with up to 3 additional neutrinos
   * kTHij -> mixing angle between states i and j
   * kDeltaij -> CP violating phase between states i and j
   * kDM1j -> mass squared difference between states 1 and j in eV^2
   */
  enum OscParams_Sterile{kTH14 = kNOscParams, kTH24, kTH34, kDM14,
                         kDelta14, kDelta24,
                         kTH15, kTH25, kTH35, kTH45, kDM15,
                         kDelta15, kDelta25, kDelta35,
                         kTH16, kTH26, kTH36, kTH46, kTH56, kDM16,
                         kDelta16, kDelta26, kDelta36, kDelta46,
                         kNOscParams_Sterile};

  /**
   * @brief Definition of extra oscillation parameters for PMNS_Decay class
   * kAlpha2 = m2/tau2, mass and lifetime of the 2nd state in the restframe, must be positive and units are eV^2
   * kAlpha3 = m3/tau3, mass and lifetime of the 3rd state in the restframe, must be positive and units are eV^2
   */
  enum OscParams_Decay{kAlpha2 = kNOscParams, kAlpha3, kNOscParams_Decay};

    /**
   * @brief Definition of extra oscillation parameters for PMNS_Deco class
   * kGamma21 decoherence parameter between mass states 1 and 2
   * kGamma31 decoherence parameter between mass states 3 and 1
   * kDecoAngle decoherence angle
   * kPower power index of decoherence energy dependence
   */
  enum OscParams_Deco{kGamma21 = kNOscParams, kGamma31, kDecoAngle, kPower, kNOscParams_Deco};

   /**
   * @brief Definition of extra oscillation parameters for PMNS_NSI class
   * kEps quantity reprensenting the intensity of NSI wrt standard matter effects between the different flavours
   * kDelta complex phases between the different flavours (non diag elements
   * kElecCoup electron coupling
   * kUpCoup u-quark coupling
   * kDownCoup d-quark coupling
   */
  enum OscParams_NSI{kEps_ee = kNOscParams, kEps_emu, kEps_etau, kEps_mumu,
                     kEps_mutau, kEps_tautau, kDelta_emu, kDelta_etau,
                     kDelta_mutau, kElecCoup, kUpCoup, kDownCoup,
                     kNOscParams_NSI};

  /**
   * @brief Definition of extra oscillation parameters for PMNS_SNSI class (scalar non standard interactions)
   * Takes class PMNS_NSI as base with all parameters associated to it
   * kLightMass mass of lightest neutrino
   * WARNING : kEps expected in units of 1/mass^2 in MeV^-2 and not dimensionless like for PMNS_NSI
   */
  enum OscParams_SNSI{kLightMass = kNOscParams_NSI, kNOscParams_SNSI};

   /**
   * @brief Definition of extra oscillation parameters for PMNS_Iter class
   * kPrec defines precision of the iterative method
   */
  enum OscParams_Iter{kPrec = kNOscParams, kNOscParams_Iter};

  /**
   * @brief Definition of extra oscillation parameters for PMNS_NUNM class (non unitary neutrino mixing)
   * kAlphaij quantify deviation from unitary wrt standard mixing between states i and j
   * kPhiij complex phases for non diagonal elements between states i and j
   * kFracVnc fraction of matter potential affecting NC
   */
  enum OscParams_NUNM{kAlpha11 = kNOscParams, kAlpha21, kAlpha31, kAlpha22, kAlpha32, kAlpha33,
                      kPhi21, kPhi31, kPhi32, kFracVnc, kNOscParams_NUNM};

  /**
   * @brief Definition of extra oscillation parameters for PMNS_LIV class (Lorentz Invariance Violation)
   * kaT_j LIV coefficient of dimension j between 2 neutrino flavors with j = 3, 5 or 7
   * kcT_j LIV coefficient of dimension j between 2 neutrino flavors with j = 4, 6 or 8
   * kDelta_j complex phase between 2 neutrino flavors for non diag element with dimension j
   */
  enum OscParams_LIV{kaT_ee_3 = kNOscParams, kaT_emu_3, kaT_etau_3, kaT_mumu_3, kaT_mutau_3,
                     kaT_tautau_3, kDelta_emu_3, kDelta_etau_3, kDelta_mutau_3,
                     kcT_ee_4, kcT_emu_4, kcT_etau_4, kcT_mumu_4, kcT_mutau_4,
                     kcT_tautau_4, kDelta_emu_4, kDelta_etau_4, kDelta_mutau_4,
                     kaT_ee_5, kaT_emu_5, kaT_etau_5, kaT_mumu_5, kaT_mutau_5,
                     kaT_tautau_5, kDelta_emu_5, kDelta_etau_5, kDelta_mutau_5,
                     kcT_ee_6, kcT_emu_6, kcT_etau_6, kcT_mumu_6, kcT_mutau_6,
                     kcT_tautau_6, kDelta_emu_6, kDelta_etau_6, kDelta_mutau_6,
                     kaT_ee_7, kaT_emu_7, kaT_etau_7, kaT_mumu_7, kaT_mutau_7,
                     kaT_tautau_7, kDelta_emu_7, kDelta_etau_7, kDelta_mutau_7,
                     kcT_ee_8, kcT_emu_8, kcT_etau_8, kcT_mumu_8, kcT_mutau_8,
                     kcT_tautau_8, kDelta_emu_8, kDelta_etau_8, kDelta_mutau_8,
                     kNOscParams_LIV};

  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nu=1,Nubar=-1};

  /**
  * @brief Different types of PMNS matrices currently supported within the analysis
  * LIV and SNSI still to be implemented at some point
  */
  enum PMNSMatrix{kFast, kPMNSSterile1, kPMNSSterile2, kPMNSSterile3,
                  kDecay, kDeco, kNSI, kSNSI, kIter, kNUNM, kLIV};

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

  /**
   * @brief Double storing the detector depth in km
   */
  double fDetDepth;

  /**
   * @brief Maximum generated flavour for ProbMatrix calculation
   */
  int fMaxGenFlavour;

  /**
   * @brief Maximum detected flavour for ProbMatrix calculation
   */
  int fMaxDetFlavour;

};

#endif
