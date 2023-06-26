#ifndef __OSCILLATOR_PROB3PPLINEAR_H__
#define __OSCILLATOR_PROB3PPLINEAR_H__

#include "OscProbCalcerBase.h"

#include "BargerPropagator.h"

/**
 * @file OscProbCalcer_Prob3ppLinear.h
 *
 * @class OscProbCalcerProb3ppLinear
 *
 * @brief Oscillation calculation engine for linear propagation in Prob3pp.
 */
class OscProbCalcerProb3ppLinear : public OscProbCalcerBase {
 public:

  /**
   * @brief Default constructor
   *
   * @param ConfigName_ Name of config used to setup the OscProbCalcerProb3ppLinear() instance
   */
  OscProbCalcerProb3ppLinear(std::string ConfigName_="");

  // ========================================================================================================================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup Prob3pp specific variables like the BargerPropagator #bNu
   */  
  void SetupPropagator();
  
  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities in Prob3pp. This links to bNu->GetProb in Prob3pp. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray.
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  void CalculateProbabilities(std::vector<FLOAT_T> OscParams);

  /**
   * @brief Return implementation specific index in the weight array for a specific combination of neutrino oscillation channel, energy and cosine zenith
   * 
   * @param NuTypeIndex The index in #fNeutrinoTypes (neutrino/antinuetrino) to return the pointer for 
   * @param OscNuIndex
   * @param EnergyIndex The index in #fEnergyArray to return the pointer for 
   * @param CosineZIndex The index in #fCosineZArray to return the pointer for 
   *
   * @return Index in #fWeightArray which corresponds to the given inputs
   */
  int ReturnWeightArrayIndex(int NuTypeIndex, int OscNuIndex, int EnergyIndex, int CosineZIndex=-1);
  
  /**
   * @brief Define the size of fWeightArray
   *
   * This is implementation specific because because Prob3pp is setup to calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  long DefineWeightArraySize();

  // ========================================================================================================================================================================
  // Functions which help setup implementation specific code

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
   * @brief Definition of oscillation parameters which are expected in this ProbGPU implementation
   */
  enum OscParams{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPATHL, kDENS, kNOscParams};

  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nubar=-1, Nu=1};

  /**
   * @brief Boolean declaring what values are being passed for the values of theta_13 (sin^2(theta) or sin^2(2*theta))
   */
  bool doubled_angle;

  /**
   * @brief BargerPropagator used within the Prob3pp calculation framework
   */
  BargerPropagator *bNu;
};

#endif
