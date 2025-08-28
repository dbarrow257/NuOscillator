
#ifndef __nuTens_H__
#define __nuTens_H__

#include "OscProbCalcerBase.h"

#include <nuTens/tensors/tensor.hpp>
#include <nuTens/propagator/const-density-solver.hpp>
#include <nuTens/propagator/DP-propagator.hpp>
#include <nuTens/propagator/pmns-matrix.hpp>

//#include "Constants/OscillatorConstants.h"

using namespace nuTens;

#include <memory>

/**
 * @file OscProbCalcer_nuTens.h
 *
 * @class OscProbCalcernuTens
 *
 * @brief Oscillation calculation engine for  propagation in nuTens.
 */
class OscProbCalcernuTens : public OscProbCalcerBase {
 public:
  /**
   * @brief Default constructor
   *
   * @param Config_ YAML::Node to setup the OscProbCalcernuTens() instance
   */
  OscProbCalcernuTens(YAML::Node Config_);

  /**
   * @brief Constructor which takes a file path, creates a YAML::Node and calls the default constructor
   *
   * @param ConfigName_ File path to config
   */
  OscProbCalcernuTens(std::string ConfigName_) : OscProbCalcernuTens(YAML::LoadFile(ConfigName_)) {}

  /**
   * @brief Destructor
   */
  virtual ~OscProbCalcernuTens();

  void SetEnergyArray(std::vector<FLOAT_T> EnergyArray) override;

 private:
  // ========================================================================================================================================================================
  // Functions which need implementation specific code

  /**
   * @brief Setup nuTens specific variables
   *
   * Setup the nuTens::Propagator instance and set all of the variables that it needs like EarthModel etc.
   */
  void SetupPropagator() override;

  inline void setParamValues(float theta12, float theta13, float theta23, float dcp, float dm12, float dm23) {

    float dm13 = dm12 + dm23;
  
    _theta12Tensor.setValue(theta12, 0, 0);
    _theta13Tensor.setValue(theta13, 0, 0);
    _theta23Tensor.setValue(theta23, 0, 0);
    _dcpTensor.setValue(dcp, 0, 0);
    _dm21Tensor.setValue(-dm12, 0, 0);
    _dm31Tensor.setValue(-dm13, 0, 0);

  }

  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator oscillation probabilities in nuTens. This links to Propagator->getProbability in nuTens. This function both calculates and stores
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
   * This is implementation specific because because nuTens is setup to calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  long DefineWeightArraySize() override;
  
  // ========================================================================================================================================================================
  //Functions which help setup implementation specific code

  // ========================================================================================================================================================================
  // Variables which are needed for implementation specific code

  /**
   * @brief Definition of oscillation parameters which are expected in this nuTens implementation
   */
  enum OscParams_PMNS{kTH12, kTH23, kTH13, kDM12, kDM23, kDCP, kPATHL, kDENS, kELECDENS, kNOscParams_PMNS};

  /**
   * @brief Define the neutrino and antineutrino values expected by this implementation
   */
  enum NuType{Nubar=-1, Nu=1};

  /**
   * @brief The mapping of the oscillation channels defined in #fOscillationChannels to the nuTens constants
   */
  std::vector<int> OscChannels;

  /**
   * @brief The number of threads being used to perform the calculation
   */
  int nThreads;

  nuTens::PMNSmatrix pmnsMatrix;
  nuTens::Tensor energiesTensor;
  nuTens::AccessedTensor<float, 2, dtypes::kCPU> _theta12Tensor = AccessedTensor<float, 2, dtypes::kCPU>::zeros({1, 3}, false);
  nuTens::AccessedTensor<float, 2, dtypes::kCPU> _theta13Tensor = AccessedTensor<float, 2, dtypes::kCPU>::zeros({1, 3}, false);
  nuTens::AccessedTensor<float, 2, dtypes::kCPU> _theta23Tensor = AccessedTensor<float, 2, dtypes::kCPU>::zeros({1, 3}, false);
  nuTens::AccessedTensor<float, 2, dtypes::kCPU> _dcpTensor = AccessedTensor<float, 2, dtypes::kCPU>::zeros({1, 3}, false);
  nuTens::AccessedTensor<float, 2, dtypes::kCPU> _dm21Tensor = AccessedTensor<float, 2, dtypes::kCPU>::zeros({1, 3}, false);
  nuTens::AccessedTensor<float, 2, dtypes::kCPU> _dm31Tensor = AccessedTensor<float, 2, dtypes::kCPU>::zeros({1, 3}, false);

  nuTens::DPpropagator tensorPropagator = DPpropagator(100.0 * units::km, false, 2.6, 4);

};

#endif
