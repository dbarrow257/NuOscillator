
#ifndef __nuTens_H__
#define __nuTens_H__

#include "OscProbCalcerBase.h"

#include <nuTens/tensors/tensor.hpp>
#include <nuTens/propagator/const-density-solver.hpp>
#include <nuTens/propagator/propagator.hpp>

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

  inline nuTens::Tensor getPMNSmatrix(float theta12, float theta23, float theta13, float dm12, float dm23, float dcp) {

    
    constexpr std::complex<float> imagUnit(0.0, 1.0);
    float m1, m2, m3;

    if(dm23 >= 0.0) {
        m3 = std::sqrt(dm23 + dm12);
        m2 = std::sqrt(dm12);
        m1 = 0.0;
    }
    else if(dm23 < 0.0) {
        m2 = std::sqrt(-dm23);
        m1 = std::sqrt(-dm23 + dm12);
        m3 = 0.0;
    }

    _masses.setValue(m1, 0, 0);
    _masses.setValue(m2, 0, 1);
    _masses.setValue(m3, 0, 2);

    _mat1.setValue({0, 0, 0}, 1.0);
    _mat1.setValue({0, 1, 1}, std::cos(theta23));
    _mat1.setValue({0, 1, 2}, std::sin(theta23));
    _mat1.setValue({0, 2, 1}, -std::sin(theta23));
    _mat1.setValue({0, 2, 2}, std::cos(theta23));

    _mat2.setValue({0, 1, 1}, 1.0);
    _mat2.setValue({0, 0, 0}, std::cos(theta13));
    _mat2.setValue({0, 0, 2}, std::sin(theta13) * std::exp(-imagUnit * dcp));
    _mat2.setValue({0, 2, 0}, -std::sin(theta13) * std::exp(imagUnit * dcp));
    _mat2.setValue({0, 2, 2}, std::cos(theta13));

    _mat3.setValue({0, 2, 2}, 1.0);
    _mat3.setValue({0, 0, 0}, std::cos(theta12));
    _mat3.setValue({0, 0, 1}, std::sin(theta12));
    _mat3.setValue({0, 1, 0}, -std::sin(theta12));
    _mat3.setValue({0, 1, 1}, std::cos(theta12));

    return nuTens::Tensor::matmul(_mat1, Tensor::matmul(_mat2, _mat3));
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


  nuTens::Tensor _mat1;
  nuTens::Tensor _mat2;
  nuTens::Tensor _mat3;

  nuTens::AccessedTensor<float, 2, dtypes::kCPU> _masses = AccessedTensor<float, 2, dtypes::kCPU>::zeros({1, 3}, false);

  nuTens::Tensor energiesTensor;

  nuTens::Propagator tensorPropagator = Propagator(3, 295 * units::km);
  std::shared_ptr<nuTens::ConstDensityMatterSolver> matterSolver = std::make_shared<nuTens::ConstDensityMatterSolver>(3, 2.7);

};

#endif
