#ifndef __OSCILLATOR_UNBINNED_BASE_H__
#define __OSCILLATOR_UNBINNED_BASE_H__

#include "OscillatorBase.h"

/**
 * @file OscillatorUnbinned.h
 *
 * @class OscillatorUnbinned
 *
 * @brief Unbinned Oscillation calculation implementation class.
 *
 * Implementation of OscillatorBase::OscillatorBase() object which uses unbinned energy and cosineZ dimensions. It is expected that the specific Energy and CosineZ
 * values to be used by the OscProbCalcerBase::OscProbCalcerBase() object are passed via SetEnergyArray() and SetCosineZArray()
 */
class OscillatorUnbinned : public OscillatorBase {
 public:

  /**
   * @brief Default constructor
   *
   * Default constructor
   */
  //DB
  OscillatorUnbinned(std::string ConfigName_);

  // ========================================================================================================================================================================
  // Public functions which are calculation implementation agnostic

    /**
   * @brief Return a pointer to the oscillation probability for the requested event attributes.
   *
   * Determine the memory address address where the calculated oscillation probability for events of the specific requested type will be stored. This will be different
   * depending on the calculation implementation. For the unbinned approach, the particular bin in which the requested energy and cosine falls will be passed directly to
   * the OscProbCalcerBase::OscProbCalcerBase() object
   *
   * @param InitNuFlav Initial neutrino flavour of the neutrino
   * @param FinalNuFlav Final neutrino flavour of the neutrino
   * @param EnergyVal True energy of the neutrino
   * @param CosineZVal True direction of the neutrino in CosineZ
   *
   * @return Pointer to the memory address where the calculated oscillation probability for events of the specific requested type will be stored
   */
  const FLOAT_T* ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal=DUMMYVAL);

  /**
   * @brief Set the Energy array in the OscProbCalcerBase::OscProbCalcerBase() object
   *
   * The unbinned approach uses true event-by-event Energy and CosineZ values. These need to be set after they are read in from MC so can't be defined a-priori. Thus, this
   * function calls OscProbCalcerBase::SetEnergyArray()
   * 
   * @param Array The energy array to set in the OscProbCalcerBase::OscProbCalcerBase() object 
   */
  void SetEnergyArray(std::vector<FLOAT_T> Array);

  /**
   * @brief Set the CosineZ array in the OscProbCalcerBase::OscProbCalcerBase() object
   *
   * The unbinned approach uses true event-by-event CosineZ and CosineZ values. These need to be set after they are read in from MC so can't be defined a-priori. Thus, this
   * function calls OscProbCalcerBase::SetCosineZArray()
   * 
   * @param Array The CosineZ array to set in the OscProbCalcerBase::OscProbCalcerBase() object 
   */
  void SetCosineZArray(std::vector<FLOAT_T> Array);
  
  // ========================================================================================================================================================================
  // Public virtual functions which need calculater specific implementations

 protected:

  // ========================================================================================================================================================================
  // Protected functions which are calculation implementation agnostic  

  // ========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic

  // ========================================================================================================================================================================
  // Basic protected variables required for oscillation probability calculation

 private:

  // ========================================================================================================================================================================
  // Basic private variables required for oscillation probability calculation

};

#endif
