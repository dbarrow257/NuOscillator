#ifndef __OSCILLATOR_BASE_H__A
#define __OSCILLATOR_BASE_H__

#include "OscProbCalcerBase.h"

/**
 * @file OscillatorBase.h
 *
 * @brief Oscillation calculation (binned, unbinned, etc.) implementation agnostic base class.
 *
 * The base class which controls the calculation of neutrino oscillation probabilities through various calculation techniques. There are
 * oscillation calculation techniques which can be implemented (e.g. binned, unbinned) which are expected to be derived classes of this
 * base function
 *
 * 
 */
class OscillatorBase {
 public:

  // ========================================================================================================================================================================
  // Public functions which are calculation implementation agnostic

  /**
   * @brief Perform a sanity check which ensures that #fOscProbCalcers have been set correctly.
   *
   * Ensure that: the #fOscProbCalcers has been initialised with atleast one OscProbCalcerBase::OscProbCalcerBase() instance, each OscProbCalcerBase::OscProbCalcerBase()
   * instance passes the sanity check, and that each OscProbCalcerBase::OscProbCalcerBase() instance has the same CosineZ treatement as the OscillatorBase()
   * instance expects.
   */
  void SanityCheck();

  /**
   * @brief Calculate the oscillation probabilities for the providied oscillation probabilities.
   *
   * For each instance of OscProbCalcerBase::OscProbCalcerBase in #fOscProbCalcers, calculate the oscillation probabilites for the provided oscillation parameters.
   *
   * @param OscParams Vector of oscillation parameters to calculate probabities at.
   */
  void CalculateProbabilities(std::vector<FLOAT_T> OscParams);

  /**
   * @brief Return number of expected oscillation parameters for a particular OscProbCalcerBase::OscProbCalcerBase() instance in #fOscProbCalcers.
   *
   * Each instance of OscProbCalcerBase::OscProbCalcerBase expects a particular number of oscillation parameters. This function returns this value for an 
   * OscProbCalcerBase::OscProbCalcerBase instance in #fOscProbCalcers.
   *
   * @param CalcerIndex Index of #fOscProbCalcers instance in which to query
   * @return The number of expected oscillation parameters from the CalcerIndex-th index in #fOscProbCalcers
   */
  int ReturnNOscParams(int CalcerIndex=0);

  /**
   * @brief Print the oscillation probabilities for a particular OscProbCalcerBase::OscProbCalcerBase() instance in #fOscProbCalcers.
   *
   * Print the calculated oscillation probability values for a particular OscProbCalcerBase::OscProbCalcerBase() instance in #fOscProbCalcers.
   * This is typically performed after a call to CalculateProbabilities()
   */
  void PrintWeights(int CalcerIndex=0);

  /**
   * @brief Setup each instance of OscProbCalcerBase::OscProbCalcerBase() instance in #fOscProbCalcers
   *
   * Setup function is public so that OscProbCalcerBase::OscProbCalcerBase() instances, when using the unbinned calculation implementation, can be setup. The unbinned
   * implementation assumes that the Energy and CosineZ arrays will be set after the OscillatorBase() object has been initialised.
   */
  void Setup();

  // ========================================================================================================================================================================
  // Public virtual functions which need calculater specific implementations

  /**
   * @brief Return a pointer to the oscillation probability for the requested event attributes
   *
   * Determine the memory address address where the calculated oscillation probability for events of the specific requested type will be stored. This will be different
   * depending on the calculation implementation. For the binned approach, the particular bin in which the requested energy and cosine falls will need to be determined.
   * For unbinned approach, the exact neutrino energy and cosine are expected to be used. Alternative approachs which smear the oscillation probability will be able to
   * use this function to apply the smear.
   *
   * @param InitNuFlav Initial neutrino flavour of the neutrino
   * @param FinalNuFlav Final neutrino flavour of the neutrino
   * @param EnergyVal True energy of the neutrino
   * @param CosineZVal True direction of the neutrino in CosineZ
   *
   * @return Pointer to the memory address where the calculated oscillation probability for events of the specific requested type will be stored
   */
  virtual const FLOAT_T* ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal=DUMMYVAL) = 0;

 protected:
  // ========================================================================================================================================================================
  // Protected functions which are calculation implementation agnostic

  /**
   * @brief Default constructor
   */
  OscillatorBase();

  /**
   * @brief Initialise an OscProbCalcerBase::OscProbCalcerBase() instance for each entry in #fOscProbCalcerImplementationToCreate
   *
   * #fOscProbCalcerImplementationToCreate is expected to be initialised within the costructor of the calculation specific code. For each entry in this vector, create an 
   * instance of OscProbCalcerBase::OscProbCalcerBase() and store it in #fOscProbCalcers. This function first parses the #fOscProbCalcerImplementationToCreate vector to
   * ensure that it is correctly filled by the calculation specific code. It then calls InitialiseOscProbCalcer(), for each entry in #fOscProbCalcerImplementationToCreate
   * and it is that function which actually initialises a specific OscProbCalcerBase::OscProbCalcerBase() object and returns it.
   */
  void InitialiseOscProbCalcers();

  /**
   * @brief Set the energy array which will be used by the OscProbCalcerBase::OscProbCalcerBase() instance stored in a particular index in #fOscProbCalcers
   *
   * Each instance of OscProbCalcerBase::OscProbCalcerBase() needs an energy array which will be used to calculate the oscillation probabilities. This function sets that
   * array
   *
   * @param Array The energy array which will be passed to the OscProbCalcerBase::OscProbCalcerBase() instance
   * @param CalcerIndex The index iin #fOscProbCalcers which will be handed the energy array
   */
  void SetEnergyArrayInCalcer(std::vector<FLOAT_T> Array, int CalcerIndex=0);

  /**
   * @brief Set the energy array which will be used by the OscProbCalcerBase::OscProbCalcerBase() instance stored in a particular index in #fOscProbCalcers
   *
   * Each instance of OscProbCalcerBase::OscProbCalcerBase() needs an energy array which will be used to calculate the oscillation probabilities. This function sets that
   * array
   *
   * @param Array The energy array which will be passed to the OscProbCalcerBase::OscProbCalcerBase() instance
   * @param CalcerIndex The index iin #fOscProbCalcers which will be handed the energy array
   */
  void SetCosineZArrayInCalcer(std::vector<FLOAT_T> Array, int CalcerIndex=0);

  /**
   *
   */
  const FLOAT_T* ReturnPointerToWeightinCalcer(int CalcerIndex, int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal=DUMMYVAL);
  
  // ========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic
  
  // ========================================================================================================================================================================
  // Basic protected variables required for oscillation probability calculation

  std::vector<std::string> fOscProbCalcerImplementationToCreate;

  bool fCosineZIgnored;

  // This is a vector object to accomodate any implementations which require multiple calculators to perform the reweight
  // For instance, this could be used to deal with the MaCh3 Event-by-Event approach by having a OscProbCalcerBase object for each oscillation channel
  int fNCalcers;

  /**
   * @brief A vector which contains all instances of OscProbCalcerBase()
   */
  std::vector<OscProbCalcerBase*> fOscProbCalcers;

  int fVerbose;
  enum Verbosity{NONE,INFO};

 private:
  OscProbCalcerBase* InitialiseOscProbCalcer(std::string OscProbCalcerImplementationToCreate);

  // ========================================================================================================================================================================
  // Basic private variables required for oscillation probability calculation
  bool fOscProbCalcerSet;
};

#endif
