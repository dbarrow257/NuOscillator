#ifndef __OSCILLATOR_BASE_H__
#define __OSCILLATOR_BASE_H__

#include "OscProbCalcer/OscProbCalcerBase.h"
#include "Constants/OscillatorConstants.h"

#include "yaml-cpp/yaml.h"

/**
 * @file OscillatorBase.h
 *
 * @class OscillatorBase
 *
 * @brief Oscillation calculation (binned, unbinned, etc.) implementation agnostic base class.
 *
 * The base class which controls the calculation of neutrino oscillation probabilities through various calculation techniques. There are
 * oscillation calculation techniques which can be implemented (e.g. binned, unbinned) which are expected to be derived classes of this
 * base function.
 */
class OscillatorBase {
 public:
   /**
    * @brief Destructor
    */
   virtual ~OscillatorBase();

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
  void CalculateProbabilities(const std::vector<FLOAT_T>& OscParams);

  void CalculateProbabilities();

  void DefineParameter(std::string ParName_, FLOAT_T* ParValue_) {
    if (!fOscProbCalcerSet) {
      std::cerr << "DefineParameter function called before OscProbCalcer set!" << std::endl;
      throw std::runtime_error("DefineParameter function called before OscProbCalcer set");
    }
    
    fOscProbCalcer->DefineParameter(ParName_,ParValue_);
  }

  /**
   * @brief Return number of expected oscillation parameters for a particular OscProbCalcerBase::OscProbCalcerBase() instance in #fOscProbCalcers.
   *
   * Each instance of OscProbCalcerBase::OscProbCalcerBase expects a particular number of oscillation parameters. This function returns this value for an 
   * OscProbCalcerBase::OscProbCalcerBase instance in #fOscProbCalcer.
   *
   * @return The number of expected oscillation parameters from #fOscProbCalcer
   */
  int ReturnNOscParams();

  /**
   * @brief Print the oscillation probabilities for a particular OscProbCalcerBase::OscProbCalcerBase() instance in #fOscProbCalcer.
   *
   * Print the calculated oscillation probability values for a particular OscProbCalcerBase::OscProbCalcerBase() instance in #fOscProbCalcers.
   * This is typically performed after a call to CalculateProbabilities()
   */
  void PrintWeights();

  /**
   * @brief Setup each instance of OscProbCalcerBase::OscProbCalcerBase() instance in #fOscProbCalcer
   *
   * Setup function is public so that OscProbCalcerBase::OscProbCalcerBase() instances, when using the unbinned calculation implementation, can be setup. The unbinned
   * implementation assumes that the Energy and CosineZ arrays will be set after the OscillatorBase() object has been initialised.
   */
  void Setup();

  /**
   * @brief Return a string which encapsulates the calculation type along with the calculation engine
   *
   * @return Returns a string which looks like CalculationType_CalculationEngine (e.g. Binned_CUDAProb3)
   */
  std::string ReturnImplementationName();

  /**
   * @brief Return the number of Energy points which are being evaluated in the OscProbCalcerBase::OscProbCalcerBase() object
   *
   * @return Returns an integer describing the number of Energy points which are being evaluated in the OscProbCalcerBase::OscProbCalcerBase() object
   */
  int ReturnNEnergyPoints();

  /**
   * @brief Return the number of CosineZ points which are being evaluated in the OscProbCalcerBase::OscProbCalcerBase() object
   *
   * @return Returns an integer describing the number of CosineZ points which are being evaluated in the OscProbCalcerBase::OscProbCalcerBase() object
   */
  int ReturnNCosineZPoints();

  /**
   * @brief Check whether a particular OscProbCalcerBase::OscProbCalcerBase() instance has a particular oscillation channel
   *
   * @param GeneratedFlavour The oscillation channel generated neutrino flavour to check for
   * @param DetectedFlavour The oscillation channel detected neutrino flavour to check for
   *
   * @return Boolean flag which describes whether the oscillation channel was found in the instance of OscProbCalcerBase::OscProbCalcerBase()
   */
  bool HasOscProbCalcerGotOscillationChannel(int GeneratedFlavour, int DetectedFlavour);

  /**
   * @brief Set the energy array which will be used by the OscProbCalcerBase::OscProbCalcerBase() instance stored in a particular index in #fOscProbCalcers
   *
   * Each instance of OscProbCalcerBase::OscProbCalcerBase() needs an energy array which will be used to calculate the oscillation probabilities. This function sets that
   * array
   *
   * @param Array The energy array which will be passed to the OscProbCalcerBase::OscProbCalcerBase() instance
   */
  void SetEnergyArrayInCalcer(std::vector<FLOAT_T> Array);

  /**
   * @brief Set the energy array which will be used by the OscProbCalcerBase::OscProbCalcerBase() instance stored in #fOscProbCalcer
   *
   * Each instance of OscProbCalcerBase::OscProbCalcerBase() needs an energy array which will be used to calculate the oscillation probabilities. This function sets that
   * array
   *
   * @param Array The energy array which will be passed to the OscProbCalcerBase::OscProbCalcerBase() instance
   */
  void SetCosineZArrayInCalcer(std::vector<FLOAT_T> Array);

  /**
   * @brief Return flag which describes whether the OscProbCalcerBase::OscProbCalcerBase() has had it's Energy and CosineZ evaluation points set in the constructor of the
   * derived OscillatorBase::OscillatorBase() object
   *
   * @return Boolean flag
   */
  bool EvalPointsSetInConstructor() {return fEvalPointsSetInConstructor;}

  /**
   * @brief Return flag which describes whether CosineZ binning is considered in the object
   *
   * @return Boolean flag
   */
  bool CosineZIgnored() {return fCosineZIgnored;}

  /**
   * @brief Return oscillation probability for a given initial and final flavour, for a given energy and cosineZ
   *
   * @param InitNuFlav Initial neutrino flavour of the neutrino
   * @param FinalNuFlav Final neutrino flavour of the neutrino
   * @param EnergyVal True energy of the neutrino
   * @param CosineZVal True direction of the neutrino in CosineZ
   *
   * @return Value of the oscillation probability for events of the specific requested type will be stored 
   */  
  FLOAT_T ReturnOscillationProbability(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal=DUMMYVAL) {
    const FLOAT_T* Pointer = ReturnWeightPointer(InitNuFlav, FinalNuFlav, EnergyVal, CosineZVal);
    return *Pointer;
  }
  
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

  /**
   * @brief Return a vector of bin edges which can be used to plot the oscillation probability
   *
   * If we want to plot the oscillation probabilties, we need to define a binning which can be used to create the histogram. This is implementation dependent
   * because the binned implementation should return the binning used to define the probability binning, but event-by-event needs to return a binning such
   * that a single evaluation point falls in a single bin
   *
   * @param ReturnEnergy Flag to return the binning for Energy or CosineZ
   *
   * @return Vector of bin edges which should be used to create a histogram for plotting purposes
   */
  virtual std::vector<FLOAT_T> ReturnBinEdgesForPlotting(bool ReturnEnergy) = 0;

 protected:
  /**
   * @brief Default constructor
   *
   * @details It is protected to prevent initialisation of base class
   *
   * @param ConfigName_ YAML config file used to set runtime constants
   */
  OscillatorBase(std::string ConfigName_);
  
  /**
   * @brief Default constructor
   *
   * @details It is protected to prevent initialisation of base class
   *
   * @param ConfigName_ YAML config file used to set runtime constants
   */
  OscillatorBase(YAML::Node Config_);

  // ========================================================================================================================================================================
  // Protected functions which are calculation implementation agnostic

  /**
   * @brief Return a pointer to the oscillation probability memory address in a particular index of #fOscProbCalcers for a particular event
   *
   * Return the memory address for the oscillation probability which is calculated for neutrinos of a particular initial and final flavour, energy and cosine value. This is
   * a memory address for a particular OscProbCalcerBase::OscProbCalcerBase() instance stored in a particular index in #fOscProbCalcers.
   *
   * @param CalcerIndex The index in #fOscProbCalcers.
   * @param InitNuFlav Initial neutrino flavour of event
   * @param FinalNuFlav Final neutrino flavour of event
   * @param EnergyVal Neutrino energy of event
   * @param CosineZVal Netrino cosine zenith direction of event
   *
   * @return Memory address associated with given event attributes for CalcerIndex-th index in #fOscProbCalcers
   */
  const FLOAT_T* ReturnPointerToWeightinCalcer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal=DUMMYVAL);
  
  // ========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic

  /**
   * @brief Do follow-up calculations with the oscillation probabilities. E.g.: averages in the SubSampling CalculationType.
   */
  virtual void PostCalculateProbabilities() {}

  /**
   * @brief Setup specific oscillator implementation
   */
  virtual void SetupOscillatorImplementation() {}

  // ========================================================================================================================================================================
  // Basic protected variables required for oscillation probability calculation

  /**
   * @brief A boolean used for declaring whether the OscillatorBase() object expects to care about the CosineZ dimension
   */
  bool fCosineZIgnored;

  /**
   * @brief Flad which describes whether the Energy and CosineZ evaluation points where set in the constructor of the OscillatorBase() derived object
   */
  bool fEvalPointsSetInConstructor;

  /**
   * @brief The instance of OscProbCalcerBase()
   */
  OscProbCalcerBase* fOscProbCalcer;

  /**
   * @brief A string describing the calculation implementation, e.g. Binned
   */
  std::string fCalculationTypeName;

  /**
   * @brief The verbosity level of console output
   */
  int fVerbose;

  /**
   * @brief YAML Config manager
   */
  YAML::Node Config;

 private:

  /**
   * @brief Return an OscProbCalcerBase::OscProbCalcerBase() object from the requested inputs
   *
   * Create and initialise #fOscProbCalcer to be an instance of OscProbCalcerBase::OscProbCalcerBase() associated with a particular implementation denoted by @param OscProbCalcerImplementationToCreate
   * and config path @param OscProbCalcerConfigname, recast it to a base object OscProbCalcerBase::OscProbCalcerBase() and returns it.
   */
  void InitialiseOscProbCalcer();

  /**
   * @brief Initialise OscillatorBase object
   *
   * Define everything which is needed for the OscillatorBase object to function
   *
   * @param Config_ YAML Node to be initialised from
   */
  void Initialise(YAML::Node Config_);

  // ========================================================================================================================================================================
  // Basic private variables required for oscillation probability calculation

  /**
   * @brief A boolean which declares whether #fOscProbCalcer has been initialised
   */
  bool fOscProbCalcerSet;

};

#endif
