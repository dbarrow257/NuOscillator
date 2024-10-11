#ifndef __OSCPROBCALCER_BASE_H__
#define __OSCPROBCALCER_BASE_H__

#include "Constants/OscillatorConstants.h"

#include <vector>
#include <string>

#include "yaml-cpp/yaml.h"

/**
 * @file OscProbCalcerBase.h
 *
 * @class OscProbCalcerBase
 *
 * @brief Oscillation calculation engine (CUDAProb3, ProbGPU, Prob3++, etc.) implementation agnostic base class.
 */
class OscProbCalcerBase {
 public:
  // ========================================================================================================================================================================
  // Public functions which are calculation implementation agnostic

   /**
    * @brief Destructor
    */
   virtual ~OscProbCalcerBase();

  /**
   * @brief Define the Energy which will be used when calculating the oscillation probabilities
   * 
   * @param EnergyArray The Energy array which will be used by the calculation engine
   */
  void SetEnergyArray(std::vector<FLOAT_T> EnergyArray);

  /**
   * @brief Define the CosineZ which will be used when calculating the oscillation probabilities
   *
   * @param CosineZArray The CosineZ array which will be used by the calculation engine
   */
  void SetCosineZArray(std::vector<FLOAT_T> CosineZArray);

  /**
   * @brief Return pointer to the weight array for a specific Energy and CosineZ
   *
   * The specific calculation engine will calculate an oscillation probability for a particular initial and final neutrino flavour, neutrino Energy and CosineZ and stores
   * it in #fWeightArray. This function returns a pointer to the correct index in #fWeightArray for the given inputs
   *
   * @param InitNuFlav Initial neutrino flavour of the neutrino
   * @param FinalNuFlav Final neutrino flavour of the neutrino
   * @param Energy True energy of the neutrino
   * @param CosineZ True direction of the neutrino in CosineZ
   * 
   * @return Pointer to the memory address where the calculated oscillation probability for events of the specific requested type will be stored
   */
  const FLOAT_T* ReturnPointerToWeight(int InitNuFlav, int FinalNuFlav, FLOAT_T Energy, FLOAT_T CosineZ=DUMMYVAL);

  /**
   * @brief General function used to call the oscillation probability calculation
   *
   * This function performs both the implementation specific CalculateProbabilities() function, along with checking whether the oscillation parameters have been
   * updated since the last call. It also calls SanitiseProbabilities().
   *
   * @param OscParams The oscillation parameters to calculate the oscillation probability at
   */
  void Reweight(const std::vector<FLOAT_T>& OscParams);

  /**
   * @brief General function used to setup all variables used within the reweighting
   *
   * Ensures that the Energy and CosineZ arrays have been set correctly, along with initialising the saved oscillation parameters and weight array (via. ResetCurrOscParams()
   * and IntialiseWeightArray() ). Then sets the propagator specific implementation SetupPropagator()
   */
  void Setup();

  /**
   * @brief Ensures that the specific implementation has been correctly initialised
   *
   * Checks that the Energy and CosineZ arrays have been passed to the OscProbCalcerBase::OscProbCalcerBase() object, the neutrino flavour mapping, weight array and 
   * propagator has also been initialised
   */
  bool SanityCheck();

  /**
   * @brief Print the calculated oscillation probabilities
   */
  void PrintWeights();

  /**
   * @brief Print the oscillation channels which a particular instance has been configured with
   */
  void PrintKnownOscillationChannels();

  /**
   * @brief Check whether a particular oscillation channel has been configured in a particular instance of OscProbCalcerBase::OscProbCalcerBase()
   *
   * @param GeneratedFlavour Generated flavour of the questioned oscillation channel
   * @param DetectedFlavour Detected flavour of the questioned oscillation channel
   *
   * @return Boolean parameter which describes if the oscillation channel is configured
   */
  bool HasOscillationChannel(int GeneratedFlavour, int DetectedFlavour);

  /**
   * @brief Return the number of oscillation parameters the specific implementation expects
   * @return Return the number of oscillation parameters the specific implementation expects
   */  
  int ReturnNOscParams() {return fNOscParams;}

  /**
   * @brief Return the oscillation parameters which were used for the last calculation
   * @return Return the oscillation parameters which were used for the last calculation
   */
  std::vector<FLOAT_T> ReturnOscParamsCurr() {return fOscParamsCurr;}

  /**
   * @brief Return the number of Energy points which are being used by the specific instance of OscProbCalcerBase::OscProbCalcerBase()
   * @return Return the number of Energy points which are being used by the specific instance of OscProbCalcerBase::OscProbCalcerBase()
   */
  int ReturnNEnergyPoints() {return fNEnergyPoints;}

  /**
   * @brief Return the number of CosineZ points which are being used by the specific instance of OscProbCalcerBase::OscProbCalcerBase()
   * @return Return the number of CosineZ points which are being used by the specific instance of OscProbCalcerBase::OscProbCalcerBase()
   */
  int ReturnNCosineZPoints() {return fNCosineZPoints;}

  /**
   * @brief Return the vector of Energy points which are being used by the specific instance of OscProbCalcerBase::OscProbCalcerBase()
   * @return Return the vector of Energy points which are being used by the specific instance of OscProbCalcerBase::OscProbCalcerBase()
   */
  std::vector<FLOAT_T> ReturnEnergyArray() {return fEnergyArray;}

  /**
   * @brief Return the vector of CosineZ points which are being used by the specific instance of OscProbCalcerBase::OscProbCalcerBase()
   * @return Return the vector of CosineZ points which are being used by the specific instance of OscProbCalcerBase::OscProbCalcerBase()
   */
  std::vector<FLOAT_T> ReturnCosineZArray() {return fCosineZArray;}

  /**
   * @brief Return the specific implementation name
   * @return Return the specific implementation name
   */
  std::string ReturnImplementationName() {return fImplementationName;}

  /**
   * @brief Return a boolean which describes whether the specific implementation cares about CosineZ
   * @return Return a boolean which describes whether the specific implementation cares about CosineZ
   */
  bool ReturnCosineZIgnored() {return fCosineZIgnored;}

  /**
   * @brief Return the number of oscillation probabilities that are being calculated
   * @return Return the number of oscillation probabilities that are being calculated
   */
  int ReturnNWeights() {return fNWeights;}

  /**
   * @brief Return the vector of oscillation probabilites which have been calculated
   * @return Return the vector of oscillation probabilites which have been calculated
   */
  std::vector<FLOAT_T> ReturnWeightArray() {return fWeightArray;}

  /**
   * @brief Return vector of oscillation probabilities with associated neutrin type, oscillation channel, Energy and CosineZ
   *
   * This uses the OscillationProbability struct
   *
   * @return Vector of OscillationProbability structs
   */
  std::vector<NuOscillator::OscillationProbability> ReturnProbabilities();

  /**
   * @brief Return whether the instance of OscProbCalcerBase() has already been assigned the Energy binning
   * @return Boolean flag
   */
  bool ReturnHasSetEnergyArray() {return fEnergyArraySet;}

  /**
   * @brief Return whether the instance of OscProbCalcerBase() has already been assigned the CosineZ binning
   * @return Boolean flag
   */
  bool ReturnHasSetCosineZArray() {return fCosineZArraySet;}

  // ========================================================================================================================================================================
  // Public virtual functions which need calculater specific implementations

 protected:

  /**
    * @brief Default constructor
    *
    * @details It is protected to prevent initialisation of base class
    *
    * @param Config_ YAML::Node used to define the object
    */
  OscProbCalcerBase(YAML::Node InputConfig_);

   /**
    * @brief Wrapper Construct
    *
    * @details If only a file path of a config has been specified, create a YAML::Node and then continue with the default constructor
    *
    * @param ConfigName_ Name of YAML config used to set runtime variables
    */
  OscProbCalcerBase(std::string ConfigName_) : OscProbCalcerBase(YAML::LoadFile(ConfigName_)) {}
  
  // ========================================================================================================================================================================
  // Protected functions which are calculation implementation agnostic  

  /**
   * @brief Some implementations don't care about the CosineZ values, so this is a method to set #fCosineZArraySet to true in those cases
   *
   * If the specific implementation does not include atmospheric implementations, then #fCosineZArraySet should be set to true such that the Sanity check passes. Also
   * sets #fCosineZIgnored to be true if it should be ignored.
   *
   * @param Ignore A boolean which declares whether the CosineZ dimension should be ignored.
   */
  void IgnoreCosineZBinning(bool Ignore);

  /**
   * @brief Check whether the oscillation parameters have changed since their previous value
   *
   * We can save a fairly expensive calculation by checking whether the oscillation parameters have been updated since the last call to CalculateProbabilities() , for the
   * price of a couple of if statements. Whilst in normal MCMC running, this will almost always return true, for sigma variations of systematics, this will almost always
   * return false
   *
   * @param OscParamsToCheck The oscillation parameters which have been requested for the next calculation
   * @return Boolean whether the OscParamsToCheck match those saved from the last calculation
   */
  bool AreOscParamsChanged(const std::vector<FLOAT_T>& OscParamsToCheck);

  /**
   * @brief Save the oscillation parameters which have been requested
   * 
   * Save the oscillation parameters which have been used in the probability calculation to #fOscParamsCurr
   *
   * @param OscParamsToCheck Parameter set to save
   */
  void SetCurrOscParams(const std::vector<FLOAT_T>& OscParamsToCheck);

  /**
   * @brief (Re-)Initialise the saved oscillation parameters in #fOscParamsCurr
   */
  void ResetCurrOscParams();

  /**
   * @brief Return the index in #fCosineZArray for a particular value of CosineZ. If it's not found, throws an error
   *
   * Determine the index in #fCosineZArray for a particular value of CosineZ
   *
   * @param CosineZVal Value to search for
   * @return Index in #fCosineZArray
   */
  int ReturnCosineZIndexFromValue(FLOAT_T CosineZVal);

  /**
   * @brief Return the index in #fEnergyArray for a particular value of Energy. If it's not found, throws an error
   *
   * Determine the index in #fEnergyArray for a particular value of Energy
   *
   * @param EnergyVal Value to search for
   * @return Index in #fEnergyArray
   */
  int ReturnEnergyIndexFromValue(FLOAT_T EnergyVal);

  /**
   * @brief Loop through #fOscillationChannels and determine which index corresponds to the requested flavours
   *
   * @param InitNuFlav Generated neutrino flavour to search for
   * @param FinalNuFlav Detected neutrino flavour to search for
   *
   * @return The index which corresponds to the requested neutrino flavours
   */
  int ReturnOscChannelIndexFromFlavours(int InitNuFlav, int FinalNuFlav);

  /**
   * @brief Return the index in #fNeutrinoTypes for a particular neutrino flavour (neutrino or antineutrino)
   *
   * Neutrino flavour mapping is defined in #fNeutrinoTypes, so for a particular events flavour, the index in the mapping is returned to aide in #fWeightArray mapping
   *
   * @param NuFlav Initial neutrino type (neutrino or antineutrino) 
   * @return Index in #fNeutrinoTypes
   */
  int ReturnNuTypeFromFlavour(int NuFlav);

  /**
   * @brief Initialise the #fNeutrinoTypes mapping array to a particular size with dummy values
   *
   * @param Size Size of array to initialise
   */
  void InitialiseNeutrinoTypesArray(int Size);
  
  /**
   * @brief Read the [OscProbCalcerSetup][Implementation][OscChannelMapping] YAML node and configure #fOscillationChannels
   */
  void InitialiseOscillationChannelMapping();

  /**
   * @brief Check that the NuType/NuFlav mapping is set correctly based on the inputs from the particular implementation
   *
   * Ensures that the mapping variables (#fNeutrinoTypes, #fOscillationChannels) are the expected size and filled with reasonable values
   */
  void CheckNuFlavourMapping();

  /**
   * @brief Initialise the array in which the oscillation probabilities will be stored with dummy values
   */
  void IntialiseWeightArray();

  /**
   * @brief Ensure that the oscillation probabilities are within [0.,1.] range, if not throw error
   */
  void SanitiseProbabilities();

  // ========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic

  /**
   * @brief Calculate some oscillation probabilities for a particular oscillation parameter set
   *
   * Calculator engine specific implementation to calculate the oscillation probabilties for the given oscillation parameters. This function both calculates and stores
   * the oscillation probabilities in #fWeightArray.
   *
   * @param OscParams The parameter set to calculate oscillation probabilities at
   */
  virtual void CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) = 0;

  /**
   * @brief Setup any implementation specific variables/functions
   *
   * Calculator engine specific implementation to setup any variables needed to calculate the oscillation probabilities
   */
  virtual void SetupPropagator() = 0;

  /**
   * @brief Return the index in the weight array for a specific combination of neutrino oscillation channel, energy and cosine zenith
   *
   * Depending upon the specific implementation, there maybe different ways to calculate the oscillation probabilites (i.e. all oscillation channels at once or seperately).
   * This function allows a mapping between the neutrino flavour to calculate for (along with Energy and CosineZ) to the index in #fWeightArray.
   * 
   * @param NuTypeIndex The index in #fNeutrinoTypes (neutrino/antinuetrino) to return the pointer for 
   * @param OscChanIndex The index in #fOscillationChannels to return the pointer for 
   * @param EnergyIndex The index in #fEnergyArray to return the pointer for 
   * @param CosineZIndex The index in #fCosineZArray to return the pointer for 
   *
   * @return Index in #fWeightArray which corresponds to the given inputs
   */
  virtual int ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex=-1) = 0;

  /**
   * @brief Define the size of fWeightArray
   *
   * This is implementation specific because some propagators calculate all oscillation channels together, whilst others calculate only a single oscillation channel.
   *
   * @return Length that #fWeightArray should be initialised to
   */
  virtual long DefineWeightArraySize() = 0;

  // ========================================================================================================================================================================
  // Basic variables required for oscillation probability calculation

  // Enums to define the mappings below. Each implementation is expected to define a mapping of which initial and neutrino flavours are considered, along with whether
  // neutrinos and antineutrinos are considered

  /**
   * @brief The number of neutrino types (neutrino and antineutrino)
   */
  int fNNeutrinoTypes;

  /**
   * @brief The mapping of neutrino types to uniquely defined integers (i.e. index)
   */
  std::vector<int> fNeutrinoTypes;

  /**
   * @brief The number of oscillation channels which have been configured
   */
  int fNOscillationChannels;

  /**
   * @brief The mapping of oscillation channels which have been requested from the YAML config
   */
  std::vector<NuOscillator::OscillationChannel> fOscillationChannels;

  /**
   * @brief The number of Energy points which are being evaluated by the oscillation probability engine
   */
  int fNEnergyPoints;

  /**
   * @brief The number of CosineZ points which are being evaluated by the oscillation probability engine
   */
  int fNCosineZPoints;

  /**
   * @brief The vector of Energy values being evaluated by the oscillation probability engine
   */
  std::vector<FLOAT_T> fEnergyArray;

  /**
   * @brief The vector of CosineZ values being evaluated by the oscillation probability engine
   */
  std::vector<FLOAT_T> fCosineZArray;

  /**
   * @brief The number of oscillation probabilities being calculated
   */
  int fNWeights;
  
  /**
   * @brief Vector which stores the oscillation probabilities
   */
  std::vector<FLOAT_T> fWeightArray;

  /**
   * @brief Store the oscillation parameter set used to calculate the previous and current probabilities
   */
  int fNOscParams;

  /**
   * @brief Vector which stores the saved oscillation probabilities
   */
  std::vector<FLOAT_T> fOscParamsCurr;

  /**
   * @brief Define the verbosity of the console output
   */
  int fVerbose;

  /**
   * @brief Define the implementation name - Can be used for recasting later (similar to TObject::InheritsFrom() etc.)
   */
  std::string fImplementationName;

  /**
   * @brief Flag to define whether the CosineZ binning has been ignored in the specific implementation
   */
  bool fCosineZIgnored;

  /**
   * @brief YAML Config object used to get runtime specific variables
   */
  YAML::Node Config;

 private:
  // ========================================================================================================================================================================
  // Sanity check variables
  /**
   * @brief Boolean declaring whether the Energy array has been set
   */
  bool fEnergyArraySet;

  /**
   * @brief Boolean declaring whether the CosineZ array has been set
   */
  bool fCosineZArraySet;
  
  /**
   * @brief Boolean declaring whether the specific propagator SetupPropagator() function has been called
   */
  bool fPropagatorSet;

  /**
   * @brief Boolean declaring whether #fWeightArray has been initialised correctly
   */
  bool fWeightArrayInit;

  /**
   * @brief Boolean declaring whether #fNeutrinoTypes and #fOscillationChannels have been initialised correctly
   */
  bool fNuMappingSet;
};

#endif
