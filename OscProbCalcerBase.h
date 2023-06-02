#ifndef __OSCILLATOR_BASE_H__
#define __OSCILLATOR_BASE_H__

#define DUMMYVAL -999

#ifdef UsingDoubles
using FLOAT_T = double;
#else
using FLOAT_T = float;
#endif

#include <vector>
#include <string>

class OscProbCalcerBase {
 public:
  OscProbCalcerBase();

  // ========================================================================================================================================================================
  // Public functions which are calculation implementation agnostic
  
  // Define the energies and cosines which will be used when calculating the oscillation probabilities
  void SetEnergyArray(std::vector<FLOAT_T> EnergyArray);
  void SetCosineZArray(std::vector<FLOAT_T> CosineZArray);

  // Return pointer to the weight array for a specific energy (and cosine)
  const FLOAT_T* ReturnPointerToWeight(int InitNuFlav, int FinalNuFlav, FLOAT_T Energy, FLOAT_T CosineZ=DUMMYVAL);

  // General function used to call the oscillation probability calculation
  void Reweight(std::vector<FLOAT_T> OscParams);

  // General function used to setup all variables/functions
  void Setup();

  // Does my instance of OscProbCalcerBase pass all the sanity check;
  void SanityCheck();

  // Print values of fWeightArray
  void PrintWeights();

  // Getters
  int ReturnExpectedNOscParams() {return fNOscParams;}
  int ReturnNEnergyPoints() {return fNEnergyPoints;}
  int ReturnNCosineZPoints() {return fNCosineZPoints;}

  // ========================================================================================================================================================================
  // Public virtual functions which need calculater specific implementations

 protected:

  // ========================================================================================================================================================================
  // Protected functions which are calculation implementation agnostic  

  // Some implementations don't care about the CosineZ values, so this is a method to set fCosineZArraySet to true in those cases
  void IgnoreCosineZBinning(bool Ignore);

  // Check whether the oscillation parameters have changed since their previous value
  bool AreOscParamsChanged(std::vector<FLOAT_T> OscParamsToCheck);

  // Save the oscillation parameters which have been requested
  void SetCurrOscParams(std::vector<FLOAT_T> OscParamsToCheck);

  // Reset the saved oscillation parameters
  void ResetCurrOscParams();

  // Return the index in fCosineZArray for a particular value of CosineZ. If it's not found, throws an error
  int ReturnCosineZIndexFromValue(FLOAT_T CosineZVal);

  // Return the index in fCosineZArray for a particular value of CosineZ. If it's not found, throws an error
  int ReturnEnergyIndexFromValue(FLOAT_T EnergyVal);

  // Calculate the index in the NuType/NuFlav mapping for a given value 
  int ReturnInitialIndexFromFlavour(int InitFlav);
  int ReturnFinalIndexFromFlavour(int FinalFlav);
  int ReturnNuTypeFromFlavour(int NuFlav);

  // Initialise the mapping arrays to a particular size with dummy values
  void InitialiseNeutrinoTypesArray(int Size);
  void InitialiseInitialFlavoursArray(int Size);
  void InitialiseFinalFlavoursArray(int Size);

  // Check that the NuType/NuFlav mapping is set correctly based on the inputs from the particular implementation
  void CheckNuFlavourMapping();

  // Initialise the array in which the oscillation probabilities will be stored.
  void IntialiseWeightArray();

  // Ensure that the oscillation probabilities are within [0.,1.] range
  void SanitiseProbabilities();

  // ========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic

  // Calculate some oscillation probabilities for a particular oscillation parameter set
  virtual void CalculateProbabilities(std::vector<FLOAT_T> OscParams) = 0;

  // Setup any implementation specific variables/functions
  virtual void SetupPropagator() = 0;

  // Return the index in the weight array for a specific combination of neutrino oscillation channel, energy and cosine zenith
  virtual int ReturnWeightArrayIndex(int NuTypeIndex, int InitNuIndex, int FinalNuIndex, int EnergyIndex, int CosineZIndex=-1) = 0;

  // Define the size of fWeightArray. This is implementation specific because some propagators calculate all
  // oscillation channels together, whilst others calculate only a single oscillation channel. 
  virtual long DefineWeightArraySize() = 0;

  // ========================================================================================================================================================================
  // Basic variables required for oscillation probability calculation

  // Enums to define the mappings below. Each implementation is expected to define a mapping of which initial and neutrino flavours are considered, along with whether
  // neutrinos and antineutrinos are considered

  //DB fXXX these parameters
  int nNeutrinoTypes;
  std::vector<int> NeutrinoTypes;
  int nInitialFlavours;
  std::vector<int> InitialFlavours;
  int nFinalFlavours;
  std::vector<int> FinalFlavours;

  // Store energy and cosine points which will be used when calculating oscillation probabilities
  int fNEnergyPoints;
  int fNCosineZPoints;
  std::vector<FLOAT_T> fEnergyArray;
  std::vector<FLOAT_T> fCosineZArray;

  // Place to store the oscillation probabilities
  int fNWeights;
  std::vector<FLOAT_T> fWeightArray;

  // Store the oscillation parameter set used to calculate the previous and current probabilities
  // These are used to check whether the oscillation parameters have been updated from the previous CalculateProbabilities(OscParams) call.
  // If they haven't don't update the already calculated probabilities
  int fNOscParams;
  std::vector<FLOAT_T> fOscParamsCurr;

  // Set some verbosity for console output
  int fVerbose;
  enum Verbosity{NONE,INFO};

  // Define the implementation name - Could be used for recasting later (similar to TObject::InheritsFrom() etc.)
  std::string fImplementationName;

 private:
  // ========================================================================================================================================================================
  // Sanity check variables
  bool fEnergyArraySet;
  bool fCosineZArraySet;
  bool fPropagatorSet;
  bool fWeightArrayInit;
  bool fNuMappingSet;
};

#endif
