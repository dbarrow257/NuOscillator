#ifndef __OSCILLATOR_BASE_H__
#define __OSCILLATOR_BASE_H__

#define DUMMYVAL -999

#ifdef UsingDoubles
using FLOAT_T = double;
#else
using FLOAT_T = float;
#endif

#include <vector>

class OscProbCalcerBase {
 public:
  OscProbCalcerBase();

  //========================================================================================================================================================================
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

  //Getters
  int ReturnExpectedNOscParams() {return fNOscParams;}
  int ReturnNEnergyPoints() {return fNEnergyPoints;}
  int ReturnNCosineZPoints() {return fNCosineZPoints;}

  //========================================================================================================================================================================
  // Public virtual functions which need calculater specific implementations

 protected:

  //========================================================================================================================================================================
  // Protected functions which are calculation implementation agnostic  

  // Check whether the oscillation parameters have changed since their previous value
  bool AreOscParamsChanged(std::vector<FLOAT_T> OscParamsToCheck);

  // Save the oscillation parameters which have been requested
  void SetCurrOscParams(std::vector<FLOAT_T> OscParamsToCheck);

  // Reet the saved oscillation parameters
  void ResetCurrOscParams();

  //========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic

  // Calculate some oscillation probabilities for a particular oscillation parameter set
  virtual void CalculateProbabilities(std::vector<FLOAT_T> OscParams) = 0;

  // Setup any implementation specific variables/functions
  virtual void SetupPropagator() = 0;

  // Return pointer to the weight array for a specific energy
  virtual const FLOAT_T* ReturnPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T Energy, FLOAT_T CosineZ=DUMMYVAL) = 0;

  // Initialise the array in which the oscillation probabilities will be stored. This is implementation specific because some propagators calculate all
  // oscillation channels together, whilst others calculate only a single oscillation channel. 
  virtual void IntiailiseWeightArray() = 0;

  //========================================================================================================================================================================
  // Sanity check variables
  bool fEnergyArraySet;
  bool fCosineZArraySet;
  bool fPropagatorSet;
  bool fWeightArrayInit;

  //========================================================================================================================================================================
  // Basic variables required for oscillation probability calculation

  // Store energy and cosine points which will be used when calculating oscillation probabilities
  int fNEnergyPoints;
  int fNCosineZPoints;
  std::vector<FLOAT_T> fEnergyArray;
  std::vector<FLOAT_T> fCosineZArray;

  // Place to store the oscillation probabilities
  std::vector<FLOAT_T> fWeightArray;

  // Store the oscillation parameter set used to calculate the previous and current probabilities
  // These are used to check whether the oscillation parameters have been updated from the previous CalculateProbabilities(OscParams) call.
  // If they haven't don't update the already calculated probabilities
  int fNOscParams;
  std::vector<FLOAT_T> fOscParamsCurr;
};

#endif
