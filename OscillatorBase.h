#ifndef __OSCILLATOR_BASE_H__
#define __OSCILLATOR_BASE_H__

#include "OscProbCalcerBase.h"

class OscillatorBase {
 public:

  // ========================================================================================================================================================================
  // Public functions which are calculation implementation agnostic
  
  //DB Info in header
  void SanityCheck();
  void PrintImplementationName(int CalcerIndex=0);
  void CalculateProbabilities(std::vector<FLOAT_T> OscParams);
  int ReturnNOscParams(int CalcerIndex=0);
  void PrintWeights(int CalcerIndex=0);

  // Setup function has to be public so Unbinned approach can do it's thing
  void Setup();

  // ========================================================================================================================================================================
  // Public virtual functions which need calculater specific implementations

  virtual const FLOAT_T* ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal=DUMMYVAL) = 0;

 protected:
  OscillatorBase();

  // ========================================================================================================================================================================
  // Protected functions which are calculation implementation agnostic  

  void InitialiseOscProbCalcers();

  void SetEnergyArrayInCalcer(std::vector<FLOAT_T> Array, int CalcerIndex=0);
  void SetCosineZArrayInCalcer(std::vector<FLOAT_T> Array, int CalcerIndex=0);

  // ========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic

  // ========================================================================================================================================================================
  // Basic protected variables required for oscillation probability calculation

  const FLOAT_T* ReturnPointerToWeightinCalcer(int CalcerIndex, int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal=DUMMYVAL);
  std::vector<std::string> fOscProbCalcerImplementationToCreate;

  bool fCosineZIgnored;

  // This is a vector object to accomodate any implementations which require multiple calculators to perform the reweight
  // For instance, this could be used to deal with the MaCh3 Event-by-Event approach by having a OscProbCalcerBase object for each oscillation channel
  int fNCalcers;
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
