#ifndef __OSCILLATOR_BASE_H__
#define __OSCILLATOR_BASE_H__

#include "OscProbCalcerBase.h"

class OscillatorBase {
 public:
  OscillatorBase(std::vector<std::string> OscProbCalcerImplementationToCreate);

  // ========================================================================================================================================================================
  // Public functions which are calculation implementation agnostic
  
  //DB
  void SanityCheck();
  void PrintImplementationName(int CalcerIndex=0);
  void SetEnergyArrayInCalcer(std::vector<FLOAT_T> Array, int CalcerIndex=0);
  void SetCosineZArrayInCalcer(std::vector<FLOAT_T> Array, int CalcerIndex=0);
  void CalculateProbabilities(std::vector<FLOAT_T> OscParams);
  void Setup();
  int ReturnNOscParams(int CalcerIndex=0);
  void PrintWeights(int CalcerIndex=0);

  // ========================================================================================================================================================================
  // Public virtual functions which need calculater specific implementations

 protected:

  // ========================================================================================================================================================================
  // Protected functions which are calculation implementation agnostic  

  // ========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic

  // ========================================================================================================================================================================
  // Basic protected variables required for oscillation probability calculation

  // This is a vector object to accomodate any implementations which require multiple calculators to perform the reweight
  // For instance, this could be used to deal with the MaCh3 Event-by-Event approach by having a OscProbCalcerBase object for each oscillation channel
  int fNCalcers;
  std::vector<OscProbCalcerBase*> OPCalcers;

  int fVerbose;
  enum Verbosity{NONE,INFO};

 private:
  OscProbCalcerBase* InitialiseOscProbCalcer(std::string OscProbCalcerImplementationToCreate);

  // ========================================================================================================================================================================
  // Basic private variables required for oscillation probability calculation
  bool fOscProbCalcerSet;
};

#endif
