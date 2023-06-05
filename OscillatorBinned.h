#ifndef __OSCILLATOR_BINNED_BASE_H__
#define __OSCILLATOR_BINNED_BASE_H__

#include "OscillatorBase.h"

class OscillatorBinned : public OscillatorBase {
 public:
  OscillatorBinned(std::vector<std::string> OscProbCalcerImplementationToCreate, bool fCosineZIgnored_=false);

  // ========================================================================================================================================================================
  // Public functions which are calculation implementation agnostic

  const FLOAT_T* ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal=DUMMYVAL);
  
  // ========================================================================================================================================================================
  // Public virtual functions which need calculater specific implementations

 protected:

  // ========================================================================================================================================================================
  // Protected functions which are calculation implementation agnostic  

  std::vector<FLOAT_T> ReadBinEdgesFromFile(std::string FileName, std::string HistogramName);
  std::vector<FLOAT_T> ReturnBinCentersFromBinEdges(std::vector<FLOAT_T> BinEdges);

  // ========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic

  // ========================================================================================================================================================================
  // Basic protected variables required for oscillation probability calculation

 private:

  // ========================================================================================================================================================================
  // Basic private variables required for oscillation probability calculation

  std::string FileName;
  std::string EnergyAxisHistName;
  std::string CosineZAxisHistName;

  std::vector<FLOAT_T> EnergyAxisBinEdges;
  std::vector<FLOAT_T> CosineZAxisBinEdges;
  std::vector<FLOAT_T> EnergyAxisBinCenters;
  std::vector<FLOAT_T> CosineZAxisBinCenters;
};

#endif
