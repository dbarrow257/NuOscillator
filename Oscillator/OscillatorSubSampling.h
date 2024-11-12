#ifndef __OSCILLATOR_SUBSAMPLING_BASE_H__
#define __OSCILLATOR_SUBSAMPLING_BASE_H__

#include "OscillatorBase.h"

/**
 * @file OscillatorSubSampling.h
 *
 * @class OscillatorSubSampling
 *
 * @brief SubSampling Oscillation calculation implementation class.
 *
 * Implementation of OscillatorBase::OscillatorBase() object which uses standard binning for the energy and cosineZ dimension, read from TFile/TH1D such that the binning in
 * each dimension is independent. 
 */
class OscillatorSubSampling : public OscillatorBase {
 public:

  /**
   * @brief Default constructor
   *
   * Default constructor
   *
   * @param ConfigName_ YAML config file used to set runtime constants
   */
  OscillatorSubSampling(std::string ConfigName_);

  OscillatorSubSampling(YAML::Node Config_);

  /**
   * @brief Destructor
   */
  virtual ~OscillatorSubSampling();

  // ========================================================================================================================================================================
  // Public functions which are calculation implementation agnostic

  /**
   * @brief Return a pointer to the oscillation probability for the requested event attributes.
   * 
   * Determine the memory address address where the calculated oscillation probability for events of the specific requested type will be stored. This will be different
   * depending on the calculation implementation. For the binned approach, the particular bin in which the requested energy and cosine falls will need to be determined.
   *
   * @param InitNuFlav Initial neutrino flavour of the neutrino
   * @param FinalNuFlav	Final neutrino flavour of the neutrino
   * @param EnergyVal True energy of the neutrino
   * @param CosineZVal True direction of the neutrino in CosineZ
   *
   * @return Pointer to the memory address where the calculated oscillation probability for events of the specific requested type will be stored
   */
  const FLOAT_T* ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal=DUMMYVAL);
  
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

  void Initialise();

  void PostCalculateProbabilities();

  void SetupOscillatorImplementation();
  
  int FindBinIndexFromEdges(FLOAT_T Val, std::vector<FLOAT_T> BinEdges);

  std::vector<FLOAT_T> AveragedOscillationProbabilities; //Vector holding averaged Probabilities [length = nBins]                                                                                                     
  std::vector< std::vector<const FLOAT_T*> > OscillationProbabilitiesToAverage; //Vector holding the pointers to the fine oscillation probabilities to average across [length = nBins][length = nFineBinsInCoarseBin]

  size_t nCoarseEnergyBins;
  size_t nCoarseCosineZBins;

  size_t nFineEnergyBins;
  size_t nFineCosineZBins;

  size_t TotalCoarseBins;

  std::vector<NuOscillator::OscillationChannel> OscillationChannels;
  size_t nOscillationChannels;

  std::vector<int> NeutrinoTypes;
  size_t nNeutrinoTypes;
  
  // ========================================================================================================================================================================
  // Basic private variables required for oscillation probability calculation

  /**
   * @brief The FileName which the binning is read from
   */
  std::string FileName;

  /**
   * @brief The name of the histogram which contains the Energy axis binning
   */
  std::string CoarseEnergyAxisHistName;
  std::string FineEnergyAxisHistName;
  /**
   * @brief The name of the histogram which contains the CosineZ axis binning
   */
  std::string CoarseCosineZAxisHistName;
  std::string FineCosineZAxisHistName;

  /**
   * @brief A vector of Energy axis bin edges
   */
  std::vector<FLOAT_T> CoarseEnergyAxisBinEdges;
  std::vector<FLOAT_T> FineEnergyAxisBinEdges;
  /**
   * @brief A vector of CosineZ axis bin edges
   */
  std::vector<FLOAT_T> CoarseCosineZAxisBinEdges;
  std::vector<FLOAT_T> FineCosineZAxisBinEdges;

  /**
   * @brief A vector of Energy axis bin centers
   */
  std::vector<FLOAT_T> FineEnergyAxisBinCenters;
  /**
   * @brief A vector of CosineZ axis bin centers
   */
  std::vector<FLOAT_T> FineCosineZAxisBinCenters;

};

#endif
