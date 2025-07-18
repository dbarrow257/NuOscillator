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

  /**
   * @brief Default constructor
   *
   * Default constructor
   *
   * @param Config_ YAML node used to set runtime constants
   */
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
  const FLOAT_T* ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal=DUMMYVAL) override;

  /**
   * @brief Return a vector of bin edges used for oscillation probability plotting
   *
   * Return the binning used to for the coarse probability binning
   *
   * @param ReturnEnergy Flag used to identify whether to return Energy or CosineZ binning
   *
   * @return Vector of bin edges which are used for plotting purposes
   *
   */
  std::vector<FLOAT_T> ReturnBinEdgesForPlotting(bool ReturnEnergy) override;
  
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

  /**
   * @brief Initialise Oscillator instance: read energies and cosineZs, set them in Calcer
   */
  void Initialise();

  /**
   * @brief Calculated coarse bin probabilities, defined as the average over the fine bin probabilities inside a particular coarse bin
   */
  void PostCalculateProbabilities() override;

  /**
   * @brief Setup the oscillator
   */
  void SetupOscillatorImplementation() override;
  
  /**
   * @brief Return the index of the bin in which a given value would be.
   *
   *
   * @param Val Value whose bin index is desired.
   * @param BinEdges Binning used
   *
   * @return Index of the bin in which the input value is located.
   *
   */
  int FindBinIndexFromEdges(FLOAT_T Val, std::vector<FLOAT_T> BinEdges);

  /**
   * @brief Vector holding averaged Probabilities [length = nBins]
   */
  std::vector<FLOAT_T> AveragedOscillationProbabilities;         

  /**
   * @brief Vector holding the pointers to the fine oscillation probabilities to average across [length = nBins][length = nFineBinsInCoarseBin]
   */                                                                                             
  std::vector< std::vector<const FLOAT_T*> > OscillationProbabilitiesToAverage; 

  /**
   * @brief Number of Coarse Energy Bins
   */
  size_t nCoarseEnergyBins;

  /**
   * @brief Number of Coarse CosineZ Bins
   */
  size_t nCoarseCosineZBins;

  /**
   * @brief Number of Fine Energy Bins
   */
  size_t nFineEnergyBins;

  /**
   * @brief Number of Fine CosineZ Bins
   */
  size_t nFineCosineZBins;

  /**
   * @brief Number of Total Coarse Bins
   */
  size_t TotalCoarseBins;

  /**
   * @brief Oscillation channels in NuOscillator::OscillationChannel format
   */
  std::vector<NuOscillator::OscillationChannel> OscillationChannels;

  /**
   * @brief Number of oscillation channels
   */
  size_t nOscillationChannels;

  /**
   * @brief Vector of indices for the neutrino types 
   */
  std::vector<int> NeutrinoTypes;

  /**
   * @brief Number of neutrino types 
   */
  size_t nNeutrinoTypes;
  
  // ========================================================================================================================================================================
  // Basic private variables required for oscillation probability calculation

  /**
   * @brief The FileName which the binning is read from
   */
  std::string FileName;

  /**
   * @brief The name of the histogram which contains the Energy axis coarse binning
   */
  std::string CoarseEnergyAxisHistName;

  /**
   * @brief The name of the histogram which contains the Energy axis fine binning
   */
  std::string FineEnergyAxisHistName;

  /**
   * @brief The name of the histogram which contains the CosineZ axis coarse binning
   */
  std::string CoarseCosineZAxisHistName;

  /**
   * @brief The name of the histogram which contains the CosineZ axis fine binning
   */
  std::string FineCosineZAxisHistName;

  /**
   * @brief A vector of Energy axis coarse bin edges
   */
  std::vector<FLOAT_T> CoarseEnergyAxisBinEdges;

  /**
   * @brief A vector of Energy axis fine bin edges
   */
  std::vector<FLOAT_T> FineEnergyAxisBinEdges;

  /**
   * @brief A vector of CosineZ axis coarse bin edges
   */
  std::vector<FLOAT_T> CoarseCosineZAxisBinEdges;

  /**
   * @brief A vector of CosineZ axis fine bin edges
   */
  std::vector<FLOAT_T> FineCosineZAxisBinEdges;

  /**
   * @brief A vector of Energy axis fine bin centers
   */
  std::vector<FLOAT_T> FineEnergyAxisBinCenters;
  
  /**
   * @brief A vector of CosineZ axis fine bin centers
   */
  std::vector<FLOAT_T> FineCosineZAxisBinCenters;

};

#endif
