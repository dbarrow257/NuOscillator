#ifndef __OSCILLATOR_BINNED_BASE_H__
#define __OSCILLATOR_BINNED_BASE_H__

#include "OscillatorBase.h"

/**
 * @file OscillatorBinned.h
 *
 * @class OscillatorBinned
 *
 * @brief Binned Oscillation calculation implementation class.
 *
 * Implementation of OscillatorBase::OscillatorBase() object which uses standard binning for the energy and cosineZ dimension, read from TFile/TH1D such that the binning in
 * each dimension is independent. 
 */
class OscillatorBinned : public OscillatorBase {
 public:

  /**
   * @brief Default constructor
   *
   * Default constructor
   *
   * @param ConfigName_ YAML config file used to set runtime constants
   */
  OscillatorBinned(std::string ConfigName_);

  /**
   * @brief Destructor
   */
  virtual ~OscillatorBinned();

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

  /**
   * @brief Read bin edges from input file and histogram
   *
   * Read the histogram (TH1) from the file (TFile). Bin edges are read for each dimension independently. If correlated binning is required, a new class implementation
   * should be added.
   *
   * @param FileName Name of TFile to read
   * @param HistogramName Name of TH1 to read from TFile
   * @param IsCosineZAxis Denotes whether the axis currently being read is the Energy or CosineZ axis
   *
   * @return Vector of bin edges
   */
  std::vector<FLOAT_T> ReadBinEdgesFromFile(std::string FileName, std::string HistogramName, bool IsCosineZAxis);

  /**
   * @brief Return the bin center values from the bin edges
   *
   * Return a vector of bin center values from the input bin edges. This is done as the OscProbCalcerBase::OscProbCalcerBase() object needs specific value to calculate
   * the oscillation probability. Hence, when using binned oscillation probabilities, a choice of what value to provide to OscProbCalcerBase::SetEnergyArray() needs to be
   * made. The choice taken here is to use the center of the bin.
   *
   * @param BinEdges A vector of bin edges
   *
   * @return A vector of bin centers
   */
  std::vector<FLOAT_T> ReturnBinCentersFromBinEdges(std::vector<FLOAT_T> BinEdges);

  // ========================================================================================================================================================================
  // Protected virtual functions which are calculation implementation agnostic

  // ========================================================================================================================================================================
  // Basic protected variables required for oscillation probability calculation

 private:

  // ========================================================================================================================================================================
  // Basic private variables required for oscillation probability calculation

  /**
   * @brief The FileName which the binning is read from
   */
  std::string FileName;

  /**
   * @brief The name of the histogram which contains the Energy axis binning
   */
  std::string EnergyAxisHistName;

  /**
   * @brief The name of the histogram which contains the CosineZ axis binning
   */
  std::string CosineZAxisHistName;

  /**
   * @brief A vector of Energy axis bin edges
   */
  std::vector<FLOAT_T> EnergyAxisBinEdges;
  /**
   * @brief A vector of CosineZ axis bin edges
   */
  std::vector<FLOAT_T> CosineZAxisBinEdges;
  /**
   * @brief A vector of Energy axis bin centers
   */
  std::vector<FLOAT_T> EnergyAxisBinCenters;
  /**
   * @brief A vector of CosineZ axis bin centers
   */
  std::vector<FLOAT_T> CosineZAxisBinCenters;
};

#endif
