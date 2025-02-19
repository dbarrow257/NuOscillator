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

  OscillatorBinned(YAML::Node Config_);

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
  const FLOAT_T* ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal=DUMMYVAL) override;

  /**
   * @brief Return a vector of bin edges used for oscillation probability plotting
   *
   * Return the binning used to group the energy/cosineZ when associating an MC event to a probability
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

  void Initialise();

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
