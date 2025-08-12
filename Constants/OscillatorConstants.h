#ifndef __OSCILLATOR_CONSTANTS__
#define __OSCILLATOR_CONSTANTS__

#define DUMMYVAL -999

#include <math.h>

#include "TFile.h"
#include "TH1.h"

#if UseDoubles==1
using FLOAT_T = double;
#else
using FLOAT_T = float;
#endif

#include <string>
#include <iostream>
#include <vector>

/**
 * @file OscillatorConstants.h
 */

namespace NuOscillator
{
  /**
   * @brief Different verbosity levels for console output
   */
  enum Verbosity{NONE=0,INFO=1,VERBOSE=2};
  
  /**
   * @brief Different neutrino flavours currently supported within the analysis
   *
   * If more need to be added, no changes should be required outside of this file
   */
  enum NeutrinoFlavours{kElectron=1,kMuon=2,kTau=3,kSterile1=4,kSterile2=5,kSterile3=6,nNeutrinoFlavours=7};
  
  /**
   * @brief Enum which fixes the ordering of the generated and detected neutrino flavours in the #OscillationChannel structure
   */
  enum {kNuFlavour_Generated=0,kNuFlavour_Detected=1,nNuFlavours=2};
  
  /**
   * @brief Structure which defines the oscillation channel generated and detected neutrino flavours
   */
  struct OscillationChannel{
    int GeneratedFlavour;
    int DetectedFlavour;
  };
  
  /**
   * @brief Structure to contain all information about the neutrino type, oscillation channel, Energy and CosineZ used to calculate a specific probability
   */ 
  struct OscillationProbability{
    int NuType;
    OscillationChannel OscChan;
    FLOAT_T Energy;
    FLOAT_T CosineZ;
    FLOAT_T Probability;
  };
  
}

/**
 * @brief Returns the basic oscillation parameters.
 * @return Vector of oscillation parameters.
 */
inline std::vector<FLOAT_T> ReturnOscParams_Basic() {
  std::vector<FLOAT_T> OscParams_Basic = {3.07e-1,5.28e-1,2.18e-2,7.53e-5,2.509e-3,-1.601};
  return OscParams_Basic;
}

/**
 * @brief Returns the oscillation parameters for atmospheric neutrinos, with production height.
 *
 * @return Vector of oscillation parameters.
 */
inline std::vector<FLOAT_T> ReturnOscParams_Atm() {
  std::vector<FLOAT_T> OscParams_Atm = {3.07e-1,5.28e-1,2.18e-2,7.53e-5,2.509e-3,-1.601,25.0};
  return OscParams_Atm;
}

/**
 * @brief Returns the oscillation parameters for beam neutrinos, with baseline and density. Without electron density
 *
 * @return Vector of oscillation parameters.
 */
inline std::vector<FLOAT_T> ReturnOscParams_Beam_woYe() {
  std::vector<FLOAT_T> OscParams_Beam_woYe = {3.07e-1,5.28e-1,2.18e-2,7.53e-5,2.509e-3,-1.601,250.0,2.6};
  return OscParams_Beam_woYe;
}

/**
 * @brief Returns the oscillation parameters for beam neutrinos, with baseline and density. With electron density
 *
 * @return Vector of oscillation parameters.
 */
inline std::vector<FLOAT_T> ReturnOscParams_Beam_wYe() {
  std::vector<FLOAT_T> OscParams_Beam_wYe = {3.07e-1,5.28e-1,2.18e-2,7.53e-5,2.509e-3,-1.601,250.0,2.6,0.5};
  return OscParams_Beam_wYe;
}

/**
 * @brief Returns the oscillation parameters for beam neutrinos, with baseline and density. With electron density and decoherence parameters.
 *
 * @return Vector of oscillation parameters.
 */
inline std::vector<FLOAT_T> ReturnOscParams_Beam_wYe_wDeco() {
  std::vector<FLOAT_T> OscParams_Beam_wYe_wDeco = {3.07e-1,5.28e-1,2.18e-2,7.53e-5,2.509e-3,-1.601,250.0,2.6,0.5,9.48e-18,2.0,1.0}; //9.48e-18
  return OscParams_Beam_wYe_wDeco;
}

/**
 * @brief Returns the oscillation parameters for beam neutrinos, with baseline and density. With electron density and Lorentz invariant parameters.
 *
 * @return Vector of oscillation parameters.
 */
inline std::vector<FLOAT_T> ReturnOscParams_Beam_wYe_wLIV() {
  std::vector<FLOAT_T> OscParams_Beam_wYe_wLIV = {3.07e-1,5.28e-1,2.18e-2,7.53e-5,2.509e-3,-1.601,250.0,2.6,0.5,1.0e-23,0.0,1.0e-23,0.0,0.0};
  return OscParams_Beam_wYe_wLIV;
}

/**
 * @brief Return vector of all config names for each oscillation engine which has been enabled
 *
 * @return Vector of paths to config files
 */
inline std::vector<std::string> ReturnKnownConfigs() {
  std::vector<std::string> ConfigNames;

#if UseCUDAProb3 == 1
  ConfigNames.push_back("./NuOscillatorConfigs/Binned_CUDAProb3.yaml");
#endif

#if UseCUDAProb3Linear == 1
  ConfigNames.push_back("./NuOscillatorConfigs/Binned_CUDAProb3Linear.yaml");
#endif

#if UseProbGPULinear == 1
  ConfigNames.push_back("./NuOscillatorConfigs/Binned_ProbGPULinear.yaml");
#endif

#if UseProb3ppLinear == 1
  ConfigNames.push_back("./NuOscillatorConfigs/Binned_Prob3ppLinear.yaml");
#endif

#if UseNuFASTLinear == 1
  ConfigNames.push_back("./NuOscillatorConfigs/Binned_NuFASTLinear.yaml");
#endif

#if UseNuSQUIDSLinear == 1
  ConfigNames.push_back("./NuOscillatorConfigs/Binned_NuSQUIDSLinear.yaml");
#endif
  
#if UseOscProb == 1
  ConfigNames.push_back("./NuOscillatorConfigs/Binned_OscProb.yaml");
#endif  

  return ConfigNames;
}

/**
 * @brief Convert a neutrino flavour string to integer
 *
 * @return Enum value in #NeutrinoFlavours
 */
inline int NeutrinoFlavour_StrToInt(const std::string& NuFlav) {
  if (NuFlav == "Electron" || NuFlav == "electron") {
    return NuOscillator::kElectron;
  } else if (NuFlav == "Muon" || NuFlav == "muon") {
    return NuOscillator::kMuon;
  } else if (NuFlav == "Tau" || NuFlav == "tau") {
    return NuOscillator::kTau;
  } else if (NuFlav == "Sterile1" || NuFlav == "sterile1") {
    return NuOscillator::kSterile1;
  } else if (NuFlav == "Sterile2" || NuFlav == "sterile2") {
    return NuOscillator::kSterile2;
  } else if (NuFlav == "Sterile3" || NuFlav == "sterile3") {
    return NuOscillator::kSterile3;
  } else {
    std::cerr << "Could not convert input string:" << NuFlav << " to known enum value in NeutrinoFlavours" << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  
  return -1;
}

/**
 * @brief Convert a neutrino flavour integer to string
 *
 * Inverse of NeutrinoFlavour_StrToInt()
 *
 * @return Name of neutrino flavour
 */
inline std::string NeutrinoFlavour_IntToStr(const int NuFlav) {
  switch (NuFlav) {
  case NuOscillator::kElectron: 
    return "Electron";
  case NuOscillator::kMuon:
    return "Muon";
  case NuOscillator::kTau:
    return "Tau";
  case NuOscillator::kSterile1:
    return "Sterile1";
  case NuOscillator::kSterile2:
    return "Sterile2";
  case NuOscillator::kSterile3:
    return "Sterile3";
  default:
    std::cerr << "Recieved unknown NeutrinoFlavour:" << NuFlav << " which is inconsistent with those in enum 'NeutrinoFlavours'" << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  return "";
}

/**
 * @brief Return the Verbosity enum value correpsonding to a particular string
 * 
 * @param Verbosity String to convert to enum value
 *
 * @return Enum value describing the verbosity level
 */
inline int Verbosity_StrToInt(const std::string& Verbosity) {
  if (Verbosity == "NONE") {
    return NuOscillator::NONE;
  } else if (Verbosity == "INFO") {
    return NuOscillator::INFO;
  } else if (Verbosity == "VERBOSE") {
    return NuOscillator::VERBOSE;
  } else {
    std::cerr << "Invalid verbosity provided:" << Verbosity << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  
  return -1;
}

/**
 * @brief Take an input string formatted as 'GeneratedNeutrinoFlavour:DetectedNeutrinoFlavour' and return an NuOscillator::OscillationChannel() structure
 *
 * @param String formatted as 'GeneratedNeutrinoFlavour:DetectedNeutrinoFlavour'
 *
 * @return NuOscillator::OscillationChannel() structure with generated and detected neutrino flavours
 */
inline NuOscillator::OscillationChannel ReturnOscillationChannel(const std::string& InputString) {
  std::string GeneratedFlavour = "";
  std::string DetectedFlavour = "";

  size_t Delimiter = InputString.find(":");
  if (Delimiter != std::string::npos) {
    GeneratedFlavour = InputString.substr(0,Delimiter);
    DetectedFlavour = InputString.substr(Delimiter+1,InputString.size());
  } else {
    std::cerr << "Expected a string formatted as: 'GeneratedNeutrinoFlavour:DetectedNeutrinoFlavour'" << std::endl;
    std::cerr << "Received:" << InputString << std::endl;
  }
  
  NuOscillator::OscillationChannel OscChannel = {NeutrinoFlavour_StrToInt(GeneratedFlavour),NeutrinoFlavour_StrToInt(DetectedFlavour)};
  return OscChannel;
}

/**
 * @brief Generate vector of logarithmically spaced points 
 *
 * @param Emin lower limit
 * @param Emax upper limit
 * @param nDiv Number of divisions
 * @return Vector of logarithmically spaced points between Emin and Emax with nDiv divisions
 */
inline std::vector<FLOAT_T> logspace(FLOAT_T Emin, FLOAT_T  Emax, int nDiv) {
  if (nDiv==0) {
    std::cerr << "Requested log spacing distribution with 0 divisions" << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  std::vector<FLOAT_T> logpoints(nDiv+1, 0.0);
  logpoints[0]=Emin;

  if (Emin == 0.) {
    Emin = 0.01;
  }

  FLOAT_T Emin_log,Emax_log;
  Emin_log = log10(Emin);
  Emax_log = log10(Emax);

  FLOAT_T step_log = (Emax_log - Emin_log)/FLOAT_T(nDiv);

  FLOAT_T EE = Emin_log+step_log;

  for (int i=1; i<nDiv; i++) {
    logpoints[i] = pow(10.,EE);
    EE += step_log;
  }

  logpoints[nDiv]=Emax;

  return logpoints;
}

/**
 * @brief Generate vector of linearly spaced points 
 *
 * @param Emin lower limit
 * @param Emax upper limit
 * @param nDiv Number of divisions
 * @return Vector of linearly spaced points between Emin and Emax with nDiv divisions
 */
inline std::vector<FLOAT_T> linspace(FLOAT_T Emin, FLOAT_T Emax, int nDiv) {
  if (nDiv==0) {
    std::cerr << "Requested linear spacing distribution with 0 divisions" << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  std::vector<FLOAT_T> linpoints(nDiv+1, 0.0);

  FLOAT_T step_lin = (Emax - Emin)/FLOAT_T(nDiv);

  FLOAT_T EE = Emin;

  for (int i=0; i<nDiv; i++) {
    if (fabs(EE)<1e-6) {EE = 0.;}

    linpoints[i] = EE;
    EE += step_lin;
  }

  linpoints[nDiv] = Emax;

  return linpoints;
}

/**
 * @brief Read bin edges from input template histogram 
 *
 * @param TFileName File name
 * @param HistogramName Histogram name
 * @param Verbose Verbosity level
 * @return Vector of bin edges
 */
inline std::vector<FLOAT_T> ReadBinEdgesFromFile(std::string TFileName, std::string HistogramName, int Verbose=NuOscillator::Verbosity::NONE) {
  std::vector<FLOAT_T> BinEdges;

  TFile* File = new TFile(TFileName.c_str());
  if (!File || File->IsZombie()) {
    std::cerr << "Could not find file:" << TFileName << std::endl;
    throw std::runtime_error("Invalid input file.");
  }

  TH1* Histogram = (TH1*)File->Get(HistogramName.c_str());
  if (!Histogram) {
    std::cerr << "Could not find Histogram:" << HistogramName << " in File:" << TFileName << std::endl;
    throw std::runtime_error("Invalid input file.");
  }

  BinEdges.resize(Histogram->GetNbinsX()+1);
  for (int iBin=0;iBin<=Histogram->GetNbinsX();iBin++) {
    BinEdges[iBin] = Histogram->GetBinLowEdge(iBin+1);
  }

  File->Close();
  delete File;

  if (Verbose >= NuOscillator::INFO) {
    std::cout << "Bin edges successfully read from File:" << TFileName << " , Histogram:" << HistogramName << " :=" << std::endl;
    for (size_t i=0;i<BinEdges.size();i++) {
      std::cout << BinEdges[i] << ", ";
    }
    std::cout << std::endl;
  }

  return BinEdges;
}

/**
 * @brief Return the bin centers given the bin edges of a template histogram
 *
 * @param BinEdges Vector of bin edges
 * @return Vector of bin centers
 */
inline std::vector<FLOAT_T> ReturnBinCentersFromBinEdges(std::vector<FLOAT_T> BinEdges) {
  int nBins = BinEdges.size()-1;
  std::vector<FLOAT_T> BinCenters = std::vector<FLOAT_T>(nBins);

  for (int iBin=0;iBin<nBins;iBin++) {
    BinCenters[iBin] = (BinEdges[iBin]+BinEdges[iBin+1])/2.0;
  }

  return BinCenters;
}

#endif
