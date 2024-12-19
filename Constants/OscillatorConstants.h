#ifndef __OSCILLATOR_CONSTANTS__
#define __OSCILLATOR_CONSTANTS__

#define DUMMYVAL -999

#include <math.h>

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
  enum Verbosity{NONE=0,INFO=1};
  
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

inline std::vector<FLOAT_T> ReturnOscParams_Atm() {
  std::vector<FLOAT_T> OscParams_Atm = {3.07e-1,5.28e-1,2.18e-2,7.53e-5,2.509e-3,-1.601,25.0};
  return OscParams_Atm;
}

inline std::vector<FLOAT_T> ReturnOscParams_Beam_woYe() {
  std::vector<FLOAT_T> OscParams_Beam_woYe = {3.07e-1,5.28e-1,2.18e-2,7.53e-5,2.509e-3,-1.601,250.0,2.6};
  return OscParams_Beam_woYe;
}

inline std::vector<FLOAT_T> ReturnOscParams_Beam_wYe() {
  std::vector<FLOAT_T> OscParams_Beam_wYe = {3.07e-1,5.28e-1,2.18e-2,7.53e-5,2.509e-3,-1.601,250.0,2.6,0.5};
  return OscParams_Beam_wYe;
}

/**
 * @brief Return vector of all config names for each oscillation engine which has been enabled
 *
 * @return Vector of paths to config files
 */
inline std::vector<std::string> ReturnKnownConfigs() {
  std::vector<std::string> ConfigNames;

#if UseCUDAProb3 == 1
  ConfigNames.push_back("./Configs/Binned_CUDAProb3.yaml");
#endif

#if UseCUDAProb3Linear == 1
  ConfigNames.push_back("./Configs/Binned_CUDAProb3Linear.yaml");
#endif

#if UseProbGPULinear == 1
  ConfigNames.push_back("./Configs/Binned_ProbGPULinear.yaml");
#endif

#if UseProb3ppLinear == 1
  ConfigNames.push_back("./Configs/Binned_Prob3ppLinear.yaml");
#endif

#if UseNuFASTLinear == 1
  ConfigNames.push_back("./Configs/Binned_NuFASTLinear.yaml");
#endif

#if UseNuSQUIDSLinear == 1
  ConfigNames.push_back("./Configs/Unbinned_NuSQUIDSLinear.yaml");
#endif  

  return ConfigNames;
}

/**
 * @brief Convert a neutrino flavour string to integer
 *
 * @return Enum value in #NeutrinoFlavours
 */
inline int NeutrinoFlavour_StrToInt(std::string NuFlav) {
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
    throw;
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
inline std::string NeutrinoFlavour_IntToStr(int NuFlav) {
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
    throw;
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
inline int Verbosity_StrToInt(std::string Verbosity) {
  if (Verbosity == "NONE") {
    return NuOscillator::NONE;
  } else if (Verbosity == "INFO") {
    return NuOscillator::INFO;
  } else {
    std::cerr << "Invalid verbosity provided:" << Verbosity << std::endl;
    throw;
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
inline NuOscillator::OscillationChannel ReturnOscillationChannel(std::string InputString) {
  std::string GeneratedFlavour = "";
  std::string DetectedFlavour = "";

  size_t Delimiter = InputString.find(":");
  if (Delimiter != std::string::npos) {
    GeneratedFlavour = InputString.substr(0,Delimiter);
    DetectedFlavour = InputString.substr(Delimiter+1,InputString.size());
  } else {
    std::cerr << "Expected a string formatted as: 'GeneratedNeutrinoFlavour:DetectedNeutrinoFlavour'" << std::endl;
    std::cerr << "Recieved:" << InputString << std::endl;
  }
  
  NuOscillator::OscillationChannel OscChannel = {NeutrinoFlavour_StrToInt(GeneratedFlavour),NeutrinoFlavour_StrToInt(DetectedFlavour)};
  return OscChannel;
}

inline std::vector<FLOAT_T> logspace(FLOAT_T Emin, FLOAT_T  Emax, int nDiv) {
  if (nDiv==0) {
    std::cerr << "Requested log spacing distribution with 0 divisions" << std::endl;
    throw;
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

inline std::vector<FLOAT_T> linspace(FLOAT_T Emin, FLOAT_T Emax, int nDiv) {
  if (nDiv==0) {
    std::cerr << "Requested linear spacing distribution with 0 divisions" << std::endl;
    throw;
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

#endif
