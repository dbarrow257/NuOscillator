#ifndef __OSCILLATOR_CONSTANTS__
#define __OSCILLATOR_CONSTANTS__

#define DUMMYVAL -999

#ifdef UseDoubles
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

/**
 * @brief Convert a neutrino flavour string to integer
 *
 * @return Enum value in #NeutrinoFlavours
 */
inline int NeutrinoFlavour_StrToInt(std::string NuFlav) {
  if (NuFlav == "Electron" || NuFlav == "electron") {
    return kElectron;
  } else if (NuFlav == "Muon" || NuFlav == "muon") {
    return kMuon;
  } else if (NuFlav == "Tau" || NuFlav == "tau") {
    return kTau;
  } else if (NuFlav == "Sterile1" || NuFlav == "sterile1") {
    return kSterile1;
  } else if (NuFlav == "Sterile2" || NuFlav == "sterile2") {
    return kSterile2;
  } else if (NuFlav == "Sterile3" || NuFlav == "sterile3") {
    return kSterile3;
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
    case kElectron: 
      return "Electron";
    case kMuon:
      return "Muon";
    case kTau:
      return "Tau";
    case kSterile1:
      return "Sterile1";
    case kSterile2:
      return "Sterile2";
    case kSterile3:
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
    return NONE;
  } else if (Verbosity == "INFO") {
    return INFO;
  } else {
    std::cerr << "Invalid verbosity provided:" << Verbosity << std::endl;
    throw;
  }
  
  return -1;
}

/**
 * @brief Take an input string formatted as 'GeneratedNeutrinoFlavour:DetectedNeutrinoFlavour' and return an OscillationChannel() structure
 *
 * @param String formatted as 'GeneratedNeutrinoFlavour:DetectedNeutrinoFlavour'
 *
 * @return OscillationChannel() structure with generated and detected neutrino flavours
 */
inline OscillationChannel ReturnOscillationChannel(std::string InputString) {
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
  
  OscillationChannel OscChannel = {NeutrinoFlavour_StrToInt(GeneratedFlavour),NeutrinoFlavour_StrToInt(DetectedFlavour)};
  return OscChannel;
}

#endif
