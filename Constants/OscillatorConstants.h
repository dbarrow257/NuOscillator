#ifndef __OSCILLATOR_CONSTANTS__
#define __OSCILLATOR_CONSTANTS__

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

//DB
enum NeutrinoFlavours{kElectron=1,kMuon=2,kTau=3,kSterile1=4,kSterile2=5,kSterile3=6,nNeutrinoFlavours=7};

//DB
enum {kNuFlavour_Generated=0,kNuFlavour_Detected=1,nNuFlavours=2};

//DB
struct OscillationChannel{
  int GeneratedFlavour;
  int DetectedFlavour;
};

//DB
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

//DB
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

//DB
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
