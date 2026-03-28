#include "Oscillator/OscillatorFactory.h"

#include "Constants/OscillatorConstants.h"

#include "TStyle.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"

#include <iostream>
#include <math.h>
#include <chrono>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

int main(int argc, char **argv) {
  if (argc != 2) {
    std::cerr << argv[0] << " InputConfig.yaml" << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  std::string ConfigName = argv[1];

  gStyle->SetOptStat(0);
  
  bool PrintWeights = true;
  std::unordered_map<std::string, FLOAT_T> OscillationParameters = ReturnOscParamsFromConfig(YAML::LoadFile(ConfigName));
 
  //Don't plot by default
  bool Plot = false;
  
  std::vector<FLOAT_T> EnergyArray = logspace(0.1,100.,1e3);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,15);

  std::cout << "========================================================" << std::endl;
  std::cout << "Starting setup in executable" << std::endl;

  OscillatorFactory* OscFactory = new OscillatorFactory();
  OscillatorBase* Oscillator;

  std::cout << "========================================================" << std::endl;
  std::cout << "Initialising " << ConfigName << std::endl;
  
  //Create OscillatorBase* object from YAML config
  Oscillator = OscFactory->CreateOscillator(ConfigName);
  
  //Check if the Energy and CosineZ evaluation points have been set in the constructor of the object (i.e. Binned where the templates have been picked up by the constructor)
  //or if we need to set them after the fact (i.e. unbinned where the points may change depending on the events etc.)
  if (!Oscillator->EvalPointsSetInConstructor()) {
    Oscillator->SetEnergyArrayInCalcer(EnergyArray);
    
    //Check if we also need to set the CosineZ binning
    if (!Oscillator->CosineZIgnored()) {
      Oscillator->SetCosineZArrayInCalcer(CosineZArray);
    }
  }

  std::cout << "========================================================" << std::endl;
  std::cout << "Setting up Oscillators" << std::endl;

  // Legacy mode should not use the DefineParameter functionality - the following will throw error
  for (auto Parameter : OscillationParameters) {
    Oscillator->DefineParameter(Parameter.first,&OscillationParameters[Parameter.first]);
  }
  
  Oscillator->Setup();
  
  std::cout << "Finished setup in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting reweight in executable" << std::endl;

  std::vector<FLOAT_T> OscParams = {3.07e-1,5.28e-1,2.18e-2,7.53e-5,2.509e-3,-1.601,250.0,2.6,0.5};
  Oscillator->CalculateProbabilities(OscParams);
  
  if (PrintWeights) {
    Oscillator->PrintWeights();
  }

  std::cout << "Finished reweight in executable" << std::endl;
  std::cout << "========================================================" << std::endl;

}
