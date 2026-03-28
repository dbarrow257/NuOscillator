#include "Oscillator/OscillatorFactory.h"

#include "Constants/OscillatorConstants.h"

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
  std::string OscillationParameterConfig = argv[1];
  std::unordered_map<std::string, FLOAT_T> OscillationParameters_Global = ReturnOscParamsFromConfig(YAML::LoadFile(OscillationParameterConfig));

  bool PrintWeights = true;
  
  std::vector<FLOAT_T> EnergyArray = logspace(0.1,100.,1e3);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,1e3);

  std::cout << "========================================================" << std::endl;
  std::cout << "Starting setup in executable" << std::endl;

  std::vector<OscillatorBase*> Oscillators;
  OscillatorFactory* OscFactory = new OscillatorFactory();
  OscillatorBase* Oscillator;

  //Get the standard set of config names
  std::vector<std::string> ConfigNames = ReturnKnownConfigs();

  std::cout << "========================================================" << std::endl;
  std::cout << "Setting up Oscillators" << std::endl;

  std::vector< std::unordered_map<std::string, FLOAT_T> > OscillationParameters_Oscillators(ConfigNames.size());
  
  for (size_t iConfig=0;iConfig<ConfigNames.size();iConfig++) {
    std::cout << "========================================================" << std::endl;
    std::cout << "Initialising " << ConfigNames[iConfig] << std::endl;
    
    //Create OscillatorBase* object from YAML config
    Oscillator = OscFactory->CreateOscillator(ConfigNames[iConfig]);

    //Check if the Energy and CosineZ evaluation points have been set in the constructor of the object (i.e. Binned where the templates have been picked up by the constructor)
    //or if we need to set them after the fact (i.e. unbinned where the points may change depending on the events etc.)
    if (!Oscillator->EvalPointsSetInConstructor()) {
      Oscillator->SetEnergyArrayInCalcer(EnergyArray);
	
      //Check if we also need to set the CosineZ binning
      if (!Oscillator->CosineZIgnored()) {
	Oscillator->SetCosineZArrayInCalcer(CosineZArray);
      }
    }

    //Overwrite the oscillation parameters in the Oscillator config from the config passed through the arguments
    OscillationParameters_Oscillators[iConfig] = ReturnOscParamsFromConfig(YAML::LoadFile(ConfigNames[iConfig]));
    for (auto Parameter : OscillationParameters_Global) {
      if (OscillationParameters_Oscillators[iConfig].count(Parameter.first)) {
	OscillationParameters_Oscillators[iConfig][Parameter.first] = Parameter.second;
      }
    }

    //Set the oscillation parameters in the Oscillator
    for (auto Parameter : OscillationParameters_Oscillators[iConfig]) {
      Oscillator->DefineParameter(Parameter.first, &OscillationParameters_Oscillators[iConfig][Parameter.first]);
    }
    
    // Setup propagators    
    Oscillator->Setup();
    
    //Append OscillatorBase* object to the vector
    Oscillators.push_back(Oscillator);
  }
  
  std::cout << "Finished setup in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting reweight in executable" << std::endl;

  // Reweight and calculate oscillation probabilities
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    std::cout << "Performing reweight in Oscillator: " << iOsc << "/" << Oscillators.size() << std::endl;
    
    Oscillators[iOsc]->CalculateProbabilities();
  
    if (PrintWeights) {
      Oscillators[iOsc]->PrintWeights();
    }
  }

  std::cout << "Finished reweight in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
}
