#include "Oscillator/OscillatorFactory.h"

#include "Constants/OscillatorConstants.h"

#include <iostream>
#include <math.h>
#include <chrono>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

int main() {
  bool PrintWeights = false;

  std::vector<FLOAT_T> EnergyArray = logspace(0.1,100.,1e3);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,1e3);

  std::cout << "========================================================" << std::endl;
  std::cout << "Starting setup in executable" << std::endl;

  std::vector<OscillatorBase*> Oscillators;
  OscillatorFactory* OscFactory = new OscillatorFactory();
  OscillatorBase* Oscillator;

  //Get the standard set of config names
  std::vector<std::string> ConfigNames = ReturnKnownConfigs();
    
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

    //Append OscillatorBase* object to the vector
    Oscillators.push_back(Oscillator);
  }

  std::cout << "========================================================" << std::endl;
  std::cout << "Setting up Oscillators" << std::endl;

  // Setup propagators
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    Oscillators[iOsc]->Setup();
  }
  
  std::cout << "Finished setup in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting reweight in executable" << std::endl;

  // Reweight and calculate oscillation probabilities
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    std::cout << "Performing reweight in Oscillator: " << iOsc << "/" << Oscillators.size() << std::endl;
    
    // These don't have to be explicilty beam or atmospheric specific, all they have to be is equal to the number of oscillation parameters expected by the implementation
    // If you have some NSO calculater, then it will work providing the length of the vector of oscillation parameters is equal to the number of expected oscillation parameters
    if (Oscillators[iOsc]->ReturnNOscParams() == (int)NuOscillator::OscParams_Beam_woYe.size()) {
      Oscillators[iOsc]->CalculateProbabilities(NuOscillator::OscParams_Beam_woYe);
    } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)NuOscillator::OscParams_Beam_wYe.size()) {
	Oscillators[iOsc]->CalculateProbabilities(NuOscillator::OscParams_Beam_wYe); 
    } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)NuOscillator::OscParams_Atm.size()) {
      Oscillators[iOsc]->CalculateProbabilities(NuOscillator::OscParams_Atm);
    } else {
      std::cerr << "Did not find viable oscillation parameters to hand to the oscillation probability calculater" << std::endl;
      std::cerr << "Oscillator->ReturnNOscParams():" << Oscillators[iOsc]->ReturnNOscParams() << std::endl;
      throw;
    }
  
    if (PrintWeights) {
      Oscillators[iOsc]->PrintWeights();
    }
  }

  std::cout << "Finished reweight in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
}
