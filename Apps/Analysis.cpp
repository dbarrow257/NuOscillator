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
  bool PrintWeights = true;

  std::vector<FLOAT_T> OscParams_Atm(7);
  OscParams_Atm[0] = 3.07e-1;
  OscParams_Atm[1] = 5.28e-1;
  OscParams_Atm[2] = 2.18e-2;
  OscParams_Atm[3] = 7.53e-5;
  OscParams_Atm[4] = 2.509e-3;
  OscParams_Atm[5] = -1.601;
  OscParams_Atm[6] = 25.0;

  std::vector<FLOAT_T> OscParams_Beam(8);
  OscParams_Beam[0] = 3.07e-1;
  OscParams_Beam[1] = 5.28e-1;
  OscParams_Beam[2] = 2.18e-2;
  OscParams_Beam[3] = 7.53e-5;
  OscParams_Beam[4] = 2.509e-3;
  OscParams_Beam[5] = -1.601;
  OscParams_Beam[6] = 250.0;
  OscParams_Beam[7] = 2.6;

  std::vector<FLOAT_T> EnergyArray = logspace(0.1,100.,1e3);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,1e3);

  std::cout << "========================================================" << std::endl;
  std::cout << "Starting setup in executable" << std::endl;

  std::vector<OscillatorBase*> Oscillators;
  std::vector<std::string> ConfigNames;

  OscillatorFactory* OscFactory = new OscillatorFactory();
  OscillatorBase* Oscillator;

#if UseCUDAProb3 == 1
  ConfigNames.push_back("./Configs/Binned_CUDAProb3.yaml");
  ConfigNames.push_back("./Configs/Unbinned_CUDAProb3.yaml");
#endif

#if UseCUDAProb3Linear == 1
  ConfigNames.push_back("./Configs/Binned_CUDAProb3Linear.yaml");
  ConfigNames.push_back("./Configs/Unbinned_CUDAProb3Linear.yaml");
#endif

#if UseProbGPULinear == 1
  ConfigNames.push_back("./Configs/Binned_ProbGPULinear.yaml");
  ConfigNames.push_back("./Configs/Unbinned_ProbGPULinear.yaml");
#endif

#if UseProb3ppLinear == 1
  ConfigNames.push_back("./Configs/Binned_Prob3ppLinear.yaml");
  ConfigNames.push_back("./Configs/Unbinned_Prob3ppLinear.yaml");
#endif

  //Alternative option to show how all information can be held in a single YAML file rather than using a preset
  //ConfigNames.push_back("./Configs/CUDAProb3_Binned-SelfContainedFile.yaml");

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
    // These don't have to be explicilty beam or atmospheric specific, all they have to be is equal to the number of oscillation parameters expected by the implementation
    // If you have some NSO calculater, then it will work providing the length of the vector of oscillation parameters is equal to the number of expected oscillation parameters
    if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Beam.size()) {
      Oscillators[iOsc]->CalculateProbabilities(OscParams_Beam); 
    } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Atm.size()) {
      Oscillators[iOsc]->CalculateProbabilities(OscParams_Atm);
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
