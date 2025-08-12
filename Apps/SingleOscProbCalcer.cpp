#include "Constants/OscillatorConstants.h"
#include "OscProbCalcer/OscProbCalcerBase.h"

#include "OscProbCalcer/OscProbCalcerFactory.h"

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
  std::string OscProbCalcerConfigname = argv[1];
  
  bool PrintWeights = true;
  std::unordered_map<std::string, FLOAT_T> OscillationParameters = ReturnOscParamsFromConfig(YAML::LoadFile(OscProbCalcerConfigname));

  std::vector<FLOAT_T> EnergyArray = logspace(0.1,100.,1e3);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,15);
  
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting setup in executable" << std::endl;

  std::cout << "========================================================" << std::endl;
  std::cout << "Initialising " << OscProbCalcerConfigname << std::endl;

  OscProbCalcerFactory* OscProbCalcFactory = new OscProbCalcerFactory();
  OscProbCalcerBase* Calcer = OscProbCalcFactory->CreateOscProbCalcer(OscProbCalcerConfigname);
  delete OscProbCalcFactory;
  std::cout << "========================================================" << std::endl;
  std::cout << "Setting up Oscillators" << std::endl;

  for (auto Parameter : OscillationParameters) {
    Calcer->DefineParameter(Parameter.first,&OscillationParameters[Parameter.first]);
  }
  
  Calcer->SetEnergyArray(EnergyArray);
  if (!Calcer->ReturnCosineZIgnored()) {
    Calcer->SetCosineZArray(CosineZArray);
  }

  Calcer->Setup();

  std::cout << "Finished setup in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting reweight in executable" << std::endl;

  // Reweight and calculate oscillation probabilities
  Calcer->Reweight();

  if (PrintWeights) {
    std::vector<NuOscillator::OscillationProbability> OscProbs = Calcer->ReturnProbabilities();
    for (int iOscProb=0;iOscProb<(int)OscProbs.size();iOscProb++) {
      std::cout << iOscProb << " " << OscProbs[iOscProb].NuType << " " << OscProbs[iOscProb].OscChan.GeneratedFlavour << " " << OscProbs[iOscProb].OscChan.DetectedFlavour << " " << OscProbs[iOscProb].Energy << " " << OscProbs[iOscProb].CosineZ << " " << OscProbs[iOscProb].Probability << std::endl;
    }
  }
  
  std::cout << "Finished reweight in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
}
