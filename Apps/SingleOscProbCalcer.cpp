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
    std::cerr << "./SingleOscProbCalcer InputConfig.yaml" << std::endl;
    throw;
  }
  std::string OscProbCalcerConfigname = argv[1];
  
  bool PrintWeights = false;

  std::vector<FLOAT_T> OscParams_Atm(8);
  OscParams_Atm[0] = 3.07e-1;
  OscParams_Atm[1] = 5.28e-1;
  OscParams_Atm[2] = 2.18e-2;
  OscParams_Atm[3] = 7.53e-5;
  OscParams_Atm[4] = 2.509e-3;
  OscParams_Atm[5] = -1.601;
  OscParams_Atm[6] = 25.0;

  std::vector<FLOAT_T> OscParams_Beam_woYe(8);
  OscParams_Beam_woYe[0] = 3.07e-1;
  OscParams_Beam_woYe[1] = 5.28e-1;
  OscParams_Beam_woYe[2] = 2.18e-2;
  OscParams_Beam_woYe[3] = 7.53e-5;
  OscParams_Beam_woYe[4] = 2.509e-3;
  OscParams_Beam_woYe[5] = -1.601;
  OscParams_Beam_woYe[6] = 250.0;
  OscParams_Beam_woYe[7] = 2.6;

  std::vector<FLOAT_T> OscParams_Beam_wYe(9);
  OscParams_Beam_wYe[0] = 3.07e-1;
  OscParams_Beam_wYe[1] = 5.28e-1;
  OscParams_Beam_wYe[2] = 2.18e-2;
  OscParams_Beam_wYe[3] = 7.53e-5;
  OscParams_Beam_wYe[4] = 2.509e-3;
  OscParams_Beam_wYe[5] = -1.601;
  OscParams_Beam_wYe[6] = 250.0;
  OscParams_Beam_wYe[7] = 2.6;
  OscParams_Beam_wYe[8] = 0.5;
  
  std::vector<FLOAT_T> EnergyArray = logspace(0.1,100.,1e3);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,1e3);

  std::cout << "========================================================" << std::endl;
  std::cout << "Starting setup in executable" << std::endl;

  std::cout << "========================================================" << std::endl;
  std::cout << "Initialising " << OscProbCalcerConfigname << std::endl;

  OscProbCalcerFactory* OscProbCalcFactory = new OscProbCalcerFactory();
  OscProbCalcerBase* Calcer = OscProbCalcFactory->CreateOscProbCalcer(OscProbCalcerConfigname);

  std::cout << "========================================================" << std::endl;
  std::cout << "Setting up Oscillators" << std::endl;

  Calcer->SetEnergyArray(EnergyArray);
  if (!Calcer->ReturnCosineZIgnored()) {
    Calcer->SetCosineZArray(CosineZArray);
  }
  Calcer->Setup();

  std::cout << "Finished setup in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting reweight in executable" << std::endl;

  // Reweight and calculate oscillation probabilities

  // These don't have to be explicilty beam or atmospheric specific, all they have to be is equal to the number of oscillation parameters expected by the implementation
  // If you have some NSO calculater, then it will work providing the length of the vector of oscillation parameters is equal to the number of expected oscillation parameters
  if (Calcer->ReturnNOscParams() == (int)OscParams_Beam_woYe.size()) {
    Calcer->Reweight(OscParams_Beam_woYe);
  } else if (Calcer->ReturnNOscParams() == (int)OscParams_Beam_wYe.size()) {
    Calcer->Reweight(OscParams_Beam_wYe); 
  } else if (Calcer->ReturnNOscParams() == (int)OscParams_Atm.size()) {
    Calcer->Reweight(OscParams_Atm);
  } else {
    std::cerr << "Did not find viable oscillation parameters to hand to the oscillation probability calculater" << std::endl;
    std::cerr << "Oscillator->ReturnNOscParams():" << Calcer->ReturnNOscParams() << std::endl;
    throw;
  }

  if (PrintWeights) {
    std::vector<NuOscillator::OscillationProbability> OscProbs = Calcer->ReturnProbabilities();
    for (int iOscProb=0;iOscProb<(int)OscProbs.size();iOscProb++) {
      std::cout << iOscProb << " " << OscProbs[iOscProb].NuType << " " << OscProbs[iOscProb].OscChan.GeneratedFlavour << " " << OscProbs[iOscProb].OscChan.DetectedFlavour << " " << OscProbs[iOscProb].Energy << " " << OscProbs[iOscProb].CosineZ << " " << OscProbs[iOscProb].Probability << std::endl;
    }
  }
  
  std::cout << "Finished reweight in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
}
