#include "Constants/OscillatorConstants.h"
#include "OscProbCalcer/OscProbCalcerBase.h"

#if UseCUDAProb3==1
#include "OscProbCalcer/OscProbCalcer_CUDAProb3.h"
#endif

#if UseCUDAProb3Linear==1
#include "OscProbCalcer/OscProbCalcer_CUDAProb3Linear.h"
#endif

#if UseProb3ppLinear==1
#include "OscProbCalcer/OscProbCalcer_Prob3ppLinear.h"
#endif

#if UseProbGPULinear==1
#include "OscProbCalcer/OscProbCalcer_ProbGPULinear.h"
#endif

#if UseNuFASTLinear==1
#include "OscProbCalcer/OscProbCalcer_NuFASTLinear.h"
#endif

#if UseOscProb==1
#include "OscProbCalcer/OscProbCalcer_OscProb.h"
#endif

#include <iostream>
#include <math.h>
#include <chrono>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "./SingleOscProbCalcer OscProbCalcerImplementationToCreate InputConfig.yaml" << std::endl;
    throw;
  }
  std::string OscProbCalcerImplementationToCreate = argv[1];
  std::string OscProbCalcerConfigname = argv[2];
  
  bool PrintWeights = true;

  std::vector<FLOAT_T> OscParams_Atm(7);
  OscParams_Atm[0] = 3.07e-1;
  OscParams_Atm[1] = 5.28e-1;
  OscParams_Atm[2] = 2.18e-2;
  OscParams_Atm[3] = 7.53e-5;
  OscParams_Atm[4] = 2.509e-3;
  OscParams_Atm[5] = -1.601;
  OscParams_Atm[6] = 25.0;

  std::vector<FLOAT_T> OscParams_Atm_OscProb(7);
  OscParams_Atm_OscProb[0] = 3.03e-1;
  OscParams_Atm_OscProb[1] = 4.51e-1;
  OscParams_Atm_OscProb[2] = 2.2225e-2;
  OscParams_Atm_OscProb[3] = 7.41e-5;
  OscParams_Atm_OscProb[4] = 2.511e-3 - 7.41e-5;
  OscParams_Atm_OscProb[5] = 232. * 3.14159265 / 180.;
  OscParams_Atm_OscProb[6] = 2;
  //OscParams_Atm_OscProb[7] = 1e-3;

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
  
  std::vector<FLOAT_T> EnergyArray = logspace(0.1,100.,2);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,1);

  std::cout << "========================================================" << std::endl;
  std::cout << "Starting setup in executable" << std::endl;

  std::cout << "========================================================" << std::endl;
  std::cout << "Initialising " << OscProbCalcerConfigname << std::endl;
  
  OscProbCalcerBase* Calcer;

  if (OscProbCalcerImplementationToCreate == "CUDAProb3") {
#if UseCUDAProb3==1
    OscProbCalcerCUDAProb3* CUDAProb3 = new OscProbCalcerCUDAProb3(OscProbCalcerConfigname);
    Calcer = (OscProbCalcerBase*)CUDAProb3;
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }

  else if (OscProbCalcerImplementationToCreate == "CUDAProb3Linear") {
#if UseCUDAProb3Linear==1
    OscProbCalcerCUDAProb3Linear* CUDAProb3Linear = new OscProbCalcerCUDAProb3Linear(OscProbCalcerConfigname);
    Calcer = (OscProbCalcerBase*)CUDAProb3Linear;
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else if (OscProbCalcerImplementationToCreate == "Prob3ppLinear") {
#if UseProb3ppLinear==1
    OscProbCalcerProb3ppLinear* Prob3ppLinear = new OscProbCalcerProb3ppLinear(OscProbCalcerConfigname);
    Calcer = (OscProbCalcerBase*)Prob3ppLinear;
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else if (OscProbCalcerImplementationToCreate == "ProbGPULinear") {
#if UseProbGPULinear==1
    OscProbCalcerProbGPULinear* ProbGPULinear = new OscProbCalcerProbGPULinear(OscProbCalcerConfigname);
    Calcer = (OscProbCalcerBase*)ProbGPULinear;
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }

  else if (OscProbCalcerImplementationToCreate == "NuFASTLinear") {
#if UseNuFASTLinear==1
    OscProbCalcerNuFASTLinear* NuFASTLinear = new OscProbCalcerNuFASTLinear(OscProbCalcerConfigname);
    Calcer = (OscProbCalcerBase*)NuFASTLinear;
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }

  else if (OscProbCalcerImplementationToCreate == "OscProb") {
#if UseOscProb==1
    OscProbCalcerOscProb* OscProb = new OscProbCalcerOscProb(OscProbCalcerConfigname);
    Calcer = (OscProbCalcerBase*)OscProb;
#else
    std::cerr << "Oscillator was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but Use" << OscProbCalcerImplementationToCreate << " is undefined. Indicates problem in setup" << std::endl;
    throw;
#endif
  }
  
  else {
    std::cerr << "Executable was requsted to create " << OscProbCalcerImplementationToCreate << " OscProbCalcer but this is not implemented within " << __FILE__ << std::endl;
    std::cerr << "Please correct the mistake or implement the new OscProbCalcer" << std::endl;
    throw;
  }

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
  if(OscProbCalcerImplementationToCreate == "OscProb") {
        Calcer->Reweight(OscParams_Atm_OscProb);
  }
  else if (Calcer->ReturnNOscParams() == (int)OscParams_Beam_woYe.size()) {
    Calcer->Reweight(OscParams_Beam_woYe);
  } else if (Calcer->ReturnNOscParams() == (int)OscParams_Beam_wYe.size()) {
    Calcer->Reweight(OscParams_Beam_wYe); 
  } else if (Calcer->ReturnNOscParams() == (int)OscParams_Atm.size()) {
    Calcer->Reweight(OscParams_Atm);
  }
   else {
    std::cerr << "Did not find viable oscillation parameters to hand to the oscillation probability calculater" << std::endl;
    std::cerr << "Oscillator->ReturnNOscParams():" << Calcer->ReturnNOscParams() << std::endl;
    throw;
  }

 

  if (PrintWeights) {
    std::vector<OscillationProbability> OscProbs = Calcer->ReturnProbabilities();
    for (int iOscProb=0;iOscProb<(int)OscProbs.size();iOscProb++) {
      std::cout << iOscProb << " " << OscProbs[iOscProb].NuType << " " << OscProbs[iOscProb].OscChan.GeneratedFlavour << " " << OscProbs[iOscProb].OscChan.DetectedFlavour << " " << OscProbs[iOscProb].Energy << " " << OscProbs[iOscProb].CosineZ << " " << OscProbs[iOscProb].Probability << std::endl;
    }
  }
  
  std::cout << "Finished reweight in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
}
