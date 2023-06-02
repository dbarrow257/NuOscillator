#include <iostream>
#include <math.h>

#include "OscProbCalcerBase.h"

#ifdef UseCUDAProb3
#include "OscProbCalcer_CUDAProb3.h"
#endif

#ifdef UseProb3ppLinear
#include "OscProbCalcer_Prob3ppLinear.h"
#endif

#ifdef UseProbGPULinear
#include "OscProbCalcer_ProbGPULinear.h"
#endif

std::vector<FLOAT_T> logspace(FLOAT_T Emin, FLOAT_T  Emax, int nDiv);
std::vector<FLOAT_T> linspace(FLOAT_T Emin, FLOAT_T Emax, int nDiv);

int main() {
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

  std::cout << "========================================================" << std::endl;
  std::cout << "Starting setup in executable" << std::endl;

  std::vector<OscProbCalcerBase*> OscProbCalcers;

#ifdef UseCUDAProb3
  OscProbCalcerCUDAProb3* CUDAProb3 = new OscProbCalcerCUDAProb3();
  OscProbCalcers.push_back((OscProbCalcerBase*)CUDAProb3);
#endif

#ifdef UseProb3ppLinear
  OscProbCalcerProb3ppLinear* Prob3ppLinear = new OscProbCalcerProb3ppLinear();
  OscProbCalcers.push_back((OscProbCalcerBase*)Prob3ppLinear);
#endif

#ifdef UseProbGPULinear
  OscProbCalcerProbGPULinear* ProbGPULinear = new OscProbCalcerProbGPULinear();
  OscProbCalcers.push_back((OscProbCalcerBase*)ProbGPULinear);
#endif

  /*
  std::vector<FLOAT_T> EnergyArray = logspace(0.01,100,1000);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,10);
  */

  std::vector<FLOAT_T> EnergyArray;
  EnergyArray.push_back(10);
  std::vector<FLOAT_T> CosineZArray;
  CosineZArray.push_back(0.);

  // Setup propagators
  for (size_t iCalcer=0;iCalcer<OscProbCalcers.size();iCalcer++) {
    OscProbCalcers[iCalcer]->SetEnergyArray(EnergyArray);
    OscProbCalcers[iCalcer]->SetCosineZArray(CosineZArray); 
    OscProbCalcers[iCalcer]->Setup();
  }
  std::cout << "Finished setup in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting reweight in executable" << std::endl;

  // Reweight and calculate oscillation probabilities
  for (size_t iCalcer=0;iCalcer<OscProbCalcers.size();iCalcer++) {

    // These don't have to be explicilty beam or atmospheric specific, all they have to be is equal to the number of oscillation parameters expected by the implementation
    // If you have some NSI calculater, then it will work providing tthe length of the vector of oscillation parameters is equal to the number of expected oscillation parameters
    if (OscProbCalcers[iCalcer]->ReturnExpectedNOscParams() == (int)OscParams_Beam.size()) {
      OscProbCalcers[iCalcer]->Reweight(OscParams_Beam); 
    } else if (OscProbCalcers[iCalcer]->ReturnExpectedNOscParams() == (int)OscParams_Atm.size()) {
      OscProbCalcers[iCalcer]->Reweight(OscParams_Atm);
    } else {
      std::cerr << "Did not find viable oscillation parameters to hand to the oscillation probability calculater" << std::endl;
      std::cerr << "OscProbCalcers[iCalcer]->ReturnExpectedNOscParams():" << OscProbCalcers[iCalcer]->ReturnExpectedNOscParams() << std::endl;
      throw;
    }

    OscProbCalcers[iCalcer]->PrintWeights();
  }

  std::cout << "Finished reweight in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
}

std::vector<FLOAT_T> logspace(FLOAT_T Emin, FLOAT_T  Emax, int nDiv) {
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

std::vector<FLOAT_T> linspace(FLOAT_T Emin, FLOAT_T Emax, int nDiv) {
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
