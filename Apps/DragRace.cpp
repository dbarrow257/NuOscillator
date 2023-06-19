#include "OscillatorBinned.h"
#include "OscillatorUnbinned.h"

#include "OscillatorConstants.h"

#include <iostream>
#include <math.h>
#include <chrono>

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

std::vector<FLOAT_T> logspace(FLOAT_T Emin, FLOAT_T  Emax, int nDiv);
std::vector<FLOAT_T> linspace(FLOAT_T Emin, FLOAT_T Emax, int nDiv);

int main() {
  enum Verbosity{NONE,INFO};
  int Verbose = NONE;
  
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

  std::vector<FLOAT_T> EnergyArray = logspace(0.1,100.,1e6);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,1);

  std::cout << "========================================================" << std::endl;
  std::cout << "Starting setup in executable" << std::endl;

  std::vector<OscillatorBase*> Oscillators;
  std::string ConfigName;

#if UseCUDAProb3 == 1

#if UseBinned == 1
  ConfigName = "";
  OscillatorBinned* Oscillator_CUDAProb3Binned = new OscillatorBinned(ConfigName);
  Oscillators.push_back((OscillatorBase*)Oscillator_CUDAProb3Binned);
#else
  ConfigName = "";
  OscillatorUnbinned* Oscillator_CUDAProb3Linear = new OscillatorUnbinned(CUDAProb3Linear_Vector,CUDAProb3Linear_Verbosity,CUDAProb3Linear_IgnoreCosineZ);
  Oscillator_CUDAProb3Linear->SetEnergyArray(EnergyArray);
  Oscillator_CUDAProb3Linear->SetCosineZArray(CosineZArray);
  Oscillators.push_back((OscillatorBase*)Oscillator_CUDAProb3Linear);
#endif

#endif

#if UseCUDAProb3Linear == 1

#if UseBinned == 1
  ConfigName = "";
  OscillatorBinned* Oscillator_CUDAProb3LinearBinned = new OscillatorBinned(ConfigName);
  Oscillators.push_back((OscillatorBase*)Oscillator_CUDAProb3LinearBinned);
#else
  ConfigName = "";
  OscillatorUnbinned* Oscillator_CUDAProb3LinearLinear = new OscillatorUnbinned(ConfigName);
  Oscillator_CUDAProb3LinearLinear->SetEnergyArray(EnergyArray);
  Oscillators.push_back((OscillatorBase*)Oscillator_CUDAProb3LinearLinear);
#endif

#endif

#if UseProbGPULinear == 1

#if UseBinned == 1
  ConfigName = "";
  OscillatorBinned* Oscillator_ProbGPUBinned = new OscillatorBinned(ConfigName);
  Oscillators.push_back((OscillatorBase*)Oscillator_ProbGPUBinned);
#else
  ConfigName = "";
  OscillatorUnbinned* Oscillator_ProbGPULinear = new OscillatorUnbinned(ConfigName);
  Oscillator_ProbGPULinear->SetEnergyArray(EnergyArray);
  Oscillators.push_back((OscillatorBase*)Oscillator_ProbGPULinear);
#endif

#endif 

#if UseProb3ppLinear == 1

#if UseBinned == 1
  ConfigName = "";
  OscillatorBinned* Oscillator_Prob3ppBinned = new OscillatorBinned(ConfigName);
  Oscillators.push_back((OscillatorBase*)Oscillator_Prob3ppBinned);
#else
  ConfigName = "";
  OscillatorUnbinned* Oscillator_Prob3ppLinear = new OscillatorUnbinned(ConfigName);
  Oscillator_Prob3ppLinear->SetEnergyArray(EnergyArray);
  Oscillators.push_back((OscillatorBase*)Oscillator_Prob3ppLinear);
#endif

#endif

  // Setup propagators
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    Oscillators[iOsc]->Setup();
  }

  std::cout << "Finished setup in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting drag race in executable" << std::endl;

  int nThrows = 1000;

  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    std::cout << Oscillators[iOsc]->ReturnImplementationName() << " starting drag race" << std::endl;

    auto t1 = high_resolution_clock::now();
    for (int iThrow=0;iThrow<nThrows;iThrow++) {      
      //Throw dcp to some new value
      FLOAT_T RandVal = rand();
      OscParams_Atm[5] = RandVal;
      OscParams_Beam[5] = RandVal;
     
      if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Beam.size()) {
	Oscillators[iOsc]->CalculateProbabilities(OscParams_Beam);
      } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Atm.size()) {
	Oscillators[iOsc]->CalculateProbabilities(OscParams_Atm);
      } else {
	std::cerr << "Did not find viable oscillation parameters to hand to the oscillation probability calculater" << std::endl;
	std::cerr << "Oscillator->ReturnNOscParams():" << Oscillators[iOsc]->ReturnNOscParams() << std::endl;
	throw;
      }
    }
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << Oscillators[iOsc]->ReturnImplementationName() << " ended drag race - took " << ms_double.count()/nThrows << " milliseconds per reweight (Using nEnergyPoints = " << Oscillators[iOsc]->ReturnNEnergyPoints() << ")" << std::endl;
  }
  
  std::cout << "Finished drag race in executable" << std::endl;
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
