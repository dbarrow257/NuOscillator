#include "Oscillator/OscillatorFactory.h"

#include "Constants/OscillatorConstants.h"

#include <iostream>
#include <math.h>
#include <chrono>

#include "TH1D.h"
#include "TCanvas.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

int main() {
  
  std::vector<FLOAT_T> EnergyArray = logspace(0.1,100.,1e5);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,1);

  std::vector<FLOAT_T> OscParams_Basic = ReturnOscParams_Basic();
  std::vector<FLOAT_T> OscParams_Atm = ReturnOscParams_Atm();
  std::vector<FLOAT_T> OscParams_Beam_woYe = ReturnOscParams_Beam_woYe();
  std::vector<FLOAT_T> OscParams_Beam_wYe = ReturnOscParams_Beam_wYe();

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
  std::cout << "Starting drag race in executable" << std::endl;

  int nThrows = 1000;
  std::vector< std::vector<double> > ReweightTimes(Oscillators.size());
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    ReweightTimes[iOsc].resize(nThrows);
  }

  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    std::cout << Oscillators[iOsc]->ReturnImplementationName() << " starting drag race" << std::endl;

    auto t1 = high_resolution_clock::now();

    for (int iThrow=0;iThrow<nThrows;iThrow++) {
      //Throw dcp to some new value
      FLOAT_T RandVal = static_cast <FLOAT_T> (rand()) / static_cast <FLOAT_T> (RAND_MAX);
      OscParams_Atm[5] = RandVal;
      OscParams_Beam_woYe[5] = RandVal;
      OscParams_Beam_wYe[5] = RandVal;

      auto t1_single = high_resolution_clock::now();
      if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Beam_woYe.size()) {
        Oscillators[iOsc]->CalculateProbabilities(OscParams_Beam_woYe);
      } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Beam_wYe.size()) {
        Oscillators[iOsc]->CalculateProbabilities(OscParams_Beam_wYe);
      } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Atm.size()) {
        Oscillators[iOsc]->CalculateProbabilities(OscParams_Atm);
      } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Basic.size()) {
        Oscillators[iOsc]->CalculateProbabilities(OscParams_Basic);
      } else {
        std::cerr << "Did not find viable oscillation parameters to hand to the oscillation probability calculater" << std::endl;
        std::cerr << "Oscillator->ReturnNOscParams():" << Oscillators[iOsc]->ReturnNOscParams() << std::endl;
        throw;
      }
      auto t2_single = high_resolution_clock::now();
      duration<double, std::milli> ms_single_double = t2_single-t1_single;
      ReweightTimes[iOsc][iThrow] = ms_single_double.count();
    }
    
    auto t2 = high_resolution_clock::now();
    duration<double, std::milli> ms_double = t2 - t1;

    std::cout << Oscillators[iOsc]->ReturnImplementationName() << " ended drag race - took " << ms_double.count()/nThrows << " milliseconds per reweight (Using nEnergyPoints = " << Oscillators[iOsc]->ReturnNEnergyPoints() << ")" << std::endl;
  }

  TCanvas* Canv = new TCanvas;

  std::vector<TH1D*> TimeDistributions(Oscillators.size());
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {

    double Min = 1e8;
    double Max = -1e8;
    for (int iThrow=0;iThrow<nThrows;iThrow++) {
      if (ReweightTimes[iOsc][iThrow] > Max) Max = ReweightTimes[iOsc][iThrow];
      if (ReweightTimes[iOsc][iThrow] < Min) Min = ReweightTimes[iOsc][iThrow];
    }
    
    TimeDistributions[iOsc] = new TH1D((Oscillators[iOsc]->ReturnImplementationName()).c_str(),";Reweight Time [ms];",30,Min-0.1*(Max-Min),Max+0.1*(Max-Min));
    for (int iThrow=0;iThrow<nThrows;iThrow++) {
      TimeDistributions[iOsc]->Fill(ReweightTimes[iOsc][iThrow]);
    }

    TimeDistributions[iOsc]->Draw();
    Canv->Print((Oscillators[iOsc]->ReturnImplementationName()+"_TimeDistribution.pdf").c_str());
  }

  std::cout << "Finished drag race in executable" << std::endl;
  std::cout << "========================================================" << std::endl;

  for (size_t iOsc = 0; iOsc < Oscillators.size(); iOsc++) {
    delete Oscillators[iOsc];
  }
  Oscillators.clear();
  delete OscFactory;
}
