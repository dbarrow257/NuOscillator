#include "Oscillator/OscillatorFactory.h"

#include "Constants/OscillatorConstants.h"

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
    std::cerr << "./SingleOscillator InputConfig.yaml" << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  std::string ConfigName = argv[1];
  
  bool PrintWeights = true;

  std::vector<FLOAT_T> EnergyArray = logspace(0.1,100.,1e3);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,1e3);

  std::vector<FLOAT_T> OscParams_Basic = ReturnOscParams_Basic();
  std::vector<FLOAT_T> OscParams_Atm = ReturnOscParams_Atm();
  std::vector<FLOAT_T> OscParams_Beam_woYe = ReturnOscParams_Beam_woYe();
  std::vector<FLOAT_T> OscParams_Beam_wYe = ReturnOscParams_Beam_wYe();

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

  Oscillator->Setup();
  
  std::cout << "Finished setup in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting reweight in executable" << std::endl;

  // These don't have to be explicilty beam or atmospheric specific, all they have to be is equal to the number of oscillation parameters expected by the implementation
  // If you have some NSO calculater, then it will work providing the length of the vector of oscillation parameters is equal to the number of expected oscillation parameters
  if (Oscillator->ReturnNOscParams() == (int)OscParams_Beam_woYe.size()) {
    Oscillator->CalculateProbabilities(OscParams_Beam_woYe);
  } else if (Oscillator->ReturnNOscParams() == (int)OscParams_Beam_wYe.size()) {
    Oscillator->CalculateProbabilities(OscParams_Beam_wYe); 
  } else if (Oscillator->ReturnNOscParams() == (int)OscParams_Atm.size()) {
    Oscillator->CalculateProbabilities(OscParams_Atm);
  } else if (Oscillator->ReturnNOscParams() == (int)OscParams_Basic.size()) {
    Oscillator->CalculateProbabilities(OscParams_Basic);
  } else {
    std::cerr << "Did not find viable oscillation parameters to hand to the oscillation probability calculater" << std::endl;
    std::cerr << "Oscillator->ReturnNOscParams():" << Oscillator->ReturnNOscParams() << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  
  if (PrintWeights) {
    Oscillator->PrintWeights();
  }

  std::cout << "Finished reweight in executable" << std::endl;
  std::cout << "========================================================" << std::endl;

  TCanvas* Canv = new TCanvas;
  Canv->SetLogx(true);

  if (Oscillator->CosineZIgnored()) {
    std::vector<FLOAT_T> EnergyBinning = Oscillator->ReturnBinEdgesForPlotting(true);
    
    TH1D* Hist = new TH1D("Probability","",EnergyBinning.size()-1,EnergyBinning.data());
    for (int xBin=1;xBin<=Hist->GetNbinsX();xBin++) {
      FLOAT_T Prob;
      
      if (Oscillator->EvalPointsSetInConstructor()) {
	Prob = Oscillator->ReturnOscillationProbability(1,1,Hist->GetXaxis()->GetBinCenter(xBin));
      } else {
	Prob = Oscillator->ReturnOscillationProbability(1,1,EnergyArray[xBin-1]);
      }
      
      Hist->SetBinContent(xBin,Prob);
    }
    
    Hist->Draw();
    Canv->Print("Probability.pdf");
  } else {
    std::vector<FLOAT_T> EnergyBinning = Oscillator->ReturnBinEdgesForPlotting(true);
    std::vector<FLOAT_T> CosineZBinning = Oscillator->ReturnBinEdgesForPlotting(false);
    
    TH2D* Hist = new TH2D("Oscillogram","",EnergyBinning.size()-1,EnergyBinning.data(),CosineZBinning.size()-1,CosineZBinning.data());
    for (int xBin=1;xBin<=Hist->GetNbinsX();xBin++) {
      for (int yBin=1;yBin<=Hist->GetNbinsY();yBin++) {
	FLOAT_T Prob;
	
	if (Oscillator->EvalPointsSetInConstructor()) {
	  Prob = Oscillator->ReturnOscillationProbability(1,1,Hist->GetXaxis()->GetBinCenter(xBin),Hist->GetYaxis()->GetBinCenter(yBin));
	} else {
	  Prob = Oscillator->ReturnOscillationProbability(1,1,EnergyArray[xBin-1],CosineZArray[yBin-1]);
	}
	
	Hist->SetBinContent(xBin,yBin,Prob);
      }
    }
    
    Hist->Draw("COLZ");
    Canv->Print("Oscillogram.pdf");
  }

}
