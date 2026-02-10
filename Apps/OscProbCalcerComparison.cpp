#include "Oscillator/OscillatorFactory.h"

#include "Constants/OscillatorConstants.h"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <chrono>

#include "TGraph.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::duration;
using std::chrono::milliseconds;

int main() {
  std::string FileExt = ".png";
  
  //============================================================================================================
  //This executable assumes comparing oscillation engines which have 9 oscillation channels (e,mu,tau)

  int nChannels = 9;
  int nInitFlav = 3;
  int nFinalFlav = 3;
  
  //============================================================================================================
  
  std::vector<FLOAT_T> EnergyArray = logspace(0.1,10.,1e3);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,1);

  std::vector<FLOAT_T> OscParams_Basic = ReturnOscParams_Basic();
  std::vector<FLOAT_T> OscParams_Atm = ReturnOscParams_Atm();
  std::vector<FLOAT_T> OscParams_Beam_woYe = ReturnOscParams_Beam_woYe();
  std::vector<FLOAT_T> OscParams_Beam_wYe = ReturnOscParams_Beam_wYe();
  std::vector<FLOAT_T> OscParams_Beam_wYe_wDeco = ReturnOscParams_Beam_wYe_wDeco();
  std::vector<FLOAT_T> OscParams_Beam_wYe_wLIV = ReturnOscParams_Beam_wYe_wLIV();
  std::vector<FLOAT_T> OscParams_Beam_wYe_wNSI = ReturnOscParams_Beam_wYe_wNSI();
  
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
  std::cout << "Starting OscProbCalcerComparison" << std::endl;

  std::vector< std::vector< std::vector<FLOAT_T> > > ProbabilityArray(Oscillators.size());
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    ProbabilityArray[iOsc].resize(nChannels);
    for (int iChan=0;iChan<nChannels;iChan++) {
      ProbabilityArray[iOsc][iChan].resize(Oscillators[iOsc]->ReturnNEnergyPoints());
    }
  }
  
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    std::cout << "Creating probability from:" << Oscillators[iOsc]->ReturnImplementationName() << std::endl;
    
    if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Basic.size()) {
      Oscillators[iOsc]->CalculateProbabilities(OscParams_Basic);
    } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Atm.size()) {
      Oscillators[iOsc]->CalculateProbabilities(OscParams_Atm);
    } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Beam_woYe.size()) {
      Oscillators[iOsc]->CalculateProbabilities(OscParams_Beam_woYe);
    } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Beam_wYe.size()) {
      Oscillators[iOsc]->CalculateProbabilities(OscParams_Beam_wYe);
    } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Beam_wYe_wDeco.size()) {
      Oscillators[iOsc]->CalculateProbabilities(OscParams_Beam_wYe_wDeco);
    } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Beam_wYe_wLIV.size()) {
      Oscillators[iOsc]->CalculateProbabilities(OscParams_Beam_wYe_wLIV);
    } else if (Oscillators[iOsc]->ReturnNOscParams() == (int)OscParams_Beam_wYe_wNSI.size()) {
      Oscillators[iOsc]->CalculateProbabilities(OscParams_Beam_wYe_wNSI);
    } else {
      std::cerr << "Did not find viable oscillation parameters to hand to the oscillation probability calculater" << std::endl;
      std::cerr << "Oscillator->ReturnNOscParams():" << Oscillators[iOsc]->ReturnNOscParams() << std::endl;
      throw std::runtime_error("Invalid setup");
    }

    for (int iInitFlav=1;iInitFlav<=nInitFlav;iInitFlav++) {
      for (int iFinalFlav=1;iFinalFlav<=nFinalFlav;iFinalFlav++) {
	for (int iEnergy=0;iEnergy<Oscillators[iOsc]->ReturnNEnergyPoints();iEnergy++) {
	  const FLOAT_T Probability = Oscillators[iOsc]->ReturnOscillationProbability(iInitFlav,iFinalFlav,EnergyArray[iEnergy]);
	  ProbabilityArray[iOsc][(iInitFlav-1)*nInitFlav+(iFinalFlav-1)][iEnergy] = Probability;
	}
      }
    }
  }

  std::cout << "Created probability array.." << std::endl;

  std::vector< std::vector<TGraph*> > Probabilities(Oscillators.size());
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    Probabilities[iOsc].resize(nChannels);
    for (int iChan=0;iChan<nChannels;iChan++) {
      Probabilities[iOsc][iChan] = new TGraph(Oscillators[iOsc]->ReturnNEnergyPoints(),EnergyArray.data(),ProbabilityArray[iOsc][iChan].data());
      Probabilities[iOsc][iChan]->SetTitle(";Neutrino Energy [GeV];Probability");
    }
  }
  
  std::vector< std::vector< std::vector<FLOAT_T> > > ProbabilityArray_Ratio(Oscillators.size());
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    ProbabilityArray_Ratio[iOsc].resize(nChannels);
    for	(int iChan=0;iChan<nChannels;iChan++) {
      ProbabilityArray_Ratio[iOsc][iChan].resize(Oscillators[iOsc]->ReturnNEnergyPoints());
      for (int iEnergy=0;iEnergy<Oscillators[iOsc]->ReturnNEnergyPoints();iEnergy++) {
	ProbabilityArray_Ratio[iOsc][iChan][iEnergy] = 2.*(ProbabilityArray[iOsc][iChan][iEnergy]-ProbabilityArray[0][iChan][iEnergy])/(ProbabilityArray[iOsc][iChan][iEnergy]+ProbabilityArray[0][iChan][iEnergy]);
      }
    }
  }
  
  std::vector< std::vector<TGraph*> > Probabilities_Ratio(Oscillators.size());
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    Probabilities_Ratio[iOsc].resize(nChannels);
    for (int iChan=0;iChan<nChannels;iChan++) {
      Probabilities_Ratio[iOsc][iChan] = new TGraph(Oscillators[iOsc]->ReturnNEnergyPoints(),EnergyArray.data(),ProbabilityArray_Ratio[iOsc][iChan].data());
      Probabilities_Ratio[iOsc][iChan]->SetTitle(";Neutrino Energy [GeV];2*(Engine[n]-Engine[0])/(Engine[n]+Engine[0]) of Probabilities");
    }
  }

  for (int iChan=0;iChan<nChannels;iChan++) {
    FLOAT_T Max = -1e8;
    FLOAT_T Min = 1e8;
    
    for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
      for (int iEnergy=0;iEnergy<Oscillators[iOsc]->ReturnNEnergyPoints();iEnergy++) {
	if (ProbabilityArray_Ratio[iOsc][iChan][iEnergy] > Max) Max = ProbabilityArray_Ratio[iOsc][iChan][iEnergy];
	if (ProbabilityArray_Ratio[iOsc][iChan][iEnergy] < Min) Min = ProbabilityArray_Ratio[iOsc][iChan][iEnergy];
      }
    }

    for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
      Probabilities_Ratio[iOsc][iChan]->SetMaximum(Max+0.1*(Max-Min));
      Probabilities_Ratio[iOsc][iChan]->SetMinimum(Min-0.1*(Max-Min));
    }
  }

  for (int iEnergy=0;iEnergy<Oscillators[0]->ReturnNEnergyPoints();iEnergy++) {
    std::cout << "Energy:" << std::setw(10) << EnergyArray[iEnergy] << " ";
    for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
      std::cout << std::setw(10) << Probabilities[iOsc][0]->GetPointY(iEnergy) << " (" << std::setw(10) << Probabilities_Ratio[iOsc][0]->GetPointY(iEnergy) << ") ";
    }
    std::cout << std::endl;
  }

  TLegend* Leg = new TLegend(0.2,0.85,0.9,0.98);
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    for (int iChan=0;iChan<nChannels;iChan++) {
      Probabilities[iOsc][iChan]->SetLineColor(40+iOsc);
      Probabilities_Ratio[iOsc][iChan]->SetLineColor(40+iOsc);
    }
    Leg->AddEntry(Probabilities[iOsc][0],(Oscillators[iOsc]->ReturnImplementationName()).c_str(),"l");
  }
  
  TCanvas* Canv1 = new TCanvas("Canv1","");
  TCanvas* Canv2 = new TCanvas("Canv2","");
  
  for (int iInitFlav=1;iInitFlav<=nInitFlav;iInitFlav++) {
    for (int iFinalFlav=1;iFinalFlav<=nFinalFlav;iFinalFlav++) {

      int iChan = (iInitFlav-1)*nInitFlav+(iFinalFlav-1);
      
      TLatex Text(0.4,0.75,(NeutrinoFlavour_IntToStr(iInitFlav)+"#rightarrow"+NeutrinoFlavour_IntToStr(iFinalFlav)).c_str());
      Text.SetTextSize(0.045);
      Text.SetNDC(true);
   
      Canv1->cd();
      Canv1->SetTopMargin(0.15);
      Canv1->SetLeftMargin(0.15);
      Canv1->SetLogx(true);
      
      for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
	if (iOsc==0) {
	  Probabilities[iOsc][iChan]->Draw("AL");
	} else {
	  Probabilities[iOsc][iChan]->Draw("L SAME");
	}
      }
      Leg->Draw("SAME");
      Text.Draw("SAME");
      Canv1->Print((Form("LinearProbabilityComparison_%i",iChan)+FileExt).c_str());
      
      Canv2->cd();
      Canv2->SetTopMargin(0.15);
      Canv2->SetLeftMargin(0.15);
      Canv2->SetLogx(true);
      
      for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
	if (iOsc==0) {
	  Probabilities_Ratio[iOsc][iChan]->Draw("AL");
	} else {
	  Probabilities_Ratio[iOsc][iChan]->Draw("L SAME");
	}
      }
      Leg->Draw("SAME");
      Text.Draw("SAME");
      Canv2->Print((Form("LinearProbabilityComparison_Ratio_%i",iChan)+FileExt).c_str());
    }
  }
  
  std::cout << "Finished executable" << std::endl;
  std::cout << "========================================================" << std::endl;
}
