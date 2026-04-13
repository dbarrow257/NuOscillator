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

int main(int argc, char **argv) {
  std::string FileExt = ".png";
  
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

  //Get the standard set of config names or use what is provided from the arguments
  std::vector<std::string> ConfigNames;
  if (argc == 1) {
    ConfigNames = ReturnKnownConfigs();
  } else {
    for (int i=1;i<argc;i++) {
      ConfigNames.push_back(argv[i]);
    }
  }

  std::cout << "========================================================" << std::endl;
  std::cout << "Setting up Oscillators" << std::endl;
  
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
      if (!Oscillator->ReturnCosineZIgnored()) {
        Oscillator->SetCosineZArrayInCalcer(CosineZArray);
      }
    }

    //Setup
    Oscillator->Setup();

    //Append OscillatorBase* object to the vector
    Oscillators.push_back(Oscillator);
  }
  delete OscFactory;
  
  if (Oscillators.size() == 0) {
    std::cerr << "Expected at least one Oscillator object to compare" << std::endl;
    throw;
  }

  std::cout << "Finished setup in executable" << std::endl;
  std::cout << "========================================================" << std::endl;
  std::cout << "Starting CompareOscillationProbabilities" << std::endl;

  //Check oscillators are all setup equally
  for (size_t iOsc=1;iOsc<Oscillators.size();iOsc++) {

    std::vector<int> NeutrinoTypes = Oscillators[iOsc]->ReturnNeutrinoTypes();
    for (size_t iNT=0;iNT<NeutrinoTypes.size();iNT++) {
      if (NeutrinoTypes[iNT] != (Oscillators[0]->ReturnNeutrinoTypes())[iNT]) {
	std::cerr << "Neutrino types are wrong!" << std::endl;
	throw;
      }
    }
    
    std::vector<NuOscillator::OscillationChannel> OscChannels = Oscillators[iOsc]->ReturnOscChannels();
    for (size_t iOscChan=0;iOscChan<OscChannels.size();iOscChan++) {
      if (OscChannels[iOscChan].GeneratedFlavour != (Oscillators[0]->ReturnOscChannels())[iOscChan].GeneratedFlavour) {
	std::cerr << "Oscillation channels are wrong!" << std::endl;
	throw;
      }
      if (OscChannels[iOscChan].DetectedFlavour != (Oscillators[0]->ReturnOscChannels())[iOscChan].DetectedFlavour) {
	std::cerr << "Oscillation channels are wrong!" << std::endl;
	throw;
      }
    }

    std::vector<FLOAT_T> EnergyArray = Oscillators[iOsc]->ReturnEnergyArray();
    for (size_t iEnergy=0;iEnergy<EnergyArray.size();iEnergy++) {
      if (EnergyArray[iEnergy] != (Oscillators[0]->ReturnEnergyArray())[iEnergy]) {
	std::cerr << "Energy are wrong!" << std::endl;
	throw;
      }
    }
    
    if (Oscillators[iOsc]->ReturnCosineZIgnored() != Oscillators[0]->ReturnCosineZIgnored()) {
      std::cerr << "CosineZIgnore wrong!" << std::endl;
      throw;
    }

    std::vector<FLOAT_T> CosineZArray = Oscillators[iOsc]->ReturnCosineZArray();
    for	(size_t iCosineZ=0;iCosineZ<CosineZArray.size();iCosineZ++) {
      if (CosineZArray[iCosineZ] != (Oscillators[0]->ReturnCosineZArray())[iCosineZ]) {
        std::cerr << "CosineZ are wrong!" << std::endl;
        throw;
      }
    }
  }

  //Build the arrays
  std::vector< std::vector< std::vector< std::vector< std::vector<FLOAT_T> > > > > ProbabilityArray(Oscillators.size());
  for (size_t iOsc=0;iOsc<Oscillators.size();iOsc++) {
    ProbabilityArray[iOsc].resize((Oscillators[iOsc]->ReturnNeutrinoTypes()).size());
      
    for (size_t iNT=0;iNT<(Oscillators[iOsc]->ReturnNeutrinoTypes()).size();iNT++) {
      ProbabilityArray[iOsc][iNT].resize((Oscillators[iOsc]->ReturnOscChannels()).size());
      
      for (int iChan=0;iChan<(Oscillators[iOsc]->ReturnOscChannels()).size();iChan++) {
	ProbabilityArray[iOsc][iNT][iChan].resize(Oscillators[iOsc]->ReturnNEnergyPoints());
	
	if (Oscillators[iOsc]->ReturnCosineZIgnored()) {
	  for (int iEnergy=0;iEnergy<Oscillators[iOsc]->ReturnNEnergyPoints();iEnergy++) {
	    ProbabilityArray[iOsc][iNT][iChan][iEnergy].resize(1);
	  }
	} else {
	  for (int iEnergy=0;iEnergy<Oscillators[iOsc]->ReturnNEnergyPoints();iEnergy++) {
	    ProbabilityArray[iOsc][iNT][iChan][iEnergy].resize(Oscillators[iOsc]->ReturnNCosineZPoints());
	  }
	}
	
      }
    }
  }

  std::cout << "========================================================" << std::endl;
  std::cout << "Calculating Probabilities..." << std::endl;
  
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

    for (size_t iNT=0;iNT<(Oscillators[iOsc]->ReturnNeutrinoTypes()).size();iNT++) {
      int NeutrinoType = (Oscillators[iOsc]->ReturnNeutrinoTypes())[iNT];
      
      for (size_t iOscChan=0;iOscChan<(Oscillators[iOsc]->ReturnOscChannels()).size();iOscChan++) {
	int InitFlav = (Oscillators[iOsc]->ReturnOscChannels())[iOscChan].GeneratedFlavour * NeutrinoType;
	int FinalFlav = (Oscillators[iOsc]->ReturnOscChannels())[iOscChan].DetectedFlavour * NeutrinoType;
	
	for (size_t iEnergy=0;iEnergy<(Oscillators[0]->ReturnEnergyArray()).size();iEnergy++) {
	  FLOAT_T Energy = (Oscillators[0]->ReturnEnergyArray())[iEnergy];
	  
	  if (Oscillators[iOsc]->ReturnCosineZIgnored()) {
	    const FLOAT_T Probability = Oscillators[iOsc]->ReturnOscillationProbability(InitFlav,FinalFlav,Energy);
	    ProbabilityArray[iOsc][iNT][iOscChan][iEnergy][0] = Probability;
	  } else {
	    for (size_t iCosineZ=0;iCosineZ<(Oscillators[iOsc]->ReturnCosineZArray()).size();iCosineZ++) {
	      FLOAT_T CosineZ = (Oscillators[iOsc]->ReturnCosineZArray())[iCosineZ];
	      const FLOAT_T Probability = Oscillators[iOsc]->ReturnOscillationProbability(InitFlav,FinalFlav,Energy,CosineZ);

	      ProbabilityArray[iOsc][iNT][iOscChan][iEnergy][iCosineZ] = Probability;
	    }
	  }
	  
	}
      }
    }
  }

  std::cout << "========================================================" << std::endl;
  std::cout << "Comparing oscillation probabilities.." << std::endl;

  std::vector<FLOAT_T> LargestDiff(Oscillators.size()-1,0.);
  
  if (Oscillators[0]->ReturnCosineZIgnored()) {
    std::cout << std::setw(10) << "Nu Type" << " : " << std::setw(12) << "Gen Flavour" << "->" << std::setw(12) << "Det Flavour" << " | " << std::setw(15) << "Energy" << " || " << std::setw(25) << Oscillators[0]->ReturnImplementationName();
    for (size_t iOsc=1;iOsc<ProbabilityArray.size();iOsc++) {
      std::cout << " | " << std::setw(25) << Oscillators[iOsc]->ReturnImplementationName() << " (" << std::setw(15) << "Differrence" << ")";
    }
    std::cout << std::endl;
  } else {
    std::cout << std::setw(10) << "Nu Type" << " : " << std::setw(12) << "Gen Flavour" << "->" << std::setw(12) << "Det Flavour" << " | " << std::setw(15) << "Energy" << " | " << std::setw(15) << "Cosine Z" << " || " << std::setw(25) << Oscillators[0]->ReturnImplementationName();
    for (size_t iOsc=1;iOsc<ProbabilityArray.size();iOsc++) {
      std::cout << " | " << std::setw(25) << Oscillators[iOsc]->ReturnImplementationName() << " (" << std::setw(15) << "Differrence" << ")";
    }
    std::cout << std::endl;
  }
  
  if (Oscillators[0]->ReturnCosineZIgnored()) {

    for (size_t iEnergy=0;iEnergy<ProbabilityArray[0][0][0].size();iEnergy++) {
      for (size_t iNT=0;iNT<ProbabilityArray[0].size();iNT++) {
	for (size_t iOscChan=0;iOscChan<ProbabilityArray[0][0].size();iOscChan++) {
	  
	  std::cout << std::setw(10) << (Oscillators[0]->ReturnNeutrinoTypes())[iNT] << " : " << std::setw(12) << (Oscillators[0]->ReturnOscChannels())[iOscChan].GeneratedFlavour << "->" << std::setw(12) << (Oscillators[0]->ReturnOscChannels())[iOscChan].DetectedFlavour << " | " << std::setw(15) << (Oscillators[0]->ReturnEnergyArray())[iEnergy] << " || " << std::setw(25) << ProbabilityArray[0][iNT][iOscChan][iEnergy][0];
	  for (size_t iOsc=1;iOsc<ProbabilityArray.size();iOsc++) {
	    std::cout << " | " << std::setw(25) << ProbabilityArray[iOsc][iNT][iOscChan][iEnergy][0] << " (" << std::setw(15) << ProbabilityArray[iOsc][iNT][iOscChan][iEnergy][0]-ProbabilityArray[0][iNT][iOscChan][iEnergy][0] << ")";
	  }
	  std::cout << std::endl;
	}
      }
    }
  } else {
    
    for (size_t iEnergy=0;iEnergy<ProbabilityArray[0][0][0].size();iEnergy++) {
      for (size_t iCosineZ=0;iCosineZ<ProbabilityArray[0][0][0][0].size();iCosineZ++) {
	for (size_t iNT=0;iNT<ProbabilityArray[0].size();iNT++) {
	  for (size_t iOscChan=0;iOscChan<ProbabilityArray[0][0].size();iOscChan++) {
	    
	    std::cout << std::setw(10) << (Oscillators[0]->ReturnNeutrinoTypes())[iNT] << " : " << std::setw(12) << (Oscillators[0]->ReturnOscChannels())[iOscChan].GeneratedFlavour << "->" << std::setw(12) << (Oscillators[0]->ReturnOscChannels())[iOscChan].DetectedFlavour << " | " << std::setw(15) << (Oscillators[0]->ReturnEnergyArray())[iEnergy] << " | " << std::setw(15) << (Oscillators[0]->ReturnCosineZArray())[iCosineZ] << " || " << std::setw(25) << ProbabilityArray[0][iNT][iOscChan][iEnergy][iCosineZ];
	    for (size_t iOsc=1;iOsc<ProbabilityArray.size();iOsc++) {
	      FLOAT_T Diff = ProbabilityArray[iOsc][iNT][iOscChan][iEnergy][iCosineZ]-ProbabilityArray[0][iNT][iOscChan][iEnergy][iCosineZ];
	      if (fabs(Diff) > LargestDiff[iOsc-1]) {LargestDiff[iOsc-1] = fabs(Diff);}
	      
	      std::cout << " | " << std::setw(25) << ProbabilityArray[iOsc][iNT][iOscChan][iEnergy][iCosineZ] << " (" << std::setw(15) << Diff << ")";
	    }
	    std::cout << std::endl;
	  }
	}
	
      }
    }
  }

  std::cout << "Largest differences found: -" << std::endl;
  for (size_t iOsc=1;iOsc<ProbabilityArray.size();iOsc++) {
    std::cout << "\t" << Oscillators[iOsc]->ReturnImplementationName() << " to " << Oscillators[0]->ReturnImplementationName() << " = " << LargestDiff[iOsc-1] << std::endl;
  }
  
  std::cout << "Finished executable" << std::endl;
  std::cout << "========================================================" << std::endl;
}
