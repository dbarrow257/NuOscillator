#include "Oscillator/OscillatorSubSampling.h"

#include <iostream>

#include "TFile.h"
#include "TH1.h"

OscillatorSubSampling::OscillatorSubSampling(std::string ConfigName_) : OscillatorBase(ConfigName_) {
  Initialise();
}

OscillatorSubSampling::OscillatorSubSampling(YAML::Node Config_) : OscillatorBase(Config_) {
  Initialise();
}

void OscillatorSubSampling::Initialise() {
  CoarseEnergyAxisBinEdges = std::vector<FLOAT_T>();
  CoarseCosineZAxisBinEdges = std::vector<FLOAT_T>();
  FineEnergyAxisBinEdges = std::vector<FLOAT_T>();
  FineCosineZAxisBinEdges = std::vector<FLOAT_T>();

  FineEnergyAxisBinCenters = std::vector<FLOAT_T>();
  FineCosineZAxisBinCenters = std::vector<FLOAT_T>();

  fCalculationTypeName = "SubSampling";

  //=======
  // Grab the following from config manager

  FileName = Config[fCalculationTypeName]["FileName"].as<std::string>();
  CoarseEnergyAxisHistName = Config[fCalculationTypeName]["CoarseEnergyAxisHistName"].as<std::string>();
  FineEnergyAxisHistName = Config[fCalculationTypeName]["FineEnergyAxisHistName"].as<std::string>();

  if (!fCosineZIgnored) {
    CoarseCosineZAxisHistName = Config[fCalculationTypeName]["CoarseCosineZAxisHistName"].as<std::string>();
    FineCosineZAxisHistName = Config[fCalculationTypeName]["FineCosineZAxisHistName"].as<std::string>();
  } else {
    CoarseCosineZAxisHistName = "Dummy";
    FineCosineZAxisHistName = "Dummy";
  }

  //=======
  // Grab the bin edges and centers for both coarse and fine binning
  
  CoarseEnergyAxisBinEdges = ReadBinEdgesFromFile(FileName,CoarseEnergyAxisHistName);
  FineEnergyAxisBinEdges = ReadBinEdgesFromFile(FileName,FineEnergyAxisHistName);
  FineEnergyAxisBinCenters = ReturnBinCentersFromBinEdges(FineEnergyAxisBinEdges);
  if (!fCosineZIgnored) {
    CoarseCosineZAxisBinEdges = ReadBinEdgesFromFile(FileName,CoarseCosineZAxisHistName);
    FineCosineZAxisBinEdges = ReadBinEdgesFromFile(FileName,FineCosineZAxisHistName);
    FineCosineZAxisBinCenters = ReturnBinCentersFromBinEdges(FineCosineZAxisBinEdges);
  } else {
    CoarseCosineZAxisBinEdges = std::vector<FLOAT_T>();
    CoarseCosineZAxisBinEdges.push_back(-std::numeric_limits<float>::max());
    CoarseCosineZAxisBinEdges.push_back(std::numeric_limits<float>::max());

    FineCosineZAxisBinEdges = std::vector<FLOAT_T>();
    FineCosineZAxisBinEdges.push_back(-std::numeric_limits<float>::max());
    FineCosineZAxisBinEdges.push_back(std::numeric_limits<float>::max());

    FineCosineZAxisBinCenters = std::vector<FLOAT_T>();
    FineCosineZAxisBinCenters.push_back(0.);
  }

  //=============
  //The energies and cosinezs which are used for probability calculation are the fine energies, so set them now

  SetEnergyArrayInCalcer(FineEnergyAxisBinCenters);
  if (!fCosineZIgnored) {
    SetCosineZArrayInCalcer(FineCosineZAxisBinCenters);
  }
  fEvalPointsSetInConstructor = true;
}

OscillatorSubSampling::~OscillatorSubSampling() {
}

void OscillatorSubSampling::SetupOscillatorImplementation() {
  NeutrinoTypes = fOscProbCalcer->ReturnNeutrinoTypes();
  nNeutrinoTypes = NeutrinoTypes.size();

  OscillationChannels = fOscProbCalcer->ReturnOscChannels();
  nOscillationChannels = OscillationChannels.size();
  
  nCoarseEnergyBins = CoarseEnergyAxisBinEdges.size()-1.;
  nCoarseCosineZBins = CoarseCosineZAxisBinEdges.size()-1.; //This will be equal to 1 for non-atmospheric OscProbCalcers because of the dummy values defined in the Intialisation function

  nFineEnergyBins = FineEnergyAxisBinCenters.size();
  nFineCosineZBins = FineCosineZAxisBinCenters.size(); //This will be equal to 1 for non-atmospheric OscProbCalcers because of the dummy values defined in the Intialisation function
  
  TotalCoarseBins = nCoarseEnergyBins*nCoarseCosineZBins*nOscillationChannels*nNeutrinoTypes;
  AveragedOscillationProbabilities.resize(TotalCoarseBins);
  OscillationProbabilitiesToAverage.resize(TotalCoarseBins);

  for (size_t iFineCosineZBin=0;iFineCosineZBin<nFineCosineZBins;iFineCosineZBin++) {
    FLOAT_T FineCosineZBinCenter = FineCosineZAxisBinCenters[iFineCosineZBin];
    int CoarseCosineZBin = FindBinIndexFromEdges(FineCosineZBinCenter,CoarseCosineZAxisBinEdges);
    
    for (size_t iFineEnergyBin=0;iFineEnergyBin<nFineEnergyBins;iFineEnergyBin++) {
      FLOAT_T FineEnergyBinCenter = FineEnergyAxisBinCenters[iFineEnergyBin];
      int CoarseEnergyBin = FindBinIndexFromEdges(FineEnergyBinCenter,CoarseEnergyAxisBinEdges);

      for (size_t iNuType=0;iNuType<nNeutrinoTypes;iNuType++) {
	int NuType = NeutrinoTypes[iNuType];
	
	for (size_t iOscChan=0;iOscChan<nOscillationChannels;iOscChan++) {
	  const FLOAT_T* OscProbPointer;
	  
	  if (!fCosineZIgnored) {
	    OscProbPointer = ReturnPointerToWeightinCalcer(NuType*OscillationChannels[iOscChan].GeneratedFlavour,NuType*OscillationChannels[iOscChan].DetectedFlavour,FineEnergyBinCenter,FineCosineZBinCenter);
	  } else {
	    OscProbPointer = ReturnPointerToWeightinCalcer(NuType*OscillationChannels[iOscChan].GeneratedFlavour,NuType*OscillationChannels[iOscChan].DetectedFlavour,FineEnergyBinCenter);
	  }
	  
	  int GlobalBin = iNuType*nOscillationChannels*nCoarseCosineZBins*nCoarseEnergyBins + iOscChan*nCoarseCosineZBins*nCoarseEnergyBins + CoarseCosineZBin*nCoarseEnergyBins + CoarseEnergyBin;
	  if ((GlobalBin < 0) || (GlobalBin >= TotalCoarseBins)) {
	    std::cerr << "Invalid global bin found:" << GlobalBin << std::endl;

	    std::cerr << "iNuType: " << iNuType << std::endl;
	    std::cerr << "iOscChan: " << iOscChan << std::endl;
	    std::cerr << "CoarseCosineZBin: " << CoarseCosineZBin << std::endl;
	    std::cerr << "CoarseEnergyBin: " << CoarseEnergyBin << std::endl;
	    
	    throw std::runtime_error("Fatal error in OscillatorSubSampling::SetupOscillatorImplementation()");
	  }
	  OscillationProbabilitiesToAverage[GlobalBin].push_back(OscProbPointer);
	}
      }
    }
  }

}

//Should do something smart like binary search but it should work for testing
int OscillatorSubSampling::FindBinIndexFromEdges(FLOAT_T Val, std::vector<FLOAT_T> BinEdges) {
  if (Val < BinEdges[0] || Val >= BinEdges[BinEdges.size()-1]) {
    std::cerr << "Value is not defined within the range of the bin edges!" << std::endl;
    std::cerr << "Val:" << Val << std::endl;
    std::cerr << "BinEdges[0]:" << BinEdges[0] << std::endl;
    std::cerr << "BinEdges[BinEdges.size()-1]:" << BinEdges[BinEdges.size()-1] << std::endl;
    throw std::runtime_error("Fatal error in OscillatorSubSampling::FindBinIndexFromEdges()") ;
  }
  for (size_t iBin=0;iBin<BinEdges.size()-1;iBin++) {
    if ((Val >= BinEdges[iBin]) && (Val < BinEdges[iBin+1])) {
      return iBin;
    }
  }

  std::cerr << "Did not find a sufficient bin for the value" << std::endl;
  throw std::runtime_error("Fatal error in OscillatorSubSampling::FindBinIndexFromEdges()");

  return -1;
}

const FLOAT_T* OscillatorSubSampling::ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal) {
  int CoarseEnergyBin = FindBinIndexFromEdges(EnergyVal,CoarseEnergyAxisBinEdges);
  int CoarseCosineZBin = FindBinIndexFromEdges(CosineZVal,CoarseCosineZAxisBinEdges);

  int OscChanIndex = -1;
  for (size_t iOscChan=0;iOscChan<OscillationChannels.size();iOscChan++) {
    if (OscillationChannels[iOscChan].GeneratedFlavour == std::abs(InitNuFlav) && OscillationChannels[iOscChan].DetectedFlavour == std::abs(FinalNuFlav)) {
      OscChanIndex = iOscChan;
      break;
    }
  }
  if (OscChanIndex == -1) {
    std::cerr << "Did not find valid oscillation channel" << std::endl;
    throw std::runtime_error("Fatal error in OscillatorSubSampling::ReturnWeightPointer()");
  }

  if (InitNuFlav*FinalNuFlav < 0) {
    std::cerr << "Invalid InitNuFlav and FinalNuFlav" << std::endl;
    std::cerr << "InitNuFlav: " << InitNuFlav << std::endl;
    std::cerr << "FinalNuFlav: " << FinalNuFlav << std::endl;
    throw std::runtime_error("Fatal error in OscillatorSubSampling::ReturnWeightPointer()");
  }
  int NuTypeIndex = fOscProbCalcer->ReturnNuTypeFromFlavour(InitNuFlav);

  int GlobalBin = NuTypeIndex*nOscillationChannels*nCoarseCosineZBins*nCoarseEnergyBins + OscChanIndex*nCoarseCosineZBins*nCoarseEnergyBins + CoarseCosineZBin*nCoarseEnergyBins + CoarseEnergyBin;
  if ((GlobalBin < 0) || (GlobalBin > AveragedOscillationProbabilities.size())) {
    std::cerr << "Invalid Global Bin index in OscillatorSubSampling::ReturnWeightPointer" << std::endl;
    std::cerr << "CoarseEnergyBin:" << CoarseEnergyBin << std::endl;
    std::cerr << "CoarseCosineZBin:" << CoarseCosineZBin << std::endl;
    std::cerr << "GlobalBin:" << GlobalBin << std::endl;
    throw std::runtime_error("Fatal error in OscillatorSubSampling::ReturnWeightPointer()");
  }
  
  return &(AveragedOscillationProbabilities[GlobalBin]);
}

void OscillatorSubSampling::PostCalculateProbabilities() {
  for (size_t iBin=0;iBin<AveragedOscillationProbabilities.size();iBin++) {
    FLOAT_T Avg = 0;
    for (size_t iPtr=0;iPtr<OscillationProbabilitiesToAverage[iBin].size();iPtr++) {
      Avg += *(OscillationProbabilitiesToAverage[iBin][iPtr]);
    }
    Avg /= OscillationProbabilitiesToAverage[iBin].size();

    if(OscillationProbabilitiesToAverage[iBin].size()==0){
      std::cerr<<"Division by zero encountered in SubSampling."<<std::endl;
      throw std::runtime_error("Fatal error in OscillatorSubSampling::PostCalculateProbabilities()");
    }

    AveragedOscillationProbabilities[iBin] = Avg;
  }
}

std::vector<FLOAT_T> OscillatorSubSampling::ReturnBinEdgesForPlotting(bool ReturnEnergy) {
  if (ReturnEnergy) {
    return CoarseEnergyAxisBinEdges;
  } else {
    return CoarseCosineZAxisBinEdges;
  }
}
