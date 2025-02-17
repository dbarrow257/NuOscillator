#include "Oscillator/OscillatorBinned.h"

#include <iostream>

#include "TFile.h"
#include "TH1.h"

OscillatorBinned::OscillatorBinned(std::string ConfigName_) : OscillatorBase(ConfigName_) {
  Initialise();
}

OscillatorBinned::OscillatorBinned(YAML::Node Config_) : OscillatorBase(Config_) {
  Initialise();
}

void OscillatorBinned::Initialise() {
  EnergyAxisBinEdges = std::vector<FLOAT_T>();
  CosineZAxisBinEdges = std::vector<FLOAT_T>();
  EnergyAxisBinCenters = std::vector<FLOAT_T>();
  CosineZAxisBinCenters = std::vector<FLOAT_T>();

  fCalculationTypeName = "Binned";

  //=======
  // Grab the following from config manager

  FileName = Config[fCalculationTypeName]["FileName"].as<std::string>();
  EnergyAxisHistName = Config[fCalculationTypeName]["EnergyAxisHistName"].as<std::string>();
  CosineZAxisHistName = Config[fCalculationTypeName]["CosineZAxisHistName"].as<std::string>();
  //=======

  EnergyAxisBinEdges = ReadBinEdgesFromFile(FileName,EnergyAxisHistName);
  EnergyAxisBinCenters = ReturnBinCentersFromBinEdges(EnergyAxisBinEdges);
  if (!fCosineZIgnored) {
    CosineZAxisBinEdges = ReadBinEdgesFromFile(FileName,CosineZAxisHistName);
    CosineZAxisBinCenters = ReturnBinCentersFromBinEdges(CosineZAxisBinEdges);
  }

  fEvalPointsSetInConstructor = true;

  SetEnergyArrayInCalcer(EnergyAxisBinCenters);
  if (!fCosineZIgnored) {
    SetCosineZArrayInCalcer(CosineZAxisBinCenters);
  }

}

OscillatorBinned::~OscillatorBinned() {
}

const FLOAT_T* OscillatorBinned::ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal) {
  FLOAT_T EnergyValBinCenter = DUMMYVAL;
  FLOAT_T CosineZValBinCenter = DUMMYVAL;

  int nEnergyBins = EnergyAxisBinCenters.size();
  if (EnergyVal < EnergyAxisBinEdges[0] || EnergyVal >= EnergyAxisBinEdges[nEnergyBins]) {
    std::cerr << "Requested Energy is not within the range of pre-defined binning (EnergyAxisBinEdges)" << std::endl;
    std::cerr << "EnergyVal:" << EnergyVal << std::endl;
    std::cerr << "EnergyAxisBinEdges[0]:" << EnergyAxisBinEdges[0] << std::endl;
    std::cerr << "nEnergyBins:" << nEnergyBins << std::endl;
    std::cerr << "EnergyAxisBinEdges[nEnergyBins+1]:" << EnergyAxisBinEdges[nEnergyBins+1] << std::endl;
  }
  int EnergyIndex = -1;
  for (int iBin=0;iBin<nEnergyBins;iBin++) {
    if (EnergyVal >= EnergyAxisBinEdges[iBin] && EnergyVal < EnergyAxisBinEdges[iBin+1]) {
      EnergyIndex = iBin;
      break;
    }
  }
  if (EnergyIndex == -1) {
    std::cerr << "Invalid bin found in OscillatorBinned::ReturnWeightPointer - Did not find the correct bin for Energy:" << EnergyVal << std::endl;
    throw std::runtime_error("Invalid setup");
  }
  EnergyValBinCenter = EnergyAxisBinCenters[EnergyIndex];

  if (!fCosineZIgnored) {
    int nCosineZBins = CosineZAxisBinCenters.size();
    if (CosineZVal < CosineZAxisBinEdges[0] || CosineZVal >= CosineZAxisBinEdges[nCosineZBins]) {
      std::cerr << "Requested CosineZ is not within the range of pre-defined binning (CosineZAxisBinEdges)" << std::endl;
      std::cerr << "CosineZVal:" << CosineZVal << std::endl;
      std::cerr << "CosineZAxisBinEdges[0]:" << CosineZAxisBinEdges[0] << std::endl;
      std::cerr << "nCosineZBins:" << nCosineZBins << std::endl;
      std::cerr << "CosineZAxisBinEdges[nCosineZBins+1]:" << CosineZAxisBinEdges[nCosineZBins+1] << std::endl;
    }
    int CosineZIndex = -1;
    for (int iBin=0;iBin<nCosineZBins;iBin++) {
      if (CosineZVal >= CosineZAxisBinEdges[iBin] && CosineZVal < CosineZAxisBinEdges[iBin+1]) {
	CosineZIndex = iBin;
	break;
      }
    }
    if (CosineZIndex == -1) {
      std::cerr << "Invalid bin found in OscillatorBinned::ReturnWeightPointer - Did not find the correct bin for CosineZ:" << CosineZVal << std::endl;
      throw std::runtime_error("Invalid setup");
    }
    CosineZValBinCenter = CosineZAxisBinCenters[CosineZIndex];
  }

  return ReturnPointerToWeightinCalcer(InitNuFlav,FinalNuFlav,EnergyValBinCenter,CosineZValBinCenter);
}

std::vector<FLOAT_T> OscillatorBinned::ReturnBinEdgesForPlotting(bool ReturnEnergy) {
  if (ReturnEnergy) {
    return EnergyAxisBinEdges;
  } else {
    return CosineZAxisBinEdges;
  }
}
