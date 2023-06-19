#include "OscillatorBinned.h"

#include <iostream>

OscillatorBinned::OscillatorBinned(std::string ConfigName_) : OscillatorBase(ConfigName_) {
  EnergyAxisBinEdges = std::vector<FLOAT_T>();
  CosineZAxisBinEdges = std::vector<FLOAT_T>();
  EnergyAxisBinCenters = std::vector<FLOAT_T>();
  CosineZAxisBinCenters = std::vector<FLOAT_T>();

  fCalculationTypeName = "Binned";

  //=======
  // Grab the following from config manager
  //DB
  //fOscProbCalcerImplementationToCreate = OscProbCalcerImplementationToCreate_;

  FileName = Config[fCalculationTypeName]["FileName"].as<std::string>();
  EnergyAxisHistName = Config[fCalculationTypeName]["EnergyAxisHistName"].as<std::string>();
  CosineZAxisHistName = Config[fCalculationTypeName]["CosineZAxisHistName"].as<std::string>();
  //=======

  InitialiseOscProbCalcers();

  EnergyAxisBinEdges = ReadBinEdgesFromFile(FileName,EnergyAxisHistName);
  EnergyAxisBinCenters = ReturnBinCentersFromBinEdges(EnergyAxisBinEdges);
  for (int CalcerIndex=0;CalcerIndex<fNCalcers;CalcerIndex++) {
    SetEnergyArrayInCalcer(EnergyAxisBinCenters, CalcerIndex);
  }

  if (!fCosineZIgnored) {
    CosineZAxisBinEdges = ReadBinEdgesFromFile(FileName,CosineZAxisHistName);
    CosineZAxisBinCenters = ReturnBinCentersFromBinEdges(CosineZAxisBinEdges);
    
    for (int CalcerIndex=0;CalcerIndex<fNCalcers;CalcerIndex++) {
      SetCosineZArrayInCalcer(CosineZAxisBinCenters, CalcerIndex);
    }
  }
}

std::vector<FLOAT_T> OscillatorBinned::ReadBinEdgesFromFile(std::string FileName, std::string HistogramName) {
  std::vector<FLOAT_T> BinEdges;

  //DB Once ROOT is linked up, use TH1 - Assumes each axis is independently binned (Could be changed with another implementation)
  /*
  TFile* File = new TFile(FileName); 
  if (!File || File->IsZombie()) {
    std::cerr << "Could not find file:" << FileName << std::endl;
    throw;
  }

  TH1* Histogram = (TH1*)File->Get(HistogramName);
  if (!Histogram) {
    std::cerr << "Could not find Histogram:" << HistogramName << " in File:" << FileNaame << std::endl;
    throw;
  }

  BinEdges.resize(Histogram->GetNBinsX()+1);
  for (size_t iBin=0;iBin<=Histogram->GetNBinsX();iBin++) {
    BinEdges[iBin] = Histogram->GetBinLowEdge(iBin+1);
  }
  */

  //===========
  //DB Once ROOT is linked up, remove the following code
  if (HistogramName=="CosineZHistName") {
    int NBins = 1;
    for (int iBin=0;iBin<NBins;iBin++) {
      BinEdges.push_back(iBin*0.002-1.0);
    }
  } else {
    int NBins = 1001;
    for (int iBin=0;iBin<NBins;iBin++) {
      BinEdges.push_back(iBin*0.1);
    }
  }
  //===========

  if (fVerbose >= INFO) {
    std::cout << "Bin edges successfully read from File:" << FileName << " , Histogram:" << HistogramName << " :=" << std::endl;
    for (size_t i=0;i<BinEdges.size();i++) {
      std::cout << BinEdges[i] << ", ";
    }
    std::cout << std::endl;
  }

  return BinEdges;
}

std::vector<FLOAT_T> OscillatorBinned::ReturnBinCentersFromBinEdges(std::vector<FLOAT_T> BinEdges) {
  int nBins = BinEdges.size();
  std::vector<FLOAT_T> BinCenters = std::vector<FLOAT_T>(nBins);

  for (int iBin=0;iBin<nBins;iBin++) {
    BinCenters[iBin] = (BinEdges[iBin]+BinEdges[iBin+1])/2.0;
  }

  return BinCenters;
}

const FLOAT_T* OscillatorBinned::ReturnWeightPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T EnergyVal, FLOAT_T CosineZVal) {
  FLOAT_T EnergyValBinCenter = DUMMYVAL;
  FLOAT_T CosineZValBinCenter = DUMMYVAL;

  int nEnergyBins = EnergyAxisBinCenters.size();
  if (EnergyVal < EnergyAxisBinEdges[0] || EnergyVal >= EnergyAxisBinEdges[nEnergyBins+1]) {
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
    throw;
  }
  EnergyValBinCenter = EnergyAxisBinCenters[EnergyIndex];

  if (fCosineZIgnored) {
    int nCosineZBins = CosineZAxisBinCenters.size();
    if (CosineZVal < CosineZAxisBinEdges[0] || CosineZVal >= CosineZAxisBinEdges[nCosineZBins+1]) {
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
      throw;
    }
    CosineZValBinCenter = CosineZAxisBinCenters[CosineZIndex];
  }

  int CalcerIndex = 0;
  return ReturnPointerToWeightinCalcer(CalcerIndex,InitNuFlav,FinalNuFlav,EnergyValBinCenter,CosineZValBinCenter);
}
