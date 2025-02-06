#include "Constants/OscillatorConstants.h"

#include <vector>

#include "TFile.h"
#include "TH1D.h"

int main(int argc, char **argv) {
  std::vector<FLOAT_T> EnergyArray = logspace(0.1,100.,50);
  std::vector<FLOAT_T> CosineZArray = linspace(-1.0,1.0,15);

  TFile* File = TFile::Open("ExampleAtmosphericBinning.root","RECREATE");
  TH1D* EnergyBinning = new TH1D("EnergyAxisBinning",";Energy [GeV]",EnergyArray.size()-1,EnergyArray.data());
  TH1D* CosineZBinning = new TH1D("CosineZAxisBinning",";CosineZ",CosineZArray.size()-1,CosineZArray.data());

  EnergyBinning->Write();
  CosineZBinning->Write();
  File->Close();
}
