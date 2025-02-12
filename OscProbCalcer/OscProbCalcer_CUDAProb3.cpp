#include "OscProbCalcer_CUDAProb3.h"

#include "constants.hpp"
#include "propagator.hpp"
#include "physics.hpp"

#if UseGPU == 1
#include "cudapropagator.cuh"
#else
#include "cpupropagator.hpp"
#endif

#include <iostream>
using namespace cudaprob3;

OscProbCalcerCUDAProb3::OscProbCalcerCUDAProb3(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config
  EarthDensityFile = Config_["OscProbCalcerSetup"]["EarthModelFileName"].as<std::string>();
  if (fVerbose >= NuOscillator::INFO){std::cout << "EarthDensityFile:" << EarthDensityFile << std::endl;}
  UseProductionHeightsAve = Config_["OscProbCalcerSetup"]["UseProductionHeightsAveraging"].as<bool>();
  if(UseProductionHeightsAve){
    ProductionHeightsFile = Config_["OscProbCalcerSetup"]["ProductionHeightsFileName"].as<std::string>();
    if (fVerbose >= NuOscillator::INFO){std::cout << "ProductionHeightsFile:" << ProductionHeightsFile << std::endl;}
    // Read histogram suffixes 
    ProductionHeightsHistFlavourSuffixes.resize(6);
    ProductionHeightsHistFlavourSuffixes[0] = Config_["OscProbCalcerSetup"]["ProductionHeightsHistFlavourSuffixes"]["Nue"].as<std::string>();
    ProductionHeightsHistFlavourSuffixes[1] = Config_["OscProbCalcerSetup"]["ProductionHeightsHistFlavourSuffixes"]["Numu"].as<std::string>();
    ProductionHeightsHistFlavourSuffixes[2] = Config_["OscProbCalcerSetup"]["ProductionHeightsHistFlavourSuffixes"]["Nutau"].as<std::string>();
    ProductionHeightsHistFlavourSuffixes[3] = Config_["OscProbCalcerSetup"]["ProductionHeightsHistFlavourSuffixes"]["Nuebar"].as<std::string>();
    ProductionHeightsHistFlavourSuffixes[4] = Config_["OscProbCalcerSetup"]["ProductionHeightsHistFlavourSuffixes"]["Numubar"].as<std::string>();
    ProductionHeightsHistFlavourSuffixes[5] = Config_["OscProbCalcerSetup"]["ProductionHeightsHistFlavourSuffixes"]["Nutaubar"].as<std::string>();
  }

  fNOscParams = kNOscParams;

  UseEarthModelSystematics = Config_["OscProbCalcerSetup"]["UseEarthModelSystematics"].as<bool>();
  if(UseEarthModelSystematics){
    if (fVerbose >= NuOscillator::INFO){std::cout<<"Using Earth Model systematics"<<std::endl;}
    nLayers = Config_["OscProbCalcerSetup"]["Layers"].as<int>();
    fNOscParams +=  2*nLayers;
  }
  //=======

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  // Implementation specific variables
  OscChannels = std::vector<int>(fNOscillationChannels,DUMMYVAL);
  for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {
    if (fOscillationChannels[iOscChannel].GeneratedFlavour == NuOscillator::kElectron) {
      if (fOscillationChannels[iOscChannel].DetectedFlavour == NuOscillator::kElectron) {OscChannels[iOscChannel] = e_e;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == NuOscillator::kMuon) {OscChannels[iOscChannel] = e_m;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == NuOscillator::kTau) {OscChannels[iOscChannel] = e_t;}
    } else if (fOscillationChannels[iOscChannel].GeneratedFlavour == NuOscillator::kMuon) {
      if (fOscillationChannels[iOscChannel].DetectedFlavour == NuOscillator::kElectron) {OscChannels[iOscChannel] = m_e;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == NuOscillator::kMuon) {OscChannels[iOscChannel] = m_m;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == NuOscillator::kTau) {OscChannels[iOscChannel] = m_t;}
    } else if (fOscillationChannels[iOscChannel].GeneratedFlavour == NuOscillator::kTau) {
      if (fOscillationChannels[iOscChannel].DetectedFlavour == NuOscillator::kElectron) {OscChannels[iOscChannel] = t_e;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == NuOscillator::kMuon) {OscChannels[iOscChannel] = t_m;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == NuOscillator::kTau) {OscChannels[iOscChannel] = t_t;}
    }
  }  

  nThreads = 1;
}

OscProbCalcerCUDAProb3::~OscProbCalcerCUDAProb3() {

}

void OscProbCalcerCUDAProb3::SetupPropagator() {

#if UseGPU == 1
  if (fVerbose >= NuOscillator::INFO) {std::cout << "Using GPU CUDAProb3 propagator" << std::endl;}
  propagator = std::unique_ptr< Propagator< FLOAT_T > > ( new CudaPropagatorSingle<FLOAT_T>(0, fNCosineZPoints, fNEnergyPoints)); // Single-GPU propagator
  fImplementationName += "-GPU";
#else

  nThreads = 1;
#if UseMultithreading == 1
#pragma omp parallel
  {
#pragma omp single
    nThreads = omp_get_num_threads();
  }
#endif

  if (fVerbose >= NuOscillator::INFO) {std::cout << "Using CPU CUDAProb3 propagator with " << nThreads << " threads" << std::endl;}
  propagator = std::unique_ptr< Propagator< FLOAT_T > > ( new CpuPropagator<FLOAT_T>(fNCosineZPoints, fNEnergyPoints, nThreads)); // MultiThread CPU propagator
  fImplementationName += "-CPU-"+std::to_string(nThreads);
#endif

  propagator->setEnergyList(fEnergyArray);
  propagator->setCosineList(fCosineZArray);
  propagator->setDensityFromFile(EarthDensityFile);

  if (fVerbose >= NuOscillator::INFO) {std::cout << "Setup CUDAProb3 oscillation probability calculater" << std::endl;}
}
 
void OscProbCalcerCUDAProb3::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
  // Oscpars, as given from MaCh3, expresses the mixing angles in sin^2(theta). This propagator expects them in theta
  for (int iOscPar=0;iOscPar<=kTH13;iOscPar++) {
    if (OscParams[iOscPar] < 0) {
      std::cerr << "Invalid oscillation parameter (Can not sqrt this value)!:" << OscParams[iOscPar] << std::endl;
      throw std::runtime_error("Invalid setup");
    }
  }

  FLOAT_T theta12 = asin(sqrt(OscParams[kTH12]));
  FLOAT_T theta23 = asin(sqrt(OscParams[kTH23]));
  FLOAT_T theta13 = asin(sqrt(OscParams[kTH13]));
  FLOAT_T dm12sq  = OscParams[kDM12];
  FLOAT_T dm23sq  = OscParams[kDM23];
  FLOAT_T dcp     = OscParams[kDCP];
  FLOAT_T prodH   = OscParams[kPRODH];

  propagator->setNeutrinoMasses(dm12sq, dm23sq);
  propagator->setProductionHeight(prodH);

  if(UseProductionHeightsAve){
    SetProductionHeightsAveraging();
  }

  if(UseEarthModelSystematics){
    ApplyEarthModelSystematics(OscParams);
  }

  // CUDAProb3 calculates oscillation probabilites for each NeutrinoType, so need to copy them from the calculator into fWeightArray
  int CopyArrSize = fNEnergyPoints * fNCosineZPoints;
  FLOAT_T* CopyArr = new FLOAT_T[CopyArrSize];

  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {

    NeutrinoType NuType;
    if (fNeutrinoTypes[iNuType]==Nubar) {
      // Haven't really thought about it, but prob3++ sets dcp->-dcp here: https://github.com/rogerwendell/Prob3plusplus/blob/fd189e232e96e2c5ebb2f7bd3a5406b288228e41/BargerPropagator.cc#L235
      // Copying that behaviour gives same behaviour as prob3++/probGPU
      propagator->setMNSMatrix(theta12, theta13, theta23, -dcp);
      NuType = cudaprob3::Antineutrino;
    } else {
      propagator->setMNSMatrix(theta12, theta13, theta23, dcp);
      NuType = cudaprob3::Neutrino;
    }

    propagator->calculateProbabilities(NuType);

    for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {
      propagator->getProbabilityArr(CopyArr,static_cast<cudaprob3::ProbType>(OscChannels[iOscChannel]));
      
      // Mapping which links the oscillation channel, neutrino type and energy/cosineZ index to the fWeightArray index
      int IndexToFill = iNuType*fNOscillationChannels*CopyArrSize + iOscChannel*CopyArrSize;
      for (int iOscProb=0;iOscProb<CopyArrSize;iOscProb++) {
	// Sometimes CUDAProb3 can return *slightly* unphysical oscillation probabilities
	CopyArr[iOscProb] = CopyArr[iOscProb] > 0.0 ? CopyArr[iOscProb] : 0.0;
	CopyArr[iOscProb] = CopyArr[iOscProb] < 1.0 ? CopyArr[iOscProb] : 1.0;
	fWeightArray[IndexToFill+iOscProb] = CopyArr[iOscProb];
      }
    }
  }

  delete[] CopyArr;
}

int OscProbCalcerCUDAProb3::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNCosineZPoints*fNEnergyPoints + OscChanIndex*fNCosineZPoints*fNEnergyPoints + EnergyIndex*fNCosineZPoints + CosineZIndex;
  return IndexToReturn;
}

long OscProbCalcerCUDAProb3::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * fNCosineZPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}

void OscProbCalcerCUDAProb3::SetProductionHeightsAveraging(){
  if (fVerbose >= NuOscillator::INFO){std::cout<<"Setting up production heights for averaging..."<<std::endl;}
  // Open prod heights file
  TFile* File = new TFile(ProductionHeightsFile.c_str());
  if (!File || File->IsZombie()) {
    std::cerr << "Could not open: " << ProductionHeightsFile << std::endl;
    throw;
  }

  // Read TH2Ds
  std::vector<std::vector<TString>> NeutrinoFlavourNames(2);
  int NNeutrinoFlavours = 3;
  NeutrinoFlavourNames[0].resize(NNeutrinoFlavours);
  NeutrinoFlavourNames[1].resize(NNeutrinoFlavours);
  NeutrinoFlavourNames[0][0] = ProductionHeightsHistFlavourSuffixes[0].c_str();
  NeutrinoFlavourNames[0][1] = ProductionHeightsHistFlavourSuffixes[1].c_str();
  NeutrinoFlavourNames[0][2] = ProductionHeightsHistFlavourSuffixes[2].c_str();
  NeutrinoFlavourNames[1][0] = ProductionHeightsHistFlavourSuffixes[3].c_str();
  NeutrinoFlavourNames[1][1] = ProductionHeightsHistFlavourSuffixes[4].c_str();
  NeutrinoFlavourNames[1][2] = ProductionHeightsHistFlavourSuffixes[5].c_str();

  std::vector<std::vector<TH3D*>> vecHist;
  vecHist.resize(fNNeutrinoTypes);
  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++){
    vecHist[iNuType].resize(NNeutrinoFlavours);
  }

  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++){
    for (int iNuFlav=0;iNuFlav<NNeutrinoFlavours;iNuFlav++){
      TString HistName = "ProductionHeight_"+NeutrinoFlavourNames[iNuType][iNuFlav];
      TH3D* Hist = (TH3D*)File->Get(HistName);

      if(!Hist){
        std::cerr << HistName << " not found in File:" << ProductionHeightsFile << std::endl;
        File->ls();
        std::cerr << __LINE__ << " : " << __FILE__ << std::endl;
        throw;
      }

      vecHist[iNuType][iNuFlav] = Hist;

      if(vecHist[iNuType][iNuFlav]->GetNbinsX()!=fNEnergyPoints){
        std::cerr << HistName << " has different number of X bins:" << vecHist[iNuType][iNuFlav]->GetNbinsX() << std::endl;
        std::cerr << "Expected:" << fNEnergyPoints << std::endl;
        std::cerr << __LINE__ << " : " << __FILE__ << std::endl;
        throw;
      }

      if(vecHist[iNuType][iNuFlav]->GetNbinsY() != fNCosineZPoints){
        std::cerr << HistName << " has different number of Y bins:" << vecHist[iNuType][iNuFlav]->GetNbinsY() << std::endl;
        std::cerr << "Expected:" << fNCosineZPoints << std::endl;
        std::cerr << __LINE__ << " : " << __FILE__ << std::endl;
        throw;
      }
    }
  }

  // Get number of height points
  int NProductionHeightAveragingBins = vecHist[0][0]->GetNbinsZ();
  if(NProductionHeightAveragingBins>cudaprob3::Constants<FLOAT_T>::MaxProdHeightBins()){
    std::cerr << "Different number of height bins:" << NProductionHeightAveragingBins << std::endl;
    std::cerr << "Expected:" << cudaprob3::Constants<FLOAT_T>::MaxProdHeightBins() << std::endl;
  }

  // Make 1D array with probabilities
  int ProductionHeightProbabilitiesListSize = NNeutrinoFlavours*fNNeutrinoTypes*fNCosineZPoints*fNEnergyPoints*NProductionHeightAveragingBins;
  std::vector<FLOAT_T> ProductionHeightProbabilitiesList(ProductionHeightProbabilitiesListSize);
  int index = 0;
  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++){
    for(int iNuFlav=0;iNuFlav<NNeutrinoFlavours;iNuFlav++){
      for (int ibin_E=0;ibin_E<fNEnergyPoints;ibin_E++){
        for (int ibin_Z=0;ibin_Z<fNCosineZPoints;ibin_Z++){
          double Total = 0.;

          for (int iProductionHeight=0;iProductionHeight<NProductionHeightAveragingBins;iProductionHeight++){
            double dP_dh = vecHist[iNuType][iNuFlav]->GetBinContent(ibin_E+1,ibin_Z+1,iProductionHeight+1);
            double dh = vecHist[iNuType][iNuFlav]->GetZaxis()->GetBinWidth(iProductionHeight+1);

            ProductionHeightProbabilitiesList[index] = dP_dh * dh;
            Total += ProductionHeightProbabilitiesList[index];

            index += 1;
          }

          if (fabs(Total-1.) > 1e-6) {
            std::cerr << "Probabilities integrated over production height do not sum to 1" << std::endl;
            std::cerr << "Total:" << Total << std::endl;
            for (int iProductionHeight=0;iProductionHeight<NProductionHeightAveragingBins;iProductionHeight++) {
              std::cout << "iProductionHeight:" << iProductionHeight << " | dP_dh:" << vecHist[iNuType][iNuFlav]->GetBinContent(ibin_E+1,ibin_Z+1,iProductionHeight+1) << std::endl;
            }
            std::cerr << __LINE__ << " : " << __FILE__ << std::endl;
            throw;
          }
        }
      }
    }
  }

  // Make list of heights
  std::vector<FLOAT_T> ProductionHeightsList(NProductionHeightAveragingBins+1);
  for(int iBinH=0; iBinH<NProductionHeightAveragingBins+1; iBinH++){
    ProductionHeightsList[iBinH] = vecHist[0][0]->GetZaxis()->GetBinLowEdge(iBinH+1);
  }
  
  // Set in propagator
  propagator->SetNumberOfProductionHeightBinsForAveraging(NProductionHeightAveragingBins);
  propagator->setProductionHeightList(ProductionHeightProbabilitiesList,ProductionHeightsList);
  
  if (fVerbose >= NuOscillator::INFO){std::cout<<"Completed SetProductionHeightsAveraging()"<<std::endl;}
}

void OscProbCalcerCUDAProb3::ApplyEarthModelSystematics(const std::vector<FLOAT_T>& OscParams){
  int n_layers = propagator->getNlayerBoundaries()-1;
  if(n_layers!=nLayers){
      std::cerr<<"Number of layers set in config differs from the one in "<<EarthDensityFile<<std::endl;
      std::cerr<<"Expected: "<<nLayers<<"    Got: "<<n_layers<<std::endl;
  }

  std::vector<FLOAT_T> EarthBoundaries(nLayers);
  std::vector<FLOAT_T> EarthWeights(nLayers);

  int kLayerBoundaries = kPRODH + 1;
  int kLayerWeights = kLayerBoundaries + nLayers;
    
  for(int iLayer = 0; iLayer<nLayers; iLayer++){
    EarthBoundaries[iLayer] = OscParams[kLayerBoundaries+iLayer];
    EarthWeights[iLayer] = OscParams[kLayerWeights+iLayer];
  }
    
  // Check if the model is a set of polynomials
  if(propagator->PolynomialDensity()){
    propagator->ModifyEarthModelPoly(EarthBoundaries, EarthWeights);
  }
  else{
    propagator->ModifyEarthModel(EarthBoundaries, EarthWeights);
  }
}