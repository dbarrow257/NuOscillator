#include "OscProbCalcer_CUDAProb3.h"

#include "constants.hpp"
#include "propagator.hpp"
#include "physics.hpp"

#include "cpupropagator.hpp"

#include <iostream>
using namespace cudaprob3;

OscProbCalcerCUDAProb3::OscProbCalcerCUDAProb3(std::string ConfigName_) : OscProbCalcerBase()
{
  ConfigName = ConfigName_; //DB Create yaml style config once built in MaCh3

  // Required variables
  fImplementationName = "CUDAProb3";
  fVerbose = INFO; //DB Get From Config

  fNOscParams = kNOscParams;

  nNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(nNeutrinoTypes);
  NeutrinoTypes[0] = Nu;
  NeutrinoTypes[1] = Nubar;

  nInitialFlavours = 2;
  InitialiseInitialFlavoursArray(nInitialFlavours);
  InitialFlavours[0] = Electron;
  InitialFlavours[1] = Muon;

  nFinalFlavours = 3;
  InitialiseFinalFlavoursArray(nFinalFlavours);
  FinalFlavours[0] = Electron;
  FinalFlavours[1] = Muon;
  FinalFlavours[2] = Tau;

  // Implementation specific variables
  OscChannels.resize(nInitialFlavours);
  for (int i=0;i<nInitialFlavours;i++) {
    OscChannels[i].resize(nFinalFlavours);
  }
  OscChannels[0][0] = e_e;
  OscChannels[0][1] = e_m;
  OscChannels[0][2] = e_t;
  OscChannels[1][0] = m_e;
  OscChannels[1][1] = m_m;
  OscChannels[1][2] = m_t;

  nThreads = 0;
  EarthDensityFile = "../CUDAProb3/models/PREM_4layer.dat"; //DB Get From Config
}

void OscProbCalcerCUDAProb3::SetupPropagator() {

#ifdef UseGPU
  if (fVerbose >= INFO) {std::cout << "Using GPU CUDAProb3 propagaator" << std::endl;}
  propagator = std::unique_ptr<Propagator<FLOAT_T>> ( new CudaPropagatorSingle<FLOAT_T>(0,nCosine, nEnergy)); // Single-GPU propagator
#else

  nThreads = 1;
#ifdef UseMultithread
#pragma omp parallel
  {
#pragma omp single
    nThreads = omp_get_num_threads();
  }
#endif

  if (fVerbose >= INFO) {std::cout << "Using CPU CUDAProb3 propagator with " << nThreads << " threads" << std::endl;}
  propagator = std::unique_ptr< Propagator< FLOAT_T > > ( new CpuPropagator<FLOAT_T>(fNCosineZPoints, fNEnergyPoints, nThreads)); // MultiThread CPU propagator
#endif

  propagator->setEnergyList(fEnergyArray);
  propagator->setCosineList(fCosineZArray);
  propagator->setDensityFromFile(EarthDensityFile);

  if (fVerbose >= INFO) {std::cout << "Setup CUDAProb3 oscillation probability calculater" << std::endl;}
}
 
void OscProbCalcerCUDAProb3::CalculateProbabilities(std::vector<FLOAT_T> OscParams) {
  // Oscpars, as given from MaCh3, expresses the mixing angles in sin^2(theta). This propagator expects them in theta
  for (int iOscPar=0;iOscPar<=kTH13;iOscPar++) {
    if (OscParams[iOscPar] < 0) {
      std::cerr << "Invalid oscillation parameter (Can not sqrt this value)!:" << OscParams[iOscPar] << std::endl;
      throw;
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

  // CUDAProb3 calculates oscillation probabilites for each NeutrinoType, so need to copy them from the calculator into fWeightArray
  int CopyArrSize = fNEnergyPoints * fNCosineZPoints;
  FLOAT_T* CopyArr = new FLOAT_T[CopyArrSize];

  for (int iNuType=0;iNuType<nNeutrinoTypes;iNuType++) {

    NeutrinoType NuType;
    if (NeutrinoTypes[iNuType]==Nubar) {
      // Haven't really thought about it, but prob3++ sets dcp->-dcp here: https://github.com/rogerwendell/Prob3plusplus/blob/fd189e232e96e2c5ebb2f7bd3a5406b288228e41/BargerPropagator.cc#L235
      // Copying that behaviour gives same behaviour as prob3++/probGPU
      propagator->setMNSMatrix(theta12, theta13, theta23, -dcp);
      NuType = cudaprob3::Antineutrino;
    } else {
      propagator->setMNSMatrix(theta12, theta13, theta23, dcp);
      NuType = cudaprob3::Neutrino;
    }

    propagator->calculateProbabilities(NuType);

    for (int iInitFlav=0;iInitFlav<nInitialFlavours;iInitFlav++) {
      for (int iFinalFlav=0;iFinalFlav<nFinalFlavours;iFinalFlav++) {
	propagator->getProbabilityArr(CopyArr,static_cast<cudaprob3::ProbType>(OscChannels[iInitFlav][iFinalFlav]));
	
	// Mapping which links the oscillation channel, neutrino type and energy/cosineZ index to the fWeightArray index
	int IndexToFill = iNuType*nInitialFlavours*nFinalFlavours*CopyArrSize + iInitFlav*nFinalFlavours*CopyArrSize + iFinalFlav*CopyArrSize;
	for (int iOscProb=0;iOscProb<CopyArrSize;iOscProb++) {

	  // Sometimes CUDAProb3 can return *slightly* unphysical oscillation probabilities
	  CopyArr[iOscProb] = CopyArr[iOscProb] > 0.0 ? CopyArr[iOscProb] : 0.0;
	  CopyArr[iOscProb] = CopyArr[iOscProb] < 1.0 ? CopyArr[iOscProb] : 1.0;

	  fWeightArray[IndexToFill+iOscProb] = CopyArr[iOscProb];
	}
      }
    }
  }

  delete[] CopyArr;
}

int OscProbCalcerCUDAProb3::ReturnWeightArrayIndex(int NuTypeIndex, int InitNuIndex, int FinalNuIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*nInitialFlavours*nFinalFlavours*fNCosineZPoints*fNEnergyPoints + InitNuIndex*nFinalFlavours*fNCosineZPoints*fNEnergyPoints + FinalNuIndex*fNCosineZPoints*fNEnergyPoints + CosineZIndex*fNEnergyPoints + EnergyIndex;

  return IndexToReturn;
}

long OscProbCalcerCUDAProb3::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNCosineZPoints * nInitialFlavours * nFinalFlavours * nNeutrinoTypes;
  return nCalculationPoints;
}
