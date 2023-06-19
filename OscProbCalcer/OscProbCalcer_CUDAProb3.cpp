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

OscProbCalcerCUDAProb3::OscProbCalcerCUDAProb3(std::string ConfigName_) : OscProbCalcerBase(ConfigName_)
{
  //=======
  //Grab information from the config
  std::string EarthDensityModelFileName = Config["CUDAProb3"]["EarthDensityModel"].as<std::string>();

  char* EnvironVal = std::getenv("CUDAProb3Source");
  if (EnvironVal == NULL) {
    std::cerr << "CUDAProb3Source environment variable is not defined!" << std::endl;
    throw;
  }
  EarthDensityFile = std::string(EnvironVal)+"/models/"+EarthDensityModelFileName;
  //=======

  fImplementationName = "CUDAProb3";
  fNOscParams = kNOscParams;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  fNInitialFlavours = 2;
  InitialiseInitialFlavoursArray(fNInitialFlavours);
  fInitialFlavours[0] = Electron;
  fInitialFlavours[1] = Muon;

  fNFinalFlavours = 3;
  InitialiseFinalFlavoursArray(fNFinalFlavours);
  fFinalFlavours[0] = Electron;
  fFinalFlavours[1] = Muon;
  fFinalFlavours[2] = Tau;

  // Implementation specific variables
  OscChannels.resize(fNInitialFlavours);
  for (int i=0;i<fNInitialFlavours;i++) {
    OscChannels[i].resize(fNFinalFlavours);
  }
  OscChannels[0][0] = e_e;
  OscChannels[0][1] = e_m;
  OscChannels[0][2] = e_t;
  OscChannels[1][0] = m_e;
  OscChannels[1][1] = m_m;
  OscChannels[1][2] = m_t;

  nThreads = 0;
}

void OscProbCalcerCUDAProb3::SetupPropagator() {

#if UseGPU == 1
  if (fVerbose >= INFO) {std::cout << "Using GPU CUDAProb3 propagator" << std::endl;}
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

  if (fVerbose >= INFO) {std::cout << "Using CPU CUDAProb3 propagator with " << nThreads << " threads" << std::endl;}
  propagator = std::unique_ptr< Propagator< FLOAT_T > > ( new CpuPropagator<FLOAT_T>(fNCosineZPoints, fNEnergyPoints, nThreads)); // MultiThread CPU propagator
  fImplementationName += "-CPU-"+std::to_string(nThreads);
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

    for (int iInitFlav=0;iInitFlav<fNInitialFlavours;iInitFlav++) {
      for (int iFinalFlav=0;iFinalFlav<fNFinalFlavours;iFinalFlav++) {
	propagator->getProbabilityArr(CopyArr,static_cast<cudaprob3::ProbType>(OscChannels[iInitFlav][iFinalFlav]));
	
	// Mapping which links the oscillation channel, neutrino type and energy/cosineZ index to the fWeightArray index
	int IndexToFill = iNuType*fNInitialFlavours*fNFinalFlavours*CopyArrSize + iInitFlav*fNFinalFlavours*CopyArrSize + iFinalFlav*CopyArrSize;
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
  int IndexToReturn = NuTypeIndex*fNInitialFlavours*fNFinalFlavours*fNCosineZPoints*fNEnergyPoints + InitNuIndex*fNFinalFlavours*fNCosineZPoints*fNEnergyPoints + FinalNuIndex*fNCosineZPoints*fNEnergyPoints + CosineZIndex*fNEnergyPoints + EnergyIndex;

  return IndexToReturn;
}

long OscProbCalcerCUDAProb3::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNCosineZPoints * fNInitialFlavours * fNFinalFlavours * fNNeutrinoTypes;
  return nCalculationPoints;
}
