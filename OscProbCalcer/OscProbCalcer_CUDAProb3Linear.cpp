#include "OscProbCalcer_CUDAProb3Linear.h"

#if UseGPU == 1
#include "beamcudapropagator.cuh"
#else
#include "beamcpupropagator.hpp"
#endif

#include <iostream>

using namespace cudaprob3;

OscProbCalcerCUDAProb3Linear::OscProbCalcerCUDAProb3Linear(std::string ConfigName_) : OscProbCalcerBase(ConfigName_)
{
  fImplementationName = "CUDAProb3Linear";
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
  OscChannels[0][0] = cudaprob3::e_e;
  OscChannels[0][1] = cudaprob3::e_m;
  OscChannels[0][2] = cudaprob3::e_t;
  OscChannels[1][0] = cudaprob3::m_e;
  OscChannels[1][1] = cudaprob3::m_m;
  OscChannels[1][2] = cudaprob3::m_t;

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);

  nThreads = 0;
}

void OscProbCalcerCUDAProb3Linear::SetupPropagator() {

#if UseGPU == 1
  if (fVerbose >= INFO) {std::cout << "Using GPU CUDAProb3Linear propagator" << std::endl;}
  //propagator = std::unique_ptr< Propagator<FLOAT_T> > ( new BeamCudaPropagatorSingle(0, fNEnergyPoints));

  cudaprob3::Propagator<double> *Oscillator = new BeamCudaPropagatorSingle(0,fNEnergyPoints);

  //cudaprob3::BeamCudaPropagatorSingle* mypropagator = new cudaprob3::BeamCudaPropagatorSingle(0, fNEnergyPoints);
  //propagator = std::unique_ptr< cudaprob3::Propagator<FLOAT_T> > (mypropagator);

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

  if (fVerbose >= INFO) {std::cout << "Using CPU CUDAProb3Linear propagator with fNEnergyPoints:" << fNEnergyPoints << " and:" << nThreads << " threads" << std::endl;}
  propagator = std::unique_ptr< Propagator< FLOAT_T > > ( new BeamCpuPropagator<FLOAT_T>(fNEnergyPoints, nThreads)); // MultiThread CPU propagator
  fImplementationName += "-CPU-"+std::to_string(nThreads);
#endif

  std::cout << "fEnergyArray.size():" << fEnergyArray.size() << std::endl;
  propagator->setEnergyList(fEnergyArray);

  if (fVerbose >= INFO) {std::cout << "Setup CUDAProb3Linear oscillation probability calculater" << std::endl;}
}
 
void OscProbCalcerCUDAProb3Linear::CalculateProbabilities(std::vector<FLOAT_T> OscParams) {
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
  FLOAT_T PathL   = OscParams[kPATHL];
  FLOAT_T Density = OscParams[kDENS];

  propagator->setNeutrinoMasses(dm12sq, dm23sq);
  propagator->setDensity(Density);
  propagator->setPathLength(PathL);

  // CUDAProb3Linear calculates oscillation probabilites for each NeutrinoType, so need to copy them from the calculator into fWeightArray
  int CopyArrSize = fNEnergyPoints;
  FLOAT_T* CopyArr = new FLOAT_T[CopyArrSize];

  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {

    cudaprob3::NeutrinoType NuType;
    if (fNeutrinoTypes[iNuType]==Nubar) {
      NuType = cudaprob3::Antineutrino;
    } else {
      NuType = cudaprob3::Neutrino;
    }
    propagator->setMNSMatrix(theta12, theta13, theta23, dcp, NuType);

    propagator->calculateProbabilities(NuType);

    for (int iInitFlav=0;iInitFlav<fNInitialFlavours;iInitFlav++) {
      for (int iFinalFlav=0;iFinalFlav<fNFinalFlavours;iFinalFlav++) {
	propagator->getProbabilityArr(CopyArr,static_cast<cudaprob3::ProbType>(OscChannels[iInitFlav][iFinalFlav]));
	
	// Mapping which links the oscillation channel, neutrino type and energy/cosineZ index to the fWeightArray index
	int IndexToFill = iNuType*fNInitialFlavours*fNFinalFlavours*CopyArrSize + iInitFlav*fNFinalFlavours*CopyArrSize + iFinalFlav*CopyArrSize;
	for (int iOscProb=0;iOscProb<CopyArrSize;iOscProb++) {

	  // Sometimes CUDAProb3Linear can return *slightly* unphysical oscillation probabilities
	  CopyArr[iOscProb] = CopyArr[iOscProb] > 0.0 ? CopyArr[iOscProb] : 0.0;
	  CopyArr[iOscProb] = CopyArr[iOscProb] < 1.0 ? CopyArr[iOscProb] : 1.0;

	  fWeightArray[IndexToFill+iOscProb] = CopyArr[iOscProb];
	}
      }
    }
  }

  delete[] CopyArr;
}

int OscProbCalcerCUDAProb3Linear::ReturnWeightArrayIndex(int NuTypeIndex, int InitNuIndex, int FinalNuIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNInitialFlavours*fNFinalFlavours*fNEnergyPoints + InitNuIndex*fNFinalFlavours*fNEnergyPoints + FinalNuIndex*fNEnergyPoints + EnergyIndex;

  return IndexToReturn;
}

long OscProbCalcerCUDAProb3Linear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNInitialFlavours * fNFinalFlavours * fNNeutrinoTypes;
  return nCalculationPoints;
}
