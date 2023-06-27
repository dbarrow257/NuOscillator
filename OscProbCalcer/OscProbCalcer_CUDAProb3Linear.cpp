#include "OscProbCalcer_CUDAProb3Linear.h"

#if UseGPU == 1
#include "beamcudapropagator.cuh"
#else
#include "beamcpupropagator.hpp"
#endif

#include <iostream>

using namespace cudaprob3;

OscProbCalcerCUDAProb3Linear::OscProbCalcerCUDAProb3Linear(std::string ConfigName_, int Instance_) : OscProbCalcerBase(ConfigName_,"CUDAProb3Linear",Instance_)
{
  //=======
  //Grab information from the config

  //=======

  fNOscParams = kNOscParams;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  // Implementation specific variables
  OscChannels = std::vector<int>(fNOscillationChannels,DUMMYVAL);
  for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {
    if (fOscillationChannels[iOscChannel].GeneratedFlavour == kElectron) {
      if (fOscillationChannels[iOscChannel].DetectedFlavour == kElectron) {OscChannels[iOscChannel] = e_e;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == kMuon) {OscChannels[iOscChannel] = e_m;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == kTau) {OscChannels[iOscChannel] = e_t;}
    } else if (fOscillationChannels[iOscChannel].GeneratedFlavour == kMuon) {
      if (fOscillationChannels[iOscChannel].DetectedFlavour == kElectron) {OscChannels[iOscChannel] = m_e;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == kMuon) {OscChannels[iOscChannel] = m_m;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == kTau) {OscChannels[iOscChannel] = m_t;}
    } else if (fOscillationChannels[iOscChannel].GeneratedFlavour == kTau) {
      if (fOscillationChannels[iOscChannel].DetectedFlavour == kElectron) {OscChannels[iOscChannel] = t_e;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == kMuon) {OscChannels[iOscChannel] = t_m;}
      if (fOscillationChannels[iOscChannel].DetectedFlavour == kTau) {OscChannels[iOscChannel] = t_t;}
    }
  }

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);

  nThreads = 1;
}

void OscProbCalcerCUDAProb3Linear::SetupPropagator() {

#if UseGPU == 1
  if (fVerbose >= INFO) {std::cout << "Using GPU CUDAProb3Linear propagator" << std::endl;}
  propagator = std::unique_ptr< Propagator<FLOAT_T> > ( new BeamCudaPropagatorSingle(0, fNEnergyPoints));
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

    for (int iOscChannel=0;iOscChannel<fNOscillationChannels;iOscChannel++) {
      propagator->getProbabilityArr(CopyArr,static_cast<cudaprob3::ProbType>(OscChannels[iOscChannel]));
	
      // Mapping which links the oscillation channel, neutrino type and energy/cosineZ index to the fWeightArray index
      int IndexToFill = iNuType*fNOscillationChannels*CopyArrSize + iOscChannel*CopyArrSize;
      for (int iOscProb=0;iOscProb<CopyArrSize;iOscProb++) {
	
	// Sometimes CUDAProb3Linear can return *slightly* unphysical oscillation probabilities
	CopyArr[iOscProb] = CopyArr[iOscProb] > 0.0 ? CopyArr[iOscProb] : 0.0;
	CopyArr[iOscProb] = CopyArr[iOscProb] < 1.0 ? CopyArr[iOscProb] : 1.0;
	
	fWeightArray[IndexToFill+iOscProb] = CopyArr[iOscProb];
      }
    }
  }

  delete[] CopyArr;
}

int OscProbCalcerCUDAProb3Linear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerCUDAProb3Linear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
