#include "OscProbCalcer_CUDAProb3.h"

#include "constants.hpp"
#include "propagator.hpp"
#include "physics.hpp"

#include "cpupropagator.hpp"

#include <iostream>
using namespace cudaprob3;

OscProbCalcerCUDAProb3::OscProbCalcerCUDAProb3() : OscProbCalcerBase()
{
  //Base variables
  fNOscParams = kNOscParams;

  //Implementation specific variables
  nNeutrinoSigns = 2;
  nInitialFlavours = 2; // =2 if excluding taus, =3 if including taus
  nFinalFlavours = 3; // =2 if excluding taus, =3 if including taus

  NeutrinoTypes.resize(nNeutrinoSigns);
  NeutrinoTypes[0] = Neutrino;
  NeutrinoTypes[1] = Antineutrino;
  
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
  EarthDensityFile = std::string();
  //DB Should grab from some kind of config manager
  EarthDensityFile = "../CUDAProb3/models/PREM_4layer.dat";
}

void OscProbCalcerCUDAProb3::SetupPropagator() {

#ifdef UseGPU
  std::cout << "Using GPU CUDAProb3 propagaator" << std::endl;
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

  std::cout << "Using CPU CUDAProb3 propagator with " << nThreads << " threads" << std::endl;
  propagator = std::unique_ptr< Propagator< FLOAT_T > > ( new CpuPropagator<FLOAT_T>(fNCosineZPoints, fNEnergyPoints, nThreads)); // MultiThread CPU propagator
#endif

  propagator->setEnergyList(fEnergyArray);
  propagator->setCosineList(fCosineZArray);
  propagator->setDensityFromFile(EarthDensityFile);

  fPropagatorSet = true;
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

  for (int iNuType=0;iNuType<nNeutrinoSigns;iNuType++) {
    
    if (NeutrinoTypes[iNuType]==Antineutrino) {
      // Haven't really thought about it, but prob3++ sets dcp->-dcp here: https://github.com/rogerwendell/Prob3plusplus/blob/fd189e232e96e2c5ebb2f7bd3a5406b288228e41/BargerPropagator.cc#L235
      // Copying that behaviour gives same behaviour as prob3++/probGPU
      propagator->setMNSMatrix(theta12, theta13, theta23, -dcp);
    } else {
      propagator->setMNSMatrix(theta12, theta13, theta23, dcp);
    }

    propagator->calculateProbabilities(static_cast<cudaprob3::NeutrinoType>(NeutrinoTypes[iNuType]));

    for (int iInitFlav=0;iInitFlav<nInitialFlavours;iInitFlav++) {
      for (int iFinalFlav=0;iFinalFlav<nFinalFlavours;iFinalFlav++) {
	propagator->getProbabilityArr(CopyArr,static_cast<cudaprob3::ProbType>(OscChannels[iInitFlav][iFinalFlav]));
	
	// Mapping which links the oscillation channel, neutrino type and energy/cosineZ index to the fWeightArray index
	int IndexToFill = NeutrinoTypes[iNuType]*nInitialFlavours*nFinalFlavours*CopyArrSize + iInitFlav*nFinalFlavours*CopyArrSize + iFinalFlav*CopyArrSize;
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

const FLOAT_T* OscProbCalcerCUDAProb3::ReturnPointer(int InitNuFlav, int FinalNuFlav, FLOAT_T Energy, FLOAT_T CosineZ) {
  //DB Lets move all of this into the base class

  int NuType = (InitNuFlav > 0) - (InitNuFlav < 0); //Calculates the sign of InitNuFlav

  NeutrinoType NuIndex = Neutrino;
  if (NuType < 0) {
    NuIndex = Antineutrino;
  }

  int InitNuIndex = -1;
  if (fabs(InitNuFlav)==1) {
    InitNuIndex = 0; //electron
  } else if (fabs(InitNuFlav)==2) {
    InitNuIndex = 1; //muon
  } else if (fabs(InitNuFlav)==3) {
    InitNuIndex = 2; //tau
  } else {
    std::cerr << "Unknown initial neutrino flavour:" << InitNuFlav << std::endl;
    throw;
  }

  int FinalNuIndex = -1;
  if (fabs(FinalNuFlav)==1) {
    FinalNuIndex = 0; //electron
  } else if (fabs(FinalNuFlav)==2) {
    FinalNuIndex = 1; //muon
  } else if (fabs(FinalNuFlav)==3) {
    FinalNuIndex = 2; //tau
  } else {
    std::cerr << "Unknown final neutrino flavour:" << FinalNuFlav << std::endl;
    throw;
  }

  if (InitNuIndex >= nInitialFlavours) {
    std::cerr << "Requested Initial Neutrino flavour is outside expected range!" << std::endl;
    std::cerr << "InitNuFlav:" << InitNuFlav << std::endl;
    std::cerr << "InitNuIndex:" << InitNuIndex << std::endl;
    std::cerr << "nInitialFlavours:" << nInitialFlavours << std::endl;
    throw;
  }

  if (FinalNuIndex >= nFinalFlavours) {
    std::cerr << "Requested Final Neutrino flavour is outside expected range!" << std::endl;
    std::cerr << "FinalNuFlav:" << FinalNuFlav << std::endl;
    std::cerr << "FinalNuIndex:" << FinalNuIndex << std::endl;
    std::cerr << "nFinalialFlavours:" << nFinalFlavours << std::endl;
    throw;
  }

  int CosineZIndex = -1;
  for (size_t iCosineZ=0;iCosineZ<fCosineZArray.size();iCosineZ++) {
    if (CosineZ == fCosineZArray[iCosineZ]) {
      CosineZIndex = iCosineZ;
      break;
    }
  }
  if (CosineZIndex == -1) {
    std::cerr << "Did not find CosineZ in the array used in calculating oscillation probabilities" << std::endl;
    std::cerr << "Requested CosineZ:" << CosineZ << std::endl;
    throw;
  }

  int EnergyIndex = -1;
  for (size_t iEnergy=0;iEnergy<fEnergyArray.size();iEnergy++) {
    if (Energy == fEnergyArray[iEnergy]) {
      EnergyIndex = iEnergy;
      break;
    }
  }
  if (EnergyIndex == -1) {
    std::cerr << "Did not find Energy in the array used in calculating oscillation probabilities" << std::endl;
    std::cerr << "Requested Energy:" << Energy << std::endl;
    throw;
  }

  int IndexToReturn = NuIndex*nInitialFlavours*nFinalFlavours*fNCosineZPoints*fNEnergyPoints + InitNuIndex*nFinalFlavours*fNCosineZPoints*fNEnergyPoints + FinalNuIndex*fNCosineZPoints*fNEnergyPoints + CosineZIndex*fNEnergyPoints + EnergyIndex;

  return &(fWeightArray[IndexToReturn]);
}

void OscProbCalcerCUDAProb3::IntiailiseWeightArray() {
  int nCalculationPoints = fNEnergyPoints * fNCosineZPoints * nInitialFlavours * nFinalFlavours * nNeutrinoSigns;
  std::cout << "Creating weight array with " << nCalculationPoints << " entries" << std::endl;

  fWeightArray = std::vector<FLOAT_T>(nCalculationPoints,DUMMYVAL);
  fWeightArrayInit = true;
}
