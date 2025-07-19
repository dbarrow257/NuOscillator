#include "OscProbCalcer_nuTens.h"

#include <iostream>

using namespace nuTens;

void OscProbCalcernuTens::SetEnergyArray(std::vector<FLOAT_T> EnergyArray) {

    OscProbCalcerBase::SetEnergyArray(EnergyArray);

    energiesTensor = Tensor::zeros({(long int)EnergyArray.size(), 1}, dtypes::kFloat, dtypes::kCPU, false);

    for (int i = 0; i < EnergyArray.size(); i++) {

        energiesTensor.setValue({i, 0}, EnergyArray[i] * units::GeV);

    }

    tensorPropagator.setEnergies(energiesTensor);
}


OscProbCalcernuTens::OscProbCalcernuTens(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  NT_PROFILE_BEGINSESSION("OscProbCalcernuTens");

  NT_PROFILE();

  //Initialise tensors

  _mat1 = Tensor::zeros({1, 3, 3}, dtypes::kComplexFloat, dtypes::kCPU).requiresGrad(false);
  _mat2 = Tensor::zeros({1, 3, 3}, dtypes::kComplexFloat, dtypes::kCPU).requiresGrad(false);
  _mat3 = Tensor::zeros({1, 3, 3}, dtypes::kComplexFloat, dtypes::kCPU).requiresGrad(false);

  // silly workaround for silly problem where setMatterSolver assumes these values are set
  // and if they aren't, errors errors errors
  Tensor pmns = getPMNSmatrix(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
  tensorPropagator.setMixingMatrix(pmns);
  tensorPropagator.setMasses(_masses);
  // I will fix this in nuTens soon :)
  
  tensorPropagator.setMatterSolver(matterSolver);
  
  //=======
  fNOscParams = kNOscParams_PMNS;

  fNNeutrinoTypes = 1;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  //fNeutrinoTypes[1] = Nubar;
  #if UseMultithreading == 1
  fImplementationName += "-CPU-"+std::to_string(omp_get_max_threads());
  #else
  fImplementationName += "-CPU-"+std::to_string(1);
  #endif
  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);
}

OscProbCalcernuTens::~OscProbCalcernuTens() {

  NT_PROFILE_ENDSESSION();

}

void OscProbCalcernuTens::SetupPropagator() {
}

void OscProbCalcernuTens::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
  // ------------------------------- //
  // Set the experimental parameters //
  // ------------------------------- //
  const double L = OscParams[kPATHL]; // km
  const double rho = OscParams[kDENS]; // g/cc
  const double Ye = OscParams[kELECDENS];

  tensorPropagator.setBaseline(L * units::km);
  matterSolver->setDensity(rho * Ye);

  // ------------------------------------------ //
  // Calculate all 9 oscillations probabilities //
  // ------------------------------------------ //

  Tensor pmnsMatrix = getPMNSmatrix(
    /*theta12=*/OscParams[kTH12], 
    /*theta23=*/OscParams[kTH23], 
    /*theta13=*/OscParams[kTH13], 
    /*dm12=*/OscParams[kDM12], 
    /*dm23=*/OscParams[kDM23], 
    /*dcp=*/OscParams[kDCP]
  );

  tensorPropagator.setMixingMatrix(pmnsMatrix);
  tensorPropagator.setMasses(_masses);

  Tensor probs = tensorPropagator.calculateProbs();
  AccessedTensor<float, 3, dtypes::kCPU> accessedProbs = AccessedTensor<float, 3, dtypes::kCPU>(probs);

  NT_PROFILE("writing probs");
  
  for (int iOscProb=0;iOscProb<fNEnergyPoints;iOscProb++) {
    
    for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {
      
      for ( int iOscChan = 0; iOscChan < fNOscillationChannels ; iOscChan++ ) {

        // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
        const int IndexToFill = iNuType * fNOscillationChannels * fNEnergyPoints + iOscChan * fNEnergyPoints;

        double Weight = accessedProbs.getValue(
          iOscProb, 
          fOscillationChannels[iOscChan].DetectedFlavour - 1,
          fOscillationChannels[iOscChan].GeneratedFlavour - 1
        );

        if(Weight>1.0) {
          Weight = 0.999999 * Weight;
        }
        
        fWeightArray[IndexToFill+iOscProb] = Weight;
        
      }
    }
  }
}

int OscProbCalcernuTens::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcernuTens::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
