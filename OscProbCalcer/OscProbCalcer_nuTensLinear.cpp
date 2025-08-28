#include "OscProbCalcer_nuTensLinear.h"

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
  
  fNOscParams = kNOscParams_PMNS;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;
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

  tensorPropagator.setParameters(_theta12Tensor, _theta23Tensor, _theta13Tensor, _dcpTensor, _dm21Tensor, _dm31Tensor);

}

void OscProbCalcernuTens::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
  // ------------------------------- //
  // Set the experimental parameters //
  // ------------------------------- //

  tensorPropagator.setBaseline(OscParams[kPATHL] * units::km);
  tensorPropagator.setDensity(OscParams[kDENS] * OscParams[kELECDENS]);

  // ------------------------------------------ //
  // Calculate all 9 oscillations probabilities //
  // ------------------------------------------ //

  setParamValues(std::asin(std::sqrt(OscParams[kTH12])), std::asin(std::sqrt(OscParams[kTH13])), std::asin(std::sqrt(OscParams[kTH23])), OscParams[kDCP], OscParams[kDM12], OscParams[kDM23]);

  Tensor probs;
  
  NT_PROFILE("writing probs");

  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {

    if(fNeutrinoTypes[iNuType] == Nu) {

      tensorPropagator.setAntiNeutrino(false);
      probs = tensorPropagator.calculateProbs();
    }
    else if(fNeutrinoTypes[iNuType] == Nubar) {
      tensorPropagator.setAntiNeutrino(true);
      probs = tensorPropagator.calculateProbs();
    }

    // get AccessedTensor for probabilities to make reading from them faster
    auto accessedProbs = AccessedTensor<float, 3, dtypes::kCPU>(probs);

    for (int iOscProb=0;iOscProb<fNEnergyPoints;iOscProb++) {
  
      for ( int iOscChan = 0; iOscChan < fNOscillationChannels ; iOscChan++ ) {

        // Mapping which links the oscillation channel, neutrino type and energy index to the fWeightArray index
        const int IndexToFill = iNuType * fNOscillationChannels * fNEnergyPoints + iOscChan * fNEnergyPoints;

        double Weight = accessedProbs.getValue(
          iOscProb, 
          fOscillationChannels[iOscChan].GeneratedFlavour - 1,
          fOscillationChannels[iOscChan].DetectedFlavour - 1
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
