#include "OscProbCalcer_NuFASTEarth.h"

#include "NuFastEarth.h"
#include "Matrix.h"

OscProbCalcerNuFASTEarth::OscProbCalcerNuFASTEarth(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config

  // detector depth in km
  DetectorDepth = Config_["OscProbCalcerSetup"]["DetectorDepth"].as<FLOAT_T>();

  //Number of Newton-Raphson iteration, set to negative for exact eigenvalues
  EigenValuePrecision = Config_["OscProbCalcerSetup"]["EigenValuePrecision"].as<int>(); 
  
  //Earth model type picked from "Prob3", "PREM4", "Full", "NUniformLayers" (see https://github.com/PeterDenton/NuFast-Earth/blob/main/src/Earth.cpp)
  EarthModel = Config_["OscProbCalcerSetup"]["EarthModel"].as<std::string>();

  if (EarthModel == "NUniformLayers") {
    //Number of layers of Earth used in the NuFAST-Earth::PREM_NUniformLayer object - Not used unless EarthModel=="NUniformLayers"
    NUniformLayers = Config_["OscProbCalcerSetup"]["NUniformLayers"].as<int>();
  }
  //=======
  
  fNOscParams = kNOscParams;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;
  
  // This implementation only considers atmopsheric propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(false);
}

OscProbCalcerNuFASTEarth::~OscProbCalcerNuFASTEarth() {
  if (EarthDensity) {delete EarthDensity;}
}

void OscProbCalcerNuFASTEarth::SetupPropagator() {
  //=============================
  //Set the Earth Model
  if (EarthModel == "Prob3") {
    EarthDensity = new PREM_Prob3();
  } else if (EarthModel == "PREM4") {
    EarthDensity = new PREM_Four();
  } else if (EarthModel == "Full") {
    EarthDensity = new PREM_Full();
  } else if (EarthModel == "NUniformLayers") {
    EarthDensity = new PREM_NUniformLayer(NUniformLayers);
  } else {
    std::cerr << "Did not find a valid EarthModel to build - given:" << EarthModel << std::endl;
    throw std::runtime_error("Invalid Earth Model in NuFASTEarth");
  }
  //=============================
  ProbEngine= new Probability_Engine();
  ProbEngine->Set_Earth(DetectorDepth, EarthDensity);
  ProbEngine->Set_Eigenvalue_Precision(EigenValuePrecision);
}

void OscProbCalcerNuFASTEarth::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
 
  const double s12sq = OscParams[kTH12];
  const double s13sq = OscParams[kTH13];
  const double s23sq = OscParams[kTH23];
  const double delta = OscParams[kDCP];
  const double Dmsq21 = OscParams[kDM12];

  //Need to convert OscParams[kDM23] to kDM31
  const double Dmsq31 = OscParams[kDM23]+OscParams[kDM12]; // eV^2

  const double ProductionHeight = OscParams[kPROD]; //km

  for (int iNuType=0;iNuType<fNNeutrinoTypes;iNuType++) {

    bool NeutrinoType = (fNeutrinoTypes[iNuType] == Nu) ? true : false;
    ProbEngine->Set_Oscillation_Parameters(s12sq, s13sq, s23sq, delta, Dmsq21, Dmsq31, NeutrinoType);
    ProbEngine->Set_Production_Height(ProductionHeight);
    ProbEngine->Set_Spectra(fEnergyArray, fCosineZArray);
    std::vector<std::vector<Matrix3r>> probabilities = ProbEngine->Get_Probabilities();
    for (int iOscChannel = 0; iOscChannel < fNOscillationChannels; iOscChannel++) {
      for (int iEnergyPoint=0;iEnergyPoint<fNEnergyPoints;iEnergyPoint++) {
        for (int iCosineZPoint=0;iCosineZPoint<fNCosineZPoints;iCosineZPoint++) {
          int Index = ReturnWeightArrayIndex(iNuType,iOscChannel,iEnergyPoint,iCosineZPoint);
          const int gflv = fOscillationChannels[iOscChannel].GeneratedFlavour-1;
          const int dflv = fOscillationChannels[iOscChannel].DetectedFlavour-1;

          fWeightArray[Index] = probabilities[iEnergyPoint][iCosineZPoint].arr[gflv][dflv];
        }
      } 
    }
  }
}

int OscProbCalcerNuFASTEarth::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNCosineZPoints*fNEnergyPoints + OscChanIndex*fNCosineZPoints*fNEnergyPoints + EnergyIndex*fNCosineZPoints + CosineZIndex;
  return IndexToReturn;
}

long OscProbCalcerNuFASTEarth::DefineWeightArraySize() {
  long nCalculationPoints = static_cast<long>(fNEnergyPoints) * fNCosineZPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
