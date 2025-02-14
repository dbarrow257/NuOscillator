#include "OscProbCalcer_NuSQUIDSLinear.h"

#include <iostream>

OscProbCalcerNuSQUIDSLinear::OscProbCalcerNuSQUIDSLinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config

  //IntegrationStep
  if (!Config_["OscProbCalcerSetup"]["IntegrationStep"]) {
    std::cerr << "Expected to find a 'IntegrationStep' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }
  integration_step = Config_["OscProbCalcerSetup"]["IntegrationStep"].as<double>();

  //NusRelativeError
  if (!Config_["OscProbCalcerSetup"]["RelativeError"]) {
    std::cerr << "Expected to find a 'RelativeError' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }
  rel_error = Config_["OscProbCalcerSetup"]["RelativeError"].as<double>();

  //NusAbsoluteError
  if (!Config_["OscProbCalcerSetup"]["AbsoluteError"]) {
    std::cerr << "Expected to find a 'AbsoluteError' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }
  abs_error = Config_["OscProbCalcerSetup"]["AbsoluteError"].as<double>();

  //OscModel
  if (!Config_["OscProbCalcerSetup"]["OscModel"]) {
    std::cerr << "Expected to find a 'OscModel' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }
  std::string osc_model = Config_["OscProbCalcerSetup"]["OscModel"].as<std::string>();
  fOscType = PMNS_StrToInt(osc_model);
  
  if (fOscType == kDecoherence) {
    if (!Config_["OscProbCalcerSetup"]["DecoherenceModel"]) {
      std::cerr << "Expected to find a 'DecoherenceModel' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
      throw;
    }
    decoherence_model = Config_["OscProbCalcerSetup"]["DecoherenceModel"].as<std::string>();
  }
  //=======

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  fNOscParams = GetNOscParams(fOscType);

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);
}

OscProbCalcerNuSQUIDSLinear::~OscProbCalcerNuSQUIDSLinear() {
}

void OscProbCalcerNuSQUIDSLinear::SetupPropagator() {

  nusquids::marray<double,1> E_range{fNEnergyPoints};
  for(int i=0; i < fNEnergyPoints; i++){
    E_range[i] = fEnergyArray[i]*units.GeV;
  }

  switch (fOscType) {
  case kDecoherence:
    nus_decoh = new nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau, nusquids::neutrino, true);
    nubars_decoh = new nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau, nusquids::antineutrino, false); // anti-neutrinos
    
    nus_base = nus_decoh;
    nubars_base = nubars_decoh;
    break;
  default:
    std::cerr << "Unknown fOscType:" << fOscType << std::endl;
    throw;
  }
  
  //Set integration step
  nus_base->Set_h_max(integration_step*units.km );
  nubars_base->Set_h_max( integration_step*units.km );

  //We set the GSL step function
  nus_base->Set_GSL_step(gsl_odeiv2_step_rk4);
  nubars_base->Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus_base->Set_rel_error(rel_error);
  nus_base->Set_abs_error(abs_error);
  nubars_base->Set_rel_error(rel_error);
  nubars_base->Set_abs_error(abs_error);

  switch (fOscType) {
  case kDecoherence:
    if (decoherence_model == "RandomizePhase") {
      nusquids_decoherence_model = nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizePhase;
    } else if (decoherence_model == "RandomizeState") {
      nusquids_decoherence_model = nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState;
    } else if (decoherence_model == "NeutrinoLoss") {
      nusquids_decoherence_model = nusquids::nuSQUIDSDecoh::DecoherenceModel::NeutrinoLoss;
    } else {
      std::cerr << "The decoherence_model requested is not implemented. Given:" << decoherence_model << std::endl;
      throw;
    }
    break;
  }
}

void OscProbCalcerNuSQUIDSLinear::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {

  // Set mixing angles and masses for neutrinos
  nus_base->Set_MixingAngle(0,1,asin(sqrt(OscParams[kTH12]))); // \theta_12
  nus_base->Set_MixingAngle(0,2,asin(sqrt(OscParams[kTH13]))); // \theta_13
  nus_base->Set_MixingAngle(1,2,asin(sqrt(OscParams[kTH23]))); // \theta_23
  nus_base->Set_SquareMassDifference(1,OscParams[kDM12]); // \Delta m_12
  nus_base->Set_SquareMassDifference(2,OscParams[kDM12] + OscParams[kDM23]); // \Delta m_13
  nus_base->Set_CPPhase(0,2,OscParams[kDCP]);

  // Set mixing angles and masses for anti-neutrinos
  nubars_base->Set_MixingAngle(0,1,asin(sqrt(OscParams[kTH12]))); // \theta_12
  nubars_base->Set_MixingAngle(0,2,asin(sqrt(OscParams[kTH13]))); // \theta_13
  nubars_base->Set_MixingAngle(1,2,asin(sqrt(OscParams[kTH23]))); // \theta_23
  nubars_base->Set_SquareMassDifference(1,OscParams[kDM12]); // \Delta m_12
  nubars_base->Set_SquareMassDifference(2,OscParams[kDM12] + OscParams[kDM23]); // \Delta m_13
  nubars_base->Set_CPPhase(0,2,OscParams[kDCP]);

  const double layer_2 = OscParams[kPATHL]*units.km;
  std::shared_ptr<nusquids::ConstantDensity> constdens_env1 = std::make_shared<nusquids::ConstantDensity>(OscParams[kDENS],OscParams[kELECDENS]); // density [gr/cm^3[, ye [dimensionless]
  std::shared_ptr<nusquids::ConstantDensity::Track> track_env1 = std::make_shared<nusquids::ConstantDensity::Track>(layer_2);

  // Set energy density for neutrinos
  nus_base->Set_Body(constdens_env1);
  nus_base->Set_Track(track_env1);

  // Set energy density for anti-neutrinos
  nubars_base->Set_Body(constdens_env1);
  nubars_base->Set_Track(track_env1);

  // Construct the initial state
  // E_range is an array that contains all the energies.
  nusquids::marray<double,1> E_range = nus_base->GetERange();
  // Array that contains the initial state of the system, fist component is energy and second every one of the flavors
  nusquids::marray<double,2> inistate{E_range.size(),static_cast<size_t>(NuOscillator::kTau)};
  
  switch (fOscType) {
  case kDecoherence:
    //Set the decoherence model and parameters for neutrinos
    nus_decoh->Set_DecoherenceGammaMatrix(nusquids_decoherence_model, OscParams[kEnergyStrength]*units.eV);
    nus_decoh->Set_DecoherenceGammaEnergyDependence(OscParams[kEnergyDep]);
    nus_decoh->Set_DecoherenceGammaEnergyScale(OscParams[kEnergyScale]*units.GeV);
    
    //Set the decoherence model and parameters for anti-neutrinos
    nubars_decoh->Set_DecoherenceGammaMatrix(nusquids_decoherence_model, OscParams[kEnergyStrength]*units.eV);
    nubars_decoh->Set_DecoherenceGammaEnergyDependence(OscParams[kEnergyDep]);
    nubars_decoh->Set_DecoherenceGammaEnergyScale(OscParams[kEnergyScale]*units.GeV);
    break;
  }
 
  // Index counter to have a handle on where neutrino oscillation probs are stored in array fWeightArray
  int index_counter = 0;

  // Loop over all neutrino flavors, set the initial state, propagate the neutrinos and store osc probs in fWeightarray in
  // order nu_e->nu_e,nu_e->nu_mu, nu_e->nu_tau,
  //       nu_mu->nu_e, nu_mu->nu_mu, nu_mu->nu_tau,
  //       nu_tau->nu_e, nu_tau->nu_mu, nu_tau->nu_tau,
  for(int nu_flavor = 0; nu_flavor < NuOscillator::kTau; nu_flavor++){
      // Set initial state for the electron neutrinos (k==0), muon neutrinos (k==1) and tau neutrinos (k==2), other flavors to 0.0
      for ( int i = 0 ; i < inistate.extent(0); i++){
          for ( int k = 0; k < inistate.extent(1); k ++){
            inistate[i][k] = (k == nu_flavor) ? 1.0 : 0.0;
          }
      }

      //Set the initial state in nuSQuIDS object
      nus_base->Set_initial_state(inistate,nusquids::flavor);

      //Propagate the neutrinos in the earth for the path defined in path
      nus_base->EvolveState();

      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nus_base->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
          index_counter++;
        }
      }
  }

  // Now the same for anti-neutrinos:
  for(int nu_flavor = 0; nu_flavor < NuOscillator::kTau; nu_flavor++){
      // Set initial state for the electron neutrinos (k==0), muon neutrinos (k==1) and tau neutrinos (k==2), other flavors to 0.0
      for ( int i = 0 ; i < inistate.extent(0); i++){
          for ( int k = 0; k < inistate.extent(1); k ++){
            inistate[i][k] = (k == nu_flavor) ? 1.0 : 0.0;
          }
      }

      //Set the initial state in nuSQuIDS object
      nubars_base->Set_initial_state(inistate,nusquids::flavor);

      //Propagate the neutrinos in the earth for the path defined in path
      nubars_base->EvolveState();

      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nubars_base->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
          index_counter++;
        }
      }
  }
}

int OscProbCalcerNuSQUIDSLinear::PMNS_StrToInt(std::string OscModel) {
  if (OscModel=="Decoherence") {
    return kDecoherence;
  }

  std::cerr << "Unknown OscModel:" << OscModel << std::endl;
  throw;
  
  return -1;
}

int OscProbCalcerNuSQUIDSLinear::GetNOscParams(int OscModel) {
  switch (OscModel) {
  case kDecoherence:
    return kNOscParams;
  default:
    std::cerr << "Unknown OscModel:" << OscModel << std::endl;
    throw;
  }

  return -1;
}

int OscProbCalcerNuSQUIDSLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerNuSQUIDSLinear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
