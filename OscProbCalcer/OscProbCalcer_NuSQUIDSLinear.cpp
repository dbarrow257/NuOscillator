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

  //Errors

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
  osc_model = Config_["OscProbCalcerSetup"]["OscModel"].as<std::string>();

  //DB ToDo: Querry enum rather than use string comparison
  if (osc_model == "Decoherence") {
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

  fOscType = PMNS_StrToInt(osc_model);
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

  nus_decoh = new nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau, nusquids::neutrino, true);
  nubars_decoh = new nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau, nusquids::antineutrino, false); // anti-neutrinos

  double eps_mutau = 1.0e-2;
  nus_NSI = new nuSQUIDSNSI(eps_mutau, E_range, NuOscillator::kTau, nusquids::neutrino, true);

  nus_LV = new nusquids::nuSQUIDSLV(E_range, NuOscillator::kTau, nusquids::neutrino, false);

  nus_pmns = new nusquids::nuSQUIDS(E_range, NuOscillator::kTau, nusquids::neutrino, false);

  nus_base = nus_decoh;

  nus_base->Set_h_max(integration_step*units.km );
  //  nus_base->Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, nus_gamma_strength*units.eV);

  //Here we set the maximum size for the integration step, important for fast or sharp variations of the density.
  nus_decoh->Set_h_max( integration_step*units.km );
  nubars_decoh->Set_h_max( integration_step*units.km );

  //We set the GSL step function
  nus_decoh->Set_GSL_step(gsl_odeiv2_step_rk4);
  nubars_decoh->Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus_decoh->Set_rel_error(rel_error);
  nus_decoh->Set_abs_error(abs_error);
  nubars_decoh->Set_rel_error(rel_error);
  nubars_decoh->Set_abs_error(abs_error);

  //Decoh specific stuff
  //DB ToDo: Querry enum rather than use string comparison
  if (osc_model == "Decoherence") {
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
  }
}

void OscProbCalcerNuSQUIDSLinear::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {

  // Set mixing angles and masses for neutrinos
  nus_decoh->Set_MixingAngle(0,1,asin(sqrt(OscParams[kTH12]))); // \theta_12
  nus_decoh->Set_MixingAngle(0,2,asin(sqrt(OscParams[kTH13]))); // \theta_13
  nus_decoh->Set_MixingAngle(1,2,asin(sqrt(OscParams[kTH23]))); // \theta_23
  nus_decoh->Set_SquareMassDifference(1,OscParams[kDM12]); // \Delta m_12
  nus_decoh->Set_SquareMassDifference(2,OscParams[kDM12] + OscParams[kDM23]); // \Delta m_13
  nus_decoh->Set_CPPhase(0,2,OscParams[kDCP]);

  // Set mixing angles and masses for anti-neutrinos
  nubars_decoh->Set_MixingAngle(0,1,asin(sqrt(OscParams[kTH12]))); // \theta_12
  nubars_decoh->Set_MixingAngle(0,2,asin(sqrt(OscParams[kTH13]))); // \theta_13
  nubars_decoh->Set_MixingAngle(1,2,asin(sqrt(OscParams[kTH23]))); // \theta_23
  nubars_decoh->Set_SquareMassDifference(1,OscParams[kDM12]); // \Delta m_12
  nubars_decoh->Set_SquareMassDifference(2,OscParams[kDM12] + OscParams[kDM23]); // \Delta m_13
  nubars_decoh->Set_CPPhase(0,2,OscParams[kDCP]);

  const double layer_2 = OscParams[kPATHL]*units.km;
  std::shared_ptr<nusquids::ConstantDensity> constdens_env1 = std::make_shared<nusquids::ConstantDensity>(OscParams[kDENS],OscParams[kELECDENS]); // density [gr/cm^3[, ye [dimensionless]
  std::shared_ptr<nusquids::ConstantDensity::Track> track_env1 = std::make_shared<nusquids::ConstantDensity::Track>(layer_2);

  // Set energy density for neutrinos
  nus_decoh->Set_Body(constdens_env1);
  nus_decoh->Set_Track(track_env1);

  // Set energy density for anti-neutrinos
  nubars_decoh->Set_Body(constdens_env1);
  nubars_decoh->Set_Track(track_env1);

  //Set the decoherence model and parameters for neutrinos
  nus_decoh->Set_DecoherenceGammaMatrix(nusquids_decoherence_model, OscParams[kEnergyStrength]*units.eV);
  nus_decoh->Set_DecoherenceGammaEnergyDependence(OscParams[kEnergyDep]);
  nus_decoh->Set_DecoherenceGammaEnergyScale(OscParams[kEnergyScale]*units.GeV);

  //Set the decoherence model and parameters for anti-neutrinos
  nubars_decoh->Set_DecoherenceGammaMatrix(nusquids_decoherence_model, OscParams[kEnergyStrength]*units.eV);
  nubars_decoh->Set_DecoherenceGammaEnergyDependence(OscParams[kEnergyDep]);
  nubars_decoh->Set_DecoherenceGammaEnergyScale(OscParams[kEnergyScale]*units.GeV);
  
  // Construct the initial state
  // E_range is an array that contains all the energies.
  nusquids::marray<double,1> E_range = nus_decoh->GetERange();
  // Array that contains the initial state of the system, fist component is energy and second every one of the flavors
  nusquids::marray<double,2> inistate{E_range.size(),static_cast<size_t>(NuOscillator::kTau)};

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
      nus_decoh->Set_initial_state(inistate,nusquids::flavor);

      //Propagate the neutrinos in the earth for the path defined in path
      nus_decoh->EvolveState();

      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nus_decoh->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
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
      nubars_decoh->Set_initial_state(inistate,nusquids::flavor);

      //Propagate the neutrinos in the earth for the path defined in path
      nubars_decoh->EvolveState();

      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nubars_decoh->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
          index_counter++;
        }
      }
  }
}

int OscProbCalcerNuSQUIDSLinear::PMNS_StrToInt(std::string PMNSType) {
  return -1;
}

int OscProbCalcerNuSQUIDSLinear::GetNOscParams(int OscType) {
  return kNOscParams;
}

int OscProbCalcerNuSQUIDSLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerNuSQUIDSLinear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
