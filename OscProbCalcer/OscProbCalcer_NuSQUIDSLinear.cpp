#include "OscProbCalcer_NuSQUIDSLinear.h"

#include <iostream>

OscProbCalcerNuSQUIDSLinear::OscProbCalcerNuSQUIDSLinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config

  //=======

  fNOscParams = kNOscParams;

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  // This implementation only considers linear propagation, thus no requirement to set cosineZ array
  IgnoreCosineZBinning(true);
}

OscProbCalcerNuSQUIDSLinear::~OscProbCalcerNuSQUIDSLinear() {
}

void OscProbCalcerNuSQUIDSLinear::SetupPropagator() {

  nusquids::marray<double,1> E_range{fNEnergyPoints};

  std::cout<< "fNEnergyPoints: " << fNEnergyPoints << std::endl;

  for(int i=0; i < fNEnergyPoints; i++){
    E_range[i] = fEnergyArray[i]*units.GeV;
  }

  nus = nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau,nusquids::neutrino,true); // neutrinos
  nubars = nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau,nusquids::antineutrino,false); // anti-neutrinos

  //Here we define the trajectory that the particle follows and the object for more examples
  // of how construct a track and object look body_track example.
  //zenith angle, neutrinos crossing the earth
  double phi = acos(0.);

  //Declaration of the body, EarthAtm is one of the predefined bodies
//  std::shared_ptr<nusquids::EarthAtm> earth_atm = std::make_shared<nusquids::EarthAtm>();
  //Definition of the track, in encodes the trajectory inside the body, here is declared with the zenith angle.
//  auto track_atm = std::make_shared<nusquids::EarthAtm::Track>(earth_atm->MakeTrack(phi));
  //We set this in the nusSQuID object.
//  nus.Set_Body(earth_atm);
//  nus.Set_Track(track_atm);
//  nubars.Set_Body(earth_atm);
//  nubars.Set_Track(track_atm);

  //Here we set the maximum size for the integration step, important for fast or sharp variations of the density.
  nus.Set_h_max( 500.0*units.km );
  nubars.Set_h_max( 500.0*units.km );

  //We set the GSL step function
  nus.Set_GSL_step(gsl_odeiv2_step_rk4);
  nubars.Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus.Set_rel_error(1.0e-5);
  nus.Set_abs_error(1.0e-5);
  nubars.Set_rel_error(1.0e-5);
  nubars.Set_abs_error(1.0e-5);

  //Set the decoherence model and parameters
//  nus.Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, 9.48e-18*units.eV);
  nus.Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, 0.0*units.eV);
  nus.Set_DecoherenceGammaEnergyDependence(2);
  nus.Set_DecoherenceGammaEnergyScale(1.0*units.TeV);
//  nubars.Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, 9.48e-18*units.eV);
  nubars.Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, 0.0*units.eV);
  nubars.Set_DecoherenceGammaEnergyDependence(2);
  nubars.Set_DecoherenceGammaEnergyScale(1.0*units.TeV);
}

void OscProbCalcerNuSQUIDSLinear::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {

//  std::cout<< "fNEnergyPoints: " << fNEnergyPoints << std::endl;

  // Set mixing angles and masses
//  nus.Set_MixingAngle(0,1,asin(sqrt(OscParams[kTH12]))); // \theta_12
//  nus.Set_MixingAngle(0,2,asin(sqrt(OscParams[kTH23]))); // \theta_13
//  nus.Set_MixingAngle(1,2,asin(sqrt(OscParams[kTH13]))); // \theta_13
  nus.Set_SquareMassDifference(1,OscParams[kDM12]); // \Delta m_12
  nus.Set_SquareMassDifference(2,OscParams[kDM23]); // \Delta m_13

  nus.Set_MixingParametersToDefault();

//  nubars.Set_MixingAngle(0,1,asin(sqrt(OscParams[kTH12]))); // \theta_12
//  nubars.Set_MixingAngle(0,2,asin(sqrt(OscParams[kTH23]))); // \theta_13
//  nubars.Set_MixingAngle(1,2,asin(sqrt(OscParams[kTH13]))); // \theta_13
  nubars.Set_SquareMassDifference(1,OscParams[kDM12]); // \Delta m_12
  nubars.Set_SquareMassDifference(2,OscParams[kDM23]); // \Delta m_13

  nubars.Set_MixingParametersToDefault();

  const double layer_2 = OscParams[kPATHL]*units.km;
  std::shared_ptr<nusquids::ConstantDensity> constdens_env1 = std::make_shared<nusquids::ConstantDensity>(OscParams[kDENS],OscParams[kELECDENS]); // density [gr/cm^3[, ye [dimensionless]
  std::shared_ptr<nusquids::ConstantDensity::Track> track_env1 = std::make_shared<nusquids::ConstantDensity::Track>(layer_2);

  nus.Set_Body(constdens_env1);
  nus.Set_Track(track_env1);

  nubars.Set_Body(constdens_env1);
  nubars.Set_Track(track_env1);

  //Construct the initial state
  //E_range is an array that contains all the energies.
  nusquids::marray<double,1> E_range = nus.GetERange();
  //Array that contains the initial state of the system, fist component is energy and second every one of the flavors
  nusquids::marray<double,2> inistate{E_range.size(),static_cast<size_t>(NuOscillator::kTau)};

//  for(int i = 0; i < fNWeights; i++){
//    fWeightArray[i] = 0.0;
//  }

  int index_counter = 0 * fNEnergyPoints;
  for(int nu_flavor = 0; nu_flavor < NuOscillator::kTau; nu_flavor++){
  // Set initial state for the electron neutrinos (k==0), muon neutrinos (k==1) and tau neutrinos (k==2), other flavors to 0.0
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        inistate[i][k] = (k == nu_flavor) ? 1.0 : 0.0;
      }
  }

  //Set the initial state in nuSQuIDS object
  nus.Set_initial_state(inistate,nusquids::flavor);

  //Propagate the neutrinos in the earth for the path defined in path
  nus.EvolveState();

  // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
  //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
  //and vacuum oscillations are solved analytically for the given energy.

  for(int fl=0; fl<NuOscillator::kTau; fl++){
    for(int i = 0; i < fNEnergyPoints; i++) {
      fWeightArray[index_counter] = nus.EvalFlavor(fl, fEnergyArray[i]*units.GeV);
      index_counter++;
    }
  }
  }

  // Now the same for anti-nus:
  for(int nu_flavor = 0; nu_flavor < NuOscillator::kTau; nu_flavor++){
  // Set initial state for the electron neutrinos (k==0), muon neutrinos (k==1) and tau neutrinos (k==2), other flavors to 0.0
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        inistate[i][k] = (k == nu_flavor) ? 1.0 : 0.0;
      }
  }

  //Set the initial state in nuSQuIDS object
  nubars.Set_initial_state(inistate,nusquids::flavor);

  //Propagate the neutrinos in the earth for the path defined in path
  nubars.EvolveState();

  // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
  //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
  //and vacuum oscillations are solved analytically for the given energy.

  for(int fl=0; fl<NuOscillator::kTau; fl++){
    for(int i = 0; i < fNEnergyPoints; i++) {
      fWeightArray[index_counter] = nubars.EvalFlavor(fl, fEnergyArray[i]*units.GeV);
      index_counter++;
    }
  }
  }
}

int OscProbCalcerNuSQUIDSLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerNuSQUIDSLinear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
