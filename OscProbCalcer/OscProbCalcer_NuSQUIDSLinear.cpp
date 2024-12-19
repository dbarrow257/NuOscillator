#include "OscProbCalcer_NuSQUIDSLinear.h"

#include <iostream>

OscProbCalcerNuSQUIDSLinear::OscProbCalcerNuSQUIDSLinear(YAML::Node Config_) : OscProbCalcerBase(Config_)
{
  //=======
  //Grab information from the config
  //ZenithAngle
  if (!Config_["OscProbCalcerSetup"]["ZenithAngle"]) {
    std::cerr << "Expected to find a 'ZenithAngle' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  zenith_angle = Config_["OscProbCalcerSetup"]["ZenithAngle"].as<double>();

  //IntegrationStep
  if (!Config_["OscProbCalcerSetup"]["IntegrationStep"]) {
    std::cerr << "Expected to find a 'IntegrationStep' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  integration_step = Config_["OscProbCalcerSetup"]["IntegrationStep"].as<double>();

  //Errors
  //NusRelativeError
  if (!Config_["OscProbCalcerSetup"]["NusRelativeError"]) {
    std::cerr << "Expected to find a 'NusRelativeError' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

//  prem_model = Config_["OscProbCalcerSetup"]["NusRelativeError"].as<std::string>();
  nus_rel_error = Config_["OscProbCalcerSetup"]["NusRelativeError"].as<double>();

  //NusAbsoluteError
  if (!Config_["OscProbCalcerSetup"]["NusAbsoluteError"]) {
    std::cerr << "Expected to find a 'NusAbsoluteError' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nus_abs_error = Config_["OscProbCalcerSetup"]["NusAbsoluteError"].as<double>();

  //NubarsRelativeError
  if (!Config_["OscProbCalcerSetup"]["NubarsRelativeError"]) {
    std::cerr << "Expected to find a 'NubarsRelativeError' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nubars_rel_error = Config_["OscProbCalcerSetup"]["NubarsRelativeError"].as<double>();

  //NubarsAbsoluteError
  if (!Config_["OscProbCalcerSetup"]["NubarsAbsoluteError"]) {
    std::cerr << "Expected to find a 'NubarsAbsoluteError' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nubars_abs_error = Config_["OscProbCalcerSetup"]["NubarsAbsoluteError"].as<double>();

  //Decoherence setup
  //NusGammaStrength
  if (!Config_["OscProbCalcerSetup"]["NusGammaStrength"]) {
    std::cerr << "Expected to find a 'NusGammaStrength' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nus_gamma_strength = Config_["OscProbCalcerSetup"]["NusGammaStrength"].as<double>();

  //NusGammaEnergyDependence
  if (!Config_["OscProbCalcerSetup"]["NusGammaEnergyDependence"]) {
    std::cerr << "Expected to find a 'NusGammaEnergyDependence' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nus_gamma_energy_dependence = Config_["OscProbCalcerSetup"]["NusGammaEnergyDependence"].as<double>();

  //NusGammaEnergyScale
  if (!Config_["OscProbCalcerSetup"]["NusGammaEnergyScale"]) {
    std::cerr << "Expected to find a 'NusGammaEnergyScale' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nus_gamma_energy_scale = Config_["OscProbCalcerSetup"]["NusGammaEnergyScale"].as<double>();

  //Nubars
  //NubarsGammaStrength
  if (!Config_["OscProbCalcerSetup"]["NubarsGammaStrength"]) {
    std::cerr << "Expected to find a 'NubarsGammaStrength' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nubars_gamma_strength = Config_["OscProbCalcerSetup"]["NusGammaStrength"].as<double>();

  //NubarsGammaEnergyDependence
  if (!Config_["OscProbCalcerSetup"]["NubarsGammaEnergyDependence"]) {
    std::cerr << "Expected to find a 'NubarsGammaEnergyDependence' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nubars_gamma_energy_dependence = Config_["OscProbCalcerSetup"]["NubarsGammaEnergyDependence"].as<double>();

  //NubarsGammaEnergyScale
  if (!Config_["OscProbCalcerSetup"]["NubarsGammaEnergyScale"]) {
    std::cerr << "Expected to find a 'NubarsGammaEnergyScale' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nubars_gamma_energy_scale = Config_["OscProbCalcerSetup"]["NubarsGammaEnergyScale"].as<double>();

  //NusBSMModel
  if (!Config_["OscProbCalcerSetup"]["NusBSMModel"]) {
    std::cerr << "Expected to find a 'NusBSMModel' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  nus_bsm_model = Config_["OscProbCalcerSetup"]["NusBSMModel"].as<std::string>();

  if (!Config_["OscProbCalcerSetup"]["PMNSType"]) {
    std::cerr << "Expected to find a 'PMNSType' Node within the 'OscProbCalcerSetup''Implementation' Node" << std::endl;
    throw;
  }

  std::string OscMatrix = Config_["OscProbCalcerSetup"]["PMNSType"].as<std::string>();

  //=======

  fNNeutrinoTypes = 2;
  InitialiseNeutrinoTypesArray(fNNeutrinoTypes);
  fNeutrinoTypes[0] = Nu;
  fNeutrinoTypes[1] = Nubar;

  fOscType = PMNS_StrToInt(OscMatrix);
//  fNOscParams = kNOscParams;
  fNOscParams = GetNOscParams(fOscType);

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

  // Initialise neutrino (nusquids::nuSQUIDSDecoh nus) and anti-neutrino nusquids::nuSQUIDSDecoh nubars) object according to BSM model set in config file

//  std::string nus_bsm_model = "Decoherence"; // Set in config file. nusquids classes defined in nusquids: "examples/Decoherence/nuSQUIDSDecoh.h"

//  nusquids::nuSQUIDSDecoh nus;
//  nusquids::nuSQUIDSDecoh nubars;

//  std::unique_ptr<void> nus;


  // Select the model based on the string nus_bsm_model
//  if( nus_bsm_model == "SM" ){
//    nus = nusquids::nuSQUIDSAtm(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  }
//  else if( nus_bsm_model == "Decoherence" ){
  // nus = nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  nusMap["SM"] = std::make_shared<nusquids::nuSQUIDSAtm(E_range, NuOscillator::kTau,nusquids::neutrino,true);
  nusMap["Decoherence"] = std::make_shared<nusquids::nuSQUIDSDecoh>(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  nusMap["LV"] = std::make_shared<nusquids::nuSQUIDSLV>(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  nusMap["NSI"] = std::make_shared<nusquids::nuSQUIDSNSI>(E_range, NuOscillator::kTau,nusquids::neutrino,true);

//  nusMap[nus_bsm_model] = std::make_shared<nusquids::nuSQUIDSDecoh>(E_range, NuOscillator::kTau,nusquids::neutrino,true);
  
//  std::shared_ptr<nusquids::nuSQUIDS> nus;

// Function to handle casting automatically
//template<typename T>
//std::shared_ptr<T> get_nus_ptr(std::map<std::string, std::shared_ptr<nusquids::nuSQUIDS>>& nusMap, const std::string& key) {
//    auto nus = nusMap[key]; // Retrieve shared pointer from the map
//    if (!nus) {
//        std::cerr << "Error: Object not found in map!" << std::endl;
//        return nullptr;
//    }
//    return std::dynamic_pointer_cast<T>(nus); // Perform the cast to the derived type
//}

//  nus = get_nus_ptr<nusquids::nuSQUIDSDecoh>(nusMap, nus_bsm_model);

  //  auto nus = std::dynamic_pointer_cast<nusquids::nuSQUIDSDecoh>(nusMap[nus_bsm_model]);
  //auto nus = std::dynamic_pointer_cast<nusquids::nuSQUIDSDecoh>(nusMap[nus_bsm_model]);
//  } 
//  else if(nus_bsm_model == "LV"){
//    nus = nusquids::nuSQUIDSLV(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  }
//  else if( nus_bsm_model == "NSI" ){
//    nus = nusquids::nuSQUIDSNSI(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  }
//  else if( nus_bsm_model == "Atmospheric Decoherence" ){
//    nus = nusquids::nuSQUIDSDecohAtm(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  }
//  else{
//    std::cout << "Unknown model!" << std::endl;
//    throw;
//  }

//  std::map<std::string, std::shared_ptr<nusquids::nuSQUIDS>> nusMap;
//  nusMap[nus_bsm_model] = std::make_shared<nusquids::nuSQUIDSDecoh>(E_range, NuOscillator::kTau,nusquids::neutrino,true);

//auto nus = std::dynamic_pointer_cast<nusquids::nuSQUIDSDecoh>(nusMap[nus_bsm_model]);

//  std::shared_ptr<nusquids::nuSQUIDS> nus;

  auto nus = std::dynamic_pointer_cast<nusquids::nuSQUIDSDecoh>(nusMap[nus_bsm_model]);

// Check if the cast succeeded
        if (nus) {
            std::cout << "Cast successful, setting mixing angle." << std::endl;
        } else {
            std::cerr << "Error: Dynamic cast to nuSQUIDSDecoh failed." << std::endl;
        }

  nubars = nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau,nusquids::antineutrino,false); // anti-neutrinos

  //Here we define the trajectory that the particle follows and the object for more examples
  // of how construct a track and object look body_track example.
  //zenith angle, neutrinos crossing the earth
  double phi = acos(zenith_angle);

  //Here we set the maximum size for the integration step, important for fast or sharp variations of the density.
//  nus.Set_h_max( integration_step*units.km );
//  nusMap[nus_bsm_model].Set_h_max( integration_step*units.km );
  nus->Set_h_max( integration_step*units.km );
  nubars.Set_h_max( integration_step*units.km );

  //We set the GSL step function
  nus->Set_GSL_step(gsl_odeiv2_step_rk4);
  nubars.Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus->Set_rel_error(nus_rel_error);
  nus->Set_abs_error(nus_abs_error);
  nubars.Set_rel_error(nubars_rel_error);
  nubars.Set_abs_error(nubars_abs_error);
  std::cout << "A" << std::endl;
  //Set the decoherence model and parameters for neutrinos
  nus->Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, nus_gamma_strength*units.eV); // reference (default) value: nus_gamma_strength = 9.48e-18 
  std::cout << "B" << std::endl;
  nus->Set_DecoherenceGammaEnergyDependence(nus_gamma_energy_dependence); // reference (default) value: nus_gamma_energy_dependence = 2
  nus->Set_DecoherenceGammaEnergyScale(nus_gamma_energy_scale*units.GeV); // reference (default) value: nus_gamma_energy_dependence = 1.0
  //Set the decoherence model and parameters for anti-neutrinos
  nubars.Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, nubars_gamma_strength*units.eV);
  nubars.Set_DecoherenceGammaEnergyDependence(nubars_gamma_energy_dependence);
  nubars.Set_DecoherenceGammaEnergyScale(nubars_gamma_energy_scale*units.GeV);
  std::cout << "C" << std::endl;
}

void OscProbCalcerNuSQUIDSLinear::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {

//  std::cout<< "fNEnergyPoints: " << fNEnergyPoints << std::endl;
  std::cout << "B" << std::endl;
//      	auto nus = std::dynamic_pointer_cast<nusquids::nuSQUIDSDecoh>(nusMap[nus_bsm_model]);
  // Select the model based on the string nus_bsm_model
//  if( nus_bsm_model == "SM" ){
//    nus = nusquids::nuSQUIDSAtm(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  }
//  else if( nus_bsm_model == "Decoherence" ){
  // nus = nusquids::nuSQUIDSDecoh(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  nusMap[nus_bsm_model] = std::make_shared<nusquids::nuSQUIDSDecoh>(E_range, NuOscillator::kTau,nusquids::neutrino,true);
  auto nus = std::dynamic_pointer_cast<nusquids::nuSQUIDSDecoh>(nusMap[nus_bsm_model]);
//  }
//  else if(nus_bsm_model == "LV"){
//    nus = nusquids::nuSQUIDSLV(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  }
//  else if( nus_bsm_model == "NSI" ){
//    nus = nusquids::nuSQUIDSNSI(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  }
//  else if( nus_bsm_model == "Atmospheric Decoherence" ){
//    nus = nusquids::nuSQUIDSDecohAtm(E_range, NuOscillator::kTau,nusquids::neutrino,true);
//  }
//  else{
//    std::cout << "Unknown model!" << std::endl;
//    throw;
//  }
  std::cout << "C" << std::endl;
	// Check if the cast succeeded
        if (nus) {
            std::cout << "Cast successful, setting mixing angle." << std::endl;
        } else {
            std::cerr << "Error: Dynamic cast to nuSQUIDSDecoh failed." << std::endl;
        }
  // Set mixing angles and masses for neutrinos
  nus->Set_MixingAngle(0,1,asin(sqrt(OscParams[kTH12]))); // \theta_12
  nus->Set_MixingAngle(0,2,asin(sqrt(OscParams[kTH13]))); // \theta_13
  nus->Set_MixingAngle(1,2,asin(sqrt(OscParams[kTH23]))); // \theta_23
  nus->Set_SquareMassDifference(1,OscParams[kDM12]); // \Delta m_12
  nus->Set_SquareMassDifference(2,OscParams[kDM12] + OscParams[kDM23]); // \Delta m_13
  nus->Set_CPPhase(0,2,OscParams[kDCP]);
//  nus->Set_MixingParametersToDefault();

  // Set mixing angles and masses for anti-neutrinos
  nubars.Set_MixingAngle(0,1,asin(sqrt(OscParams[kTH12]))); // \theta_12
  nubars.Set_MixingAngle(0,2,asin(sqrt(OscParams[kTH13]))); // \theta_13
  nubars.Set_MixingAngle(1,2,asin(sqrt(OscParams[kTH23]))); // \theta_23
  nubars.Set_SquareMassDifference(1,OscParams[kDM12]); // \Delta m_12
  nubars.Set_SquareMassDifference(2,OscParams[kDM12] + OscParams[kDM23]); // \Delta m_13
  nubars.Set_CPPhase(0,2,OscParams[kDCP]);
//  nubars.Set_MixingParametersToDefault();

  // Declaration of the body, EarthAtm is one of the predefined bodies
//  std::shared_ptr<nusquids::EarthAtm> earth_atm = std::make_shared<nusquids::EarthAtm>();

  // Definition of the track, in encodes the trajectory inside the body, here is declared with the zenith angle.
//  auto track_atm = std::make_shared<nusquids::EarthAtm::Track>(earth_atm->MakeTrack(phi));
  // We set this in the nusSQuID object.

  const double layer_2 = OscParams[kPATHL]*units.km;
  std::shared_ptr<nusquids::ConstantDensity> constdens_env1 = std::make_shared<nusquids::ConstantDensity>(OscParams[kDENS],OscParams[kELECDENS]); // density [gr/cm^3[, ye [dimensionless]
  std::shared_ptr<nusquids::ConstantDensity::Track> track_env1 = std::make_shared<nusquids::ConstantDensity::Track>(layer_2);

  // Set energy density for neutrinos
//  nus.Set_Body(earth_atm);
//  nus.Set_Track(track_atm);
  nus->Set_Body(constdens_env1);
  nus->Set_Track(track_env1);

  // Set energy density for anti-neutrinos
//  nubars.Set_Body(earth_atm);
//  nubars.Set_Track(track_atm);
  nubars.Set_Body(constdens_env1);
  nubars.Set_Track(track_env1);

  // Construct the initial state
  // E_range is an array that contains all the energies.
  nusquids::marray<double,1> E_range = nus->GetERange();
  // Array that contains the initial state of the system, fist component is energy and second every one of the flavors
  nusquids::marray<double,2> inistate{E_range.size(),static_cast<size_t>(NuOscillator::kTau)};

//  for(int i = 0; i < fNWeights; i++){
//    fWeightArray[i] = 0.0;
//  }

  // Index counter to have a handle on where neutrino oscillation probs are stored in array fWeightArray
  int index_counter = 0 * fNEnergyPoints;

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
      nus->Set_initial_state(inistate,nusquids::flavor);

      //Propagate the neutrinos in the earth for the path defined in path
      nus->EvolveState();

      // Number of energies we want the result, notice that this can be larger than the number of the internal grid of 
      //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
      //and vacuum oscillations are solved analytically for the given energy.
      for(int fl=0; fl<NuOscillator::kTau; fl++){
        for(int i = 0; i < fNEnergyPoints; i++) {
          fWeightArray[index_counter] = nus->EvalFlavor(fl, fEnergyArray[i]*units.GeV);
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

int OscProbCalcerNuSQUIDSLinear::PMNS_StrToInt(std::string PMNSType) {
  if (PMNSType == "Fast" || PMNSType == "fast") {
    return kFast;
  } 
  else if (PMNSType == "SM") {
    return kPMNSSM;
  }
  else if (PMNSType == "Sterile+1" || PMNSType == "Sterile1") {
    return kPMNSSterile1;
  } 
  else if (PMNSType == "Sterile+2" || PMNSType == "Sterile2") {
    return kPMNSSterile2;
  }
  else if (PMNSType == "Sterile+3" || PMNSType == "Sterile3") {
    return kPMNSSterile3;
  }
  else if (PMNSType == "Decay" || PMNSType == "decay") {
    return kDecay;
  }
  else if (PMNSType == "Deco" || PMNSType == "deco") {
    return kDeco;
  }
  else if (PMNSType == "NSI" || PMNSType == "nsi") {
    return kNSI;
  }
  else if (PMNSType == "LIV" || PMNSType == "liv") {
    return kLIV;
  }
  else if (PMNSType == "NUNM" || PMNSType == "nunm") {
    return kNUNM;
  }
  else if (PMNSType == "SNSI" || PMNSType == "snsi") {
    return kSNSI;
  }
  else if (PMNSType == "Iter" || PMNSType == "iter") {
    return kIter;
  }
  else {
    std::cerr << "Invalid PMNS matrix type provided:" << PMNSType << std::endl;
    throw;
  }
  
  return -1;
}

int OscProbCalcerNuSQUIDSLinear::GetNOscParams(int OscType) {
  if (OscType == kFast) {
    return kNOscParams;
  }
  else if (OscType == kPMNSSM) {
    return kNOscParams;
  }
  else if (OscType == kPMNSSterile1) {
    return kNOscParams+6;
  }
  else if (OscType == kPMNSSterile2) {
    return kNOscParams+14;
  }
  else if (OscType == kPMNSSterile3) {
    return kNOscParams+24;
  }
  else if (OscType == kDecay) {
    return kNOscParams+2;
  }
  else if (OscType == kDeco) {
    return kNOscParams+4;
  }
  else if (OscType == kNSI) {
    return kNOscParams+12;
  }
  else if (OscType == kIter) {
    return kNOscParams+1;
  }
  else if (OscType == kNUNM) {
    return kNOscParams+10;
  }
  else if (OscType == kLIV) {
    std::cerr << "LIV PMNS matrix not implemented yet" << std::endl;
    throw;
  }
  else if (OscType == kSNSI) {
    std::cerr << "SNSI PMNS matrix not implemented yet" << std::endl;
    throw;
  }
  else {
    std::cerr << "Invalid PMNS matrix type provided:" << OscType << std::endl;
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
