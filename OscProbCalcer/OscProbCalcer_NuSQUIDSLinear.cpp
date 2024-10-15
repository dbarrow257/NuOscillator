#include "OscProbCalcer_NuSQUIDSLinear.h"

#include <iostream>

#include "nuSQuIDS/nuSQuIDS.h"
#include "examples/Decoherence/nuSQUIDSDecoh.h"

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

  squids::Const units;
  unsigned int numneu = 3;
  nusquids::nuSQUIDSDecoh nus(nusquids::logspace(1.e2*units.GeV,1.e6*units.GeV,200),numneu,nusquids::neutrino,true);

  //Here we define the trajectory that the particle follows and the object for more examples
  // of how construct a track and object look body_track example.
  //zenith angle, neutrinos crossing the earth
  double phi = acos(0.);
  //Declaration of the body, EarthAtm is one of the predefined bodies
  std::shared_ptr<nusquids::EarthAtm> earth_atm = std::make_shared<nusquids::EarthAtm>();
  //Definition of the track, in encodes the trajectory inside the body, here is declared with the zenith angle.
  auto track_atm = std::make_shared<nusquids::EarthAtm::Track>(earth_atm->MakeTrack(phi));
  //We set this in the nusSQuID object.
  nus.Set_Body(earth_atm);
  nus.Set_Track(track_atm);

  // set mixing angles and masses
  nus.Set_MixingAngle(0,1,0.563942);
  nus.Set_MixingAngle(0,2,0.154085);
  nus.Set_MixingAngle(1,2,0.785398);
  nus.Set_SquareMassDifference(1,7.65e-05);
  nus.Set_SquareMassDifference(2,0.00247);

  //Here we set the maximum size for the integration step, important for fast or sharp variations of the density.
  nus.Set_h_max( 500.0*units.km );

  //We set the GSL step function
  nus.Set_GSL_step(gsl_odeiv2_step_rk4);

  //Setting the numerical precision of gsl integrator.
  nus.Set_rel_error(1.0e-5);
  nus.Set_abs_error(1.0e-5);

    //Construct the initial state
  //E_range is an array that contains all the energies.
  nusquids::marray<double,1> E_range = nus.GetERange();
  //Array that contains the initial state of the system, fist component is energy and second every one of the flavors
  nusquids::marray<double,2> inistate{E_range.size(),numneu};
  double N0 = 1.0e18;
  //Se a power low spectra for the muon neutrinos (k==1), other flavors to 0.0
  for ( int i = 0 ; i < inistate.extent(0); i++){
      for ( int k = 0; k < inistate.extent(1); k ++){
        inistate[i][k] = (k == 1) ? N0*pow(E_range[i],-2) : 0.0;
      }
  }

  //Set the initial state in nuSQuIDS object
  nus.Set_initial_state(inistate,nusquids::flavor);

  //Set the decoherence model and parameters
  nus.Set_DecoherenceGammaMatrix(nusquids::nuSQUIDSDecoh::DecoherenceModel::RandomizeState, 9.48e-18*units.eV);
  nus.Set_DecoherenceGammaEnergyDependence(2);
  nus.Set_DecoherenceGammaEnergyScale(1.0*units.TeV);

  //Propagate the neutrinos in the earth for the path defined in path
  nus.EvolveState();

  std::cout << std::endl << "Writing the outputs..." << std::endl;

  //number of energies we want the result, notice that this can be larger than the number of the internal grid of 
  //the nuSQuIDS object, a linear interpolation between the quantum density matrices in the interaction picture is used
  //and vacuum oscillations are solved analytically for the given energy.
  int Nen =1000;
  double lEmin=2;
  double lEmax=6;
  
  std::cout << "# log10(E) E flux_i fluxRatio_i . . . ." << std::endl;
  for(double lE=lEmin; lE<lEmax; lE+=(lEmax-lEmin)/(double)Nen){
    double E=pow(10.0,lE)*units.GeV;
    std::cout << lE << " " << E << " ";
    for(int fl=0; fl<numneu; fl++){
      std::cout << " " <<  nus.EvalFlavor(fl, E) << " " <<  nus.EvalFlavor(fl, E)/(N0*pow(E,-2));
    }
    std::cout << std::endl;
  }
  
  throw;
}

OscProbCalcerNuSQUIDSLinear::~OscProbCalcerNuSQUIDSLinear() {
}

void OscProbCalcerNuSQUIDSLinear::SetupPropagator() {
}

void OscProbCalcerNuSQUIDSLinear::CalculateProbabilities(const std::vector<FLOAT_T>& OscParams) {
}

int OscProbCalcerNuSQUIDSLinear::ReturnWeightArrayIndex(int NuTypeIndex, int OscChanIndex, int EnergyIndex, int CosineZIndex) {
  int IndexToReturn = NuTypeIndex*fNOscillationChannels*fNEnergyPoints + OscChanIndex*fNEnergyPoints + EnergyIndex;
  return IndexToReturn;
}

long OscProbCalcerNuSQUIDSLinear::DefineWeightArraySize() {
  long nCalculationPoints = fNEnergyPoints * fNOscillationChannels * fNNeutrinoTypes;
  return nCalculationPoints;
}
