#include "Oscillator/OscillatorFactory.h"

#include "OscillatorBinned.h"
#include "OscillatorUnbinned.h"
#include "OscillatorSubSampling.h"


#include <iostream>

OscillatorFactory::OscillatorFactory() {
}

OscillatorFactory::~OscillatorFactory() {

}

OscillatorBase* OscillatorFactory::CreateOscillator(std::string ConfigName_) {
  std::cout << "OscillatorFactory creating OscillatorBase object from config: " << ConfigName_ << std::endl;
  YAML::Node Config = YAML::LoadFile(ConfigName_);
  return CreateOscillator(Config);
}

OscillatorBase* OscillatorFactory::CreateOscillator(YAML::Node Config) {
  // Create config manager
  std::string OscillatorType = Config["General"]["CalculationType"].as<std::string>();
  OscillatorBase* Oscillator = NULL;

  if (OscillatorType == "Unbinned") {
    OscillatorUnbinned* ImpOscillator = new OscillatorUnbinned(Config);
    Oscillator = (OscillatorBase*)ImpOscillator;
  } else if (OscillatorType == "Binned") {
    OscillatorBinned* ImpOscillator = new OscillatorBinned(Config);
    Oscillator = (OscillatorBase*)ImpOscillator;
  } else if (OscillatorType == "SubSampling") {
    OscillatorSubSampling* ImpOscillator = new OscillatorSubSampling(Config);
    Oscillator = (OscillatorBase*)ImpOscillator;
  } else {
    std::cerr << "OscillatorFactory was provided with unknown calculation type:" << OscillatorType << std::endl;
    std::cerr << "Please fix any mistakes or implement the calculator type at:" << __LINE__ << " : " << __FILE__ << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  if (Oscillator == NULL) {
    std::cerr << "Did not successfully setup the correct return variable 'Oscillator'" << std::endl;
    throw std::runtime_error("Invalid setup");
  }

  return Oscillator;
}
