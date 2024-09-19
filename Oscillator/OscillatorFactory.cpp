#include "OscillatorFactory.h"

#include <iostream>

OscillatorFactory::OscillatorFactory() {
}

OscillatorFactory::~OscillatorFactory() {

}

OscillatorBase* OscillatorFactory::CreateOscillator(std::string ConfigName_) {
  // Create config manager
  std::cout << "OscillatorFactory creating OscillatorBase object from config: " << ConfigName_ << std::endl;
  YAML::Node Config = YAML::LoadFile(ConfigName_);

  std::string OscillatorType = Config["General"]["CalculationType"].as<std::string>();
  OscillatorBase* Oscillator = NULL;

  if (OscillatorType == "Unbinned") {
    OscillatorUnbinned* ImpOscillator = new OscillatorUnbinned(ConfigName_);
    Oscillator = (OscillatorBase*)ImpOscillator;
  } else if (OscillatorType == "Binned") {
    OscillatorBinned* ImpOscillator = new OscillatorBinned(ConfigName_);
    Oscillator = (OscillatorBase*)ImpOscillator;
  } else {
    std::cerr << "OscillatorFactory was provided with unknown calculation type:" << OscillatorType << std::endl;
    std::cerr << "Please fix any mistakes or implement the calculator type at:" << __LINE__ << " : " << __FILE__ << std::endl;
    throw;
  }

  if (Oscillator == NULL) {
    std::cerr << "Did not successfully setup the correct return variable 'Oscillator'" << std::endl;
    throw;
  }

  return Oscillator;
}
