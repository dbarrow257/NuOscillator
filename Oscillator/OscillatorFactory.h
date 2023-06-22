#ifndef __OSCILLATORFACTORY_H__
#define __OSCILLATORFACTORY_H__

#include "OscillatorBase.h"
#include "OscillatorBinned.h"
#include "OscillatorUnbinned.h"

#include "yaml-cpp/yaml.h"

class OscillatorFactory {
 public:
  OscillatorFactory();

  OscillatorBase* CreateOscillator(std::string ConfigName_);

 protected:

 private:
};

#endif
