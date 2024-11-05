#ifndef __OSCILLATORFACTORY_H__
#define __OSCILLATORFACTORY_H__

#include "OscillatorBase.h"
#include "OscillatorBinned.h"
#include "OscillatorUnbinned.h"

#include "yaml-cpp/yaml.h"

/**
 * @file OscillatorFactory.h
 * @class OscillatorFactory
 * @brief Wrapper class which creates instances of OscillatorBase::OscillatorBase() objects
 */
class OscillatorFactory {
 public:

  /**
   * @brief Default constructor
   */
  OscillatorFactory();

  /**
   * @brief Destructor
   */
  virtual ~OscillatorFactory();

  /**
   * @brief Create an instance of OscillatorBase::OscillatorBase() objects from a YAML config. This currently includes OscillatorBinned::OscillatorBinned() and 
   * OscillatorUnbinned::OscillatorUnbinned() objects
   *
   * @param ConfigName_ Path to YAML config file
   *
   * @return OscillatorBase::OscillatorBase() object
   */
  OscillatorBase* CreateOscillator(std::string ConfigName_);

  /**
   * @brief Create an instance of OscillatorBase::OscillatorBase() objects from a YAML config. This currently includes OscillatorBinned::OscillatorBinned() and
   * OscillatorUnbinned::OscillatorUnbinned() objects
   *
   * @param Config YAML config node
   *
   * @return OscillatorBase::OscillatorBase() object
   */
  OscillatorBase* CreateOscillator(YAML::Node Config);

 protected:

 private:
};

#endif
