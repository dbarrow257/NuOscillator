#ifndef __OSCPROBCALCERFACTORY_H__
#define __OSCPROBCALCERFACTORY_H__

#include "OscProbCalcerBase.h"

/**
 * @file OscProbCalcerFactory.h
 *
 * @class OscProbCalcerFactory
 *
 * @brief Factory method for creating instances of OscProbCalcerBase::OscProbCalcerBase() obejcts
 */
class OscProbCalcerFactory {
 public:

  /**
   * @brief Default constructor
   */
  OscProbCalcerFactory();

  /**
   * @brief Destructor
   */
  virtual ~OscProbCalcerFactory();

  /**
   * @brief Create an instance of OscProbCalcerBase::OscProbCalcerBase() objects from a YAML config. 
   *
   * @param ConfigName_ Path to YAML config file
   *
   * @return OscillatorBase::OscillatorBase() object
   */
  OscProbCalcerBase* CreateOscProbCalcer(std::string OscProbCalcerConfigName_);

  OscProbCalcerBase* CreateOscProbCalcer(YAML::Node OscProbCalcerConfig);

 protected:

 private:
};

#endif
