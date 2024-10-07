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
   * @param OscProbCalcerConfigName_ Path to YAML config file
   *
   * @return Intialised OscProbCalcerBase::OscProbCalcerBase() object
   */
  OscProbCalcerBase* CreateOscProbCalcer(std::string OscProbCalcerConfigName_);

  /**
   * @brief Create an instance of OscProbCalcerBase::OscProbCalcerBase() objects from a YAML Node
   *
   * @param OscProbCalcerConfig Instance of YAML node
   *
   * @return Intialised OscProbCalcerBase::OscProbCalcerBase() object
   */
  OscProbCalcerBase* CreateOscProbCalcer(YAML::Node OscProbCalcerConfig);

 protected:

 private:
};

#endif
