#ifndef BSG_OPTIONCONTAINER
#define BSG_OPTIONCONTAINER

#include <fstream>
#include <iostream>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

/**
 * Macro to easily get the options from the OptionContainer object
 *
 * @param a variable type
 * @param b variable name
 */
#define GetBSGOpt(a, b) bsg::BSGOptionContainer::GetInstance().GetBSGOption<a>(#b)
/**
 * Macro to check whether a certain option was present
 *
 * @param a variable name
 */
#define OptExists(a) bsg::BSGOptionContainer::GetInstance().Exists(#a)

namespace bsg {

namespace po = boost::program_options;

/**
 * Class that combines all options from commandline, configuration files and
 * environment variables.
 * Implemented as a Singleton
 */
class BSGOptionContainer {
 public:
  /**
   * Single Constructor
   */
  static BSGOptionContainer& GetInstance(int argc = 0, char** argv = NULL) {
    static BSGOptionContainer instance(argc, argv);
    return instance;
  }
  /**
   * Get the option from the container
   *
   * @template T variable type
   * @param name vriable name
   */
  template <typename T>
  T GetBSGOption(std::string name) {
    return vm[name].as<T>();
  }
  /**
   * Check whether an options was given
   *
   * @param name variable name
   */
  bool Exists(std::string name) { return (bool)vm.count(name); }
  inline static po::options_description GetGenericOptions() {
    return genericOptions;
  };
  inline static po::options_description GetSpectrumOptions() {
    return spectrumOptions;
  };
  inline static po::options_description GetConfigOptions() {
    return configOptions;
  };
  inline static po::options_description GetTransitionOptions() {
    return transitionOptions;
  };

 private:
  static po::variables_map vm;
  static po::options_description genericOptions;
  static po::options_description spectrumOptions;
  static po::options_description configOptions;
  static po::options_description transitionOptions;
  BSGOptionContainer(int, char**);
  BSGOptionContainer(BSGOptionContainer const& copy);
  BSGOptionContainer& operator=(BSGOptionContainer const& copy);
};

}

#endif
