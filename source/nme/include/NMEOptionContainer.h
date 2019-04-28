#ifndef NME_OPTIONCONTAINER
#define NME_OPTIONCONTAINER

#include <fstream>
#include <iostream>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;

/**
 * Macro to easily get the options from the OptionContainer object
 *
 * @param a variable type
 * @param b variable name
 */
#define GetNMEOpt(a, b) nme::NMEOptionContainer::GetInstance().GetNMEOption<a>(#b)
/**
 * Macro to check whether a certain option was present
 *
 * @param a variable name
 */
#define NMEOptExists(a) nme::NMEOptionContainer::GetInstance().Exists(#a)

namespace nme {

/**
 * Class that combines all options from commandline, configuration files and
 * environment variables.
 * Implemented as a Singleton
 */
class NMEOptionContainer {
 public:
  /**
   * Single Constructor
   */
  static NMEOptionContainer& GetInstance(int argc = 0, char** argv = NULL,
  bool parseCmdLine = true, std::string configName = "", std::string inputName = "") {
    static NMEOptionContainer instance(argc, argv, parseCmdLine, configName, inputName);
    return instance;
  }
  /**
   * Get the option from the container
   *
   * @template T variable type
   * @param name vriable name
   */
  template <typename T>
  T GetNMEOption(std::string name) {
    try {
      return vm[name].as<T>();
    } catch (boost::bad_any_cast& e) {
      std::cerr << "ERROR: Option \"" << name << "\" not defined. " << std::endl;
      throw e;
    }
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
  inline static po::options_description GetEnvOptions() { return envOptions; };

 private:
  static po::variables_map vm;
  static po::options_description genericOptions;
  static po::options_description spectrumOptions;
  static po::options_description configOptions;
  static po::options_description envOptions;
  static po::options_description transitionOptions;
  NMEOptionContainer(int, char**, bool, std::string, std::string);
  NMEOptionContainer(NMEOptionContainer const& copy);
  NMEOptionContainer& operator=(NMEOptionContainer const& copy);
};

}

#endif
