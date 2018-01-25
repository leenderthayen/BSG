#ifndef NMEOPTIONS
#define NMEOPTIONS

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
#define GetNMEOpt(a, b) NMEOptions::GetInstance().GetNMEOption<a>(#b)
/**
 * Macro to check whether a certain option was present
 * 
 * @param a variable name
 */
#define NMOptExists(a) NMEOptions::GetInstance().Exists(#a)

namespace po = boost::program_options;

/**
 * Class that combines all options from commandline, configuration files and 
 * environment variables.
 * Implemented as a Singleton
 */
class NMEOptions {
 public:
  /**
   * Single Constructor
   */
  static NMEOptions& GetInstance(int argc = 0, char** argv = NULL) {
    static NMEOptions instance(argc, argv);
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
    return vm[name].as<T>();
  }
  /**
   * Check whether an options was given
   * 
   * @param name variable name
   */
  bool Exists(std::string name) { return (bool)vm.count(name); }
  inline static po::options_description GetGenericOptions() { return genericOptions; };
  inline static po::options_description GetSpectrumOptions() { return spectrumOptions; };
  inline static po::options_description GetConfigOptions() {
    return configOptions;
  };
  inline static po::options_description GetTransitionOptions() { return transitionOptions; };
  inline static po::options_description GetEnvOptions() { return envOptions; };

 private:
  static po::variables_map vm;
  static po::options_description genericOptions;
  static po::options_description spectrumOptions;
  static po::options_description configOptions;
  static po::options_description envOptions;
  static po::options_description transitionOptions;
  NMEOptions(int, char**);
  NMEOptions(NMEOptions const& copy);
  NMEOptions& operator=(NMEOptions const& copy);
};

#endif
