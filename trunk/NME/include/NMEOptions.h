#ifndef NMEOPTIONS
#define NMEOPTIONS

#include <fstream>
#include <iostream>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

#define GetOpt(a, b) NMEOptions::GetInstance().GetOption< a >( #b )
#define OptExists(a) NMEOptions::GetInstance().Exists( #a )

namespace po = boost::program_options;

class NMEOptions {
 public:
  static NMEOptions& GetInstance(int argc = 0, char** argv = NULL) {
    static NMEOptions instance(argc, argv);
    return instance;
  }
  template <typename T>
  T GetOption(std::string name) {
    return vm[name].as<T>();
  }
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