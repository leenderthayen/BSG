#ifndef OPTIONCONTAINER
#define OPTIONCONTAINER

#include <fstream>
#include <iostream>
#include <string>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/parsers.hpp>
#include <boost/program_options/variables_map.hpp>

namespace po = boost::program_options;

class OptionContainer {
 public:
  static OptionContainer& GetInstance(int argc = 0, char** argv = NULL) {
    static OptionContainer instance(argc, argv);
    return instance;
  }
  template <typename T>
  T GetOption(std::string name) {
    return vm[name].as<T>();
  }
  bool Exists(std::string name) { return (bool)vm.count(name); }
  inline static po::options_description GetCmdOptions() { return cmdOptions; };
  inline static po::options_description GetConfigOptions() {
    return configOptions;
  };
  inline static po::options_description GetEnvOptions() { return envOptions; };

 private:
  static po::variables_map vm;
  static po::options_description cmdOptions;
  static po::options_description configOptions;
  static po::options_description envOptions;
  OptionContainer(int, char**);
  OptionContainer(OptionContainer const& copy);
  OptionContainer& operator=(OptionContainer const& copy);
};

#endif
