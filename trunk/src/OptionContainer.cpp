#include "OptionContainer.hh"
#include <iostream>

using std::cout;
using std::endl;

po::options_description OptionContainer::genericOptions("Generic options");
po::options_description OptionContainer::spectrumOptions("Spectrum shape options");
po::options_description OptionContainer::configOptions("Configuration file options");
po::options_description OptionContainer::envOptions("Environment options");
po::variables_map OptionContainer::vm;

OptionContainer::OptionContainer(int argc, char** argv) {
  std::string configName = "config.txt";
  genericOptions.add_options()
      ("help,h", "Produce help message")
      ("verbosity,v", po::value<int>()->default_value(1),
      "Set verbosity (0 = silent)")
      ("config,c", po::value<std::string>(&configName),
       "Change the configuration file.")
      ("input,i", po::value<std::string>(),
       "Specify input file containing transition and nuclear data")
      ("output,o", po::value<std::string>(),
       "Specify the output file name.");

  spectrumOptions.add_options()
      ("fermi,f", "Turn off the Fermi Function.")
      ("phase,p", "Turn off the phase space factor. ")
      ("l0,l", "Turn off L0 correction.")
      ("C", "Turn off C correction.")
      ("relativistic,R", "Turn off relativistic corrections.")
      ("deformation,D", "Turn off deformation corrections.")
      ("U", "Turn off U correction.")
      ("Q", "Turn off Q correction.")
      ("radiative,r", "Turn off radiative corrections.")
      ("recoil,n", "Turn off kinematic recoil corrections.")
      ("screening,s", "Turn off atomic screening.")
      ("exchange,x", "Turn off atomic exchange.")
      ("mismatch,m", "Turn off the atomic mismatch correction.")
      ("weakmagnetism,b", po::value<double>()->default_value(0.),
       "Specify the b/Ac1 value in the Holstein formalism.")
      ("inducedtensor,d", po::value<double>()->default_value(0.),
       "Specify the d/Ac1 value in the Holstein formalism.")
      ("X", po::value<double>()->default_value(0.),
       "Specify the ratio of matrix elements M121/M101 in the "
       "Behrens-Buhring formalism.")
      ("step,S", po::value<double>()->default_value(0.1),
       "Specify the stepsize in keV.");

  configOptions.add(spectrumOptions);
  configOptions.add_options()
      ("General.Folder", po::value<std::string>()->default_value("."),
      "Set the folder name for the results to be placed in.");

  po::store(po::parse_command_line(argc, argv, genericOptions), vm);
  po::store(po::parse_command_line(argc, argv, spectrumOptions), vm);
  po::store(po::parse_config_file(configName.c_str(), configOptions), vm);

  envOptions.add_options()(
      "ExchangeData",
      po::value<std::string>()->default_value("data/ExchangeData.dat"),
      "File location of the atomic exchange fit parameters.");

  po::store(po::parse_environment(envOptions, "BDSM_"), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << genericOptions << endl;
    cout << spectrumOptions << endl;
    cout << configOptions << endl;
    cout << envOptions << endl;
  }
}
