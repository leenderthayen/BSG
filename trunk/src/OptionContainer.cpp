#include "OptionContainer.h"
#include <iostream>

using std::cout;
using std::endl;

po::options_description OptionContainer::genericOptions("Generic options");
po::options_description OptionContainer::spectrumOptions("Spectrum shape options");
po::options_description OptionContainer::configOptions("Configuration file options");
po::options_description OptionContainer::transitionOptions("Transition information");
po::options_description OptionContainer::envOptions("Environment options");
po::variables_map OptionContainer::vm;

OptionContainer::OptionContainer(int argc, char** argv) {
  transitionOptions.add_options()
      ("Transition.Process", po::value<std::string>(),
       "Set the decay process: B+, B-")
      ("Transition.Type", po::value<std::string>(),
       "Set the decay type: fermi, gamow-teller or mixed")
      ("Transition.MixingRatio", po::value<double>()->default_value(0.),
       "Set the mixing ratio defined as g_A M_{GT} / M_{F}. Defaults to zero.")
      ("Transition.QValue", po::value<double>(),
       "Set the Q value of the decay in keV.")
      ("NuclearProperties.DaughterZ", po::value<int>(),
       "Set the proton number of the daughter nucleus.")
      ("NuclearProperties.DaughterA", po::value<int>(),
       "Set the total number of nucleons of the daughter nucleus.")
      ("NuclearProperties.DaughterRadius", po::value<double>(),
       "Set the nuclear radius of the daughter nucleus in fm.")
      ("NuclearProperties.MotherBeta2", po::value<double>()->default_value(0.),
       "Set the quadrupole deformation beta2 parameter for the mother nucleus.")
      ("NuclearProperties.MotherBeta4", po::value<double>()->default_value(0.),
       "Set the octupole deformation beta4 parameter for the mother nucleus.")
      ("NuclearProperties.DaughterBeta2", po::value<double>()->default_value(0.),
       "Set the quadrupole deformation beta2 parameter for the daughter nucleus.")
      ("NuclearProperties.DaughterBeta4", po::value<double>()->default_value(0.),
       "Set the octupole deformation beta4 parameter for the daughter nucleus.")
      ("NuclearProperties.MotherSpinParity", po::value<int>(),
       "Set the spin times 2 and parity of the mother nucleus: [+/-]2Ji")
      ("NuclearProperties.DaughterSpinParity", po::value<int>(),
       "Set the spin times 2 and parity of the daughter nucleus: [+/-]2Jf");

  std::string configName = "config.txt";
  std::string inputName = "test.ini";
  genericOptions.add_options()
      ("help,h", "Produce help message")
      ("verbosity,v", po::value<int>()->default_value(1),
      "Set verbosity (0 = silent)")
      ("config,c", po::value<std::string>(&configName),
       "Change the configuration file.")
      ("input,i", po::value<std::string>(&inputName),
       "Specify input file containing transition and nuclear data")
      ("output,o", po::value<std::string>()->default_value("output.txt"),
       "Specify the output file name.");

  spectrumOptions.add_options()
      ("fermi,f", "Turn off the Fermi Function.")
      ("phase,p", "Turn off the phase space factor. ")
      ("l0,l", "Turn off L0 correction.")
      ("C,C", "Turn off C correction.")
      ("relativistic,R", "Turn off relativistic corrections.")
      ("deformation,D", "Turn off deformation corrections.")
      ("U,U", "Turn off U correction.")
      ("Q,Q", "Turn off Q correction.")
      ("radiative,r", "Turn off radiative corrections.")
      ("recoil,n", "Turn off kinematic recoil corrections.")
      ("screening,s", "Turn off atomic screening.")
      ("exchange,x", "Turn off atomic exchange.")
      ("mismatch,m", "Turn off the atomic mismatch correction.")
      ("weakmagnetism,b", po::value<double>()->default_value(0.),
       "Specify the b/Ac1 value in the Holstein formalism.")
      ("inducedtensor,d", po::value<double>()->default_value(0.),
       "Specify the d/Ac1 value in the Holstein formalism.")
      ("ratioM121,X", po::value<double>()->default_value(0.),
       "Specify the ratio of matrix elements M121/M101 in the "
       "Behrens-Buhring formalism.")
      ("begin,B", po::value<double>()->default_value(0.),
       "Specify the starting energy in keV from which to start the spectrum calculation.")
      ("end,E", po::value<double>()->default_value(0.),
       "Specify the last energy in keV for which to calculate the spectrum.")
      ("step,S", po::value<double>()->default_value(0.1),
       "Specify the stepsize in keV.")
      ("neutrino,v", "Turn off the generation of the neutrino spectrum.");

  configOptions.add(spectrumOptions);
  configOptions.add_options()
      ("General.Folder", po::value<std::string>()->default_value("."),
      "Set the folder name for the results to be placed in.")
      ("Computational.Potential", po::value<std::string>()->default_value("SHO"),
       "Set the potential used for the calculation of the matrix elements. SHO:"
       " Spherical Harmonic Oscillator; WS: Woods-Saxon; DWS: Deformed Woods-Saxon")
      ("Constants.gA", po::value<double>()->default_value(1.2723),
       "Set the weak coupling constant.")
      ("Constants.gP", po::value<double>()->default_value(0.),
       "Set the induced pseudoscalar coupling constant.")
      ("Constants.gM", po::value<double>()->default_value(4.706),
       "Set the weak magnetism coupling constant.");

  envOptions.add_options()
      ("ExchangeData",
       po::value<std::string>()->default_value("data/ExchangeData.dat"),
       "File location of the atomic exchange fit parameters.");

  po::options_description cmdOptions;
  cmdOptions.add(genericOptions).add(spectrumOptions);
  po::store(po::parse_command_line(argc, argv, cmdOptions), vm);
  po::notify(vm);

  std::ifstream configStream(configName.c_str());
  if (!configStream.is_open()) {
    std::cerr << "ERROR: " << configName << " cannot be found.\n\n"
              << std::endl;
  }
  else {
    po::store(po::parse_config_file(configStream, configOptions), vm);
  }
  std::ifstream inputStream(inputName.c_str());
  if (!inputStream.is_open()) {
    std::cerr << "ERROR: " << inputName << " cannot be found.\n\n"
              << std::endl;
  }
  else {
    po::store(po::parse_config_file(inputStream, transitionOptions), vm);
  }
  po::store(po::parse_environment(envOptions, "BDSM_"), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << transitionOptions << endl;
    cout << "\n\n***************************************************************\n\n" << endl;
    cout << genericOptions << endl;
    cout << spectrumOptions << endl;
    cout << configOptions << endl;
    cout << envOptions << endl;
  }
}
