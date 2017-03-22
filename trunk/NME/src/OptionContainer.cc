#include "OptionContainer.h"
#include <iostream>

using std::cout;
using std::endl;

po::options_description OptionContainer::genericOptions("Generic options");
po::options_description OptionContainer::configOptions("Configuration file options");
po::options_description OptionContainer::transitionOptions("Transition information");
po::options_description OptionContainer::envOptions("Environment options");
po::variables_map OptionContainer::vm;

OptionContainer::OptionContainer(int argc, char** argv) {
  transitionOptions.add_options()
      ("Transition.Process", po::value<std::string>(),
       "Set the decay process: B+, B-")
      ("NuclearProperties.DaughterZ", po::value<int>(),
       "Set the proton number of the daughter nucleus.")
      ("NuclearProperties.DaughterA", po::value<int>(),
       "Set the total number of nucleons of the daughter nucleus.")
      ("NuclearProperties.DaughterRadius", po::value<double>(),
       "Set the nuclear radius of the daughter nucleus in fm.")
      ("NuclearProperties.MotherBeta2", po::value<double>()->default_value(0.),
       "Set the quadrupole deformation beta2 parameter for the mother nucleus.")
      ("NuclearProperties.MotherBeta4", po::value<double>()->default_value(0.),
       "Set the hexadecupole deformation beta4 parameter for the mother nucleus.")
      ("NuclearProperties.DaughterBeta2", po::value<double>()->default_value(0.),
       "Set the quadrupole deformation beta2 parameter for the daughter nucleus.")
      ("NuclearProperties.DaughterBeta4", po::value<double>()->default_value(0.),
       "Set the hexadecupole deformation beta4 parameter for the daughter nucleus.")
      ("NuclearProperties.MotherSpinParity", po::value<int>(),
       "Set the spin times 2 and parity of the mother nucleus: [+/-]2Ji")
      ("NuclearProperties.DaughterSpinParity", po::value<int>(),
       "Set the spin times 2 and parity of the daughter nucleus: [+/-]2Jf")
      ("NuclearProperties.MotherIsospin", po::value<int>(),
       "Set the isospin times 2 and parity of the mother nucleus: [+/-]2Ti")
      ("NuclearProperties.DaughterIsospin", po::value<int>(),
       "Set the isospin times 2 and parity of the daughter nucleus: [+/-]2Tf")
      ("NuclearProperties.MotherExcitationEnergy", po::value<double>()->default_value(0.),
       "Set the excitation energy of the mother nucleus in MeV")
      ("NuclearProperties.DaughterExcitationEnergy", po::value<double>()->default_value(0.),
       "Set the excitation energy of the daughter nucleus in MeV");

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
      ("weakmagnetism,b", "Calculate the weak magnetism form factor b/Ac")
      ("inducedtensor,d", "Calculate the induced tensor form factor d/Ac")
      ("matrixelement,M", po::value<std::string>(),
       "Calculate the matrix element ^XM_{yyy} written as Xyyy");

  configOptions.add_options()
      ("Computational.Method", po::value<std::string>()->default_value("ESP"),
       "Set the method to use when calculating matrix elements. "
       "Defaults to Extreme Single Particle (ESP).")
      ("Computational.Potential", po::value<std::string>()->default_value("SHO"),
       "Set the potential used for the calculation of the matrix elements. SHO:"
       " Spherical Harmonic Oscillator; WS: Woods-Saxon; DWS: Deformed Woods-Saxon")
      ("Computational.SurfaceThickness", po::value<double>()->default_value(0.650),
       "Surface thickness of the Woods-Saxon potential in fm.")
      ("Computational.V0neutron", po::value<double>()->default_value(53.0),
       "Set the depth of the Woods-Saxon potential for neutrons in MeV.")
      ("Computational.V0proton", po::value<double>()->default_value(48.0),
       "Set the depth of the Woods-Saxon potential for protons in MeV.")
      ("Computational.V0Sneutron", po::value<double>()->default_value(7.2),
       "Set the magnitude of the spin-orbit potential for neutrons in MeV.")
      ("Computational.V0Sproton", po::value<double>()->default_value(7.2),
       "Set the magnitude of the spin-orbit potential for protons in MeV.")
      ("Constants.gA", po::value<double>()->default_value(1.2723),
       "Set the weak coupling constant.")
      ("Constants.gP", po::value<double>()->default_value(0.),
       "Set the induced pseudoscalar coupling constant.")
      ("Constants.gM", po::value<double>()->default_value(4.706),
       "Set the weak magnetism coupling constant.");

  po::store(po::parse_command_line(argc, argv, genericOptions), vm);
  po::notify(vm);

  std::ifstream configStream(configName.c_str());
  if (!configStream.is_open()) {
    std::cerr << "ERROR: " << configName << " cannot be found.\n\n"
              << std::endl;
  }
  else {
    po::store(po::parse_config_file(configStream, configOptions, true), vm);
  }
  std::ifstream inputStream(inputName.c_str());
  if (!inputStream.is_open()) {
    std::cerr << "ERROR: " << inputName << " cannot be found.\n\n"
              << std::endl;
  }
  else {
    po::store(po::parse_config_file(inputStream, transitionOptions, true), vm);
  }
  po::store(po::parse_environment(envOptions, "BSG_"), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << transitionOptions << endl;
    cout << "\n\n***************************************************************\n\n" << endl;
    cout << genericOptions << endl;
    cout << configOptions << endl;
  }
}
