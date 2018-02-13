#include "NMEOptions.h"
#include <iostream>

using std::cout;
using std::endl;

po::options_description NMEOptions::genericOptions("Generic options");
po::options_description NMEOptions::configOptions("Configuration file options");
po::options_description NMEOptions::transitionOptions("Transition information");
po::options_description NMEOptions::envOptions("Environment options");
po::variables_map NMEOptions::vm;

NMEOptions::NMEOptions(int argc, char** argv) {
  transitionOptions.add_options()("Transition.Process",
                                  po::value<std::string>(),
                                  "Set the decay process: B+, B-")(
      "Transition.DensityMatrixFile", po::value<std::string>(),
      "Set the file name containing the one-body density matrix elements.")(
      "Daughter.Z", po::value<int>(),
      "Set the proton number of the daughter nucleus.")(
      "Mother.Z", po::value<int>(),
      "Set the proton number of the mother nucleus.")(
      "Daughter.A", po::value<int>(),
      "Set the total number of nucleons of the daughter nucleus.")(
      "Mother.A", po::value<int>(),
      "Set the total number of nucleons of the mother nucleus.")(
      "Daughter.Radius", po::value<double>(),
      "Set the nuclear radius of the daughter nucleus in fm.")(
      "Mother.Radius", po::value<double>(),
      "Set the nucleus radius of the mother nucleus in fm.")(
      "Mother.Beta2", po::value<double>()->default_value(0.),
      "Set the quadrupole deformation beta2 parameter for the mother nucleus.")(
      "Mother.Beta4", po::value<double>()->default_value(0.),
      "Set the hexadecupole deformation beta4 parameter for the mother "
      "nucleus.")("Mother.Beta6", po::value<double>()->default_value(0.),
                  "Set the hexadecupole deformation beta6 parameter for the "
                  "mother nucleus.")("Daughter.Beta2",
                                     po::value<double>()->default_value(0.),
                                     "Set the quadrupole deformation beta2 "
                                     "parameter for the daughter nucleus.")(
      "Daughter.Beta4", po::value<double>()->default_value(0.),
      "Set the hexadecupole deformation beta4 parameter for the daughter "
      "nucleus.")("Daughter.Beta6", po::value<double>()->default_value(0.),
                  "Set the hexadecupole deformation beta6 parameter for the "
                  "daughter nucleus.")(
      "Mother.SpinParity", po::value<int>(),
      "Set the spin times 2 and parity of the mother nucleus: [+/-]2Ji")(
      "Daughter.SpinParity", po::value<int>(),
      "Set the spin times 2 and parity of the daughter nucleus: [+/-]2Jf")(
      "Mother.Isospin", po::value<int>(),
      "Set the isospin times 2 and parity of the mother nucleus: [+/-]2Ti")(
      "Daughter.Isospin", po::value<int>(),
      "Set the isospin times 2 and parity of the daughter nucleus: [+/-]2Tf")(
      "Mother.ExcitationEnergy", po::value<double>()->default_value(0.),
      "Set the excitation energy of the mother nucleus in MeV")(
      "Daughter.ExcitationEnergy", po::value<double>()->default_value(0.),
      "Set the excitation energy of the daughter nucleus in MeV")(
      "Mother.ForcedSPSpin", po::value<int>()->default_value(0),
      "Set the spin of the transforming nucleon when Computational.ForceSpin "
      "is turned on.")("Daughter.ForcedSPSpin",
                       po::value<int>()->default_value(0),
                       "Set the spin of the transforming nucleon when "
                       "Computational.ForceSpin is turned on.");

  std::string configName = "config.txt";
  std::string inputName = "test.ini";
  genericOptions.add_options()("help,h", "Produce help message")(
      "config,c", po::value<std::string>(&configName),
      "Change the configuration file.")(
      "input,i", po::value<std::string>(&inputName),
      "Specify input file containing transition and nuclear data")(
      "output,o", po::value<std::string>()->default_value("output"),
      "Specify the output file name.")(
      "weakmagnetism,b", "Calculate the weak magnetism form factor b/Ac")(
      "inducedtensor,d", "Calculate the induced tensor form factor d/Ac")(
      "matrixelement,M", po::value<std::string>(),
      "Calculate the matrix element ^XM_{yyy} written as Xyyy");

  configOptions.add_options()(
      "Computational.Method", po::value<std::string>()->default_value("ESP"),
      "Set the method to use when calculating matrix elements. "
      "Defaults to Extreme Single Particle (ESP).")(
      "Computational.Potential", po::value<std::string>()->default_value("SHO"),
      "Set the potential used for the calculation of the matrix elements. SHO:"
      " Spherical Harmonic Oscillator; WS: Woods-Saxon; DWS: Deformed "
      "Woods-Saxon")(
      "Computational.ForceSpin", po::value<bool>()->default_value(true),
      "Choose the spin state with the given spin within EnergyMargin")(
      "Computational.EnergyMargin", po::value<double>()->default_value(0.5),
      "Set the energy margin in choosing the correct spin state in the ESP "
      "method.")("Computational.ReversedGhallagher",
                 po::value<bool>()->default_value(false),
                 "Reverse the Ghallagher coupling rules")(
      "Computational.OverrideSPCoupling",
      po::value<bool>()->default_value(false),
      "Override the single particle coupling rules and set the coupled spin to "
      "the final state.")(
      "Computational.SurfaceThickness",
      po::value<double>()->default_value(0.650),
      "Surface thickness of the Woods-Saxon potential in fm.")(
      "Computational.Vneutron", po::value<double>()->default_value(49.6),
      "Set the depth of the Woods-Saxon potential for neutrons in MeV.")(
      "Computational.Vproton", po::value<double>()->default_value(49.6),
      "Set the depth of the Woods-Saxon potential for protons in MeV.")(
      "Computational.Xneutron", po::value<double>()->default_value(0.86),
      "Set the neutron asymmetry parameter in the Woods-Saxon potential")(
      "Computational.Xproton", po::value<double>()->default_value(0.86),
      "Set the proton asymmetry parameter in the Woods-Saxon potential")(
      "Computational.V0Sneutron", po::value<double>()->default_value(7.2),
      "Set the magnitude of the spin-orbit potential for neutrons in MeV.")(
      "Computational.V0Sproton", po::value<double>()->default_value(7.2),
      "Set the magnitude of the spin-orbit potential for protons in MeV.")(
      "Constants.gA", po::value<double>()->default_value(1.2723),
      "Set the weak coupling constant.")(
      "Constants.gAeff", po::value<double>()->default_value(1.1),
      "Set the effective value for gA in nuclei.")(
      "Constants.gP", po::value<double>()->default_value(0.),
      "Set the induced pseudoscalar coupling constant.")(
      "Constants.gM", po::value<double>()->default_value(4.706),
      "Set the weak magnetism coupling constant.");

  po::store(po::command_line_parser(argc, argv)
                .options(genericOptions)
                .allow_unregistered()
                .run(),
            vm);
  // po::store(po::parse_command_line(argc, argv, genericOptions), vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << transitionOptions << endl;
    cout << "\n\n**************************************************************"
            "*\n\n" << endl;
    cout << genericOptions << endl;
    cout << configOptions << endl;
  } else {
    std::ifstream configStream(configName.c_str());
    if (!configStream.is_open()) {
      std::cerr << "ERROR: " << configName << " cannot be found.\n\n"
                << std::endl;
    } else {
      po::store(po::parse_config_file(configStream, configOptions, true), vm);
    }
    std::ifstream inputStream(inputName.c_str());
    if (!inputStream.is_open()) {
      std::cerr << "ERROR: " << inputName << " cannot be found.\n\n"
                << std::endl;
    } else {
      po::store(po::parse_config_file(inputStream, transitionOptions, true),
                vm);
    }
    po::notify(vm);
  }
}
