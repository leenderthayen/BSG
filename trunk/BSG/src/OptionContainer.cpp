#include "OptionContainer.h"
#include "NMEOptions.h"
#include <iostream>

using std::cout;
using std::endl;

po::options_description OptionContainer::genericOptions("Generic options");
po::options_description OptionContainer::spectrumOptions(
    "Spectrum shape options");
po::options_description OptionContainer::configOptions(
    "Configuration file options");
po::options_description OptionContainer::transitionOptions(
    "Transition information");
po::options_description OptionContainer::envOptions("Environment options");
po::variables_map OptionContainer::vm;

OptionContainer::OptionContainer(int argc, char** argv) {
  transitionOptions.add_options()("Transition.Process",
                                  po::value<std::string>(),
                                  "Set the decay process: B+, B-")(
      "Transition.Type", po::value<std::string>(),
      "Set the decay type: fermi, gamow-teller or mixed")(
      "Transition.MixingRatio", po::value<double>()->default_value(0.),
      "Set the mixing ratio defined as g_A M_{GT} / M_{F}. Defaults to zero.")(
      "Transition.QValue", po::value<double>(),
      "Set the Q value of the decay in keV.")(
      "Transition.AtomicEnergyDeficit", po::value<double>()->default_value(0.),
      "Set the endpoint energy deficit due to atomic excitations. Turning this on will turn off the atomic mismatch correction")(
      "Transition.PartialHalflife", po::value<double>(),
      "Set the lifetime in seconds of the beta branch.")(
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
      "Set the excitation energy of the daughter nucleus in MeV");

  std::string configName = "config.txt";
  std::string inputName = "";
  genericOptions.add_options()("help,h", "Produce help message")
      ("config,c", po::value<std::string>(&configName),
       "Change the configuration file.")(
      "input,i", po::value<std::string>(&inputName),
      "Specify input file containing transition and nuclear data")(
      "output,o", po::value<std::string>()->default_value("output"),
      "Specify the output file name.");

  spectrumOptions.add_options()("fermi,f", "Turn off the Fermi Function.")(
      "phase,p", "Turn off the phase space factor. ")(
      "l0,l", "Turn off L0 correction.")("C,C", "Turn off C correction.")(
      "relativistic,R", "Turn off relativistic corrections.")(
      "deformation,D", "Turn off deformation corrections.")(
      "U,U", "Turn off U correction.")("Q,Q", "Turn off Q correction.")(
      "radiative,r", "Turn off radiative corrections.")(
      "recoil,n", "Turn off kinematic recoil corrections.")(
      "screening,s", "Turn off atomic screening.")("exchange,x",
                                                   "Turn off atomic exchange.")(
      "mismatch,m", "Turn off the atomic mismatch correction.")(
      "weakmagnetism,b", po::value<double>(),
      "Specify the b/Ac1 value in the Holstein formalism.")(
      "inducedtensor,d", po::value<double>(),
      "Specify the d/Ac1 value in the Holstein formalism.")(
      "ratioM121,X", po::value<double>(),
      "Specify the ratio of matrix elements M121/M101 in the "
      "Behrens-Buhring formalism.")("begin,B",
                                    po::value<double>()->default_value(0.1),
                                    "Specify the starting energy in keV from "
                                    "which to start the spectrum calculation.")(
      "end,E", po::value<double>()->default_value(0.),
      "Specify the last energy in keV for which to calculate the spectrum.")(
      "step,S", po::value<double>()->default_value(0.1),
      "Specify the stepsize in keV.")(
      "neutrino,v", "Turn off the generation of the neutrino spectrum.")(
      "connect", "Turn on the connection between BSG and NME for the calculation the C_I correction, thereby using the single particle states from the latter")(
      "shape", po::value<std::string>()->default_value("Fermi"),
      "Set the base shape name for the electrostatic finite size correction.")(
      "vnew", po::value<std::vector<double> >()->multitoken(),
      "Set the first three coefficients of the radial power expansion of the new nuclear electrostatic potential")(
      "vold", po::value<std::vector<double> >()->multitoken(),
      "Set the first three coefficients of the radial power expansion of the old nuclear electrostatic potential");

  configOptions.add(spectrumOptions);
  configOptions.add_options()(
      "General.Folder", po::value<std::string>()->default_value("."),
      "Set the folder name for the results to be placed in.")(
      "Constants.gA", po::value<double>(),
      "Specify the GT coupling constant, gA")(
      "Constants.gM", po::value<double>(),
      "Specify the weak magnetism coupling constant, gM")(
      "Constants.gP", po::value<double>(),
      "Specify the induced pseudoscalar coupling constant, gP");

  envOptions.add_options()(
      "ExchangeData",
      po::value<std::string>()->default_value("data/ExchangeData.dat"),
      "File location of the atomic exchange fit parameters.");

  /** Parse command line options
   * Included: generic options & spectrum shape options
   */
  po::options_description cmdOptions;
  cmdOptions.add(genericOptions).add(spectrumOptions);
  po::store(po::command_line_parser(argc, argv)
                .options(cmdOptions)
                .allow_unregistered()
                .run(),
            vm);
  po::notify(vm);

  if (vm.count("help")) {
    cout << transitionOptions << endl;
    cout << "\n\n**************************************************************"
            "*\n\n" << endl;
    cout << genericOptions << endl;
    cout << spectrumOptions << endl;
    cout << configOptions << endl;
    cout << envOptions << endl;
  } else {
    /** Parse configuration file
     * Included: configOptions & spectrumOptions
     */
    std::ifstream configStream(configName.c_str());
    if (!configStream.is_open()) {
      std::cerr << "ERROR: " << configName << " cannot be found.\n\n"
                << std::endl;
    } else {
      po::store(po::parse_config_file(configStream, configOptions, true), vm);
    }
    /** Parse .ini file
     */
    std::ifstream inputStream(inputName.c_str());
    if (!inputStream.is_open()) {
      std::cerr << "ERROR: " << inputName << " cannot be found.\n\n"
                << std::endl;
    } else {
      po::store(po::parse_config_file(inputStream, transitionOptions, true),
                vm);
    }
    /** Parse environment options
     */
    po::store(po::parse_environment(envOptions, "BSG_"), vm);
    po::notify(vm);
  }

  NMEOptions::GetInstance(argc, argv);
}
