#include "BSGOptionContainer.h"
#include "NMEOptionContainer.h"
#include "spdlog/spdlog.h"
#include <iostream>

using std::cout;
using std::endl;

po::options_description bsg::BSGOptionContainer::genericOptions("Generic options");
po::options_description bsg::BSGOptionContainer::spectrumOptions(
    "Spectrum shape options");
po::options_description bsg::BSGOptionContainer::configOptions(
    "Spectral configuration file options");
po::options_description bsg::BSGOptionContainer::transitionOptions(
    "Transition information");
po::variables_map bsg::BSGOptionContainer::vm;

bsg::BSGOptionContainer::BSGOptionContainer(int argc, char** argv) {
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
      "Set the endpoint energy deficit due to atomic excitations. Turning this "
      "on will turn off the atomic mismatch correction")(
      "Transition.PartialHalflife", po::value<double>(),
      "Set the lifetime in seconds of the beta branch.")(
      "Transition.LogFt", po::value<double>(),
      "Set the externally calculated log ft value.")(
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
  genericOptions.add_options()("help,h", "Produce help message")(
      "config,c", po::value<std::string>(&configName),
      "Change the configuration file.")(
      "exchangedata,e",
      po::value<std::string>()->default_value("ExchangeData.dat"),
      "Set the location of the atomic exchange parameters file.")(
      "input,i", po::value<std::string>(&inputName),
      "Specify input file containing transition and nuclear data")(
      "output,o", po::value<std::string>()->default_value("output"),
      "Specify the output file name.");

  spectrumOptions.add_options()("Spectrum.Fermi,f",
                                po::value<bool>()->default_value(true),
                                "Turn off the Fermi Function.")(
      "Spectrum.Phasespace,p", po::value<bool>()->default_value(true),
      "Turn off the phase space factor. ")(
      "Spectrum.ESFiniteSize,l", po::value<bool>()->default_value(true),
      "Turn off L0 correction.")(
      "Spectrum.C,C", po::value<bool>()->default_value(true),
      "Turn off C correction.")(
      "Spectrum.Isovector,I", po::value<bool>()->default_value(true),
      "Turn off the isovector correction to C")(
      "Spectrum.Relativistic,R", po::value<bool>()->default_value(true),
      "Turn off relativistic corrections.")(
      "Spectrum.ESDeformation,D", po::value<bool>()->default_value(true),
      "Turn off deformation corrections.")(
      "Spectrum.U,U", po::value<bool>()->default_value(true),
      "Turn off U correction.")(
      "Spectrum.CoulombRecoil,Q", po::value<bool>()->default_value(true),
      "Turn off Q correction.")(
      "Spectrum.Radiative,r", po::value<bool>()->default_value(true),
      "Turn off radiative corrections.")(
      "Spectrum.Recoil,n", po::value<bool>()->default_value(true),
      "Turn off kinematic recoil corrections.")(
      "Spectrum.Screening,s", po::value<bool>()->default_value(true),
      "Turn off atomic screening.")(
      "Spectrum.Exchange,x", po::value<bool>()->default_value(true),
      "Turn off atomic exchange.")(
      "Spectrum.ModGaussFit,M", po::value<double>(),
      "Use a specific 'A' value for the Modified Gaussian fit")(
      "Spectrum.AtomicMismatch,m", po::value<bool>()->default_value(true),
      "Turn off the atomic mismatch correction.")(
      "Spectrum.WeakMagnetism,b", po::value<double>(),
      "Specify the b/Ac1 value in the Holstein formalism.")(
      "Spectrum.InducedTensor,d", po::value<double>(),
      "Specify the d/Ac1 value in the Holstein formalism.")(
      "Spectrum.Lambda,L", po::value<double>(),
      "Specify the ratio of matrix elements M121/M101 in the "
      "Behrens-Buehring formalism.")(
      "Spectrum.Begin,B", po::value<double>()->default_value(0.1),
      "Specify the starting energy in keV from which to start the spectrum calculation.")(
      "Spectrum.End,E", po::value<double>()->default_value(0.),
      "Specify the last energy in keV for which to calculate the spectrum.")(
      "Spectrum.StepSize,S", po::value<double>()->default_value(0.1),
      "Specify the stepsize in keV.")(
      "Spectrum.Steps,N", po::value<int>(),
      "Specify the number of steps in the total spectrum")(
      "Spectrum.Neutrino,v", po::value<bool>()->default_value(true),
      "Turn off the generation of the neutrino spectrum.")(
      "Spectrum.Connect", po::value<bool>()->default_value(false),
      "Turn on the connection between BSG and NME for the calculation the C_I "
      "correction, thereby using the single particle states from the latter")(
      "Spectrum.ESShape", po::value<std::string>()->default_value("Fermi"),
      "Set the base shape name for the electrostatic finite size correction.")(
      "Spectrum.vnew", po::value<std::vector<double> >()->multitoken(),
      "Set the first three coefficients of the radial power expansion of the "
      "new nuclear electrostatic potential")(
      "Spectrum.vold", po::value<std::vector<double> >()->multitoken(),
      "Set the first three coefficients of the radial power expansion of the "
      "old nuclear electrostatic potential")(
      "Spectrum.NSShape", po::value<std::string>()->default_value("ModGauss"),
      "Set the shape of the weak charge distribution used in the convolutional C correction.");

  configOptions.add(spectrumOptions);
  configOptions.add_options()(
      "General.Folder", po::value<std::string>()->default_value("."),
      "Set the folder name for the results to be placed in.")(
      "Constants.gA", po::value<double>()->default_value(1.2723),
      "Specify the GT coupling constant, gA")(
      "Constants.gM", po::value<double>()->default_value(4.706),
      "Specify the weak magnetism coupling constant, gM")(
      "Constants.gP", po::value<double>()->default_value(0.),
      "Specify the induced pseudoscalar coupling constant, gP");

  /** Parse command line options
   * Included: generic options & spectrum shape options
   */
  po::options_description cmdOptions;
  cmdOptions.add(genericOptions).add(configOptions);
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
    cout << cmdOptions << endl;
  } else {
    /** Parse configuration file
     * Included: configOptions & spectrumOptions
     */
    std::ifstream configStream(configName.c_str());
    if (!configStream.is_open()) {
      std::cout << "WARNING: " << configName << " cannot be found.\n\n"
                << std::endl;
    } else {
      po::store(po::parse_config_file(configStream, configOptions, true), vm);
    }
    std::ifstream inputStream(inputName.c_str());
    if (!inputStream.is_open()) {
      std::cerr << "ERROR: " << inputName << " cannot be found.\n\n"
                << std::endl;
    } else {
      po::store(po::parse_config_file(inputStream, transitionOptions, true), vm);
    }
    po::notify(vm);
  }

  nme::NMEOptionContainer::GetInstance(argc, argv, true);
}
