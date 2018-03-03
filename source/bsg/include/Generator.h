#ifndef GENERATOR
#define GENERATOR

#include <vector>
#include <string>
#include <tuple>
#include "NuclearStructureManager.h"
#include "NuclearUtilities.h"
#include "spdlog/spdlog.h"

class Generator {
 private:
  /**
   * Enum to distinguish beta+/-
   */
  enum BetaType { BETA_PLUS = -1, BETA_MINUS = 1 };
  /**
   * Enum to distinguish the type of beta decay
   */
  enum DecayType { FERMI, GAMOW_TELLER, MIXED };
  double aNeg[7]; /**< array containing Wilkinson's fit coefficients for the L0 correction for beta- decay */
  double aPos[7]; /**< array containing Wilkinson's fit coefficients for the L0 correction for beta+ decay */
  double exPars[9]; /**< array for the fit coefficients of the atomic exchange correction */
  double QValue; /**< Q value of the particular decay */
  double W0;  /**< the total electron energy in units of its rest mass */
  double R;   /**< the nuclear radius in natural units (HBAR=c=m_e=1) */
  double A;   /**< Mass number */
  double Z;   /**< the proton number of the daughter nucleus */
  double mixingRatio; /**< the mixing ratio of Fermi vs Gamow-Teller decay */
  double hoFit; /**< the fit value obtained after fitting the nuclear charge distribution with a Modified Gaussian */
  double daughterBeta2; /**< the quadrupole deformation of the daughter nucleus */
  double motherBeta2; /**< the quadrupole deofrmation of the mother nucleus */
  double motherExcitationEn; /**< the excitation energy in keV of the mother state */
  double daughterExcitationEn; /**< the excitation energy in keV of the daughter state */
  double atomicEnergyDeficit; /**< explicit energy difference between Q value and endpoint energy due to incomplete atomic overlap (atomic excitations) */
  int motherSpinParity; /**< the double of the spin parity of the mother state */
  int daughterSpinParity; /**< the double of the spin parity of the daughter state */
  BetaType betaType;  /**< internal state of the beta type */
  DecayType decayType; /**< internal state of the decay type */
  std::vector<double> vOld; /**< power expansion in r for old electrostatic potential */
  std::vector<double> vNew; /**< power expansion in r for new electrostatic potential */
  std::string ESShape; /**< name denoting the base shape to use for the U correction */
  std::string NSShape; /**< name denoting the base shape to use for the C correction */

  NuclearStructure::SingleParticleState spsi, spsf; /**< single particle states calculated from the NME library and used in the C_I correction when turned on */

  NuclearStructure::NuclearStructureManager* nsm; /**< pointer to the nuclear structure manager */

  std::vector<std::vector<double> > spectrum; /**< vector of vectors containing the calculated spectrum */

  /// recoil correction form factors
  double fb, fc1, fd, ratioM121;
  double bAc, dAc;
  double gA, gP, gM; /**< coupling constants of the weak Hamiltonian */

  std::shared_ptr<spdlog::logger> consoleLogger;
  std::shared_ptr<spdlog::logger> debugFileLogger;
  std::shared_ptr<spdlog::logger> rawSpectrumLogger;
  std::shared_ptr<spdlog::logger> resultsFileLogger;

  /**
   * Calculate the required nuclear matrix elements if they are not given from the commandline
   */
  void GetMatrixElements();

  /**
   * Initialize the hard-coded parameters aNeq and aPos for the Wilkinson L0 correction
   */
  void InitializeL0Constants();
  /**
   * Load the fit parameters of the atomic exchange correction from a file
   */
  void LoadExchangeParameters();

  /**
   * Initialize all constants such as Z, R, etc, taken from config files
   */
  void InitializeConstants();

  /**
   * Initialize all parameters related to the electrostatic shape
   */
  void InitializeShapeParameters();

  /**
   * Initialize all loggers
   */
  void InitializeLoggers();

  /**
   * Initialize all Nuclear structure manager-related stuff, like single particle states when C_I and NME are connected, and calculating the required matrix elements
   */
  void InitializeNSMInfo();

  /**
   * Construct the output file
   */
  void PrepareOutputFile();

  /**
   * Calculate the properly normalized ft value
   * @param partialHalflife the halflife of the transition
   */
  double CalculateLogFtValue(double partialHalflife);

  double CalculateMeanEnergy();

 public:
  /**
   * Constructor for Generator.
   * Performs the full initialization of the Generator by fetching arguments
   * from the commandline or config file, performing the L0 initialization and charge distribution fitting
   */
  Generator();
  /**
   * Destructor for Generator.
   * Deletes the reference to the nuclear structure manager
   */
  ~Generator();

  /**
   * Calculates the beta spectrum by filling the spectrum variable.
   * 
   * @returns spectrum variable
   */
  std::vector<std::vector<double> > CalculateSpectrum();
  /**
   * Calculate the decay rate at energy W
   * 
   * @param W the total electron energy in units of its rest mass
   * @returns the decay rate at energy W
   */
  std::tuple<double, double> CalculateDecayRate(double W);
};

#endif
