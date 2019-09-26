#ifndef GENERATOR
#define GENERATOR

#include <vector>
#include <string>
#include <tuple>
#include "NuclearStructureManager.h"
#include "NuclearUtilities.h"
#include "spdlog/spdlog.h"

#include "particle_data_structure/core/particle_definition.h"

namespace BSG {

  /**
   * Enum to distinguish beta+/-
   */
  enum BetaType { BETA_PLUS = -1, BETA_MINUS = 1 };
  /**
   * Enum to distinguish the type of beta decay
   */
  enum DecayType { FERMI, GAMOW_TELLER, MIXED };

  struct BetaParams {
    int Zi; /**< Proton number of initial state */
    int Zf; /**< Proton number of final state */
    int A; /**< Mass number of initial and final states */
    double R; /**< Nuclear radius of the final state */
    double W0;  /**< the total electron energy in units of its rest mass */
    double mixingRatio; /**< the mixing ratio of Fermi vs Gamow-Teller decay */
    BetaType betaType;  /**< internal state of the beta type */
    DecayType decayType; /**< internal state of the decay type */
    double exPars[9]; /**< array for the fit coefficients of the atomic exchange correction */
  };

  struct NuclearParams {
    ParticleDefinition* initialParticle, finalParticle;
    int initState, finalState;
    double daughterBeta2; /**< the quadrupole deformation of the daughter nucleus */
    double motherBeta2; /**< the quadrupole deofrmation of the mother nucleus */
    nme::NuclearStructure::SingleParticleState spsi, spsf; /**< single particle states calculated from the NME library and used in the C_I correction when turned on */
  };

class Generator {
 private:
  static double aNeg[7]; /**< array containing Wilkinson's fit coefficients for the L0 correction for beta- decay */
  static double aPos[7]; /**< array containing Wilkinson's fit coefficients for the L0 correction for beta+ decay */

  double QValue; /**< Q value of the particular decay */
  double hoFit; /**< the fit value obtained after fitting the nuclear charge distribution with a Modified Gaussian */
  double atomicEnergyDeficit; /**< explicit energy difference between Q value and endpoint energy due to incomplete atomic overlap (atomic excitations) */

  nme::NuclearStructure::NuclearStructureManager* nsm; /**< pointer to the nuclear structure manager */

  std::vector<std::vector<double> >* spectrum; /**< vector of vectors containing the calculated spectrum */

  /// recoil correction form factors
  double fb, fc1, fd, ratioM121;
  double bAc, dAc;
  double gA, gP, gM; /**< coupling constants of the weak Hamiltonian */

  std::shared_ptr<spdlog::logger> consoleLogger;
  std::shared_ptr<spdlog::logger> debugFileLogger;
  std::shared_ptr<spdlog::logger> rawSpectrumLogger;
  std::shared_ptr<spdlog::logger> resultsFileLogger;

  std::string outputName;

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

  void SetInitialState(Nucleus*, int);

  void SetFinalState(Nucleus*, int);

  void SetDecayKinematics(ReactionChannel*);

  /**
   * Calculates the beta spectrum by filling the spectrum variable.
   *
   * @returns spectrum variable
   */
  std::vector<std::vector<double> >* CalculateSpectrum();

  /**
   * Construct the output file
   */
  void PrepareOutputFile();

  /**
   * Calculate the decay rate at energy W
   *
   * @param W the total electron energy in units of its rest mass
   * @returns the decay rate at energy W
   */
  std::tuple<double, double> CalculateDecayRate(double W);

  inline void SetOutputName(std::string _output) { outputName = _output; };
};

}

#endif
