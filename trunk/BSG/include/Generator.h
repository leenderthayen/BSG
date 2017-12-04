#ifndef GENERATOR
#define GENERATOR

#include <vector>
#include "NuclearStructureManager.h"

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
  double W0;  /**< the total electron energy in units of its rest mass */
  double R;   /**< the nuclear radius in natural units (hbar=c=m_e=1) */
  double A;   /**< Mass number */
  double Z;   /**< the proton number of the daughter nucleus */
  double mixingRatio; /**< the mixing ratio of Fermi vs Gamow-Teller decay */
  double hoFit; /**< the fit value obtained after fitting the nuclear charge distribution with a Modified Gaussian */
  double daughterBeta2; /**< the quadrupole deformation of the daughter nucleus */
  double motherBeta2; /**< the quadrupole deofrmation of the mother nucleus */
  double motherExcitationEn; /**< the excitation energy in keV of the mother state */
  double daughterExcitationEn; /**< the excitation energy in keV of the daughter state */
  int motherSpinParity; /**< the double of the spin parity of the mother state */
  int daughterSpinParity; /**< the double of the spin parity of the daughter state */
  BetaType betaType;  /**< internal state of the beta type */
  DecayType decayType; /**< internal state of the decay type */

  NuclearStructure::NuclearStructureManager* nsm; /**< pointer to the nuclear structure manager */

  std::vector<std::vector<double> > spectrum; /**< vector of vectors containing the calculated spectrum */

  /// recoil correction form factors
  double fb, fc1, fd, ratioM121;
  double gA, gP, gM; /**< coupling constants of the weak Hamiltonian */

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
  double* CalculateDecayRate(double W);
  /**
   * Write the spectrum to file using the file from the config options
   */
  void WriteSpectrumToFile();
};

#endif
