#ifndef GENERATOR
#define GENERATOR

#include <string>
#include <iostream>
#include <vector>

/**
 *@brief Generator class.
 *
 * A class to be used for spectrum generation.
 * The ReadFromFile() or the alternative constructor has to be called before
 *Initialize().
 * The spectrum will be written to a stream by WriteSpectrumToStream().
 */
class Generator {
 private:
  /// type of beta decay, negative or positive
  enum BetaType { BETA_PLUS = -1, BETA_MINUS = 1 };
  enum DecayType { FERMI, GAMOW_TELLER, MIXED };
  double aNeg[7];
  double aPos[7];
  double exPars[9];
  double W0;     ///< total endpoint energy in electron mass
  double gamma;  ///< \f$ \gamma = \sqrt{1-\alpha^2Z^2} \f$ where \f$ \alpha \f$
  /// is the fine structure constant
  double R;   ///< nuclear radius in natural units
  double An;  ///< Mass number \f$ A \f$
  double Zd;  ///< number of protons of the daughter nucleus \f$ Z \f$
  double EndPointEnergy;  ///< \a kinetic beta endpoint energy, in keV
  double MixingRatio;
  double Afit;
  double beta2, beta4;
  BetaType fBetaType;  ///< type of beta decay, negative or positive
  DecayType fDecayType;
  std::string exParamFile;

  /// recoil correction form factors
  double fb, fc1, fd;
  double gA, gP, gM;

  double bAc;

  // options to control what goes into the Spectrum Shape
  // if false then they are NOT applied
  bool fPhaseSpace;              ///< apply the phase space factors?
  bool fFermiFunction;           ///< apply the Fermi function?
  bool fCCorrection;            ///< apply the C0 corrections?
  bool fCICorrection;            ///< apply the CI corrections?
  //bool fQCDInducedCorrection;    ///< apply the QCD induced recoil corrections?
  bool fRelativisticCorrection;  ///< apply the relativistic corrections?
  bool fL0Correction;            ///< apply the L0 corrections?
  bool fUCorrection;             ///< apply the U correction?
  bool fQCorrection;             ///< apply the Q correction?
  bool fRadiativeCorrection;     ///< apply the radiative corrections?
  bool fRecoilCorrection;  ///< apply the corrections for the recoiling nucleus?
  bool fAtomicScreeningCorrection;  ///< apply the screening corrections?
  bool fAtomicExchangeCorrection;   ///< apply the atomic exchange corrections?
  bool fAtomicMismatchCorrection;   ///< apply the atomic mismatch corrections?
  bool fDeformationCorrection;

  bool fCalcNeutrinoSpectrum;

 public:
  /// default constructor, will initialise to 60Co
  Generator();
  /// nothing is created in the heap, therefore we have an empty destructor
  ~Generator(){};

  /** Constructor, read decay details from FileName
   *
   * The file should contain at least 8 lines (the rest are ignored).
   *-# Explanation, skipped
   *-# An
   *-# Explanation, skipped
   *-# Zd
   *-# Explanation, skipped
   *-# EndPointEnergy (keV)
   *-# Explanation, skipped
   *-# Decay type (0=electron,1=positron)
   */
  Generator(char* FileName);

  /// Load parameters from FileName e.g. Mass, EndpointEnergy
  void ReadFromFile(char* FileName);

  /// set mass number A of the nucleus
  void SetMass(double mass);
  /// set number of protons Z of the daughter nucleus
  void SetDaughterZ(double Z);
  /// set beta endpoint energy (keV)
  void SetEndPointEnergy(double e0);

  void SetPhaseSpace(bool ps) { fPhaseSpace = ps; };
  void SetFermiFunction(bool ff) { fFermiFunction = ff; };
  void SetCCorrection(bool cc) { fCCorrection = cc; };
  void SetCICorrection(bool cic) { fCICorrection = cic; };
  void SetRelativisticCorrection(bool rc) { fRelativisticCorrection = rc; };
  void SetDeformationCorrection(bool dc) { fDeformationCorrection = dc; };
  void SetL0Correction(bool l0c) { fL0Correction = l0c; };
  void SetUCorrection(bool uc) { fUCorrection = uc; };
  void SetQCorrection(bool qc) { fQCorrection = qc; };
  void SetRadiativeCorrection(bool rc) { fRadiativeCorrection = rc; };
  void SetRecoilCorrection(bool rc) { fRecoilCorrection = rc; };
  void SetAtomicScreeningCorrection(bool as) {
    fAtomicScreeningCorrection = as;
  };
  void SetAtomicExchangeCorrection(bool ec) { fAtomicExchangeCorrection = ec; };
  void SetAtomicMismatchCorrection(bool amc) {
    fAtomicMismatchCorrection = amc;
  };
  void SetCalculateNeutrinoSpectrum(bool vs) { fCalcNeutrinoSpectrum = vs; };
  void SetWeakMagnetism(double b) { fc1 = 1.; bAc = b; };
  
  double CalculateWeakMagnetism();
  double CalculateInducedTensor();
  double CalculateInducedPseudoscalar();
  void CalculateMatrixElements();

  /**
   * Initialize the spectrum generator.
   * Must be called after the constructor, or any change to the decay parameters
   * (A,Z,E0).
   */
  void Initialize();

  void LoadExchangeParameters();

  /// The beta spectrum shape at the given energy (keV)
  std::vector<double> GetSpectrumShape(double Energy);

  /** Writes the beta spectrum shape to a <tab> separated file, suitable for
   *plotting with gnuplot.
   *@param ofile pointer to a stream where the information is written
   *@param Step energy steps, in keV
   */
  void WriteSpectrumToStream(std::ostream* ofile, double Step);

  /** Print information about the SpectrumGenerator such as Mass,
   * EndpointEnergy, decay type, other options
   * @param ofile pointer to a stream where the information is written
   */
  void PrintStatus(std::ostream* ofile);
};

#endif  // SPECTRUM_H
