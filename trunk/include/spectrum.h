#ifndef SPECTRUM_H
#define SPECTRUM_H

#include <string>
#include <iostream>
#include <vector>
#include <gsl/gsl_complex.h>

#include "constants.h"

/**
 *@brief Spectrum generator class.
 *
 * A class to be used for spectrum generation.
 * The ReadFromFile() or the alternative constructor has to be called before
 *Initialize().
 * The spectrum will be written to a stream by WriteSpectrumToStream().
 */
class Spectrum {
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
  double deformation;
  BetaType fBetaType;  ///< type of beta decay, negative or positive
  DecayType fDecayType;
  std::string exParamFile;

  /// recoil correction form factors
  double fb, fc1, fd;
  double gA, gP;

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

  // the different corrections
  /**
   * @brief Fermi function
   * @param W total energy in units of electron mass
   * @return the value
   *
   * Taken from Wilkinson, the "traditional form", whatever that means.
   * \f[ F(Z,W) = 2(\gamma+1) \Gamma(2\gamma+1)^{-2} (2pR)^{2(\gamma-1)}
   *e^{\pi\alpha Z W /p} |\Gamma(\gamma + i \alpha Z W /p) |^2\f]
   * the last factor \f$ |\Gamma(\gamma + i \alpha Z W /p) |^2 \f$ is calculated
   *using GSL. The GSL function returns the magnitude and phase, and here, as we
   *are interested in the absolute value of the result we just use the
   *magnitude.
   * The calculation of the 2nd factor \f$ \Gamma(2\gamma+1)^{-2} \f$ and the
   *last factor \f$ |\Gamma(\gamma + i \alpha Z W /p) |^2 \f$ is combined, i.e.
   * \f[ \ln\frac{a^2}{b^2} = 2\ln\frac{a}{b} = 2(\ln a - \ln b) \f]
   *
   */
  double FermiFunction(double W);

  /**
   * @brief C corrections
   * @param W total energy in units of electron mass
   * @return the value
   *
   * Long description.
   */
  double CCorrection(double W);

  double CICorrection(double W);

  //double QCDInducedCorrection(double W);

  double RelativisticCorrection(double W);

  double DeformationCorrection(double W);

  /**
   * @brief L0Correction
   * @param W total energy in units of electron mass
   * @return
   *
   *
   */
  double L0Correction(double W, double r);

  double UCorrection(double W);

  double QCorrection(double W);

  /**
   * @brief Radiative corrections
   * @param W total energy in units of electron mass
   * @return first and second order radiative corrections
   *
   * The first order correction is \f$ \delta_1=\frac{\alpha}{2\pi}g(W_0,W) \f$,
   *where \f$ g \f$ is defined by Sirlin (1967).
   * \f{eqnarray*} g(W_0,W) = & 3 \ln{\frac{m_p}{m_e}} - \frac{3}{4} +
   *\frac{4}{\beta}\rm{Spence}{\frac{2\beta}{1+\beta}} + 4 \left(
   *\frac{\tanh^{-1}{\beta}}{\beta}-1 \right) \\
   * & \times \left[ \frac{W_0-W}{3W} - \frac{3}{2} + \ln{2(W_0-W)} \right]  \\
   * & + \frac{\tanh^{-1}{\beta}}{\beta} \left[ 2(1+\beta^2) +
   *\frac{(W_0-W)^2}{6W^2} -4\tanh^{-1}{\beta}\right]
   * \f}
   *
   * where \f$ \tanh^{-1} \f$ is the inverse hyperbolic tangent function, \f$
   *\beta = \sqrt{W^2-1} \f$ and the Spence function is defined elsewhere.
   * @see Spence()
   */
  double RadiativeCorrection(double W);

  double NeutrinoRadiativeCorrection(double W);

  /**
   *@brief corrections from the recoiling nucleus
   *@param W total energy in units of electron mass
   *
   * Extend this!
   */
  double RecoilCorrection(double W);

  /**
   *@brief screening corrections as per Huber.
   *@param W total energy in units of electron mass
   *
   * The value of Ntilde is calculated during initialization.
   */
  double AtomicScreeningCorrection(double W);

  double AtomicExchangeCorrection(double W);

  double AtomicMismatchCorrection(double W);

  // helper functions
  /**
   * @brief Spence function
   * @param x
   * @return -gs_sf_dilog(x)
   *
   * The Spence function is included here because the Wilkinson article refers
   *to it. It's value is actually equal (with a sign change) to the DiLogarithm,
   *which is used from GSL.
   * \f[ \rm{Spence} = \int_0^x\frac{\ln{(1-t)}}{t}\mathrm{d}t \equiv -
   *\sum_{k=1}^{k=\infty}\frac{x^k}{k^2} \equiv -\mathrm{Li}_2(x) \f]
   */
  double Spence(double x);

  // helper functions
  /**
   * @brief Spence function
   * @param x
   * @return -gs_sf_dilog(x)
   *
   * The Spence function is included here because the Wilkinson article refers
   *to it. It's value is actually equal (with a sign change) to the DiLogarithm,
   *which is used from GSL.
   * \f[ \rm{Spence} = \int_0^x\frac{\ln{(1-t)}}{t}\mathrm{d}t \equiv -
   *\sum_{k=1}^{k=\infty}\frac{x^k}{k^2} \equiv -\mathrm{Li}_2(x) \f]
   */

 public:
  /// default constructor, will initialise to 60Co
  Spectrum();
  /// nothing is created in the heap, therefore we have an empty destructor
  ~Spectrum(){};

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
  Spectrum(char* FileName);

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
  //void SetQCDInducedCorrection(bool qic) { fQCDInducedCorrection = qic; };
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

  /**
   * @brief initialize the recoil corrections
   * @param FileName
   *
   * Parameters are read from file FileName
   */
  //void InitRecoil(char* FileName);

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
