#ifndef SPECTRALFUNCTIONS
#define SPECTRALFUNCTIONS

#include <string>
#include <iostream>
#include <vector>
#include <gsl/gsl_complex.h>

#include "Constants.h"

namespace SpectralFunctions {
/// type of beta decay, negative or positive
enum BetaType { BETA_PLUS = -1, BETA_MINUS = 1 };
enum DecayType { FERMI, GAMOW_TELLER, MIXED };

// the different corrections
double PhaseSpace(double W, double W0, int motherSpinParity,
                  int daughterSpinParity);
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
double FermiFunction(double W, int Z, double R, int betaType);

/**
 * @brief C corrections
 * @param W total energy in units of electron mass
 * @return the value
 *
 * Long description.
 */
double CCorrection(double W, double W0, int Z, int A, double R, int betaType,
                   double hoFit, int decayType, double gA, double gP,
                   double fc1, double fb, double fd, double ratioM121);

double CICorrection(double W, double W0, int Z, int A, double R, int betaType,
                    int decayType);

double RelativisticCorrection(double W, double W0, int Z, int A, double R,
                              int betaType, int decayType);

double DeformationCorrection(double W, double W0, int Z, double R, double beta2,
                             int betaType);

/**
 * @brief L0Correction
 * @param W total energy in units of electron mass
 * @return
 *
 *
 */
double L0Correction(double W, int Z, double r, int betaType, double aPos[],
                    double aNeg[]);

double UCorrection(double W, int Z, int betaType);

double QCorrection(double W, double W0, int Z, int A, int betaType,
                   int decayType, double mixingRatio);

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
double RadiativeCorrection(double W, double W0, int Z, double R, int betaType,
                           double gA, double gM);

double NeutrinoRadiativeCorrection(double Wv);

/**
 *@brief corrections from the recoiling nucleus
 *@param W total energy in units of electron mass
 *
 * Extend this!
 */
double RecoilCorrection(double W, double W0, int A, int decayType,
                        double mixingRatio);

/**
 *@brief screening corrections as per Huber.
 *@param W total energy in units of electron mass
 *
 * The value of Ntilde is calculated during initialization.
 */
double AtomicScreeningCorrection(double W, int Z, int betaType);

double AtomicExchangeCorrection(double W, double exPars[9]);

double AtomicMismatchCorrection(double W, double W0, int Z, int A,
                                int betaType);

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
}

#endif  // SPECTRALFUNCTIONS
