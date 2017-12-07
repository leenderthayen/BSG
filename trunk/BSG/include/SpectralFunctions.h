#ifndef SPECTRALFUNCTIONS
#define SPECTRALFUNCTIONS

#include <string>
#include <iostream>
#include <vector>
#include <gsl/gsl_complex.h>

#include "Constants.h"
#include "NuclearUtilities.h"

/**
 * Namespace containing all correction factors to the allowed beta spectrum shape
 * as described by Hayen et al., Rev. Mod. Phys. XXXXXX
 */
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
 * @param Z the proton number of the daughter
 * @param R the nuclear radius of the daughter
 * @param betaType the BetaType of the transition
 * @return the value
 *
 * The point charge Fermi function @f$ F_0(Z, W) @f$
 * \f[ F_0(Z,W) = 2(\gamma+1) \Gamma(2\gamma+1)^{-2} (2pR)^{2(\gamma-1)}
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
 * @brief C correction
 * @param W total energy in units of electron mass
 * @param W0 the total endpoint energy in units of the electron rest mass
 * @param Z proton number
 * @param A mass number
 * @param R nuclear radius in natural units
 * @param betaType the BetaType of the transition
 * @param hoFit the fitted A value for the Modified Gaussian distribution
 * @param decayType the DecayType of the transition
 * @param gA the axial vector coupling constant
 * @param gP the induced pseudoscalar coupling constant
 * @param fc1 the c1 form factor as per Holstein (@f$ g_A M_{GT} @f$)
 * @param fb the b form factor as per Holstein
 * @param fd the d form factor as per Holstein
 * @param ratioM121 the ratio of the @f$^AM_{121}@f$ and @f$ ^AM_{101}@f$ matrix elements in the Behrens-Buehring formalism
 * @return the value
 *
 * The C correction describing effects of finite nuclear size and induced currents.
 */
double CCorrection(double W, double W0, int Z, int A, double R, int betaType,
                   double hoFit, int decayType, double gA, double gP,
                   double fc1, double fb, double fd, double ratioM121);

/**
 * Isovector correction to the charge density-calculated C correction
 * 
 * @param W electron total energy in units of its rest mass
 * @param W0 total endpoint energy in units of the electron rest mass
 * @param Z proton number
 * @param A mass number
 * @param R nuclear radius in natural units
 * @param betaType BetaType of the transition
 * @param decayType DecayType of the transition
 */
double CICorrection(double W, double W0, int Z, int A, double R, int betaType,
                    int decayType);

/**
 * Isovector correction to the charge density-calculated C correction
 * Gets called when connection with NME is turned on, thereby using actual
 * single particle wave functions from a (deformed) Woods-Saxon potential
 * 
 * @param W electron total energy in units of its rest mass
 * @param W0 total endpoint energy i units of the electron rest mass
 * @param Z proton number
 * @param R nuclear radius in natural units
 * @betaType beta type of the transition
 * @param spsi NuclearStructure::SingleParticleState object denoting the initial state
 * @param spsf NuclearStructure::SingleParticleState object denoting the final state
 */
double CICorrection(double W, double W0, double Z, double R, int betaType,
    NuclearStructure::SingleParticleState spsi, NuclearStructure::SingleParticleState spsf);

/**
 * Relativistic matrix element correction to the vector part of the C correction
 * as per Wilkinson. 
 * 
 * @param W electron total energy in units of its rest mass
 * @param W0 total endpoint energy in units of the electron rest mass
 * @param Z proton number
 * @param A mass number
 * @param R nuclear radius in natural units
 * @param betaType BetaType of the transition
 * @param decayType DecayType of the transition
 */
double RelativisticCorrection(double W, double W0, int Z, int A, double R,
                              int betaType, int decayType);

/**
 * Deformation correction to L0
 * 
 * @param W electron total energy in units of its rest mass
 * @param W0 total endpoint energy in units of the electron rest mass
 * @param Z proton number
 * @param A mass number
 * @param R nuclear radius in natural units
 * @param beta2 quadrupole deformation
 * @param betaType BetaType of the transition
 */
double DeformationCorrection(double W, double W0, int Z, double R, double beta2,
                             int betaType);

/**
 * Electrostatic finite size correction to the point charge Fermi function
 * written as @f$L_0(Z, W) @f$
 *
 * @param W electron total energy in units of its rest mass
 * @param W0 total endpoint energy in units of the electron rest mass
 * @param Z proton number
 * @param r radius in natural units
 * @param betaType BetaType of the transition
 * @param aPos array of fitted constants for beta+ decay
 * @param aNeg array of fitted constants for beta- decay
 */
double L0Correction(double W, int Z, double r, int betaType, double aPos[],
                    double aNeg[]);

/**
 * Correction to L0 by moving from a uniformly charged sphere to a more elaborate nuclear shape
 * 
 * @param W electron total energy in units of its rest mass
 * @param Z proton number
 * @param R nuclear radius in natural units
 * @param betaType BetaType of the transition
 * @param baseShape string denoting the name of the base shape. Currently only Fermi is implemented
 * @param v vector representing the first 3 terms in an even-r power expansion of the base shape
 * @param vp vector representing the first 3 terms in an even-r power expansion of the new shape
 */
double UCorrection(double W, int Z, double R, int betaType, std::string baseShape, std::vector<double>& v, std::vector<double>& vp);

/**
 * Correction to L0 by calculating @f$ \frac{L_0'}{L_0} @f$ using a power expansion of the potentials
 * 
 * @param W electron total energy in units of its rest mass
 * @param Z proton number
 * @param R nuclear radius in natural units
 * @param betaType BetaType of the transition
 * @param v vector representing the first 3 terms in an even-r power expansion of the base shape
 * @param vp vector representing the first 3 terms in an even-r power expansion of the new shape
 * @see UCorrection
 */
double UCorrection(double W, int Z, double R, int betaType, std::vector<double>& v, std::vector<double>& vp);

/**
 * Electromagnetic correction to the Fermi function due to the change in the 
 * electromagnetic field of the recoiling nucleus compared to it standing still
 * 
 * @param W electron total energy in units of its rest mass
 * @param W0 total endpoint energy in units of the electron rest mass
 * @param Z proton number
 * @param A mass number
 * @param betaType BetaType of the transition
 * @param decayType decsay type of the transition
 * @double mixingRatio mixing ratio of Fermi versus Gamow-Teller
 */
double QCorrection(double W, double W0, int Z, int A, int betaType,
                   int decayType, double mixingRatio);

/**
 * The radiative correction up to order @f$ \alpha^3Z^2 @f$ as per Sirlin
 * 
 * @param W total electron energy in units of its rest mass
 * @param W0 total endpoint energy in units of the electron rest mass
 * @param Z proton number
 * @param R nuclear radius in natural units
 * @param betaType BetaType of the transition
 * @param gA axial vector coupling constant @f$ g_A @f$
 * @param gM nucleon isovector moment @f$ g_M = 4.706 @f$
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

/**
 * Radiative correction to the neutrino spectrum by Sirlin and Marciano
 * 
 * @param Wv neutrino total energy in units of the electron rest mass
 */
double NeutrinoRadiativeCorrection(double Wv);

/**
 * Kinematic recoil correction to the beta spectrum due to the enlargement of the phase space
 * 
 * @param W electron total energy in units of its rest mass
 * @param W0 total endpoint energy in units of the electron rest mass
 * @param A mass number
 * @param decayType decsay type of the transition
 * @double mixingRatio mixing ratio of Fermi versus Gamow-Teller
 */
double RecoilCorrection(double W, double W0, int A, int decayType,
                        double mixingRatio);

/**
 * Correction due to atomic screening calculated using the Salvat potential
 * 
 * @param W electron total energy in units of its rest mass
 * @param Z proton number
 * @param betaType the BetaType of the transition
 */
double AtomicScreeningCorrection(double W, int Z, int betaType);

/**
 * The atomic exchange correction where an electron decays into a bound state of the daughter atom
 * and its corresponding interference with the direct process
 * 
 * @param W electron total energy in units of its rest mass
 * @param exPars array containing the 9 fit parameters required for the analytical parametrisation
 */
double AtomicExchangeCorrection(double W, double exPars[9]);

/**
 * Correction due to the mismatch between initial and final atomic states, causing the
 * endpoint of the transition to get smaller
 * 
 * @param W electron total energy in units of its rest mass
 * @param W0 total endpoint energy in units of the electron rest mass
 * @param Z proton number
 * @param A mass number
 * @param betaType BetaType of the transition
 */
double AtomicMismatchCorrection(double W, double W0, int Z, int A,
                                int betaType);

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
