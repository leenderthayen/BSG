#ifndef NUCLEARUTILITIES
#define NUCLEARUTILITIES

#include <vector>
#include <cmath>

namespace NuclearStructure {

/**
 * Return the sign of a double
 *
 * @param x double
 */
inline double sign(double x) { return (x == 0) ? 0 : x / std::abs(x); }

/**
 * The Kronecker delta for two doubles
 * 
 * @param x first double
 * @param y second double
 */
inline int delta(double x, double y) { return (x == y) ? 1 : 0; }

/**
 * Enum denoting the beta type
 */
enum BetaType { BETA_PLUS = -1, BETA_MINUS = 1 };

/**
 * Struct combining the properties of a harmonic oscillator component of a wave function
 * s is +-1 depending on whether it is j=l+-1/2
 */
struct WFComp {
  double C; /**< Amplitude of the state */
  int n; /**< principal quantum number */
  int l; /**< orbital quantum number */
  int s; /**< spin projection along z, +- 1 */
};

/**
 * Struct representing a single particle state as a collection of harmonic oscillator states
 * and containing all information about (iso)spin, projections a possible symmetry axis
 * in the case of deformation, excitation energy, and Nilsson quantum numbers
 */
struct SingleParticleState {
  int dO; /**< double of @f$ \Omega @f$ */
  int dK; /**< double of K */
  int parity; /**< parity of selected state, +- 1 */
  int lambda; /**< @f$ \Lambda @f$ value */
  int nDom; /**< dominant radial quantum number */
  int nZ; /**< nZ quantum number in Nilsson state */
  int isospin; /**< isospin of state. Proton has isospin -1, neutron has isospin +1 */
  double energy; /**< energy of the state */
  std::vector<WFComp> componentsHO; /**< vector containing the individual WFComp components */
};

/**
 * Struct representing a one body transition between two single particle states with
 * a certain strength given by external input
 */
struct OneBodyTransition {
  double obdme; /**< one body density matrix element */
  int dKi; /**< double of @f$ K_i @f$ */
  int dKf; /**< double of @f$ K_f @f$ */
  SingleParticleState spsi; /**< initial single particle state */
  SingleParticleState spsf; /**< final single particle state */
};

/***
 * Struct representing a nuclear state
 */
struct Nucleus {
  int Z; /**< proton number */
  int A; /**< mass number */
  int dJ; /**< double of nuclear spin */
  double R; /**< nuclear radius in natural units */
  double excitationEnergy; /**< excitation energy in keV */
  double beta2; /**< @f$ \beta_2 @f$ quadrupole deformation */
  double beta4; /**< @f$ \beta_4 @f$ hexadecupole deformation */
  double beta6; /**< @f$ \beta_6 @f$ deformation */
};

}
#endif
