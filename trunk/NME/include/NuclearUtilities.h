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
  double C;
  int n, l, s;
};

/**
 * Struct representing a single particle state as a collection of harmonic oscillator states
 * and containing all information about (iso)spin, projections a possible symmetry axis
 * in the case of deformation, excitation energy, and Nilsson quantum numbers
 */
struct SingleParticleState {
  int dO;
  int dK;
  int parity;
  int lambda;
  int nDom;
  int nZ;
  //Proton has isospin -1, neutron has isospin +1
  int isospin;
  double energy;
  std::vector<WFComp> componentsHO;
};

/**
 * Struct representing a one body transition between two single particle states with
 * a certain strength given by external input
 */
struct OneBodyTransition {
  double obdme;
  int dKi, dKf;
  SingleParticleState spsi, spsf;
};

/***
 * Struct representing a nuclear state
 */
struct Nucleus {
  int Z, A;
  int dJ;
  double R;
  double excitationEnergy;
  double beta2, beta4, beta6;
};

}
#endif
