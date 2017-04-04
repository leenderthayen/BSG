#ifndef NUCLEARUTILITIES
#define NUCLEARUTILITIES

#include <vector>
#include <cmath>

namespace NuclearStructure {

inline double sign(double x) { return (x == 0) ? 0 : x / std::abs(x); }

inline int delta(double x, double y) { return (x == y) ? 1 : 0; }

enum BetaType { BETA_PLUS = -1, BETA_MINUS = 1 };

// s is +-1 depending of whether it is j=l+-1/2
struct WFComp {
  double C;
  int n, l, s;
};

struct SingleParticleState {
  int dO;
  int dK;
  int parity;
  //Proton has isospin -1, neutron has isospin +1
  int isospin;
  double energy;
  std::vector<WFComp> componentsHO;
};

struct OneBodyTransition {
  double obdme;
  int dKi, dKf;
  SingleParticleState spsi, spsf;
};

struct Nucleus {
  int Z, A;
  int dJ;
  double R;
  double excitationEnergy;
  double beta2, beta4, beta6;
};

}
#endif
