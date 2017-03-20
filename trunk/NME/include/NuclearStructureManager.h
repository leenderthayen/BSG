#ifndef NUCLEAR_STRUCTURE_MANAGER
#define NUCLEAR_STRUCTURE_MANAGER

#include <string>
#include <vector>

namespace NuclearStructure {

// s is +-1 depending of whether it is j=l+-1/2
struct WFComp {
  double C;
  int n, l, s;
};

struct SingleParticleState {
  int dO, dK;
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
  double beta2, beta4;
};

class NuclearStructureManager {
  public:
    NuclearStructureManager() {};
    NuclearStructureManager(Nucleus, Nucleus);
    ~NuclearStructureManager();

    void SetDaughterNucleus(int, int, int, int, double, double, double, double);
    void SetMotherNucleus(int, int, int, int, double, double, double, double);
    void Initialize(std::string);

    double CalculateMatrixElement(bool, int, int, int);
    double CalculateWeakMagnetism();
    double CalculateInducedTensor();

    void AddOneBodyTransition(double, int, int, SingleParticleState, SingleParticleState);

  private:
    Nucleus mother, daughter;
    enum BetaType { BETA_PLUS = -1, BETA_MINUS = 1 };
    BetaType betaType;
    std::vector<OneBodyTransition> oneBodyTransitions;
};
}
#endif
