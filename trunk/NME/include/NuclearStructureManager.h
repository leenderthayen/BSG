#ifndef NUCLEAR_STRUCTURE_MANAGER
#define NUCLEAR_STRUCTURE_MANAGER

namespace NuclearStructure {

// s is +-1 depending of whether it is j=l+-1/2
struct WFComp {
  double C;
  int n, l, s;
};

struct SingleParticleState {
  int dO, dK;
  int parity;
  double energy;
  std::vector<WFComp> componentsHO;
};

struct OneBodyTransition {
  double obdme;
  SingleParticleState spsi, spsf;
};

class NuclearStructureManager {

  public:
    NuclearStructureManager();
    ~NuclearStructureManager();

    double CalculateMatrixElement(bool, int, int, int);
    double CalculateWeakMagnetism();
    double CalculateInducedTensor();
    double CalculateRatioM121();

  private:
    int Z, N, A;
    double beta2, beta4;
    int dJi, dJf;
    std::vector<OneBodyTransition> oneBodyTransitions;
}

}
#endif
