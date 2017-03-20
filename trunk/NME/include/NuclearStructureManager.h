#ifndef NUCLEAR_STRUCTURE_MANAGER
#define NUCLEAR_STRUCTURE_MANAGER

#include <string>
#include <vector>

#include "NuclearUtilities.h"

namespace NuclearStructure {

class NuclearStructureManager {
  public:
    NuclearStructureManager() {};
    NuclearStructureManager(BetaType, Nucleus, Nucleus);
    ~NuclearStructureManager();

    void SetDaughterNucleus(int, int, int, double, double, double, double);
    void SetMotherNucleus(int, int, int, double, double, double, double);
    void Initialize(std::string, std::string);

    double CalculateMatrixElement(bool, int, int, int);
    double CalculateWeakMagnetism();
    double CalculateInducedTensor();

    void AddOneBodyTransition(double, int, int, SingleParticleState, SingleParticleState);

  private:
    Nucleus mother, daughter;
    BetaType betaType;
    std::vector<OneBodyTransition> oneBodyTransitions;
    std::string method, potential;

    void GetESPStates(SingleParticleState&, SingleParticleState&, std::string, int&, int&, double&);
    void GetESPOrbitalNumbers(int&, int&, int&, int&, int&, int&);
};

}
#endif
