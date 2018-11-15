#ifndef NUCLEAR_STRUCTURE_MANAGER
#define NUCLEAR_STRUCTURE_MANAGER

#include <string>
#include <vector>
#include <map>

#include "NuclearUtilities.h"
#include "spdlog/spdlog.h"

namespace NuclearStructure {

/**
 * Manager class dealing with the calculation of single particle states and the
 * calculation of matrix elements
 */
class NuclearStructureManager {
 public:
  /**
   * Constructor
   */
  NuclearStructureManager();
  /**
   * Overloaded constructor
   *
   * @param bt beta type of the transition
   * @param init Nucleus object for the initial state
   * @param fin Nucleus object for the final state
   */
  NuclearStructureManager(BetaType bt, Nucleus init, Nucleus fin);
  /**
   * Destructor
   */
  ~NuclearStructureManager(){};

  /**
   * Creates a daughter Nucleus object
   *
   * @param Z proton number
   * @param A mass number
   * @param dJ double of spin
   * @param R nuclear radius in natural units
   * @param excitationEnergy excitationEnergy in keV
   * @param beta2 quadrupole deformation
   * @param beta4 hexadecupole deformation
   * @param beta6 beta6 deformation
   */
  void SetDaughterNucleus(int Z, int A, int dJ, double R,
                          double excitationEnergy, double beta2, double beta4,
                          double beta6);
  /**
   * Creates a mother Nucleus object
   *
   * @param Z proton number
   * @param A mass number
   * @param dJ double of spin
   * @param R nuclear radius in natural units
   * @param excitationEnergy excitationEnergy in keV
   * @param beta2 quadrupole deformation
   * @param beta4 hexadecupole deformation
   * @param beta6 beta6 deformation
   */
  void SetMotherNucleus(int Z, int A, int dJ, double R, double excitationEnergy,
                        double beta2, double beta4, double beta6);
  /**
   * Creates all one body transitions and accompanying single particle states
   *
   * @param m method used
   * @param p potential to choose
   */
  void Initialize(std::string m, std::string p);
  /**
   * Set the beta type of the transition
   */
  inline void SetBetaType(BetaType bt) { betaType = bt; };

  /**
   * Calculate the @f[ ^{V/A}M_{KLs} @f] matrix element in the Behrens-Buehring
   *formalism
   *
   * @param V boolean for vector (true) or axial vector (false) matrix element
   *type
   * @param K spherical tensor rank of the operator
   * @param L orbital angular momentum of the operator
   * @param s index specifying simple or vector spherical harmonics
   */
  double CalculateReducedMatrixElement(bool V, int K, int L, int s);
  /**
   * Calculate b/Ac in the Holstein formalism
   */
  double CalculateWeakMagnetism();
  /**
   * Calculate d/Ac in the Holstein formalism
   */
  double CalculateInducedTensor();

  /**
   * Add a one body transition to the list
   *
   * @param K angular momentum coupling of [a+a]_K
   * @param obdme one body density matrix element
   * @param dKi double the K value of the initial state
   * @param dKf double the K value of the final state
   * @param spsi initial single particle state
   * @param spsf final single particle state
   */
  void AddReducedOneBodyTransitionDensity(int K, double robdme, int dKi, int dKf,
                            SingleParticleState spsi, SingleParticleState spsf);

  void GetESPStates(SingleParticleState&, SingleParticleState&, int&, int&);

 private:
  Nucleus mother, daughter;
  BetaType betaType;
  std::map<int, std::vector<ReducedOneBodyTransitionDensity> > reducedOneBodyTransitionDensities;
  std::string method, potential;

  std::shared_ptr<spdlog::logger> consoleLogger;
  std::shared_ptr<spdlog::logger> debugFileLogger;
  std::shared_ptr<spdlog::logger> nmeResultsLogger;

  bool initialized = false;

  void InitializeLoggers();
  void InitializeConstants();

  void GetESPOrbitalNumbers(int&, int&, int&, int&, int&, int&);
  double GetESPManyParticleCoupling(int, ReducedOneBodyTransitionDensity&);
  bool BuildDensityMatrixFromFile(std::string);
  void ReadNuShellXOBD(std::string);
};
}
#endif
