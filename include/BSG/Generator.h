#ifndef BSG_GENERATOR_H
#define BSG_GENERATOR_H

//Own libraries
#include "PDS/Core/ParticleDefinition.h"
#include "NHL/Containers.h"

#include <vector>
#include <string>
#include <tuple>

#include "spdlog/spdlog.h"


namespace BSG {

  //Container for all relevant transition info
  struct TransitionInfo {
    PDS::core::Particle initNucleus, finalNucleus;
    NHL::BetaType betaType;
    double Q;
    double atomicEnergyDeficit;
    double partialLifetime;
    double logft;
  };

  //Struct that is passed to the Spectral Corrections
  struct BetaParams {
    int Zi; /**< Proton number of initial state */
    int Zf; /**< Proton number of final state */
    int A; /**< Mass number of initial and final states */
    double R; /**< Nuclear radius of the final state */
    double W0;  /**< the total electron energy in units of its rest mass */
    double mixingRatio; /**< the mixing ratio of Fermi vs Gamow-Teller decay */
    BetaType betaType;  /**< internal state of the beta type */
    DecayType decayType; /**< internal state of the decay type */
    double exPars[9]; /**< array for the fit coefficients of the atomic exchange correction */
  };

  class Generator {
   public:
    /**
     * Constructor for Generator.
     * Performs the full initialization of the Generator by fetching arguments
     * from the commandline or config file, performing the L0 initialization and charge distribution fitting
     */
    Generator(std::string _ouput = "output");
    /**
     * Destructor for Generator.
     * Deletes the reference to the nuclear structure manager
     */
    ~Generator() { delete nsm; };

    bool Initialize(std::string, std::string);

    inline void SetInitialState(Particle _p) { transitionInfo.initNucleus = _p; }
    inline void SetFinalState(Particle _p) { transitionInfo.finalNucleus = _p; }

    inline void SetNSM(NME::NuclearStructureManager* _nsm) { delete nsm; nsm = _nsm; }

    void SetDecayKinematics(ReactionChannel*);

    /**
     * Calculates the beta spectrum by filling the spectrum variable.
     *
     * @returns spectrum variable
     */
    std::vector<std::vector<double> >* CalculateSpectrum();

    /**
     * Construct the output file
     */
    void PrepareOutputFile();

    /**
     * Calculate the decay rate at energy W
     *
     * @param W the total electron energy in units of its rest mass
     * @returns the decay rate at energy W
     */
    std::tuple<double, double> CalculateDecayRate(double W);

    inline void SetOutputName(std::string _output) { outputName = _output; };

  private:
    ConfigOptions configOptions;
    static double aNeg[7]; /**< array containing Wilkinson's fit coefficients for the L0 correction for beta- decay */
    static double aPos[7]; /**< array containing Wilkinson's fit coefficients for the L0 correction for beta+ decay */

    NME::NuclearStructureManager* nsm; /**< pointer to the nuclear structure manager */

    std::vector<std::vector<double> >* spectrum; /**< vector of vectors containing the calculated spectrum */

    std::shared_ptr<spdlog::logger> consoleLogger;
    std::shared_ptr<spdlog::logger> debugFileLogger;
    std::shared_ptr<spdlog::logger> rawSpectrumLogger;
    std::shared_ptr<spdlog::logger> resultsFileLogger;

    std::string outputName;

    bool InitializeParticlesFromFile(std::string);
    bool InitializeOptionsFromConfigFile(std::string);

    /**
    * Calculate the properly normalized ft value
    * @param partialHalflife the halflife of the transition
    */
    double CalculateLogFtValue(double partialHalflife);

    double CalculateMeanEnergy();

    void InitializeBetaParams();

    /**
    * Calculate the required nuclear matrix elements if they are not given from the commandline
    */
    void GetMatrixElements();

    /**
    * Initialize the hard-coded parameters aNeq and aPos for the Wilkinson L0 correction
    */
    void InitializeL0Constants();
    /**
    * Load the fit parameters of the atomic exchange correction from a file
    */
    void LoadExchangeParameters();

    /**
    * Initialize all constants such as Z, R, etc, taken from config files
    */
    void InitializeConstants();

    /**
    * Initialize all parameters related to the electrostatic shape
    */
    void InitializeShapeParameters();

    /**
    * Initialize all loggers
    */
    void InitializeLoggers();

    /**
    * Initialize all Nuclear structure manager-related stuff, like single particle states when C_I and NME are connected, and calculating the required matrix elements
    */
    void InitializeNSMInfo();
    bool InitializeParticlesFromFile(std::string);
    bool InitializeOptionsFromConfigFile(std::string);
  };

}

#endif
