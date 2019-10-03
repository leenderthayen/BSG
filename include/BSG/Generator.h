#ifndef BSG_GENERATOR_H
#define BSG_GENERATOR_H

//Same package
#include "BSG/ConfigParser.h"

//Own libraries
#include "PDS/Core/ParticleDefinition.h"
#include "NHL/Containers.h"

#include "spdlog/spdlog.h"

#include <vector>
#include <string>
#include <tuple>
#include <array>
#include <map>

namespace BSG {

  //Struct that is passed to the Spectral Corrections
  struct BetaParams {
    int Zi = 0; /**< Proton number of initial state */
    int Zf = 0; /**< Proton number of final state */
    int A = 0; /**< Mass number of initial and final states */
    double R = 0.; /**< Nuclear radius of the final state */
    double W0 = 1.;  /**< the total electron energy in units of its rest mass */
    double mixingRatio = 0.; /**< the mixing ratio of Fermi vs Gamow-Teller decay */
    NHL::BetaType betaType = NHL::BETA_MINUS;  /**< internal state of the beta type */
    NHL::BetaDecayType decayType = NHL::BetaDecayType::GAMOWTELLER; /**< internal state of the decay type */
    std::array<double, 9> exPars; /**< array for the fit coefficients of the atomic exchange correction */
    std::array<double, 7> aNeg;
    std::array<double, 7> aPos;
  };

  class Generator {
   public:
    /**
     * Constructor for Generator.
     * Performs the full initialization of the Generator by fetching arguments
     * from the commandline or config file, performing the L0 initialization and charge distribution fitting
     */
    Generator(std::string _ouput = "bsg_output");
    /**
     * Destructor for Generator.
     * Deletes the reference to the nuclear structure manager
     */
    ~Generator() {};

    bool Initialize(std::string, std::string, int argc = 0, const char** argv = nullptr);
    bool InitializeTransitionFromFile(std::string);
    void InitializeOptionsFromConfigFile(std::string, int argc = 0, const char** = nullptr);

    inline void SetInitialState(PDS::core::Particle _p) { transitionOptions.initNucleus = _p; InitializeBetaParams(); }
    inline void SetFinalState(PDS::core::Particle _p) { transitionOptions.finalNucleus = _p; InitializeBetaParams(); }
    inline void SetQValue(double q) { transitionOptions.QValue = q; InitializeBetaParams(); }

    //inline void SetNSM(NME::NuclearStructureManager* _nsm) { delete nsm; nsm = _nsm; }

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

    inline const BetaParams GetBetaParams() const { return betaParams; }

  private:
    TransitionOptions transitionOptions;
    ConfigOptions configOptions;
    BetaParams betaParams;
    static constexpr int bDim1 = 7;
    static constexpr int bDim2 = 6;
    double bNeg[bDim1][bDim2]; /**< array containing Wilkinson's fit coefficients for the L0 correction for beta- decay */
    double bPos[bDim1][bDim2]; /**< array containing Wilkinson's fit coefficients for the L0 correction for beta+ decay */
    std::map<int, std::array<double, 9> > exchangeCoefficients;

    //NME::NuclearStructureManager* nsm; /**< pointer to the nuclear structure manager */

    std::vector<std::vector<double> >* spectrum; /**< vector of vectors containing the calculated spectrum */

    std::shared_ptr<spdlog::logger> consoleLogger;
    std::shared_ptr<spdlog::logger> debugFileLogger;
    std::shared_ptr<spdlog::logger> rawSpectrumLogger;
    std::shared_ptr<spdlog::logger> resultsFileLogger;

    std::string outputName;

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
  };

}

#endif
