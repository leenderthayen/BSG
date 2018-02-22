#ifndef NUCLEARUTILITIES
#define NUCLEARUTILITIES

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iterator>
#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>

namespace NuclearStructure {

static const char* atoms[118] = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti",
"V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh",
"Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er",
"Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa",
"U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Uut", "Fl",
"Uup", "Lv", "Uus", "Uuo"};

/**
 * Return the sign of a double
 *
 * @param x double
 */
inline int sign(double x) { return (x == 0) ? 0 : (int)(x / std::abs(x)); }

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

namespace GeneralUtilities {

inline std::vector<std::vector<std::string> > getCSVData(std::string filename, std::string delimeter) {
  std::ifstream file(filename);

  std::vector<std::vector<std::string> > dataList;

  std::string line = "";
  // Iterate through each line and split the content using delimeter
  while (getline(file, line))
  {
    std::vector<std::string> vec;
    boost::algorithm::split(vec, line, boost::is_any_of(delimeter));
    dataList.push_back(vec);
  }
  // Close the File
  file.close();
 
  return dataList;
}
}
#endif
