#include <iostream>
#include <vector>
#include <string>

#include "Generator.h"
#include "OptionContainer.h"
#include "NilssonOrbits.h"
#include "ChargeDistributions.h"

#include "BSGConfig.h"

using std::cout;
using std::endl;
using std::cerr;

std::string lastUpdate = "February 24, 2017";
std::string author = "L. Hayen (leendert.hayen@kuleuven.be)";

namespace NO = NuclearStructure::nilsson;
namespace CD = ChargeDistributions;

void ShowIntro() {
  cout << "**************************************************\n";
  cout << "                  BSG v." << BSG_VERSION << "\n\n";
  cout << "        Last update: " << lastUpdate << endl;
  cout << "      " << author << endl;
  cout << "**************************************************\n\n";
}

void TestNilssonMethods() {
  cout << "Testing Nilsson methods" << endl;

  /*int dim = 2;
  double M[dim*(dim+1)/2] = {1.0, 1.0, 2.0};

  std::vector<double> eVals;
  double eVecs[dim*dim] = {};

  nilsson::Eigen(M, dim, eVecs, eVals);

  cout << "Eigenvalues: " << endl;
  for (int i = 0; i < eVals.size(); i++) {
    cout << eVals[i] << " ";
  }
  cout << endl;
  cout << "Eigenvectors" << endl;
  for (int j = 0; j < dim; j++) {
    for (int i = 0; i < dim; i++) {
      cout << eVecs[j*dim+i] << " ";
    }
    cout << endl;
  }*/


  double V0 = 53.0;
  double R = 3.655;
  double A0 = 0.650;
  double VS = 7.2;
  double Z = 0.0;
  double A = 25.0;
  int nMax = 4;
  double SW[2][84] = {};
  double SDW[462] = {};
  double spin = 1.5;
  double beta2 = 0.4;
  double beta4 = 0.05;
  bool report = true;

  NO::Calculate(spin, beta2, beta4, V0, R, A0, VS, A, Z, nMax, true);

  double nu = 1.;
  cout << "Adv: " << CD::GetRadialMEHO(1, 0, 2, 1, 0, nu) << endl;;
  cout << "Normal: " << std::pow(CD::GetRMSHO(1, 0, nu), 2.) << endl;
}

int main(int argc, char** argv) {
  ShowIntro();

  //TestNilssonMethods();
  OptionContainer::GetInstance(argc, argv);

  if (OptExists(input)) {
    Generator* gen = new Generator();
    gen->CalculateSpectrum();
    gen->WriteSpectrumToFile();
    delete gen;
  }

  return 0;
}
