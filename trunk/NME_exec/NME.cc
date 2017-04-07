#include <iostream>
#include <string>

#include "NuclearStructureManager.h"
#include "NMEOptions.h"

#include "boost/algorithm/string.hpp"

using std::cout;
using std::endl;
using std::cerr;

std::string version = "0.10";
std::string lastUpdate = "February 24, 2017";
std::string author = "L. Hayen (leendert.hayen@kuleuven.be)";

void ShowIntro() {
  cout << "**************************************************\n";
  cout << "                  NME v." << version << "\n\n";
  cout << "        Last update: " << lastUpdate << endl;
  cout << "      " << author << endl;
  cout << "**************************************************\n\n";
}

int main(int argc, char** argv) {
  ShowIntro();

  NMEOptions::GetInstance(argc, argv);

  if (OptExists(input)) {
    NuclearStructure::NuclearStructureManager* nsm = new NuclearStructure::NuclearStructureManager();
    if (OptExists(weakmagnetism)) {
      cout << "b/Ac: " << nsm->CalculateWeakMagnetism() << endl;
    }
    if (OptExists(inducedtensor)) {
      cout << "d/Ac: " << nsm->CalculateInducedTensor() << endl;
    }
    if (OptExists(matrixelement)) {
      std::string me = GetOpt(std::string, matrixelement);
      bool V = (me[0] == 'V');
      cout << me[0] << "M" << me.substr(1, 3) << ": " << nsm->CalculateMatrixElement(V, (int)(me[1]-'0'), (int)(me[2]-'0'), (int)(me[3]-'0')) << endl;
    }
  }

  return 0;
}
