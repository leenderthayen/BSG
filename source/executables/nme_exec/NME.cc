#include <iostream>
#include <string>

#include "NuclearStructureManager.h"
#include "NMEOptions.h"

#include "boost/algorithm/string.hpp"

using std::cout;
using std::endl;

int main(int argc, char** argv) {
  NMEOptions::GetInstance(argc, argv);

  if (NMEOptExists(input)) {
    NuclearStructure::NuclearStructureManager* nsm = new NuclearStructure::NuclearStructureManager();
    if (NMEOptExists(weakmagnetism)) {
      cout << "b/Ac: " << nsm->CalculateWeakMagnetism() << endl;
    }
    if (NMEOptExists(inducedtensor)) {
      cout << "d/Ac: " << nsm->CalculateInducedTensor() << endl;
    }
    if (NMEOptExists(matrixelement)) {
      std::string me = GetNMEOpt(std::string, matrixelement);
      bool V = (me[0] == 'V');
      cout << me[0] << "M" << me.substr(1, 3) << ": " << nsm->CalculateReducedMatrixElement(V, (int)(me[1]-'0'), (int)(me[2]-'0'), (int)(me[3]-'0')) << endl;
    }
  }

  return 0;
}
