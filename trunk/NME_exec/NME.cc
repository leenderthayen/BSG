#include <iostream>
#include <string>

#include "NuclearStructureManager.h"
#include "NMEOptions.h"

#include "boost/algorithm/string.hpp"

using std::cout;
using std::endl;

int main(int argc, char** argv) {
  NMEOptions::GetInstance(argc, argv);

  if (NMOptExists(input)) {
    NuclearStructure::NuclearStructureManager* nsm = new NuclearStructure::NuclearStructureManager();
    if (NMOptExists(weakmagnetism)) {
      cout << "b/Ac: " << nsm->CalculateWeakMagnetism() << endl;
    }
    if (NMOptExists(inducedtensor)) {
      cout << "d/Ac: " << nsm->CalculateInducedTensor() << endl;
    }
    if (NMOptExists(matrixelement)) {
      std::string me = GetNMOpt(std::string, matrixelement);
      bool V = (me[0] == 'V');
      cout << me[0] << "M" << me.substr(1, 3) << ": " << nsm->CalculateMatrixElement(V, (int)(me[1]-'0'), (int)(me[2]-'0'), (int)(me[3]-'0')) << endl;
    }
  }

  return 0;
}
