#include <iostream>

#include "Generator.h"
#include "OptionContainer.h"
#include "NilssonOrbits.h"

using std::cout;
using std::endl;
using std::cerr;

double version = 0.9;
char* lastUpdate = "February 16, 2017";
char* author = "L. Hayen (leendert.hayen@kuleuven.be)";

void ShowIntro() {
  cout << "**************************************************\n";
  cout << "                  BDSM v." << version << "\n\n";
  cout << "        Last update: " << lastUpdate << endl;
  cout << "      " << author << endl;
  cout << "**************************************************\n\n";
}

void TestNilssonMethods() {
  cout << "Testing Nilsson methods" << endl;
  double V0 = 53.0;
  double R = 3.655;
  double A0 = 0.650;
  double VS = 7.2;
  double Z = 0.0;
  double A = 25.0;
  int nMax = 4;
  double SW[2][84] = {};
  double SDW[462] = {};

  nilsson::WoodsSaxon(V0, R, A0, VS, A, Z, nMax, SW, SDW);
}

int main(int argc, char** argv) {
  ShowIntro();

  TestNilssonMethods();
  OptionContainer::GetInstance(argc, argv);

  if (OptExists(input)) {
    Generator* gen = new Generator();
    gen->CalculateSpectrum();
    gen->WriteSpectrumToFile();
    delete gen;
  }

  return 0;
}
