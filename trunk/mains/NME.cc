#include <iostream>

#include "NuclearStructureManager.h"
#include "OptionContainer.h"

using std::cout;
using std::endl;
using std::cerr;

char* version = "0.10";
char* lastUpdate = "February 24, 2017";
char* author = "L. Hayen (leendert.hayen@kuleuven.be)";

void ShowIntro() {
  cout << "**************************************************\n";
  cout << "                  NME v." << version << "\n\n";
  cout << "        Last update: " << lastUpdate << endl;
  cout << "      " << author << endl;
  cout << "**************************************************\n\n";
}


int main(int argc, char** argv) {
  ShowIntro();

  /*OptionContainer::GetInstance(argc, argv);

  if (OptExists(input)) {
    Generator* gen = new Generator();
    gen->CalculateSpectrum();
    gen->WriteSpectrumToFile();
    delete gen;
  }*/

  return 0;
}
