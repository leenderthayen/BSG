#include <iostream>

#include "Generator.h"
#include "OptionContainer.h"

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

int main(int argc, char** argv) {
  ShowIntro();
  OptionContainer::GetInstance(argc, argv);

  if (OptExists(input)) {
    Generator* gen = new Generator();
    gen->CalculateSpectrum();
    gen->WriteSpectrumToFile();
    delete gen;
  }

  return 0;
}
