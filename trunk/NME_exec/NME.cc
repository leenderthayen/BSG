#include <iostream>
#include <string>

#include "NuclearStructureManager.h"
#include "OptionContainer.h"

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

  OptionContainer::GetInstance(argc, argv);

  if (OptExists(input)) {
    NuclearStructureManager* nsm = new NuclearStructureManager();
  }

  return 0;
}
