#include <iostream>

#include "Generator.h"
#include "OptionContainer.h"

using std::cout;
using std::endl;
using std::cerr;

double version = 0.9;

void ShowIntro() {
  cout << "*************************************\n";
  cout << "            BDSM v." << version << "\n";
  cout << "*************************************\n";
}

int main(int argc, char** argv) {
  ShowIntro();
  OptionContainer::GetInstance(argc, argv);

  return 0;
}
