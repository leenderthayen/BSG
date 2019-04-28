#include "Generator.h"
#include "BSGOptionContainer.h"
#include <iostream>
#include <chrono>
#include <string>

int main(int argc, char** argv) {
  bsg::BSGOptionContainer::GetInstance(argc, argv);

  if (OptExists(input)) {
    bsg::Generator* gen = new bsg::Generator();
    gen->CalculateSpectrum();
    delete gen;
  }

  return 0;
}
