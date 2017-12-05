#include "Generator.h"
#include "OptionContainer.h"

int main(int argc, char** argv) {
  OptionContainer::GetInstance(argc, argv);

  if (OptExists(input)) {
    Generator* gen = new Generator();
    gen->CalculateSpectrum();
    delete gen;
  }

  return 0;
}
