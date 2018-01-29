#include "Generator.h"
#include "OptionContainer.h"
#include <iostream>
#include <chrono>

int main(int argc, char** argv) {
  auto start = std::chrono::steady_clock::now();

  OptionContainer::GetInstance(argc, argv);

  auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);

  std::cout << "milliseconds since start: " << elapsed.count() << "\n";

  if (OptExists(input)) {
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "milliseconds since start: " << elapsed.count() << "\n";
    Generator* gen = new Generator();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "milliseconds since start: " << elapsed.count() << "\n";
    gen->CalculateSpectrum();
    elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - start);
    std::cout << "milliseconds since start: " << elapsed.count() << "\n";
    delete gen;
  }

  return 0;
}
