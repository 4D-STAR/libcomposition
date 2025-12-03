#include "fourdst/composition/composition.h"
#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/atomic/species.h"

#include <unordered_map>

#include <iostream>

int main() {

  std::unordered_map<fourdst::atomic::Species, double> abundances({{fourdst::atomic::H_1, 0.7}, {fourdst::atomic::He_4, 0.28}, {fourdst::atomic::C_12, 0.02}});

  fourdst::composition::Composition comp(abundances);

  std::cout << comp << std::endl;



}
