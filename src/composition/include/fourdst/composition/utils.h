#pragma once

#include "fourdst/composition/composition.h"
#include "fourdst/atomic/atomicSpecies.h"

#include <vector>

namespace fourdst::composition {
    Composition buildCompositionFromMassFractions(
        const std::vector<std::string>& symbols,
        const std::vector<double>& massFractions
    );

    Composition buildCompositionFromMassFractions(
        const std::vector<atomic::Species>& species,
        const std::vector<double>& massFractions
    );

    Composition buildCompositionFromMassFractions(
        const std::set<atomic::Species>& species,
        const std::vector<double>& massFractions
    );
}