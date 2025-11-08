#pragma once

#include "fourdst/composition/composition.h"
#include "fourdst/atomic/atomicSpecies.h"

#include <vector>

namespace fourdst::composition {
    /**
     * @brief Build a Composition object from symbols and their corresponding mass fractions.
     * @param symbols The symbols to register.
     * @param massFractions The corresponding mass fractions for each symbol.
     * @return A Composition object constructed from the provided symbols and mass fractions.
     * @throws exceptions::UnknownSymbolError if any symbol is invalid. Symbols are invalid if they are not registered at compile time in the atomic species database (`fourdst/atomic/species.h`).
     * @throws exceptions::InvalidCompositionError if the provided mass fractions do not sum to within one part in 10^10 of 1.0.
     * @throws exceptions::InvalidCompositionError if the number of symbols does not match the number of mass fractions.
     */
    Composition buildCompositionFromMassFractions(
        const std::vector<std::string>& symbols,
        const std::vector<double>& massFractions
    );

    /**
     * @brief Build a Composition object from species and their corresponding mass fractions.
     * @param species The species to register.
     * @param massFractions The corresponding mass fractions for each species.
     * @return A Composition object constructed from the provided species and mass fractions.
     * @throws exceptions::InvalidCompositionError if the provided mass fractions do not sum to within one part in 10^10 of 1.0.
     * @throws exceptions::InvalidCompositionError if the number of species does not match the number of mass fractions.
     */
    Composition buildCompositionFromMassFractions(
        const std::vector<atomic::Species>& species,
        const std::vector<double>& massFractions
    );

    /**
     * @brief Build a Composition object from species in a set and their corresponding mass fractions.
     * @param species The species to register.
     * @param massFractions The corresponding mass fractions for each species.
     * @return A Composition object constructed from the provided species and mass fractions.
     * @throws exceptions::InvalidCompositionError if the provided mass fractions do not sum to within one part in 10^10 of 1.0.
     * @throws exceptions::InvalidCompositionError if the number of species does not match the number of mass fractions.
     *
     * @note This is the version of the function which the other overloads ultimately call.
     */
    Composition buildCompositionFromMassFractions(
        const std::set<atomic::Species>& species,
        const std::vector<double>& massFractions
    );
}