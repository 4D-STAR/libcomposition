#include "fourdst/composition/composition.h"
#include "fourdst/composition/exceptions/exceptions_composition.h"
#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils.h"

#include <ranges>
#include <vector>
#include <set>
#include <string>

namespace {
    std::optional<fourdst::atomic::Species> getSpecies(const std::string& symbol) {
        if (!fourdst::atomic::species.contains(symbol)) {
            return std::nullopt;
        }
        return fourdst::atomic::species.at(symbol);
    }

    void throw_unknown_symbol(quill::Logger* logger, const std::string& symbol) {
        throw fourdst::composition::exceptions::UnknownSymbolError("Symbol " + symbol + " is not a valid species symbol (not in the species database)");
    }
}

namespace fourdst::composition {
    Composition buildCompositionFromMassFractions(
        const std::set<atomic::Species> &species,
        const std::vector<double> &massFractions
    ) {
        Composition composition;

        for (const auto& [sp, xi] : std::views::zip(species, massFractions)) {
            composition.registerSpecies(sp);
            composition.setMolarAbundance(sp, xi/sp.mass());
        }

        return composition;
    }

    Composition buildCompositionFromMassFractions(const std::vector<atomic::Species> &species, const std::vector<double> &massFractions) {
        return buildCompositionFromMassFractions(std::set<atomic::Species>(species.begin(), species.end()), massFractions);
    }

    Composition buildCompositionFromMassFractions(const std::vector<std::string> &symbols, const std::vector<double> &massFractions) {
        std::set<atomic::Species> species;
        for (const auto& symbol : symbols) {
            auto result = getSpecies(symbol);
            if (!result) {
                throw_unknown_symbol(nullptr, symbol);
            }
            species.insert(result.value());
        }
        return buildCompositionFromMassFractions(species, massFractions);
    }



}