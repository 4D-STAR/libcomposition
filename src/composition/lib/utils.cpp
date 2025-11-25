#include "fourdst/composition/composition.h"
#include "fourdst/composition/exceptions/exceptions_composition.h"
#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/atomic/species.h"
#include "fourdst/composition/utils.h"
#include "fourdst/logging/logging.h"

#include <numeric>
#include <ranges>
#include <vector>
#include <set>
#include <string>

#include "quill/LogMacros.h"

namespace {
    std::optional<fourdst::atomic::Species> getSpecies(const std::string& symbol) {
        if (!fourdst::atomic::species.contains(symbol)) {
            return std::nullopt;
        }
        return fourdst::atomic::species.at(symbol);
    }

    quill::Logger* getLogger() {
        static quill::Logger* logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        return logger;
    }

    void throw_unknown_symbol(const std::string& symbol) {
        LOG_ERROR(getLogger(), "Symbol {} is not a valid species symbol (not in the species database)", symbol);
        throw fourdst::composition::exceptions::UnknownSymbolError("Symbol " + symbol + " is not a valid species symbol (not in the species database)");
    }
}

namespace fourdst::composition {
    Composition buildCompositionFromMassFractions(
        const std::set<atomic::Species> &species,
        const std::vector<double> &massFractions
    ) {
        const double sum = std::accumulate(
            massFractions.begin(),
            massFractions.end(),
            0.0
        );

        if (std::abs(sum - 1.0) > 1e-10) {
            throw exceptions::InvalidCompositionError(
                "Mass fractions must sum to 1.0, got " +  std::to_string(sum)
            );
        }

        if (species.size() != massFractions.size()) {
            throw exceptions::InvalidCompositionError(
                "The number of species and mass fractions must be equal. Got " +
                std::to_string(species.size()) + " species and " +
                std::to_string(massFractions.size()) + " mass fractions."
            );
        }

        Composition composition;

        for (const auto& [sp, xi] : std::views::zip(species, massFractions)) {
            composition.registerSpecies(sp);
            composition.setMolarAbundance(sp, xi/sp.mass());
        }

        return composition;
    }

    Composition buildCompositionFromMassFractions(const std::vector<atomic::Species> &species, const std::vector<double> &massFractions) {
        std::set<atomic::Species> speciesSet(species.begin(), species.end());
        std::vector<double> sortedMassFractions;

        sortedMassFractions.resize(massFractions.size());
        for (const auto& [s, xi] : std::views::zip(species, massFractions)) {
            const size_t index = std::distance(speciesSet.begin(), speciesSet.find(s));
            assert (index < sortedMassFractions.size());
            sortedMassFractions[index] = xi;
        }

        return buildCompositionFromMassFractions(speciesSet, sortedMassFractions);
    }

    Composition buildCompositionFromMassFractions(const std::vector<std::string> &symbols, const std::vector<double> &massFractions) {
        std::set<atomic::Species> species;
        for (const auto& symbol : symbols) {
            auto result = getSpecies(symbol);
            if (!result) {
                throw_unknown_symbol(symbol);
            }
            species.insert(result.value());
        }

        std::vector<double> sortedMassFractions(massFractions.size());
        for (const auto& [symbol, xi] : std::views::zip(symbols, massFractions)) {
            auto result = getSpecies(symbol);
            if (!result) {
                throw_unknown_symbol(symbol);
            }
            const size_t index = std::distance(species.begin(), species.find(result.value()));
            assert (index < sortedMassFractions.size());
            sortedMassFractions[index] = xi;
        }
        return buildCompositionFromMassFractions(species, sortedMassFractions);
    }

    Composition buildCompositionFromMassFractions(const std::unordered_map<atomic::Species, double>& massFractionsMap) {
        std::set<atomic::Species> species;
        std::vector<double> massFractions;

        massFractions.reserve(massFractionsMap.size());

        for (const auto &sp: massFractionsMap | std::views::keys) {
            species.insert(sp);
        }

        massFractions.resize(massFractionsMap.size());
        for (const auto& [sp, xi] : massFractionsMap) {
            const size_t index = std::distance(species.begin(), species.find(sp));
            assert (index < massFractions.size());
            massFractions[index] = xi;
        }

        return buildCompositionFromMassFractions(species, massFractions);
    }

    Composition buildCompositionFromMassFractions(std::map<atomic::Species, double> massFractions) {
        std::set<atomic::Species> species;
        std::vector<double> massFractionVector;

        massFractionVector.reserve(massFractions.size());

        for (const auto& [sp, xi] : massFractions) {
            species.insert(sp);
            massFractionVector.push_back(xi);
        }

        return buildCompositionFromMassFractions(species, massFractionVector);
    }

    Composition buildCompositionFromMassFractions(std::map<std::string, double> massFractions) {
        std::set<atomic::Species> species;
        std::vector<double> massFractionVector;


        for (const auto &symbol: massFractions | std::views::keys) {
            auto result = getSpecies(symbol);
            if (!result) {
                throw_unknown_symbol(symbol);
            }
            species.insert(result.value());
        }

        massFractionVector.resize(massFractions.size());

        for (const auto& [symbol, xi] : massFractions) {
            auto result = getSpecies(symbol);
            if (!result) {
                throw_unknown_symbol(symbol);
            }
            const size_t index = std::distance(species.begin(), species.find(result.value()));
            assert (index < massFractionVector.size());
            massFractionVector[index] = xi;
        }


        return buildCompositionFromMassFractions(species, massFractionVector);
    }

    Composition buildCompositionFromMassFractions(const std::unordered_map<std::string, double>& massFractions) {
        std::set<atomic::Species> species;
        std::vector<double> massFractionVector;

        for (const auto &symbol: massFractions | std::views::keys) {
            auto result = getSpecies(symbol);
            if (!result) {
                throw_unknown_symbol(symbol);
            }
            species.insert(result.value());
        }

        massFractionVector.resize(massFractions.size());
        for (const auto& [sp, xi] : massFractions) {
            auto result = getSpecies(sp);
            if (!result) {
                throw_unknown_symbol(sp);
            }
            const size_t index = std::distance(species.begin(), species.find(result.value()));
            assert (index < massFractionVector.size());
            massFractionVector[index] = xi;
        }

        return buildCompositionFromMassFractions(species, massFractionVector);
    }

}