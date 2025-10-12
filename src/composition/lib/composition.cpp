/* ***********************************************************************
//
//   Copyright (C) 2025 -- The 4D-STAR Collaboration
//   File Author: Emily Boudreaux
//   Last Modified: October 6, 2025
//
//   4DSSE is free software; you can use it and/or modify
//   it under the terms and restrictions the GNU General Library Public
//   License version 3 (GPLv3) as published by the Free Software Foundation.
//
//   4DSSE is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
//   See the GNU Library General Public License for more details.
//
//   You should have received a copy of the GNU Library General Public License
//   along with this software; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
// *********************************************************************** */
#include "quill/LogMacros.h"

#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <array>
#include <ranges>
#include <algorithm>


#include <utility>

#include "fourdst/composition/atomicSpecies.h"
#include "fourdst/composition/species.h"
#include "fourdst/composition/composition.h"

#include <numeric>

#include "fourdst/composition/exceptions/exceptions_composition.h"

namespace {
    template<typename A, typename B>
    std::vector<A> sortVectorBy(
        std::vector<A> toSort,
        const std::vector<B>& by
    ) {
        std::vector<std::size_t> indices(by.size());
        for (size_t i = 0; i < indices.size(); i++) {
            indices[i] = i;
        }

        std::ranges::sort(indices, [&](size_t a, size_t b) {
            return by[a] < by[b];
        });

        std::vector<A> sorted;
        sorted.reserve(indices.size());

        for (const auto idx: indices) {
            sorted.push_back(toSort[idx]);
        }

        return sorted;
    }
}

namespace fourdst::composition {

    CompositionEntry::CompositionEntry() :
    m_symbol("H-1"),
    m_isotope(fourdst::atomic::species.at("H-1")),
    m_initialized(false) {}

    CompositionEntry::CompositionEntry(
        const std::string& symbol,
        const bool massFracMode
    ) :
    m_symbol(symbol),
    m_isotope(atomic::species.at(symbol)),
    m_massFracMode(massFracMode) {
        setSpecies(symbol);
    }

    CompositionEntry::CompositionEntry(const CompositionEntry& entry) = default;

    void CompositionEntry::setSpecies(const std::string& symbol) {
        if (m_initialized) {
            throw exceptions::EntryAlreadyInitializedError("Composition entry is already initialized.");
        }
        if (!fourdst::atomic::species.contains(symbol)) {
            throw exceptions::InvalidSpeciesSymbolError("Invalid symbol.");
        }
        m_symbol = symbol;
        m_isotope = atomic::species.at(symbol);
        m_initialized = true;
    }

    std::string CompositionEntry::symbol() const {
        return m_symbol;
    }

    double CompositionEntry::mass_fraction() const {
        if (!m_massFracMode) {
            throw exceptions::CompositionModeError("Composition entry is in number fraction mode.");
        }
        // X_i = (moles_i / mass_total) * (mass_i / moles_i) = m_molesPerMass * A_i
        return m_molesPerMass * m_isotope.mass();
    }

    double CompositionEntry::number_fraction() const {
        if (m_massFracMode) {
            throw exceptions::CompositionModeError("Composition entry is in mass fraction mode.");
        }
        // In number fraction mode, the value is cached during the mode switch.
        return m_cachedNumberFraction;
    }

    double CompositionEntry::number_fraction(
        const double totalMolesPerMass
    ) const {
        // n_i = (moles_i / mass_total) / (moles_total / mass_total)
        if (totalMolesPerMass == 0.0) return 0.0;
        return m_molesPerMass / totalMolesPerMass;
    }

    double CompositionEntry::rel_abundance() const {
        return m_molesPerMass;
    }

    atomic::Species CompositionEntry::isotope() const {
        return m_isotope;
    }

    void CompositionEntry::setMassFraction(
        const double mass_fraction
    ) {
        if (!m_massFracMode) {
            throw exceptions::CompositionModeError("Composition entry is in number fraction mode.");
        }
        // Set the invariant from the given mass fraction
        if (m_isotope.mass() == 0.0) {
             m_molesPerMass = 0.0;
        } else {
             m_molesPerMass = mass_fraction / m_isotope.mass();
        }
    }

    void CompositionEntry::setNumberFraction(
        const double number_fraction
    ) {
        if (m_massFracMode) {
            throw exceptions::CompositionModeError("Composition entry is in mass fraction mode.");
        }
        // In number fraction mode, we only cache the value. The invariant
        // m_molesPerMass cannot be calculated until finalize() provides global context.
        m_cachedNumberFraction = number_fraction;
    }

    bool CompositionEntry::setMassFracMode(
        [[maybe_unused]] const double meanMolarMass
    ) {
        if (m_massFracMode) {
            return false;
        }
        m_massFracMode = true;
        // The invariant m_molesPerMass does not change when switching mode.
        // The cached number fraction is now stale, but that's okay.
        return true;
    }

    bool CompositionEntry::setNumberFracMode(
        const double totalMolesPerMass
    ) {
        if (!m_massFracMode) {
            return false;
        }
        m_massFracMode = false;
        // Calculate and cache the number fraction for the new mode.
        m_cachedNumberFraction = number_fraction(totalMolesPerMass);
        return true;
    }

    bool CompositionEntry::getMassFracMode() const {
        return m_massFracMode;
    }

        Composition::Composition(
        const std::vector<std::string>& symbols
    ) {
        for (const auto& symbol : symbols) {
            registerSymbol(symbol);
        }
    }

    Composition::Composition(
        const std::set<std::string>& symbols
    ) {
        for (const auto& symbol : symbols) {
            registerSymbol(symbol);
        }
    }

    Composition::Composition(
        const std::vector<std::string>& symbols,
        const std::vector<double>& fractions,
        const bool massFracMode
    ) : m_massFracMode(massFracMode) {
        if (symbols.size() != fractions.size()) {
            LOG_CRITICAL(m_logger, "The number of symbols and fractions must be equal (got {} symbols and {} fractions).", symbols.size(), fractions.size());
            throw exceptions::InvalidCompositionError("The number of symbols and fractions must be equal. Got " + std::to_string(symbols.size()) + " symbols and " + std::to_string(fractions.size()) + " fractions.");
        }

        validateComposition(fractions);

        for (const auto &symbol : symbols) {
            registerSymbol(symbol, m_massFracMode);
        }

        for (size_t i = 0; i < symbols.size(); ++i) {
            if (m_massFracMode) {
                setMassFraction(symbols[i], fractions[i]);
            } else {
                setNumberFraction(symbols[i], fractions[i]);
            }
        }
        if (const bool didFinalize = finalize(); !didFinalize) {
            std::string msg = "Failed to finalize composition on construction. ";
            msg += "Construction of a composition object requires that the sum of the fractions vector be 1.\n";
            LOG_CRITICAL(m_logger, "{}", msg);
            throw exceptions::InvalidCompositionError(msg);
        }
    }

    Composition::Composition(const Composition &composition) {
        m_finalized = composition.m_finalized;
        m_specificNumberDensity = composition.m_specificNumberDensity;
        m_meanParticleMass = composition.m_meanParticleMass;
        m_massFracMode = composition.m_massFracMode;
        m_registeredSymbols = composition.m_registeredSymbols;
        m_compositions = composition.m_compositions;
    }

    Composition& Composition::operator=(const Composition &other) {
        if (this != &other) {
            m_finalized               = other.m_finalized;
            m_specificNumberDensity   = other.m_specificNumberDensity;
            m_meanParticleMass        = other.m_meanParticleMass;
            m_massFracMode            = other.m_massFracMode;
            m_registeredSymbols       = other.m_registeredSymbols;
            m_compositions            = other.m_compositions;
        }
        return *this;
    }

    void Composition::registerSymbol(
        const std::string& symbol,
        const bool massFracMode
    ) {
        if (!isValidSymbol(symbol)) {
            LOG_ERROR(m_logger, "Invalid symbol: {}", symbol);
            throw exceptions::InvalidSymbolError("Invalid symbol: " + symbol);
        }

        if (m_registeredSymbols.empty()) {
            m_massFracMode = massFracMode;
        } else {
            if (m_massFracMode != massFracMode) {
                LOG_ERROR(m_logger, "Composition is in {} fraction mode. Cannot register symbol ({}) in {} fraction mode.", m_massFracMode ? "mass" : "number", symbol, massFracMode ? "mass" : "number");
                throw exceptions::CompositionModeError("Composition mode mismatch.");
            }
        }

        if (m_registeredSymbols.contains(symbol)) {
            LOG_WARNING(m_logger, "Symbol {} is already registered.", symbol);
            return;
        }

        m_registeredSymbols.insert(symbol);
        m_compositions[symbol] = CompositionEntry(symbol, m_massFracMode);
        m_finalized = false;
        LOG_TRACE_L3(m_logger, "Registered symbol: {}", symbol);
    }

    void Composition::registerSymbol(
        const std::vector<std::string>& symbols,
        const bool massFracMode
    ) {
        for (const auto& symbol : symbols) {
            registerSymbol(symbol, massFracMode);
        }
    }

    void Composition::registerSpecies(
        const atomic::Species &species,
        const bool massFracMode
    ) {
        registerSymbol(std::string(species.name()), massFracMode);
    }

    void Composition::registerSpecies(
        const std::vector<atomic::Species> &species,
        const bool massFracMode
    ) {
        for (const auto& s : species) {
            registerSpecies(s, massFracMode);
        }
    }

    std::set<std::string> Composition::getRegisteredSymbols() const {
        return m_registeredSymbols;
    }

    std::set<atomic::Species> Composition::getRegisteredSpecies() const {
        std::set<atomic::Species> result;
        for (const auto& entry : m_compositions | std::views::values) {
            result.insert(entry.isotope());
        }
        return result;
    }

    bool Composition::isValidSymbol(
        const std::string& symbol
    ) {
        return atomic::species.contains(symbol);
    }

    void Composition::validateComposition(const std::vector<double>& fractions) const {
        if (!isValidComposition(fractions)) {
            LOG_ERROR(m_logger, "Invalid composition.");
            throw exceptions::InvalidCompositionError("Invalid composition.");
        }
    }

    bool Composition::isValidComposition(const std::vector<double>& fractions) const {
        const double sum = std::accumulate(fractions.begin(), fractions.end(), 0.0);
        if (sum < 0.999999 || sum > 1.000001) {
            LOG_ERROR(m_logger, "The sum of fractions must be equal to 1 (expected 1, got {}).", sum);
            return false;
        }
        return true;
    }

    double Composition::setMassFraction(const std::string& symbol, const double& mass_fraction) {
        if (!m_registeredSymbols.contains(symbol)) {
            LOG_ERROR(m_logger, "Symbol {} is not registered.", symbol);
            throw exceptions::UnregisteredSymbolError("Symbol (" + symbol + ") is not registered.");
        }
        if (!m_massFracMode) {
            LOG_ERROR(m_logger, "Composition is in number fraction mode.");
            throw exceptions::CompositionModeError("Composition is in number fraction mode.");
        }
        if (mass_fraction < 0.0 || mass_fraction > 1.0) {
            LOG_ERROR(m_logger, "Mass fraction must be between 0 and 1 for symbol {}. Currently it is {}.", symbol, mass_fraction);
            throw exceptions::InvalidCompositionError("Mass fraction must be between 0 and 1.");
        }
        m_finalized = false;
        const double old_mass_fraction = m_compositions.at(symbol).mass_fraction();
        m_compositions.at(symbol).setMassFraction(mass_fraction);
        return old_mass_fraction;
    }

    std::vector<double> Composition::setMassFraction(const std::vector<std::string>& symbols, const std::vector<double>& mass_fractions) {
        if (symbols.size() != mass_fractions.size()) {
            throw exceptions::InvalidCompositionError("The number of symbols and mass fractions must be equal.");
        }
        std::vector<double> old_mass_fractions;
        old_mass_fractions.reserve(symbols.size());
        for (size_t i = 0; i < symbols.size(); ++i) {
            old_mass_fractions.push_back(setMassFraction(symbols[i], mass_fractions[i]));
        }
        return old_mass_fractions;
    }

    double Composition::setNumberFraction(
        const std::string& symbol,
        const double& number_fraction
    ) {
        if (!m_registeredSymbols.contains(symbol)) {
            LOG_ERROR(m_logger, "Symbol {} is not registered.", symbol);
            throw exceptions::UnregisteredSymbolError("Symbol (" + symbol + ") is not registered.");
        }
        if (m_massFracMode) {
            LOG_ERROR(m_logger, "Composition is in mass fraction mode.");
            throw exceptions::CompositionModeError("Composition is in mass fraction mode.");
        }
        if (number_fraction < 0.0 || number_fraction > 1.0) {
            LOG_ERROR(m_logger, "Number fraction must be between 0 and 1 for symbol {}. Currently it is {}.", symbol, number_fraction);
            throw exceptions::InvalidCompositionError("Number fraction must be between 0 and 1.");
        }
        m_finalized = false;
        const double old_number_fraction = m_compositions.at(symbol).number_fraction();
        m_compositions.at(symbol).setNumberFraction(number_fraction);
        return old_number_fraction;
    }

    std::vector<double> Composition::setNumberFraction(
        const std::vector<std::string>& symbols,
        const std::vector<double>& number_fractions
    ) {
        if (symbols.size() != number_fractions.size()) {
            throw exceptions::InvalidCompositionError("The number of symbols and number fractions must be equal.");
        }
        std::vector<double> old_number_fractions;
        old_number_fractions.reserve(symbols.size());
        for (size_t i = 0; i < symbols.size(); ++i) {
            old_number_fractions.push_back(setNumberFraction(symbols[i], number_fractions[i]));
        }
        return old_number_fractions;
    }

    double Composition::setMassFraction(
        const atomic::Species &species,
        const double &mass_fraction
    ) {
        return setMassFraction(std::string(species.name()), mass_fraction);
    }

    std::vector<double> Composition::setMassFraction(
        const std::vector<atomic::Species> &species,
        const std::vector<double> &mass_fractions
    ) {
        std::vector<std::string> symbols;
        symbols.reserve(species.size());
        for(const auto& s : species) symbols.emplace_back(s.name());
        return setMassFraction(symbols, mass_fractions);
    }

    double Composition::setNumberFraction(
        const atomic::Species &species,
        const double &number_fraction
    ) {
        return setNumberFraction(std::string(species.name()), number_fraction);
    }

    std::vector<double> Composition::setNumberFraction(
        const std::vector<atomic::Species> &species,
        const std::vector<double> &number_fractions
    ) {
        std::vector<std::string> symbols;
        symbols.reserve(species.size());
        for(const auto& s : species) symbols.push_back(std::string(s.name()));
        return setNumberFraction(symbols, number_fractions);
    }

    bool Composition::finalize(const bool norm) {
        m_specificNumberDensity = 0.0;
        m_meanParticleMass = 0.0;
        m_finalized = m_massFracMode ? finalizeMassFracMode(norm) : finalizeNumberFracMode(norm);
        m_cache.clear();
        return m_finalized;
    }

    bool Composition::finalizeMassFracMode(const bool norm) {
        std::vector<double> mass_fractions;
        mass_fractions.reserve(m_compositions.size());
        for (const auto &entry: m_compositions | std::views::values) {
            mass_fractions.push_back(entry.mass_fraction());
        }

        double sum = std::accumulate(mass_fractions.begin(), mass_fractions.end(), 0.0);
        if (norm && sum > 0) {
            for (auto& [symbol, entry] : m_compositions) {
                setMassFraction(symbol, entry.mass_fraction() / sum);
            }
            // Recalculate fractions vector after normalization for validation
            mass_fractions.clear();
            for (const auto &entry: m_compositions | std::views::values) {
                mass_fractions.push_back(entry.mass_fraction());
            }
        }

        try {
            validateComposition(mass_fractions);
        } catch ([[maybe_unused]] const exceptions::InvalidCompositionError& e) {
            LOG_ERROR(m_logger, "Composition is invalid after mass frac finalization (Total mass {}).", sum);
            return false;
        }

        for (const auto &entry: m_compositions | std::views::values) {
            m_specificNumberDensity += entry.rel_abundance(); // rel_abundance is now consistently moles/mass
        }

        if (m_specificNumberDensity > 0) {
            m_meanParticleMass = 1.0 / m_specificNumberDensity;
        }
        return true;
    }

    bool Composition::finalizeNumberFracMode(const bool norm) {
        std::vector<double> number_fractions;
        number_fractions.reserve(m_compositions.size());
        for (const auto &entry: m_compositions | std::views::values) {
            number_fractions.push_back(entry.number_fraction());
        }

        double sum = std::accumulate(number_fractions.begin(), number_fractions.end(), 0.0);
        if (norm && sum > 0) {
            for (auto& [symbol, entry] : m_compositions) {
                setNumberFraction(symbol, entry.number_fraction() / sum);
            }
            // Recalculate fractions vector after normalization for validation
            number_fractions.clear();
            for (const auto &entry: m_compositions | std::views::values) {
                number_fractions.push_back(entry.number_fraction());
            }
        }

        try {
            validateComposition(number_fractions);
        } catch ([[maybe_unused]] const exceptions::InvalidCompositionError& e) {
            LOG_ERROR(m_logger, "Composition is invalid after number frac finalization (Total number frac {}).", sum);
            return false;
        }

        // Calculate mean particle mass <A> = sum(n_i * A_i)
        for (const auto &entry: m_compositions | std::views::values) {
            m_meanParticleMass += entry.number_fraction() * entry.isotope().mass();
        }

        for (auto &entry: m_compositions | std::views::values) {
            const double X_i = (m_meanParticleMass > 0) ? (entry.number_fraction() * entry.isotope().mass() / m_meanParticleMass) : 0.0;
            entry.m_massFracMode = true;
            entry.setMassFraction(X_i);
            entry.m_massFracMode = false;
        }

        if (m_meanParticleMass > 0) {
            m_specificNumberDensity = 1.0 / m_meanParticleMass;
        }
        return true;
    }

    Composition Composition::mix(const Composition& other, const double fraction) const {
        if (!m_finalized || !other.m_finalized) {
            LOG_ERROR(m_logger, "Compositions have not both been finalized. Hint: Consider running .finalize() on both compositions before mixing.");
            throw exceptions::CompositionNotFinalizedError("Compositions have not been finalized (Hint: Consider running .finalize() on both compositions before mixing).");
        }

        if (fraction < 0.0 || fraction > 1.0) {
            LOG_ERROR(m_logger, "Mixing fraction must be between 0 and 1. Currently it is {}.", fraction);
            throw exceptions::InvalidCompositionError("Mixing fraction must be between 0 and 1. Currently it is " + std::to_string(fraction) + ".");
        }

        std::set<std::string> mixedSymbols = other.getRegisteredSymbols();
        // Get the union of the two sets of symbols to ensure all species are included in the new composition.
        mixedSymbols.insert(m_registeredSymbols.begin(), m_registeredSymbols.end());

        Composition mixedComposition(mixedSymbols);
        for (const auto& symbol : mixedSymbols) {
            double otherMassFrac = 0.0;

            const double thisMassFrac = hasSymbol(symbol) ? getMassFraction(symbol) : 0.0;
            otherMassFrac = other.hasSymbol(symbol) ? other.getMassFraction(symbol) : 0.0;

            // The mixing formula is a linear interpolation of mass fractions.
            double massFraction = fraction * thisMassFrac + otherMassFrac * (1-fraction);
            mixedComposition.setMassFraction(symbol, massFraction);
        }
        if (const bool didFinalize = mixedComposition.finalize(); !didFinalize) {
            std::string msg = "Failed to finalize mixed composition. ";
            msg += "This likely indicates an issue with the input compositions not summing to 1.\n";
            LOG_CRITICAL(m_logger, "{}", msg);
            throw exceptions::InvalidCompositionError(msg);
        }
        return mixedComposition;
    }

    double Composition::getMassFraction(const std::string& symbol) const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (!m_compositions.contains(symbol)) {
            LOG_ERROR(m_logger, "Symbol {} is not in the composition.", symbol);
            std::string currentSymbols;
            size_t count = 0;
            for (const auto& sym : m_compositions | std::views::keys) {
                currentSymbols += sym;
                if (count < m_compositions.size() - 2) {
                    currentSymbols += ", ";
                } else if (count == m_compositions.size() - 2) {
                    currentSymbols += ", and ";
                }
                count++;
            }
            throw exceptions::UnregisteredSymbolError("Symbol(" + symbol + ") is not in the current composition. Current composition has symbols: " + currentSymbols + ".");
        }
        if (m_massFracMode) {
            return m_compositions.at(symbol).mass_fraction();
        }

        return m_compositions.at(symbol).mass_fraction();
    }

    double Composition::getMassFraction(
        const atomic::Species &species
    ) const {
        return getMassFraction(std::string(species.name()));
    }

    std::unordered_map<std::string, double> Composition::getMassFraction() const {
        std::unordered_map<std::string, double> mass_fractions;
        for (const auto &symbol: m_compositions | std::views::keys) {
            mass_fractions[symbol] = getMassFraction(symbol);
        }
        return mass_fractions;
    }


    double Composition::getNumberFraction(
        const std::string& symbol
    ) const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (!m_compositions.contains(symbol)) {
            LOG_ERROR(m_logger, "Symbol {} is not in the composition.", symbol);
            throw exceptions::CompositionNotFinalizedError("Symbol " + symbol + " is not in the composition.");
        }
        if (!m_massFracMode) {
            return m_compositions.at(symbol).number_fraction();
        }
        return m_compositions.at(symbol).number_fraction(m_specificNumberDensity);
    }

    double Composition::getNumberFraction(
        const atomic::Species &species
    ) const {
        return getNumberFraction(std::string(species.name()));
    }

    std::unordered_map<std::string, double> Composition::getNumberFraction() const {
        std::unordered_map<std::string, double> number_fractions;
        for (const auto &symbol: m_compositions | std::views::keys) {
            number_fractions[symbol] = getNumberFraction(symbol);
        }
        return number_fractions;
    }

    double Composition::getMolarAbundance(
        const std::string &symbol
    ) const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (!m_compositions.contains(symbol)) {
            LOG_ERROR(m_logger, "Symbol {} is not in the composition.", symbol);
            throw exceptions::UnregisteredSymbolError("Symbol " + symbol + " is not in the composition.");
        }
        return getMassFraction(symbol) / m_compositions.at(symbol).isotope().mass();

    }

    double Composition::getMolarAbundance(
        const atomic::Species &species
    ) const {
        return getMolarAbundance(std::string(species.name()));
    }

    std::pair<CompositionEntry, GlobalComposition> Composition::getComposition(
        const std::string& symbol
    ) const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (!m_compositions.contains(symbol)) {
            LOG_ERROR(m_logger, "Symbol {} is not in the composition.", symbol);
            throw exceptions::UnregisteredSymbolError("Symbol " + symbol + " is not in the composition.");
        }
        return {m_compositions.at(symbol), {m_specificNumberDensity, m_meanParticleMass}};
    }

    std::pair<CompositionEntry, GlobalComposition> Composition::getComposition(
        const atomic::Species &species
    ) const {
        return getComposition(std::string(species.name()));
    }

    std::pair<std::unordered_map<std::string, CompositionEntry>, GlobalComposition> Composition::getComposition() const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        return {m_compositions, {m_specificNumberDensity, m_meanParticleMass}};
    }

    double Composition::getMeanParticleMass() const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        return m_meanParticleMass;
    }

    double Composition::getMeanAtomicNumber() const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition must be finalized before getting the mean atomic mass number. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition not finalized. Cannot retrieve mean atomic mass number. Hint: Consider running .finalize().");
        }

        double zSum = 0.0;

        for (const auto &val: m_compositions | std::views::values) {
            // Sum of (X_i * Z_i / A_i)
            zSum += (val.mass_fraction() * val.m_isotope.z())/val.m_isotope.a();
        }

        // <Z> = <A> * sum(X_i * Z_i / A_i)
        const double mean_A = m_meanParticleMass * zSum;
        return mean_A;
    }

    double Composition::getElectronAbundance() const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition must be finalized before getting the electron abundance. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition not finalized. Cannot retrieve electron abundance. Hint: Consider running .finalize().");
        }

        if (m_cache.Ye.has_value()) {
            return m_cache.Ye.value();
        }

        double Ye = 0.0;
        for (const auto &val: m_compositions | std::views::values) {
            Ye += (val.mass_fraction() * val.m_isotope.z())/val.m_isotope.a();
        }
        m_cache.Ye = Ye;
        return Ye;
    }

    Composition Composition::subset(
        const std::vector<std::string>& symbols,
        const std::string& method
    ) const {
        if (const std::array<std::string, 2> methods = {"norm", "none"}; std::ranges::find(methods, method) == methods.end()) {
            const std::string errorMessage = "Invalid method: " + method + ". Valid methods are 'norm' and 'none'.";
            LOG_ERROR(m_logger, "Invalid method: {}. Valid methods are norm and none.", method);
            throw exceptions::InvalidMixingMode(errorMessage);
        }

        Composition subsetComposition;
        for (const auto& symbol : symbols) {
            if (!m_compositions.contains(symbol)) {
                LOG_ERROR(m_logger, "Symbol {} is not in the composition.", symbol);
                throw exceptions::UnregisteredSymbolError("Symbol " + symbol + " is not in the composition.");
            }
            subsetComposition.registerSymbol(symbol);
            subsetComposition.setMassFraction(symbol, m_compositions.at(symbol).mass_fraction());
        }
        if (method == "norm") {
            if (const bool isNorm = subsetComposition.finalize(true); !isNorm) {
                LOG_ERROR(m_logger, "Subset composition is invalid. (Unable to finalize with normalization).");
                throw exceptions::FailedToFinalizeCompositionError("Subset composition is invalid. (Unable to finalize with normalization).");
            }
        }
        return subsetComposition;
    }

    void Composition::setCompositionMode(
        const bool massFracMode
    ) {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Mode cannot be set unless composition is finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Mode cannot be set unless composition is finalized. Hint: Consider running .finalize().");
        }

        bool okay;
        for (auto &entry: m_compositions | std::views::values) {
            if (massFracMode) {
                okay = entry.setMassFracMode(m_meanParticleMass);
            } else {
                okay = entry.setNumberFracMode(m_specificNumberDensity);
            }
            if (!okay) {
                LOG_ERROR(m_logger, "Composition mode could not be set due to some unknown error.");
                throw std::runtime_error("Composition mode could not be set due to an unknown error.");
            }
        }
        m_massFracMode = massFracMode;
    }

    CanonicalComposition Composition::getCanonicalComposition(
        const bool harsh
    ) const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (m_cache.canonicalComp.has_value()) {
            return m_cache.canonicalComp.value(); // Short circuit if we have cached the canonical composition
        }
        CanonicalComposition canonicalComposition;
        const std::array<std::string, 7> canonicalH = {
            "H-1", "H-2", "H-3", "H-4", "H-5", "H-6", "H-7"
        };
        const std::array<std::string, 8> canonicalHe = {
            "He-3", "He-4", "He-5", "He-6", "He-7", "He-8", "He-9", "He-10"
        };
        for (const auto& symbol : canonicalH) {
            if (hasSymbol(symbol)) {
                canonicalComposition.X += getMassFraction(symbol);
            }
        }
        for (const auto& symbol : canonicalHe) {
            if (hasSymbol(symbol)) {
                canonicalComposition.Y += getMassFraction(symbol);
            }
        }

        for (const auto& symbol : getRegisteredSymbols()) {
            const bool isHSymbol = std::ranges::find(canonicalH, symbol) != std::end(canonicalH);
            // ReSharper disable once CppTooWideScopeInitStatement
            const bool isHeSymbol = std::ranges::find(canonicalHe, symbol) != std::end(canonicalHe);

            if (isHSymbol || isHeSymbol) {
                continue; // Skip canonical H and He symbols
            }

            canonicalComposition.Z += getMassFraction(symbol);
        }

        // ReSharper disable once CppTooWideScopeInitStatement
        const double Z = 1.0 - (canonicalComposition.X + canonicalComposition.Y);
        if (std::abs(Z - canonicalComposition.Z) > 1e-6) {
            if (!harsh) {
                LOG_WARNING(m_logger, "Validation composition Z (X-Y = {}) is different than canonical composition Z ({}) (∑a_i where a_i != H/He).", Z, canonicalComposition.Z);
            }
            else {
                LOG_ERROR(m_logger, "Validation composition Z (X-Y = {}) is different than canonical composition Z ({}) (∑a_i where a_i != H/He).", Z, canonicalComposition.Z);
                throw std::runtime_error("Validation composition Z (X-Y = " + std::to_string(Z) + ") is different than canonical composition Z (" + std::to_string(canonicalComposition.Z) + ") (∑a_i where a_i != H/He).");
            }
        }
        m_cache.canonicalComp = canonicalComposition;
        return canonicalComposition;
    }

    std::vector<double> Composition::getMassFractionVector() const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (m_cache.massFractions.has_value()) {
            return m_cache.massFractions.value(); // Short circuit if we have cached the mass fractions
        }

        std::vector<double> massFractionVector;
        std::vector<double> speciesMass;

        massFractionVector.reserve(m_compositions.size());
        speciesMass.reserve(m_compositions.size());

        for (const auto &entry: m_compositions | std::views::values) {
            massFractionVector.push_back(entry.mass_fraction());
            speciesMass.push_back(entry.isotope().mass());
        }

        std::vector<double> massFractions = sortVectorBy(massFractionVector, speciesMass);
        m_cache.massFractions = massFractions; // Cache the result
        return massFractions;

    }

    std::vector<double> Composition::getNumberFractionVector() const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (m_cache.numberFractions.has_value()) {
            return m_cache.numberFractions.value(); // Short circuit if we have cached the number fractions
        }

        std::vector<double> numberFractionVector;
        std::vector<double> speciesMass;

        numberFractionVector.reserve(m_compositions.size());
        speciesMass.reserve(m_compositions.size());

        for (const auto &entry: m_compositions | std::views::values) {
            numberFractionVector.push_back(entry.number_fraction());
            speciesMass.push_back(entry.isotope().mass());
        }

        std::vector<double> numberFractions = sortVectorBy(numberFractionVector, speciesMass);
        m_cache.numberFractions = numberFractions; // Cache the result
        return numberFractions;
    }

    std::vector<double> Composition::getMolarAbundanceVector() const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (m_cache.molarAbundances.has_value()) {
            return m_cache.molarAbundances.value(); // Short circuit if we have cached the molar abundances
        }

        std::vector<double> molarAbundanceVector;
        std::vector<double> speciesMass;

        molarAbundanceVector.reserve(m_compositions.size());
        speciesMass.reserve(m_compositions.size());

        for (const auto &entry: m_compositions | std::views::values) {
            molarAbundanceVector.push_back(getMolarAbundance(entry.isotope()));
            speciesMass.push_back(entry.isotope().mass());
        }

        std::vector<double> molarAbundances = sortVectorBy(molarAbundanceVector, speciesMass);
        m_cache.molarAbundances = molarAbundances; // Cache the result
        return molarAbundances;

    }

    size_t Composition::getSpeciesIndex(
        const std::string &symbol
    ) const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (!m_compositions.contains(symbol)) {
            LOG_ERROR(m_logger, "Symbol {} is not in the composition.", symbol);
            throw exceptions::UnregisteredSymbolError("Symbol " + symbol + " is not in the composition.");
        }
        if (m_cache.sortedSymbols.has_value()) {
            return std::distance(
                m_cache.sortedSymbols->begin(),
                std::ranges::find(
                    m_cache.sortedSymbols.value().begin(),
                    m_cache.sortedSymbols.value().end(),
                    symbol
                )
            );
        }

        std::vector<std::string> symbols;
        std::vector<double> speciesMass;

        symbols.reserve(m_compositions.size());
        speciesMass.reserve(m_compositions.size());

        for (const auto &entry: m_compositions | std::views::values) {
            symbols.emplace_back(entry.isotope().name());
            speciesMass.push_back(entry.isotope().mass());
        }

        std::vector<std::string> sortedSymbols = sortVectorBy(symbols, speciesMass);
        m_cache.sortedSymbols = sortedSymbols;
        return std::distance(sortedSymbols.begin(), std::ranges::find(sortedSymbols, symbol));
    }

    size_t Composition::getSpeciesIndex(
        const atomic::Species &species
    ) const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (!m_compositions.contains(static_cast<std::string>(species.name()))) {
            LOG_ERROR(m_logger, "Species {} is not in the composition.", species.name());
            throw exceptions::UnregisteredSymbolError("Species " + std::string(species.name()) + " is not in the composition.");
        }
        if (m_cache.sortedSpecies.has_value()) {
            return std::distance(
                m_cache.sortedSpecies->begin(),
                std::ranges::find(
                    m_cache.sortedSpecies.value().begin(),
                    m_cache.sortedSpecies.value().end(),
                    species
                )
            );
        }

        std::vector<atomic::Species> speciesVector;
        std::vector<double> speciesMass;

        speciesVector.reserve(m_compositions.size());
        speciesMass.reserve(m_compositions.size());

        for (const auto &entry: m_compositions | std::views::values) {
            speciesVector.emplace_back(entry.isotope());
            speciesMass.push_back(entry.isotope().mass());
        }

        std::vector<atomic::Species> sortedSpecies = sortVectorBy(speciesVector, speciesMass);
        m_cache.sortedSpecies = sortedSpecies;
        return std::distance(sortedSpecies.begin(), std::ranges::find(sortedSpecies, species));
    }

    atomic::Species Composition::getSpeciesAtIndex(
        size_t index
    ) const {
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (index >= m_compositions.size()) {
            LOG_ERROR(m_logger, "Index {} is out of bounds for composition of size {}.", index, m_compositions.size());
            throw std::out_of_range("Index " + std::to_string(index) + " is out of bounds for composition of size " + std::to_string(m_compositions.size()) + ".");
        }
        if (m_cache.sortedSpecies.has_value()) {
            return m_cache.sortedSpecies.value().at(index);
        }

        std::vector<atomic::Species> speciesVector;
        std::vector<double> speciesMass;

        speciesVector.reserve(m_compositions.size());
        speciesMass.reserve(m_compositions.size());

        for (const auto &entry: m_compositions | std::views::values) {
            speciesVector.emplace_back(entry.isotope());
            speciesMass.push_back(entry.isotope().mass());
        }

        std::vector<atomic::Species> sortedSymbols = sortVectorBy(speciesVector, speciesMass);
        return sortedSymbols.at(index);
    }

    bool Composition::hasSymbol(
        const std::string& symbol
    ) const {
        return m_compositions.contains(symbol);
    }

    bool Composition::hasSpecies(const fourdst::atomic::Species &species) const {
        for (const auto &entry: m_compositions | std::views::values) {
            if (entry.isotope() == species) {
                return true;
            }
        }
        return false;
    }

    bool Composition::contains(
        const atomic::Species &isotope
    ) const {
        // Check if the isotope's symbol is in the composition
        if (!m_finalized) {
            LOG_ERROR(m_logger, "Composition has not been finalized. Hint: Consider running .finalize().");
            throw exceptions::CompositionNotFinalizedError("Composition has not been finalized. Hint: Consider running .finalize().");
        }
        if (const auto symbol = static_cast<std::string>(isotope.name()); m_compositions.contains(symbol)) {
            return true;
        }
        return false;
    }

    /// OVERLOADS

    Composition Composition::operator+(
        const Composition& other
    ) const {
        return mix(other, 0.5);
    }

    std::ostream& operator<<(
        std::ostream& os,
        const GlobalComposition& comp
    ) {
        os << "Global Composition: \n";
        os << "\tSpecific Number Density: " << comp.specificNumberDensity << "\n";
        os << "\tMean Particle Mass: " << comp.meanParticleMass << "\n";
        return os;
    }

    std::ostream& operator<<(
        std::ostream& os,
        const CompositionEntry& entry
    ) {
        os << "<" << entry.m_symbol << " : m_frac = " << entry.mass_fraction() << ">";
        return os;
    }

    std::ostream& operator<<(
        std::ostream& os,
        const Composition& composition
    ) {
        os << "Composition(finalized: " << (composition.m_finalized ? "true" : "false") << ", " ;
        size_t count = 0;
        for (const auto &entry: composition.m_compositions | std::views::values) {
            os << entry;
            if (count < composition.m_compositions.size() - 1) {
                os << ", ";
            }
            count++;
        }
        os << ")";
        return os;
    }

} // namespace fourdst::composition
