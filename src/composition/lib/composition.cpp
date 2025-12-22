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
#include <ranges>
#include <algorithm>
#include <set>
#include <string>


#include <utility>

#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/atomic/species.h"
#include "fourdst/composition/composition.h"

#include <numeric>

#include "fourdst/composition/utils/composition_hash.h"
#include "fourdst/composition/utils.h"

#include "fourdst/composition/exceptions/exceptions_composition.h"

namespace {
    void throw_unknown_symbol(quill::Logger* logger, const std::string& symbol) {
        LOG_ERROR(logger, "Symbol {} is not a valid species symbol (not in the species database)", symbol);
        throw fourdst::composition::exceptions::UnknownSymbolError("Symbol " + symbol + " is not a valid species symbol (not in the species database)");
    }

    void throw_unregistered_symbol(quill::Logger* logger, const std::string& symbol) {
        LOG_ERROR(logger, "Symbol {} is not registered in the composition.", symbol);
        throw fourdst::composition::exceptions::UnregisteredSymbolError("Symbol " + symbol + " is not registered in the composition.");
    }
}

namespace fourdst::composition {

    /////////////////////////////////////////////
    /// Constructors without molar abundances ///
    /// These all delegate to the ctor        ///
    /// vector<Species>                       ///
    /////////////////////////////////////////////

    Composition::Composition(
        const std::set<std::string>& symbols
    ) : Composition(symbols | std::ranges::to<std::vector>()) {}

    Composition::Composition(
        const std::set<atomic::Species> &species
    ) : Composition(species | std::ranges::to<std::vector>()) {}

    Composition::Composition(
        const std::unordered_set<std::string> &symbols
    ) : Composition(symbols | std::ranges::to<std::vector>()) {}

    Composition::Composition(
        const std::unordered_set<atomic::Species> &species
    ) : Composition(species | std::ranges::to<std::vector>()) {}

    Composition::Composition(
        const std::vector<std::string>& symbols
    ) : Composition(symbolVectorToSpeciesVector(symbols)) {}

    Composition::Composition(
        const std::vector<atomic::Species> &species
    ) {
        m_species = species;
        std::ranges::sort(m_species, [&](const atomic::Species& a, const atomic::Species& b) {
            return a < b;
        });

        const auto last = std::ranges::unique(m_species).begin();
        m_species.erase(last, m_species.end());

        m_molarAbundances.resize(m_species.size(), 0.0);
    }

    //////////////////////////////////////////
    /// Constructors with molar abundances ///
    /// These all delegate to the ctor     ///
    /// vector<Species, vector<double>>    ///
    //////////////////////////////////////////

    Composition::Composition(
        const std::vector<std::string>& symbols,
        const std::vector<double>& molarAbundances
    ) : Composition(symbolVectorToSpeciesVector(symbols), molarAbundances) {}

    Composition::Composition(
        const std::set<std::string> &symbols,
        const std::vector<double> &molarAbundances
    ) : Composition(symbolVectorToSpeciesVector(symbols | std::ranges::to<std::vector>()), molarAbundances) {}


    Composition::Composition(
        const std::unordered_map<std::string, double> &symbolMolarAbundances
    ) : Composition(
        symbolMolarAbundances | std::views::keys | std::ranges::to<std::vector>(),
        symbolMolarAbundances | std::views::values | std::ranges::to<std::vector>()
        ) {}

    Composition::Composition(
        const std::map<std::string, double> &symbolMolarAbundances
    ) : Composition(
        symbolMolarAbundances | std::views::keys | std::ranges::to<std::vector>(),
        symbolMolarAbundances | std::views::values | std::ranges::to<std::vector>()
        ) {}

    Composition::Composition(
        const std::unordered_map<atomic::Species, double> &speciesMolarAbundances
    ) : Composition(
        speciesMolarAbundances | std::views::keys | std::ranges::to<std::vector>(),
        speciesMolarAbundances | std::views::values | std::ranges::to<std::vector>()
        ) {}

    Composition::Composition(
        const std::map<atomic::Species, double> &speciesMolarAbundances
    ) : Composition(
        speciesMolarAbundances | std::views::keys | std::ranges::to<std::vector>(),
        speciesMolarAbundances | std::views::values | std::ranges::to<std::vector>()
        ) {}

    Composition::Composition(
        const std::vector<atomic::Species> &species,
        const std::vector<double> &molarAbundances
    ) {
        if (__builtin_expect(species.size() != molarAbundances.size(), 0)) {
            LOG_CRITICAL(getLogger(), "The number of species and molarAbundances must be equal (got {} species and {} molarAbundances).", species.size(), molarAbundances.size());
            throw exceptions::InvalidCompositionError("The number of species and fractions must be equal. Got " + std::to_string(species.size()) + " species and " + std::to_string(molarAbundances.size()) + " fractions.");
        }

        const size_t numSpecies = species.size();
        m_species.reserve(numSpecies);
        m_molarAbundances.reserve(numSpecies);

        for (size_t i = 0; i < numSpecies; ++i) {
            m_species.push_back(species[i]);
            if (__builtin_expect(molarAbundances[i] < 0.0, 0)) {
                LOG_CRITICAL(getLogger(), "Molar abundance for species {} is negative (y = {}). Molar abundances must be non-negative.", species[i].name(), molarAbundances[i]);
                throw exceptions::InvalidCompositionError("Molar abundance for species " + std::string(species[i].name()) + " is negative (y = " + std::to_string(molarAbundances[i]) + "). Molar abundances must be non-negative.");
            }
            m_molarAbundances.push_back(molarAbundances[i]);
        }

        auto combined = std::views::zip(m_species, m_molarAbundances);

        std::ranges::sort(combined, [](const auto& a, const auto& b) -> bool {
            const auto& spA = std::get<0>(a);
            const auto& spB = std::get<0>(b);

            if (spA != spB) {
                return spA < spB;
            }

            return std::get<1>(a) > std::get<1>(b);
        });

        auto [first, last] = std::ranges::unique(combined, [](const auto& a, const auto& b) {
            return std::get<0>(a) == std::get<0>(b);
        });

        const auto newEndIndex = std::distance(combined.begin(), first);
        m_species.erase(m_species.begin() + newEndIndex, m_species.end());
        m_molarAbundances.erase(m_molarAbundances.begin() + newEndIndex, m_molarAbundances.end());
    }

    ////////////////////////////////////////////
    /// Copy and conversion constructors     ///
    ////////////////////////////////////////////

    Composition::Composition(const Composition &composition) {
        m_species = composition.m_species;
        m_molarAbundances = composition.m_molarAbundances;
    }

    Composition::Composition(const CompositionAbstract &composition) {
        for (const auto& species : composition.getRegisteredSpecies()) {
            registerSpecies(species);
            setMolarAbundance(species, composition.getMolarAbundance(species));
        }
    }

    Composition& Composition::operator=(
        const Composition &other
    ) {
        if (this != &other) {
            m_species = other.m_species;
            m_molarAbundances   = other.m_molarAbundances;
        }
        m_cache.clear();
        return *this;
    }

    Composition & Composition::operator=(const CompositionAbstract &other) {
        m_species.clear();
        m_molarAbundances.clear();
        m_cache.clear();
        for (const auto& species : other.getRegisteredSpecies()) {
            registerSpecies(species);
            setMolarAbundance(species, other.getMolarAbundance(species));
        }
        return *this;
    }

    std::unique_ptr<CompositionAbstract> Composition::clone() const {
        return std::make_unique<Composition>(*this);
    }

    //------------------------------------------
    // Registration methods
    //------------------------------------------

    void Composition::registerSymbol(
        const std::string& symbol
    ) {
        const auto result = getSpecies(symbol);
        if (!result) {
            throw_unknown_symbol(getLogger(), symbol);
        }

        registerSpecies(result.value());
    }

    void Composition::registerSymbol(
        const std::vector<std::string>& symbols
    ) {
        registerSpecies(symbolVectorToSpeciesVector(symbols));
    }

    void Composition::registerSpecies(
        const atomic::Species &species
    ) noexcept {
        if (const auto it = std::ranges::lower_bound(m_species, species); it == m_species.end() || *it != species) {
            const auto index = std::distance(m_species.begin(), it);
            m_species.insert(it, species);
            m_molarAbundances.insert(m_molarAbundances.begin() + index, 0.0);
            m_cache.clear();
        }
    }

    void Composition::registerSpecies(
        const std::vector<atomic::Species> &species
    ) noexcept {
        // We do not simply call registerSpecies(species) here as that would have a complexity of O(n^2) due to constantly
        // reinserting into the vector. Rather we build the vector once and then sort it

        if (species.empty()) return;

        const size_t total_size = m_species.size() + species.size();
        m_species.reserve(total_size);
        m_molarAbundances.reserve(total_size);

        for (const auto& sp : species) {
            m_species.push_back(sp);
            m_molarAbundances.push_back(0.0);
        }

        auto combined = std::views::zip(m_species, m_molarAbundances);

        std::ranges::sort(combined, [](const auto& a, const auto& b) {
            const auto& speciesA = std::get<0>(a);
            const auto& speciesB = std::get<0>(b);

            if (speciesA != speciesB) {
                return speciesA < speciesB;
            }

            return std::get<1>(a) > std::get<1>(b);
        });

        auto [first, last] = std::ranges::unique(combined, [](const auto& a, const auto& b) {
            return std::get<0>(a) == std::get<0>(b);
        });

        const auto newEndIndex = std::distance(combined.begin(), first);

        m_species.erase(m_species.begin() + newEndIndex, m_species.end());
        m_molarAbundances.erase(m_molarAbundances.begin() + newEndIndex, m_molarAbundances.end());

        m_cache.clear();
    }

    std::set<std::string> Composition::getRegisteredSymbols() const noexcept {
        std::set<std::string> symbols;
        for (const auto& species : m_species) {
            symbols.insert(std::string(species.name()));
        }
        return symbols;
    }

    const std::vector<atomic::Species> &Composition::getRegisteredSpecies() const noexcept {
        return m_species;
    }


    //------------------------------------------
    // Molar abundance setters
    //------------------------------------------

    void Composition::setMolarAbundance(
        const std::string &symbol,
        const double &molar_abundance
    ) {
        const auto species = getSpecies(symbol);
        if (__builtin_expect(!species, 0)) {
            throw_unknown_symbol(getLogger(), symbol);
        }

        setMolarAbundance(species.value(), molar_abundance);
    }

    void Composition::setMolarAbundance(
        const atomic::Species &species,
        const double &molar_abundance
    ) {
        if (__builtin_expect(molar_abundance < 0.0, 0)) {
            LOG_ERROR(getLogger(), "Molar abundance must be non-negative for symbol {}. Currently it is {}.", species.name(), molar_abundance);
            throw exceptions::InvalidCompositionError("Molar abundance must be non-negative, got " + std::to_string(molar_abundance) + " for symbol " + std::string(species.name()) + ".");
        }

        const std::expected<std::ptrdiff_t, SpeciesIndexLookupError> speciesIndexResult = findSpeciesIndex(species);
        if (__builtin_expect(!speciesIndexResult, 0)) {
            throw_unregistered_symbol(getLogger(), std::string(species.name()));
        }

        assert(static_cast<size_t>(speciesIndexResult.value()) < m_molarAbundances.size());

        m_molarAbundances[speciesIndexResult.value()] = molar_abundance;
        m_cache.clear();
    }


    ////----------------------------------------------
    ///   Methods which set multiple molar abundances
    ///   delegate to vector<Species>, vector<double>
    ///-----------------------------------------------

    void Composition::setMolarAbundance(
        const std::vector<std::string> &symbols,
        const std::vector<double> &molar_abundances
    ) {
        setMolarAbundance(symbolVectorToSpeciesVector(symbols), molar_abundances);
    }

    void Composition::setMolarAbundance(
        const std::set<std::string> &symbols,
        const std::vector<double> &molar_abundances
    ) {
        setMolarAbundance(symbolVectorToSpeciesVector(symbols | std::ranges::to<std::vector>()), molar_abundances);
    }

    void Composition::setMolarAbundance(
        const std::set<atomic::Species> &species,
        const std::vector<double> &molar_abundances
    ) {
        setMolarAbundance(species | std::ranges::to<std::vector>(), molar_abundances);
    }

    void Composition::setMolarAbundance(
        const std::vector<atomic::Species> &species,
        const std::vector<double> &molar_abundances
    ) {
        if (__builtin_expect(species.size() != molar_abundances.size(), 0)) {
            LOG_CRITICAL(getLogger(), "The number of species and molar_abundances must be equal (got {} species and {} molar_abundances).", species.size(), molar_abundances.size());
            throw exceptions::InvalidCompositionError("The number of species and fractions must be equal. Got " + std::to_string(species.size()) + " species and " + std::to_string(molar_abundances.size()) + " fractions.");
        }

        if (species.empty()) return;

        if (species.size() == m_species.size()) {
            if (species == m_species) {
                for (const auto& [sp, y] : std::views::zip(species, molar_abundances)) {
                    if (__builtin_expect(y < 0.0, 0)) {
                        LOG_ERROR(getLogger(), "Molar abundance must be non-negative. Instead got {} for species {}.", y, sp.name());
                        throw exceptions::InvalidCompositionError("Molar abundance must be non-negative. Instead got " + std::to_string(y) + " for species " + std::string(sp.name()) + ".");
                    }
                }

                m_molarAbundances = molar_abundances;
                m_cache.clear();
                return;
            }
        }

        for (size_t i  = 0; i < species.size(); ++i) {
            const double y = molar_abundances[i];
            const auto& sp = species[i];
            if (__builtin_expect(y < 0.0, 0)) {
                LOG_CRITICAL(getLogger(), "Molar abundance must be non-negative. Instead got {} for species {}.", y, sp.name());
                throw exceptions::InvalidCompositionError("Molar abundance must be non-negative. Instead got " + std::to_string(y) + " for species " + std::string(sp.name()) + ".");
            }

            const std::expected<std::ptrdiff_t, SpeciesIndexLookupError> speciesIndexResult = findSpeciesIndex(sp);
            if (__builtin_expect(!speciesIndexResult, 0)) {
                throw_unregistered_symbol(getLogger(), std::string(sp.name()));
            }

            const std::ptrdiff_t speciesIndex = speciesIndexResult.value();

            m_molarAbundances[speciesIndex] = y;
        }

        m_cache.clear();
    }


    //------------------------------------------
    // Fraction and abundance getters
    //------------------------------------------

    double Composition::getMassFraction(const std::string& symbol) const {
        const auto species = getSpecies(symbol);
        if (!species) {
            throw_unknown_symbol(getLogger(), symbol);
        }
        return getMassFraction(species.value());
    }

    double Composition::getMassFraction(
        const atomic::Species &species
    ) const {
        const std::expected<std::ptrdiff_t, SpeciesIndexLookupError> speciesIndexResult = findSpeciesIndex(species);
        if (!speciesIndexResult) {
            throw_unregistered_symbol(getLogger(), std::string(species.name()));
        }

        double totalMass = 0;
        double speciesMass = 0;
        for (const auto& [sp, y] : *this) {
            const double contrib = y * sp.mass();
            totalMass += contrib;
            if (sp == species) {
                speciesMass = contrib;
            }
        }
        return speciesMass / totalMass;
    }

    std::unordered_map<atomic::Species, double> Composition::getMassFraction() const noexcept {
        std::unordered_map<atomic::Species, double> mass_fractions;
        for (const auto &species: *this | std::views::keys) {
            mass_fractions.emplace(species, getMassFraction(species));
        }
        return mass_fractions;
    }


    double Composition::getNumberFraction(
        const std::string& symbol
    ) const {
        const auto species = getSpecies(symbol);
        if (!species) {
            throw_unknown_symbol(getLogger(), symbol);
        }
        return getNumberFraction(species.value());
    }

    double Composition::getNumberFraction(
        const atomic::Species &species
    ) const {
        const std::expected<std::ptrdiff_t, SpeciesIndexLookupError> speciesIndexResult = findSpeciesIndex(species);
        if (!speciesIndexResult) {
            throw_unregistered_symbol(getLogger(), std::string(species.name()));
        }
        const std::ptrdiff_t speciesIndex = speciesIndexResult.value();

        const double total_moles_per_gram = std::accumulate(
            m_molarAbundances.begin(),
            m_molarAbundances.end(),
            0.0
        );
        return m_molarAbundances[speciesIndex] / total_moles_per_gram;
    }

    std::unordered_map<atomic::Species, double> Composition::getNumberFraction() const noexcept {
        std::unordered_map<atomic::Species, double> number_fractions;
        for (const auto &species: m_species) {
            number_fractions.emplace(species, getNumberFraction(species));
        }
        return number_fractions;
    }

    double Composition::getMolarAbundance(
        const std::string &symbol
    ) const {
        const auto species = getSpecies(symbol);
        if (!species) {
            throw_unknown_symbol(getLogger(), symbol);
        }
        return getMolarAbundance(species.value());

    }

    double Composition::getMolarAbundance(
        const atomic::Species &species
    ) const {
        const std::expected<std::ptrdiff_t, SpeciesIndexLookupError> speciesIndexResult = findSpeciesIndex(species);
        if (!speciesIndexResult) {
            throw_unregistered_symbol(getLogger(), std::string(species.name()));
        }
        const std::ptrdiff_t speciesIndex = speciesIndexResult.value();
        return m_molarAbundances[speciesIndex];
    }

    //------------------------------------------
    // Derived property getters
    //------------------------------------------

    double Composition::getMeanParticleMass() const noexcept {
        double totalMass = 0.0;
        double totalMoles = 0.0;

        for (size_t i = 0; i < m_species.size(); ++i) {
            totalMoles += m_molarAbundances[i];
            totalMass  += m_molarAbundances[i] * m_species[i].mass();
        }

        return totalMass / totalMoles;
    }


    double Composition::getElectronAbundance() const noexcept {
        double Ye = 0.0;
        for (const auto& [species, y] : *this) {
            Ye += species.z() * y;
        }
        return Ye;
    }


    CanonicalComposition Composition::getCanonicalComposition(
    ) const {
        using namespace fourdst::atomic;

        if (m_cache.canonicalComp.has_value()) {
            return m_cache.canonicalComp.value(); // Short circuit if we have cached the canonical composition
        }
        CanonicalComposition canonicalComposition;
        static const std::unordered_set<Species> canonicalH = {H_1, H_2, H_3, H_4, H_5, H_6, H_7};
        static const std::unordered_set<Species> canonicalHe = {He_3, He_4, He_5, He_6, He_7, He_8, He_9, He_10};

        for (const auto& symbol : canonicalH) {
            if (contains(symbol)) {
                canonicalComposition.X += getMassFraction(symbol);
            }
        }
        for (const auto& symbol : canonicalHe) {
            if (contains(symbol)) {
                canonicalComposition.Y += getMassFraction(symbol);
            }
        }

        for (const auto& species : m_species) {
            if (canonicalH.contains(species) || canonicalHe.contains(species)) {
                continue; // Skip canonical H and He symbols
            }

            canonicalComposition.Z += getMassFraction(species);
        }

        // ReSharper disable once CppTooWideScopeInitStatement
        const double Z = 1.0 - (canonicalComposition.X + canonicalComposition.Y);
        if (std::abs(Z - canonicalComposition.Z) > 1e-16) {
            LOG_ERROR(getLogger(), "Validation composition Z (X-Y = {}) is different than canonical composition Z ({}) (∑a_i where a_i != H/He).", Z, canonicalComposition.Z);
            throw exceptions::InvalidCompositionError("Validation composition Z (X-Y = " + std::to_string(Z) + ") is different than canonical composition Z (" + std::to_string(canonicalComposition.Z) + ") (∑a_i where a_i != H/He).");
        }
        m_cache.canonicalComp = canonicalComposition;
        return canonicalComposition;
    }

    //------------------------------------------
    // Vector getters
    //------------------------------------------

    std::vector<double> Composition::getMassFractionVector() const noexcept {
        if (m_cache.massFractions.has_value()) {
            return m_cache.massFractions.value(); // Short circuit if we have cached the mass fractions
        }

        std::vector<double> massFractionVector;

        massFractionVector.reserve(m_molarAbundances.size());

        for (const auto &species: m_species) {
            massFractionVector.push_back(getMassFraction(species));
        }

        m_cache.massFractions = massFractionVector; // Cache the result
        return massFractionVector;

    }

    std::vector<double> Composition::getNumberFractionVector() const noexcept {
        if (m_cache.numberFractions.has_value()) {
            return m_cache.numberFractions.value(); // Short circuit if we have cached the number fractions
        }

        std::vector<double> numberFractionVector;

        numberFractionVector.reserve(m_molarAbundances.size());

        for (const auto &species: m_species) {
            numberFractionVector.push_back(getNumberFraction(species));
        }

        m_cache.numberFractions = numberFractionVector; // Cache the result
        return numberFractionVector;
    }

    std::vector<double> Composition::getMolarAbundanceVector() const noexcept {
        return m_molarAbundances;
    }

    //------------------------------------------
    // Species index getters and lookups
    //------------------------------------------

    size_t Composition::getSpeciesIndex(
        const std::string &symbol
    ) const {
        const auto species = getSpecies(symbol);
        if (!species) {
            throw_unknown_symbol(getLogger(), symbol);
        }

        return getSpeciesIndex(species.value());
    }

    size_t Composition::getSpeciesIndex(
        const atomic::Species &species
    ) const {
        std::expected<std::ptrdiff_t, SpeciesIndexLookupError> speciesIndexResult = findSpeciesIndex(species);
        if (!speciesIndexResult) {
            switch (speciesIndexResult.error()) {
                case SpeciesIndexLookupError::NO_REGISTERED_SPECIES:
                    [[fallthrough]];
                case SpeciesIndexLookupError::SPECIES_NOT_FOUND:
                    throw_unregistered_symbol(getLogger(), std::string(species.name()));
                default:
                    throw std::logic_error("Unhandled SpeciesIndexLookupError in Composition::getSpeciesIndex");
            }
        }

        return static_cast<size_t>(speciesIndexResult.value());
    }

    atomic::Species Composition::getSpeciesAtIndex(
        const size_t index
    ) const {
        if (index >= m_species.size()) {
            LOG_ERROR(getLogger(), "Index {} is out of bounds for registered species (size {}).", index, m_species.size());
            throw std::out_of_range("Index " + std::to_string(index) + " is out of bounds for registered species (size " + std::to_string(m_species.size()) + ").");
        }

        return m_species[index];
    }


    //------------------------------------------
    // Utility methods
    //------------------------------------------

    std::size_t Composition::hash() const {
        if (m_cache.hash.has_value()) {
            return m_cache.hash.value();
        }
        std::size_t hash = utils::CompositionHash::hash_exact(*this);
        m_cache.hash = hash;
        return hash;
    }

    bool Composition::contains(
        const atomic::Species &species
    ) const noexcept {
        return std::ranges::binary_search(m_species, species);
    }

    bool Composition::contains(
        const std::string &symbol
    ) const {
        const auto species = getSpecies(symbol);
        if (!species) {
            throw_unknown_symbol(getLogger(), symbol);
        }
        return contains(species.value());
    }

    size_t Composition::size() const noexcept {
        return m_species.size();
    }

    std::expected<std::ptrdiff_t, Composition::SpeciesIndexLookupError> Composition::findSpeciesIndex(const atomic::Species &species) const noexcept {
        if (m_species.empty()) return std::unexpected(SpeciesIndexLookupError::NO_REGISTERED_SPECIES);

        const auto it = std::ranges::lower_bound(m_species, species);

        if (it == m_species.end() || *it != species) {
            return std::unexpected(SpeciesIndexLookupError::SPECIES_NOT_FOUND);
        }

        return std::distance(m_species.begin(), it);
    }

    std::vector<atomic::Species> Composition::symbolVectorToSpeciesVector(const std::vector<std::string> &symbols) {
        std::vector<atomic::Species> species;
        species.reserve(symbols.size());


        for (const auto& symbol : symbols) {
            const auto speciesResult = getSpecies(symbol);
            if (!speciesResult) {
                throw_unknown_symbol(getLogger(), symbol);
            }
            species.push_back(speciesResult.value());
        }

        return species;
    }


    //------------------------------------------
    // Stream operator
    //------------------------------------------

    std::ostream& operator<<(
        std::ostream& os,
        const Composition& composition
    ) {
        os << "Composition(Mass Fractions => [";
        size_t count = 0;
        for (const auto &species : composition.m_species) {
            os << species << ": " << composition.getMassFraction(species);
            if (count < composition.size() - 1) {
                os << ", ";
            }
            count++;
        }
        os << "])";
        return os;
    }

} // namespace fourdst::composition
