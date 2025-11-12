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

    std::optional<fourdst::atomic::Species> getSpecies(const std::string& symbol) {
        if (!fourdst::atomic::species.contains(symbol)) {
            return std::nullopt;
        }
        return fourdst::atomic::species.at(symbol);
    }

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
        const std::vector<atomic::Species> &species
    ) {
        for (const auto& s : species) {
            registerSpecies(s);
        }
    }

    Composition::Composition(
        const std::set<atomic::Species> &species
    ) {
        for (const auto& s : species) {
            registerSpecies(s);
        }
    }

    Composition::Composition(
        const std::vector<std::string>& symbols,
        const std::vector<double>& molarAbundances
    ) {
        if (symbols.size() != molarAbundances.size()) {
            LOG_CRITICAL(getLogger(), "The number of symbols and molarAbundances must be equal (got {} symbols and {} molarAbundances).", symbols.size(), molarAbundances.size());
            throw exceptions::InvalidCompositionError("The number of symbols and fractions must be equal. Got " + std::to_string(symbols.size()) + " symbols and " + std::to_string(molarAbundances.size()) + " fractions.");
        }

        for (const auto &[symbol, y] : std::views::zip(symbols, molarAbundances)) {
            registerSymbol(symbol);
            setMolarAbundance(symbol, y);
        }
    }

    Composition::Composition(
        const std::vector<atomic::Species> &species,
        const std::vector<double> &molarAbundances
    ) {
        if (species.size() != molarAbundances.size()) {
            LOG_CRITICAL(getLogger(), "The number of species and molarAbundances must be equal (got {} species and {} molarAbundances).", species.size(), molarAbundances.size());
            throw exceptions::InvalidCompositionError("The number of species and fractions must be equal. Got " + std::to_string(species.size()) + " species and " + std::to_string(molarAbundances.size()) + " fractions.");
        }

        for (const auto& [s, y] : std::views::zip(species, molarAbundances)) {
            registerSpecies(s);
            setMolarAbundance(s, y);
        }
    }

    Composition::Composition(
        const std::set<std::string> &symbols,
        const std::vector<double> &molarAbundances
    ) {
        if (symbols.size() != molarAbundances.size()) {
            LOG_CRITICAL(getLogger(), "The number of symbols and molarAbundances must be equal (got {} symbols and {} molarAbundances).", symbols.size(), molarAbundances.size());
            throw exceptions::InvalidCompositionError("The number of symbols and fractions must be equal. Got " + std::to_string(symbols.size()) + " symbols and " + std::to_string(molarAbundances.size()) + " fractions.");
        }

        for (const auto& [symbol, y] : std::views::zip(sortVectorBy<std::string>(std::vector<std::string>(symbols.begin(), symbols.end()), molarAbundances), molarAbundances)) {
            registerSymbol(symbol);
            setMolarAbundance(symbol, y);
        }
    }

    Composition::Composition(
        const Composition &composition
    ) {
        m_registeredSpecies = composition.m_registeredSpecies;
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
            m_registeredSpecies = other.m_registeredSpecies;
            m_molarAbundances   = other.m_molarAbundances;
        }
        return *this;
    }

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
        for (const auto& symbol : symbols) {
            registerSymbol(symbol);
        }
    }

    void Composition::registerSpecies(
        const atomic::Species &species
    ) noexcept {
        m_registeredSpecies.insert(species);
        if (!m_molarAbundances.contains(species)) {
            m_molarAbundances.emplace(species, 0.0);
        }
    }

    void Composition::registerSpecies(
        const std::vector<atomic::Species> &species
    ) noexcept {
        for (const auto& s : species) {
            registerSpecies(s);
        }
    }

    std::set<std::string> Composition::getRegisteredSymbols() const noexcept {
        std::set<std::string> symbols;
        for (const auto& species : m_registeredSpecies) {
            symbols.insert(std::string(species.name()));
        }
        return symbols;
    }

    const std::set<atomic::Species> &Composition::getRegisteredSpecies() const noexcept {
        return m_registeredSpecies;
    }


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
        if (!m_molarAbundances.contains(species)) {
            throw_unregistered_symbol(getLogger(), std::string(species.name()));
        }
        std::map<atomic::Species, double> raw_mass;
        double totalMass = 0;
        for (const auto& [sp, y] : m_molarAbundances) {
            const double contrib = y * sp.mass();
            totalMass += contrib;
            raw_mass.emplace(sp, contrib);
        }
        return raw_mass.at(species) / totalMass;
    }

    std::unordered_map<atomic::Species, double> Composition::getMassFraction() const noexcept {
        std::unordered_map<atomic::Species, double> mass_fractions;
        for (const auto &species: m_molarAbundances | std::views::keys) {
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
        if (!m_molarAbundances.contains(species)) {
            throw_unregistered_symbol(getLogger(), std::string(species.name()));
        }
        double total_moles_per_gram = 0.0;
        for (const auto &y: m_molarAbundances | std::views::values) {
            total_moles_per_gram += y;
        }
        return m_molarAbundances.at(species) / total_moles_per_gram;
    }

    std::unordered_map<atomic::Species, double> Composition::getNumberFraction() const noexcept {
        std::unordered_map<atomic::Species, double> number_fractions;
        for (const auto &species: m_molarAbundances | std::views::keys) {
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
        if (!m_molarAbundances.contains(species)) {
            throw_unregistered_symbol(getLogger(), std::string(species.name()));
        }
        return m_molarAbundances.at(species);
    }

    double Composition::getMeanParticleMass() const noexcept {
        std::vector<double> X = getMassFractionVector();
        double sum = 0.0;
        for (const auto& [species, x] : std::views::zip(m_registeredSpecies, X)) {
            sum += x/species.mass();
        }

        return 1.0 / sum;
    }

    double Composition::getElectronAbundance() const noexcept {
        double Ye = 0.0;
        for (const auto& [species, y] : m_molarAbundances) {
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
        const std::set<Species> canonicalH = {H_1, H_2, H_3, H_4, H_5, H_6, H_7};
        const std::set<Species> canonicalHe = {He_3, He_4, He_5, He_6, He_7, He_8, He_9, He_10};

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

        for (const auto& species : m_molarAbundances | std::views::keys) {
            const bool isHIsotope = canonicalH.contains(species);
            const bool isHeIsotope = canonicalHe.contains(species);

            if (isHIsotope || isHeIsotope) {
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

    std::vector<double> Composition::getMassFractionVector() const noexcept {
        if (m_cache.massFractions.has_value()) {
            return m_cache.massFractions.value(); // Short circuit if we have cached the mass fractions
        }

        std::vector<double> massFractionVector;
        std::vector<double> speciesMass;

        massFractionVector.reserve(m_molarAbundances.size());
        speciesMass.reserve(m_molarAbundances.size());

        for (const auto &species: m_molarAbundances | std::views::keys) {
            massFractionVector.push_back(getMassFraction(species));
            speciesMass.push_back(species.mass());
        }

        std::vector<double> massFractions = sortVectorBy(massFractionVector, speciesMass);
        m_cache.massFractions = massFractions; // Cache the result
        return massFractions;

    }

    std::vector<double> Composition::getNumberFractionVector() const noexcept {
        if (m_cache.numberFractions.has_value()) {
            return m_cache.numberFractions.value(); // Short circuit if we have cached the number fractions
        }

        std::vector<double> numberFractionVector;
        std::vector<double> speciesMass;

        numberFractionVector.reserve(m_molarAbundances.size());
        speciesMass.reserve(m_molarAbundances.size());

        for (const auto &species: m_molarAbundances | std::views::keys) {
            numberFractionVector.push_back(getNumberFraction(species));
            speciesMass.push_back(species.mass());
        }

        std::vector<double> numberFractions = sortVectorBy(numberFractionVector, speciesMass);
        m_cache.numberFractions = numberFractions; // Cache the result
        return numberFractions;
    }

    std::vector<double> Composition::getMolarAbundanceVector() const noexcept {
        if (m_cache.molarAbundances.has_value()) {
            return m_cache.molarAbundances.value(); // Short circuit if we have cached the molar abundances
        }

        std::vector<double> molarAbundanceVector;
        std::vector<double> speciesMass;

        molarAbundanceVector.reserve(m_molarAbundances.size());
        speciesMass.reserve(m_molarAbundances.size());

        for (const auto &[species, y]: m_molarAbundances) {
            molarAbundanceVector.push_back(y);
            speciesMass.push_back(species.mass());
        }

        std::vector<double> molarAbundances = sortVectorBy(molarAbundanceVector, speciesMass);
        m_cache.molarAbundances = molarAbundances; // Cache the result
        return molarAbundances;

    }

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
        if (!m_registeredSpecies.contains(species)) {
            LOG_ERROR(getLogger(), "Species {} is not in the composition.", species.name());
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

        speciesVector.reserve(m_molarAbundances.size());
        speciesMass.reserve(m_molarAbundances.size());

        for (const auto &s: m_registeredSpecies) {
            speciesVector.emplace_back(s);
            speciesMass.push_back(s.mass());
        }

        std::vector<atomic::Species> sortedSpecies = sortVectorBy(speciesVector, speciesMass);
        m_cache.sortedSpecies = sortedSpecies;
        return std::distance(sortedSpecies.begin(), std::ranges::find(sortedSpecies, species));
    }

    atomic::Species Composition::getSpeciesAtIndex(
        const size_t index
    ) const {
        if (m_cache.sortedSpecies.has_value()) {
            return m_cache.sortedSpecies.value().at(index);
        }

        std::vector<atomic::Species> speciesVector;
        std::vector<double> speciesMass;

        speciesVector.reserve(m_molarAbundances.size());
        speciesMass.reserve(m_molarAbundances.size());

        for (const auto &species: m_registeredSpecies) {
            speciesVector.emplace_back(species);
            speciesMass.push_back(species.mass());
        }

        std::vector<atomic::Species> sortedSymbols = sortVectorBy(speciesVector, speciesMass);
        if (index >= sortedSymbols.size()) {
            LOG_ERROR(getLogger(), "Index {} is out of range for composition of size {}.", index, sortedSymbols.size());
            throw std::out_of_range("Index " + std::to_string(index) + " is out of range for composition of size " + std::to_string(sortedSymbols.size()) + ".");
        }
        return sortedSymbols.at(index);
    }

    std::unique_ptr<CompositionAbstract> Composition::clone() const {
        return std::make_unique<Composition>(*this);
    }

    bool Composition::contains(
        const atomic::Species &species
    ) const noexcept {
        return m_registeredSpecies.contains(species);
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
        return m_registeredSpecies.size();
    }

    void Composition::setMolarAbundance(
        const std::string &symbol,
        const double &molar_abundance
    ) {
        const auto species = getSpecies(symbol);
        if (!species) {
            throw_unknown_symbol(getLogger(), symbol);
        }

        setMolarAbundance(species.value(), molar_abundance);
    }

    void Composition::setMolarAbundance(
        const atomic::Species &species,
        const double &molar_abundance
    ) {
        if (!m_registeredSpecies.contains(species)) {
            throw_unregistered_symbol(getLogger(), std::string(species.name()));
        }
        if (molar_abundance < 0.0) {
            LOG_ERROR(getLogger(), "Molar abundance must be non-negative for symbol {}. Currently it is {}.", species.name(), molar_abundance);
            throw exceptions::InvalidCompositionError("Molar abundance must be non-negative, got " + std::to_string(molar_abundance) + " for symbol " + std::string(species.name()) + ".");
        }
        m_molarAbundances.at(species) = molar_abundance;
    }

    void Composition::setMolarAbundance(
        const std::vector<std::string> &symbols,
        const std::vector<double> &molar_abundances
    ) {
        for (const auto& [symbol, y] : std::views::zip(symbols, molar_abundances)) {
            setMolarAbundance(symbol, y);
        }
    }

    void Composition::setMolarAbundance(
        const std::vector<atomic::Species> &species,
        const std::vector<double> &molar_abundances
    ) {
        for (const auto& [s, y] : std::views::zip(species, molar_abundances)) {
            setMolarAbundance(s, y);
        }
    }

    void Composition::setMolarAbundance(
        const std::set<std::string> &symbols,
        const std::vector<double> &molar_abundances
    ) {
        for (const auto& [symbol, y] : std::views::zip(symbols, molar_abundances)) {
            setMolarAbundance(symbol, y);
        }
    }

    void Composition::setMolarAbundance(
        const std::set<atomic::Species> &species,
        const std::vector<double> &molar_abundances
    ) {
        for (const auto& [s, y] : std::views::zip(species, molar_abundances)) {
            setMolarAbundance(s, y);
        }
    }

    /// OVERLOADS

    std::ostream& operator<<(
        std::ostream& os,
        const Composition& composition
    ) {
        os << "Composition(Mass Fractions => [";
        size_t count = 0;
        for (const auto &species : composition.m_registeredSpecies) {
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
