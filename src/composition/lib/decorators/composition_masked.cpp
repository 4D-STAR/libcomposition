#include "fourdst/composition/decorators/composition_masked.h"

#include "fourdst/atomic/species.h"
#include <memory>

namespace fourdst::composition {
MaskedComposition::MaskedComposition(
        const CompositionAbstract& baseComposition,
        const std::set<atomic::Species>& activeSpecies
    ) :
    CompositionDecorator(baseComposition.clone()),
    m_activeSpecies(activeSpecies) {
        for (const auto& species : m_activeSpecies) {
            if (CompositionDecorator::contains(species)) {
                m_masked_composition.emplace(species, CompositionDecorator::getMolarAbundance(species));
            } else {
                m_masked_composition.emplace(species, 0.0);
            }
        }
    }

    bool MaskedComposition::contains(const atomic::Species &species) const noexcept{
        if (m_activeSpecies.contains(species)) {
            return true;
        }
        return false;
    }

    bool MaskedComposition::contains(const std::string &symbol) const {
        if (!atomic::species.contains(symbol)) {
            throw exceptions::UnknownSymbolError("Cannot find species '" + symbol + "' in base composition");
        }
        const atomic::Species& species = atomic::species.at(symbol);
        if (m_activeSpecies.contains(species)) {
            return true;
        }

        return false;
    }

    const std::set<atomic::Species>& MaskedComposition::getRegisteredSpecies() const noexcept {
        return m_activeSpecies;
    }

    std::set<std::string> MaskedComposition::getRegisteredSymbols() const noexcept {
        std::set<std::string> symbols;
        for (const auto& species : m_activeSpecies) {
            symbols.insert(std::string(species.name()));
        }
        return symbols;
    }

    size_t MaskedComposition::size() const noexcept {
        return m_activeSpecies.size();
    }

    std::unordered_map<atomic::Species, double> MaskedComposition::getMassFraction() const noexcept {
        std::unordered_map<atomic::Species, double> massFractions;
        for (const auto& species : m_activeSpecies) {
            if (CompositionDecorator::contains(species)) {
                massFractions[species] = CompositionDecorator::getMassFraction(species);
            } else {
                massFractions[species] = 0.0;
            }
        }
        return massFractions;
    }

    std::unordered_map<atomic::Species, double> MaskedComposition::getNumberFraction() const noexcept {
        std::unordered_map<atomic::Species, double> numberFractions;
        for (const auto& species : m_activeSpecies) {
            if (CompositionDecorator::contains(species)) {
                numberFractions[species] = CompositionDecorator::getNumberFraction(species);
            } else {
                numberFractions[species] = 0.0;
            }
        }
        return numberFractions;
    }

    double MaskedComposition::getMassFraction(const std::string &symbol) const  {
        if (!contains(symbol)) {
            throw exceptions::UnregisteredSymbolError("Species '" + symbol + "' is not part of the active species in the MaskedComposition.");
        }
        if (CompositionDecorator::contains(symbol)) {
            return CompositionDecorator::getMassFraction(symbol);
        } return 0.0;
    }
    double MaskedComposition::getMassFraction(const atomic::Species &species) const {
        if (!contains(species)) {
            throw exceptions::UnregisteredSymbolError("Species '" + std::string(species.name()) + "' is not part of the active species in the MaskedComposition.");
        }
        if (CompositionDecorator::contains(species)) {
            return CompositionDecorator::getMassFraction(species);
        } return 0.0;
    }
    double MaskedComposition::getNumberFraction(const std::string &symbol) const {
        if (!contains(symbol)) {
            throw exceptions::UnregisteredSymbolError("Species '" + symbol + "' is not part of the active species in the MaskedComposition.");
        }
        if (CompositionDecorator::contains(symbol)) {
            return CompositionDecorator::getNumberFraction(symbol);
        } return 0.0;
    }
    double MaskedComposition::getNumberFraction(const atomic::Species &species) const {
        if (!contains(species)) {
            throw exceptions::UnregisteredSymbolError("Species '" + std::string(species.name()) + "' is not part of the active species in the MaskedComposition.");
        }
        if (CompositionDecorator::contains(species)) {
            return CompositionDecorator::getNumberFraction(species);
        } return 0.0;
    }
    double MaskedComposition::getMolarAbundance(const std::string &symbol) const {
        if (!contains(symbol)) {
            throw exceptions::UnregisteredSymbolError("Species '" + symbol + "' is not part of the active species in the MaskedComposition.");
        }
        if (CompositionDecorator::contains(symbol)) {
            return CompositionDecorator::getMolarAbundance(symbol);
        } return 0.0;
    }
    double MaskedComposition::getMolarAbundance(const atomic::Species &species) const {
        if (!contains(species)) {
            throw exceptions::UnregisteredSymbolError("Species '" + std::string(species.name()) + "' is not part of the active species in the MaskedComposition.");
        }
        if (CompositionDecorator::contains(species)) {
            return CompositionDecorator::getMolarAbundance(species);
        } return 0.0;
    }
    double MaskedComposition::getMeanParticleMass() const noexcept {
        double meanParticleMass = 0.0;
        for (const auto& species : m_activeSpecies) {
            if (CompositionDecorator::contains(species)) {
                const double numberFraction = CompositionDecorator::getNumberFraction(species);
                const double atomicMass = species.mass();
                meanParticleMass += numberFraction * atomicMass;
            }
        }
        return meanParticleMass;
    }

    double MaskedComposition::getElectronAbundance() const noexcept {
        double Ye = 0.0;
        for (const auto& species : m_activeSpecies) {
            if (CompositionDecorator::contains(species)) {
                Ye += CompositionDecorator::getMolarAbundance(species) * species.z();
            }
        }
        return Ye;
    }

    std::vector<double> MaskedComposition::getMassFractionVector() const noexcept {
        std::vector<double> massFractions;
        massFractions.reserve(m_activeSpecies.size());
        for (const auto& species : m_activeSpecies) {
            massFractions.push_back(getMassFraction(species));
        }
        return massFractions;
    }

    std::vector<double> MaskedComposition::getNumberFractionVector() const noexcept {
        std::vector<double> numberFractions;
        numberFractions.reserve(m_activeSpecies.size());
        for (const auto& species : m_activeSpecies) {
            numberFractions.push_back(getNumberFraction(species));
        }
        return numberFractions;
    }

    std::vector<double> MaskedComposition::getMolarAbundanceVector() const noexcept {
        std::vector<double> molarAbundances;
        molarAbundances.reserve(m_activeSpecies.size());
        for (const auto& species : m_activeSpecies) {
            molarAbundances.push_back(getMolarAbundance(species));
        }
        return molarAbundances;
    }

    size_t MaskedComposition::getSpeciesIndex(const std::string &symbol) const {
        if (!contains(symbol)) {
            throw exceptions::UnregisteredSymbolError("Species '" + symbol + "' is not part of the active species in the MaskedComposition.");
        }
        return std::distance(m_activeSpecies.begin(), m_activeSpecies.find(atomic::species.at(symbol)));
    }

    size_t MaskedComposition::getSpeciesIndex(const atomic::Species &species) const {
        return std::distance(m_activeSpecies.begin(), m_activeSpecies.find(species));
    }

    atomic::Species MaskedComposition::getSpeciesAtIndex(const size_t index) const {
        if (index >= m_activeSpecies.size()) {
            throw std::out_of_range("Index " + std::to_string(index) + " is out of bounds for active species of size " + std::to_string(m_activeSpecies.size()) + ".");
        }
        auto it = m_activeSpecies.begin();
        std::advance(it, index);
        return *it;
    }

    std::unique_ptr<CompositionAbstract> MaskedComposition::clone() const {
        return std::make_unique<MaskedComposition>(*m_base_composition, m_activeSpecies);
    }

    std::map<atomic::Species, double>::iterator MaskedComposition::begin() {
        return m_masked_composition.begin();
    }

    std::map<atomic::Species, double>::iterator MaskedComposition::end() {
        return m_masked_composition.end();
    }

    std::map<atomic::Species, double>::const_iterator MaskedComposition::begin() const {
        return m_masked_composition.cbegin();
    }

    std::map<atomic::Species, double>::const_iterator MaskedComposition::end() const {
        return m_masked_composition.cend();
    }
};
