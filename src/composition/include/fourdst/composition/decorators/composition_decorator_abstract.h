#pragma once

#include "fourdst/atomic/atomicSpecies.h"

#include "fourdst/composition/composition_abstract.h"

#include <utility>
#include <set>
#include <unordered_map>
#include <map>
#include <vector>

namespace fourdst::composition {

    class CompositionDecorator: public CompositionAbstract {
    public:
        explicit CompositionDecorator(std::unique_ptr<CompositionAbstract> decorator) : m_base_composition(std::move(decorator)) {};
        [[nodiscard]] bool contains(const atomic::Species &species) const noexcept override { return m_base_composition->contains(species); };
        [[nodiscard]] bool contains(const std::string& symbol) const override { return m_base_composition->contains(symbol); };
        [[nodiscard]] size_t size() const noexcept override { return m_base_composition->size(); };
        [[nodiscard]] std::set<std::string> getRegisteredSymbols() const noexcept override { return m_base_composition->getRegisteredSymbols(); };
        [[nodiscard]] const std::set<atomic::Species> &getRegisteredSpecies() const noexcept override { return m_base_composition->getRegisteredSpecies(); };
        [[nodiscard]] std::unordered_map<atomic::Species, double> getMassFraction() const noexcept override { return m_base_composition->getMassFraction(); };
        [[nodiscard]] std::unordered_map<atomic::Species, double> getNumberFraction() const noexcept override { return m_base_composition->getNumberFraction(); };
        [[nodiscard]] double getMassFraction(const std::string& symbol) const override { return m_base_composition->getMassFraction(symbol); };
        [[nodiscard]] double getMassFraction(const atomic::Species& species) const override { return m_base_composition->getMassFraction(species); };
        [[nodiscard]] double getNumberFraction(const std::string& symbol) const override { return m_base_composition->getNumberFraction(symbol); };
        [[nodiscard]] double getNumberFraction(const atomic::Species& species) const override { return m_base_composition->getNumberFraction(species); };
        [[nodiscard]] double getMolarAbundance(const std::string& symbol) const override { return m_base_composition->getMolarAbundance(symbol); };
        [[nodiscard]] double getMolarAbundance(const atomic::Species& species) const override { return m_base_composition->getMolarAbundance(species); };
        [[nodiscard]] double getMeanParticleMass() const noexcept override { return m_base_composition->getMeanParticleMass(); };
        [[nodiscard]] double getElectronAbundance() const noexcept override { return m_base_composition->getElectronAbundance(); };
        [[nodiscard]] std::vector<double> getMassFractionVector() const noexcept override { return m_base_composition->getMassFractionVector(); };
        [[nodiscard]] std::vector<double> getNumberFractionVector() const noexcept override { return m_base_composition->getNumberFractionVector(); };
        [[nodiscard]] std::vector<double> getMolarAbundanceVector() const noexcept override { return m_base_composition->getMolarAbundanceVector(); };
        [[nodiscard]] size_t getSpeciesIndex(const std::string& symbol) const override { return m_base_composition->getSpeciesIndex(symbol); };
        [[nodiscard]] size_t getSpeciesIndex(const atomic::Species& species) const override { return m_base_composition->getSpeciesIndex(species); };
        [[nodiscard]] atomic::Species getSpeciesAtIndex(const size_t index) const override { return m_base_composition->getSpeciesAtIndex(index); }
        [[nodiscard]] size_t hash() const override { return m_base_composition->hash(); };

        [[nodiscard]] std::map<atomic::Species, double>::iterator begin() override { return m_base_composition->begin(); };
        [[nodiscard]] std::map<atomic::Species, double>::iterator end() override { return m_base_composition->end(); };

        [[nodiscard]] std::map<atomic::Species, double>::const_iterator begin() const override { return std::as_const(*m_base_composition).begin(); };
        [[nodiscard]] std::map<atomic::Species, double>::const_iterator end() const override { return std::as_const(*m_base_composition).end(); };
    protected:
        std::unique_ptr<CompositionAbstract> m_base_composition;
    };
}