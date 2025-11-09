#pragma once

#include "fourdst/atomic/atomicSpecies.h"

#include <string>
#include <unordered_map>
#include <set>
#include <vector>

namespace fourdst::composition {
    /**
     * @brief Abstract base class for chemical composition representations.
     *
     * The purpose of this class is to define a standard interface for all composition types.
     * Children of this class are responsible for implementing the setter methods, but any
     * object that is a child of CompositionAbstract will always have these getter methods.
     *
     * This ensures that all derived composition classes provide a consistent API for querying
     * composition properties, regardless of how the data is set or stored.
     *
     * @par Example
     * @code
     * class MyComposition : public CompositionAbstract {
     *     // ...implement all pure virtual methods...
     * };
     *
     * MyComposition comp;
     * if (comp.hasSymbol("H")) {
     *     double mf = comp.getMassFraction("H");
     * }
     * std::set<std::string> symbols = comp.getRegisteredSymbols();
     * @endcode
     */
    class CompositionAbstract {
    public:
        /**
         * @brief Virtual destructor.
         */
        virtual ~CompositionAbstract() = default;

        /**
         * @brief Check if the composition contains the given species.
         * @param species The atomic species to check.
         * @return True if the species is contained, false otherwise.
         */
        [[nodiscard]] virtual bool contains(const fourdst::atomic::Species& species) const noexcept = 0;

        /**
         * @brief Check if the composition contains the given species.
         * @param symbol The symbol of the atomic species to check.
         * @return True if the species is contained, false otherwise.
         */
        [[nodiscard]] virtual bool contains(const std::string& symbol) const = 0;

        [[nodiscard]] virtual size_t size() const noexcept = 0;

        /**
         * @brief Get all registered chemical symbols in the composition.
         * @return A set of registered chemical symbols.
         */
        [[nodiscard]] virtual std::set<std::string> getRegisteredSymbols() const noexcept = 0;

        /**
         * @brief Get all registered atomic species in the composition.
         * @return A set of registered atomic species.
         */
        [[nodiscard]] virtual const std::set<fourdst::atomic::Species> &getRegisteredSpecies() const noexcept = 0;

        /**
         * @brief Get the mass fraction for all registered symbols.
         * @return An unordered map from symbol to mass fraction.
         */
        [[nodiscard]] virtual std::unordered_map<fourdst::atomic::Species, double> getMassFraction() const noexcept = 0;

        /**
         * @brief Get the number fraction for all registered symbols.
         * @return An unordered map from symbol to number fraction.
         */
        [[nodiscard]] virtual std::unordered_map<fourdst::atomic::Species, double> getNumberFraction() const noexcept = 0;

        /**
         * @brief Get the mass fraction for a given symbol.
         * @param symbol The chemical symbol.
         * @return The mass fraction for the symbol.
         */
        [[nodiscard]] virtual double getMassFraction(const std::string& symbol) const = 0;

        /**
         * @brief Get the mass fraction for a given species.
         * @param species The atomic species.
         * @return The mass fraction for the species.
         */
        [[nodiscard]] virtual double getMassFraction(const fourdst::atomic::Species& species) const = 0;

        /**
         * @brief Get the number fraction for a given symbol.
         * @param symbol The chemical symbol.
         * @return The number fraction for the symbol.
         */
        [[nodiscard]] virtual double getNumberFraction(const std::string& symbol) const = 0;

        /**
         * @brief Get the number fraction for a given species.
         * @param species The atomic species.
         * @return The number fraction for the species.
         */
        [[nodiscard]] virtual double getNumberFraction(const fourdst::atomic::Species& species) const = 0;

        /**
         * @brief Get the molar abundance for a given symbol.
         * @param symbol The chemical symbol.
         * @return The molar abundance for the symbol.
         */
        [[nodiscard]] virtual double getMolarAbundance(const std::string& symbol) const = 0;

        /**
         * @brief Get the molar abundance for a given species.
         * @param species The atomic species.
         * @return The molar abundance for the species.
         */
        [[nodiscard]] virtual double getMolarAbundance(const fourdst::atomic::Species& species) const = 0;

        /**
         * @brief Get the mean particle mass of the composition.
         * @return The mean particle mass.
         */
        [[nodiscard]] virtual double getMeanParticleMass() const noexcept = 0;

        /**
         * @brief Get the electron abundance of the composition.
         * @return The electron abundance.
         */
        [[nodiscard]] virtual double getElectronAbundance() const noexcept = 0;

        /**
         * @brief Get the mass fraction as a vector.
         * @return A vector of mass fractions for all species.
         */
        [[nodiscard]] virtual std::vector<double> getMassFractionVector() const noexcept = 0;

        /**
         * @brief Get the number fraction as a vector.
         * @return A vector of number fractions for all species.
         */
        [[nodiscard]] virtual std::vector<double> getNumberFractionVector() const noexcept = 0;

        /**
         * @brief Get the molar abundance as a vector.
         * @return A vector of molar abundances for all species.
         */
        [[nodiscard]] virtual std::vector<double> getMolarAbundanceVector() const noexcept = 0;

        /**
         * @brief Get the index of a species by symbol.
         * @param symbol The chemical symbol.
         * @return The index of the species.
         */
        [[nodiscard]] virtual size_t getSpeciesIndex(const std::string& symbol) const = 0;

        /**
         * @brief Get the index of a species.
         * @param species The atomic species.
         * @return The index of the species.
         */
        [[nodiscard]] virtual size_t getSpeciesIndex(const fourdst::atomic::Species& species) const = 0;

        /**
         * @brief Get the species at a given index.
         * @param index The index of the species.
         * @return The atomic species at the specified index.
         */
        [[nodiscard]] virtual atomic::Species getSpeciesAtIndex(size_t index) const = 0;

        [[nodiscard]] virtual std::unique_ptr<CompositionAbstract> clone() const = 0;

        [[nodiscard]] virtual std::map<atomic::Species, double>::iterator begin() = 0;
        [[nodiscard]] virtual std::map<atomic::Species, double>::iterator end() = 0;
        [[nodiscard]] virtual std::map<atomic::Species, double>::const_iterator begin() const = 0;
        [[nodiscard]] virtual std::map<atomic::Species, double>::const_iterator end() const = 0;
    };

    // ReSharper disable once CppClassCanBeFinal
    class CompositionDecorator : public CompositionAbstract {
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
        [[nodiscard]] atomic::Species getSpeciesAtIndex(const size_t index) const override { return m_base_composition->getSpeciesAtIndex(index); };
        [[nodiscard]] std::unique_ptr<CompositionAbstract> clone() const override {
            return std::make_unique<CompositionDecorator>(m_base_composition->clone());
        }
        [[nodiscard]] std::map<atomic::Species, double>::iterator begin() override {
            return m_base_composition->begin();
        }
        [[nodiscard]] std::map<atomic::Species, double>::iterator end() override {
            return m_base_composition->end();
        }
        [[nodiscard]] std::map<atomic::Species, double>::const_iterator begin() const override {
            return std::as_const(*m_base_composition).begin();
        }
        [[nodiscard]] std::map<atomic::Species, double>::const_iterator end() const override {
            return std::as_const(*m_base_composition).end();
        }
    private:
        std::unique_ptr<CompositionAbstract> m_base_composition;
    };
}