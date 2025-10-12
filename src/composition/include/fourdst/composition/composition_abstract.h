#pragma once

#include "fourdst/composition/atomicSpecies.h"

#include <string>
#include <unordered_map>
#include <set>

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
     * @brief Check if a chemical symbol is registered in the composition.
     * @param symbol The chemical symbol to check (e.g., "H", "He").
     * @return True if the symbol is present, false otherwise.
     */
    [[nodiscard]] virtual bool hasSymbol(const std::string& symbol) const = 0;

    /**
     * @brief Check if a species is registered in the composition.
     * @param species The atomic species to check.
     * @return True if the species is present, false otherwise.
     */
    [[nodiscard]] virtual bool hasSpecies(const fourdst::atomic::Species& species) const = 0;

    /**
     * @brief Check if the composition contains the given species.
     * @param species The atomic species to check.
     * @return True if the species is contained, false otherwise.
     */
    [[nodiscard]] virtual bool contains(const fourdst::atomic::Species& species) const = 0;

    /**
     * @brief Get all registered chemical symbols in the composition.
     * @return A set of registered chemical symbols.
     */
    [[nodiscard]] virtual std::set<std::string> getRegisteredSymbols() const = 0;

    /**
     * @brief Get all registered atomic species in the composition.
     * @return A set of registered atomic species.
     */
    [[nodiscard]] virtual std::set<fourdst::atomic::Species> getRegisteredSpecies() const = 0;

    /**
     * @brief Get the mass fraction for all registered symbols.
     * @return An unordered map from symbol to mass fraction.
     */
    [[nodiscard]] virtual std::unordered_map<std::string, double> getMassFraction() const = 0;

    /**
     * @brief Get the number fraction for all registered symbols.
     * @return An unordered map from symbol to number fraction.
     */
    [[nodiscard]] virtual std::unordered_map<std::string, double> getNumberFraction() const = 0;

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
    [[nodiscard]] virtual double getMeanParticleMass() const = 0;

    /**
     * @brief Get the mean atomic number of the composition.
     * @return The mean atomic number.
     */
    [[nodiscard]] virtual double getMeanAtomicNumber() const = 0;

    /**
     * @brief Get the electron abundance of the composition.
     * @return The electron abundance.
     */
    [[nodiscard]] virtual double getElectronAbundance() const = 0;

    /**
     * @brief Get the mass fraction as a vector.
     * @return A vector of mass fractions for all species.
     */
    [[nodiscard]] virtual std::vector<double> getMassFractionVector() const = 0;

    /**
     * @brief Get the number fraction as a vector.
     * @return A vector of number fractions for all species.
     */
    [[nodiscard]] virtual std::vector<double> getNumberFractionVector() const = 0;

    /**
     * @brief Get the molar abundance as a vector.
     * @return A vector of molar abundances for all species.
     */
    [[nodiscard]] virtual std::vector<double> getMolarAbundanceVector() const = 0;

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
    [[nodiscard]] virtual fourdst::atomic::Species getSpeciesAtIndex(size_t index) const = 0;
};