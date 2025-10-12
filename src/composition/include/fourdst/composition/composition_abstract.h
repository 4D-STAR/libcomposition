#pragma once

#include "fourdst/composition/atomicSpecies.h"

#include <string>
#include <unordered_map>
#include <set>

class CompositionAbstract {
public:
    virtual ~CompositionAbstract() = default;
    [[nodiscard]] virtual bool hasSymbol(const std::string& symbol) const = 0;
    [[nodiscard]] virtual bool hasSpecies(const fourdst::atomic::Species& species) const = 0;

    [[nodiscard]] virtual bool contains(const fourdst::atomic::Species& species) const = 0;

    [[nodiscard]] virtual std::set<std::string> getRegisteredSymbols() const = 0;
    [[nodiscard]] virtual std::set<fourdst::atomic::Species> getRegisteredSpecies() const = 0;

    [[nodiscard]] virtual std::unordered_map<std::string, double> getMassFraction() const = 0;
    [[nodiscard]] virtual std::unordered_map<std::string, double> getNumberFraction() const = 0;

    [[nodiscard]] virtual double getMassFraction(const std::string& symbol) const = 0;
    [[nodiscard]] virtual double getMassFraction(const fourdst::atomic::Species& species) const = 0;

    [[nodiscard]] virtual double getNumberFraction(const std::string& symbol) const = 0;
    [[nodiscard]] virtual double getNumberFraction(const fourdst::atomic::Species& species) const = 0;

    [[nodiscard]] virtual double getMolarAbundance(const std::string& symbol) const = 0;
    [[nodiscard]] virtual double getMolarAbundance(const fourdst::atomic::Species& species) const = 0;

    [[nodiscard]] virtual double getMeanParticleMass() const = 0;
    [[nodiscard]] virtual double getMeanAtomicNumber() const = 0;
    [[nodiscard]] virtual double getElectronAbundance() const = 0;

    [[nodiscard]] virtual std::vector<double> getMassFractionVector() const = 0;
    [[nodiscard]] virtual std::vector<double> getNumberFractionVector() const = 0;
    [[nodiscard]] virtual std::vector<double> getMolarAbundanceVector() const = 0;

    [[nodiscard]] virtual size_t getSpeciesIndex(const std::string& symbol) const = 0;
    [[nodiscard]] virtual size_t getSpeciesIndex(const fourdst::atomic::Species& species) const = 0;

    [[nodiscard]] virtual fourdst::atomic::Species getSpeciesAtIndex(size_t index) const = 0;
};