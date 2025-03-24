#include "composition.h"
#include "quill/LogMacros.h"

#include <stdexcept>

#include "atomicSpecies.h"

using namespace composition;

Composition::Composition(const std::vector<std::string>& symbols) {
    for (const auto& symbol : symbols) {
        registerSymbol(symbol);
    }
}

void Composition::validateComposition(const std::vector<double>& mass_fractions) const {
    if (!isValidComposition(mass_fractions)) {
        LOG_ERROR(m_logger, "Invalid composition.");
        throw std::runtime_error("Invalid composition.");
    }
}

bool Composition::isValidComposition(const std::vector<double>& mass_fractions) const {
    double sum = 0.0;
    for (const auto& mass_fraction : mass_fractions) {
        sum += mass_fraction;
    }
    if (sum < 0.99999999 || sum > 1.000000001) {
        LOG_ERROR(m_logger, "The sum of mass fractions must be equal to 1.");
        return false;
    }

    return true;
}

Composition::Composition(const std::vector<std::string>& symbols, const std::vector<double>& mass_fractions) {
    if (symbols.size() != mass_fractions.size()) {
        LOG_ERROR(m_logger, "The number of symbols and mass fractions must be equal.");
        throw std::runtime_error("The number of symbols and mass fractions must be equal.");
    }

    validateComposition(mass_fractions);
    finalize();

    for (const auto &symbol : symbols) {
        registerSymbol(symbol);
    }

    for (size_t i = 0; i < symbols.size(); ++i) {
        setComposition(symbols[i], mass_fractions[i]);
    }
}

void Composition::registerSymbol(const std::string& symbol) {
    if (!isValidSymbol(symbol)) {
        LOG_ERROR(m_logger, "Invalid symbol: {}", symbol);
        throw std::runtime_error("Invalid symbol.");
    }

    if (m_registeredSymbols.find(symbol) != m_registeredSymbols.end()) {
        LOG_WARNING(m_logger, "Symbol {} is already registered.", symbol);
        return;
    }

    m_registeredSymbols.insert(symbol);
    LOG_INFO(m_logger, "Registered symbol: {}", symbol);
}

bool Composition::isValidSymbol(const std::string& symbol) const {
    return chemSpecies::species.count(symbol) > 0;
}

double Composition::setComposition(const std::string& symbol, const double& mass_fraction) {
    m_finalized = false;
    if (m_registeredSymbols.find(symbol) == m_registeredSymbols.end()) {
        LOG_ERROR(m_logger, "Symbol {} is not registered.", symbol);
        throw std::runtime_error("Symbol is not registered.");
    }

    if (mass_fraction < 0.0 || mass_fraction > 1.0) {
        LOG_ERROR(m_logger, "Mass fraction must be between 0 and 1.");
        throw std::runtime_error("Mass fraction must be between 0 and 1.");
    }

    double old_mass_fraction = 0.0;
    if (m_compositions.find(symbol) != m_compositions.end()) {
        old_mass_fraction = m_compositions[symbol].mass_fraction;
    }
    m_compositions[symbol] = {symbol, mass_fraction};


    return old_mass_fraction;
}

std::vector<double> Composition::setComposition(const std::vector<std::string>& symbols, const std::vector<double>& mass_fractions) {
    m_finalized = false;
    if (symbols.size() != mass_fractions.size()) {
        LOG_ERROR(m_logger, "The number of symbols and mass fractions must be equal.");
        throw std::runtime_error("The number of symbols and mass fractions must be equal.");
    }

    std::vector<double> old_mass_fractions;
    old_mass_fractions.reserve(symbols.size());
    for (size_t i = 0; i < symbols.size(); ++i) {
        old_mass_fractions.push_back(setComposition(symbols[i], mass_fractions[i]));
    }
    return old_mass_fractions;
}

std::set<std::string> Composition::getRegisteredSymbols() const {
    return m_registeredSymbols;
}

std::unordered_map<std::string, CompositionEntry> Composition::getCompositions() const {
    if (!m_finalized) {
        LOG_ERROR(m_logger, "Composition has not been finalized.");
        throw std::runtime_error("Composition has not been finalized (Consider running .finalize()).");
    }
    return m_compositions;
}

bool Composition::finalize() {
    m_finalized = true;
    std::vector<double> mass_fractions;
    mass_fractions.reserve(m_compositions.size());
    for (const auto& [_, entry] : m_compositions) {
        mass_fractions.push_back(entry.mass_fraction);
    }
    try {
        validateComposition(mass_fractions);
    } catch (const std::runtime_error& e) {
        double massSum = 0.0;
        for (const auto& [_, entry] : m_compositions) {
            massSum += entry.mass_fraction;
        }
        LOG_ERROR(m_logger, "Composition is invalid (Total mass {}).", massSum);
        m_finalized = false;
        return false;
    }
    return true;
}

CompositionEntry Composition::getComposition(const std::string& symbol) const {
    if (!m_finalized) {
        LOG_ERROR(m_logger, "Composition has not been finalized.");
        throw std::runtime_error("Composition has not been finalized (Consider running .finalize()).");
    }
    if (m_compositions.count(symbol) == 0) {
        LOG_ERROR(m_logger, "Symbol {} is not in the composition.", symbol);
        throw std::runtime_error("Symbol is not in the composition.");
    }
    return m_compositions.at(symbol);
}