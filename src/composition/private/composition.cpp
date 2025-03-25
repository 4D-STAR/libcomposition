#include "composition.h"
#include "quill/LogMacros.h"

#include <stdexcept>
#include <unordered_map>
#include <vector>
#include <array>

#include <utility>

#include "atomicSpecies.h"

using namespace composition;

CompositionEntry::CompositionEntry() : m_symbol("H-1"), m_isotope(chemSpecies::species.at("H-1")), m_initialized(false) {}

CompositionEntry::CompositionEntry(const std::string& symbol, bool massFracMode) : m_symbol(symbol), m_isotope(chemSpecies::species.at(symbol)), m_massFracMode(massFracMode) {
    setSpecies(symbol);
}

CompositionEntry::CompositionEntry(const CompositionEntry& entry) : 
    m_symbol(entry.m_symbol), 
    m_isotope(entry.m_isotope), 
    m_massFracMode(entry.m_massFracMode), 
    m_massFraction(entry.m_massFraction), 
    m_numberFraction(entry.m_numberFraction), 
    m_relAbundance(entry.m_relAbundance),
    m_initialized(entry.m_initialized) {}

void CompositionEntry::setSpecies(const std::string& symbol) {
    if (m_initialized) {
        throw std::runtime_error("Composition entry is already initialized.");
    }
    if (chemSpecies::species.count(symbol) == 0) {
        throw std::runtime_error("Invalid symbol.");
    }
    m_symbol = symbol;
    m_isotope = chemSpecies::species.at(symbol);
    m_initialized = true;
}

std::string CompositionEntry::symbol() const {
    return m_symbol;
}

double CompositionEntry::mass_fraction() const {
    if (!m_massFracMode) {
        throw std::runtime_error("Composition entry is in number fraction mode.");
    }
    return m_massFraction;
}

double CompositionEntry::mass_fraction(double meanMolarMass) const {
    if (m_massFracMode) {
        return m_massFraction;
    }
    return m_relAbundance / meanMolarMass;
}


double CompositionEntry::number_fraction() const {
    if (m_massFracMode) {
        throw std::runtime_error("Composition entry is in mass fraction mode.");
    }
    return m_numberFraction;
}

double CompositionEntry::number_fraction(double totalMoles) const {
    if (m_massFracMode) {
        return m_relAbundance / totalMoles;
    }
    return m_numberFraction;
}

double CompositionEntry::rel_abundance() const {
    return m_relAbundance;
}

chemSpecies::Species CompositionEntry::isotope() const {
    return m_isotope;
}

void CompositionEntry::setMassFraction(double mass_fraction) {
    if (!m_massFracMode) {
        throw std::runtime_error("Composition entry is in number fraction mode.");
    }
    m_massFraction = mass_fraction;
    m_relAbundance = m_massFraction / m_isotope.mass();
}

void CompositionEntry::setNumberFraction(double number_fraction) {
    if (m_massFracMode) {
        throw std::runtime_error("Composition entry is in mass fraction mode.");
    }
    m_numberFraction = number_fraction;
    m_relAbundance = m_numberFraction * m_isotope.mass();
}

bool CompositionEntry::setMassFracMode(double meanParticleMass) {
    if (m_massFracMode) {
        return false;
    }
    m_massFracMode = true;
    m_massFraction = m_relAbundance / meanParticleMass;
    return true;
}

bool CompositionEntry::setNumberFracMode(double specificNumberDensity) {
    if (!m_massFracMode) {
        return false;
    }
    m_massFracMode = false;
    m_numberFraction = m_relAbundance / specificNumberDensity;
    return true;
}

bool CompositionEntry::getMassFracMode() const {
    return m_massFracMode;
}

Composition::Composition(const std::vector<std::string>& symbols) {
    for (const auto& symbol : symbols) {
        registerSymbol(symbol);
    }
}

Composition::Composition(const std::vector<std::string>& symbols, const std::vector<double>& fractions, bool massFracMode) : m_massFracMode(massFracMode) {
    if (symbols.size() != fractions.size()) {
        LOG_ERROR(m_logger, "The number of symbols and fractions must be equal.");
        throw std::runtime_error("The number of symbols and fractions must be equal.");
    }

    validateComposition(fractions);

    for (const auto &symbol : symbols) {
        registerSymbol(symbol);
    }

    for (size_t i = 0; i < symbols.size(); ++i) {
        if (m_massFracMode) {
            setMassFraction(symbols[i], fractions[i]);
        } else {
            setNumberFraction(symbols[i], fractions[i]);
        }
    }
    finalize();
}

void Composition::registerSymbol(const std::string& symbol, bool massFracMode) {
    if (!isValidSymbol(symbol)) {
        LOG_ERROR(m_logger, "Invalid symbol: {}", symbol);
        throw std::runtime_error("Invalid symbol.");
    }

    // If no symbols have been registered allow mode to be set
    if (m_registeredSymbols.size() == 0) {
        m_massFracMode = massFracMode;
    } else {
        if (m_massFracMode != massFracMode) {
            LOG_ERROR(m_logger, "Composition is in mass fraction mode. Cannot register symbol in number fraction mode.");
            throw std::runtime_error("Composition is in mass fraction mode. Cannot register symbol in number fraction mode.");
        }
    }

    if (m_registeredSymbols.find(symbol) != m_registeredSymbols.end()) {
        LOG_WARNING(m_logger, "Symbol {} is already registered.", symbol);
        return;
    }

    m_registeredSymbols.insert(symbol);
    CompositionEntry entry(symbol, m_massFracMode);
    m_compositions[symbol] = entry;
    LOG_INFO(m_logger, "Registered symbol: {}", symbol);
}

void Composition::registerSymbol(const std::vector<std::string>& symbols, bool massFracMode) {
    for (const auto& symbol : symbols) {
        registerSymbol(symbol, massFracMode);
    }
}

std::set<std::string> Composition::getRegisteredSymbols() const {
    return m_registeredSymbols;
}

void Composition::validateComposition(const std::vector<double>& fractions) const {
    if (!isValidComposition(fractions)) {
        LOG_ERROR(m_logger, "Invalid composition.");
        throw std::runtime_error("Invalid composition.");
    }
}

bool Composition::isValidComposition(const std::vector<double>& fractions) const {
    double sum = 0.0;
    for (const auto& fraction : fractions) {
        sum += fraction;
    }
    if (sum < 0.999999 || sum > 1.000001) {
        LOG_ERROR(m_logger, "The sum of fractions must be equal to 1.");
        return false;
    }

    return true;
}

bool Composition::isValidSymbol(const std::string& symbol) const {
    return chemSpecies::species.count(symbol) > 0;
}

double Composition::setMassFraction(const std::string& symbol, const double& mass_fraction) {
    if (m_registeredSymbols.find(symbol) == m_registeredSymbols.end()) {
        LOG_ERROR(m_logger, "Symbol {} is not registered.", symbol);
        throw std::runtime_error("Symbol is not registered.");
    }

    if (!m_massFracMode) {
        LOG_ERROR(m_logger, "Composition is in number fraction mode.");
        throw std::runtime_error("Composition is in number fraction mode.");
    }

    if (mass_fraction < 0.0 || mass_fraction > 1.0) {
        LOG_ERROR(m_logger, "Mass fraction must be between 0 and 1 for symbol {}. Currently it is {}.", symbol, mass_fraction);
        throw std::runtime_error("Mass fraction must be between 0 and 1.");
    }

    m_finalized = false;
    double old_mass_fraction = m_compositions.at(symbol).mass_fraction();
    m_compositions.at(symbol).setMassFraction(mass_fraction);

    return old_mass_fraction;
}

std::vector<double> Composition::setMassFraction(const std::vector<std::string>& symbols, const std::vector<double>& mass_fractions) {
    if (symbols.size() != mass_fractions.size()) {
        LOG_ERROR(m_logger, "The number of symbols and mass fractions must be equal.");
        throw std::runtime_error("The number of symbols and mass fractions must be equal.");
    }

    std::vector<double> old_mass_fractions;
    old_mass_fractions.reserve(symbols.size());
    for (size_t i = 0; i < symbols.size(); ++i) {
        old_mass_fractions.push_back(setMassFraction(symbols[i], mass_fractions[i]));
    }
    return old_mass_fractions;
}

double Composition::setNumberFraction(const std::string& symbol, const double& number_fraction) {
    if (m_registeredSymbols.find(symbol) == m_registeredSymbols.end()) {
        LOG_ERROR(m_logger, "Symbol {} is not registered.", symbol);
        throw std::runtime_error("Symbol is not registered.");
    }

    if (m_massFracMode) {
        LOG_ERROR(m_logger, "Composition is in mass fraction mode.");
        throw std::runtime_error("Composition is in mass fraction mode.");
    }

    if (number_fraction < 0.0 || number_fraction > 1.0) {
        LOG_ERROR(m_logger, "Number fraction must be between 0 and 1 for symbol {}. Currently it is {}.", symbol, number_fraction);
        throw std::runtime_error("Number fraction must be between 0 and 1.");
    }

    m_finalized = false;
    double old_number_fraction = m_compositions.at(symbol).number_fraction();
    m_compositions.at(symbol).setNumberFraction(number_fraction);

    return old_number_fraction;
}

std::vector<double> Composition::setNumberFraction(const std::vector<std::string>& symbols, const std::vector<double>& number_fractions) {
    if (symbols.size() != number_fractions.size()) {
        LOG_ERROR(m_logger, "The number of symbols and number fractions must be equal.");
        throw std::runtime_error("The number of symbols and number fractions must be equal.");
    }

    std::vector<double> old_number_fractions;
    old_number_fractions.reserve(symbols.size());
    for (size_t i = 0; i < symbols.size(); ++i) {
        old_number_fractions.push_back(setNumberFraction(symbols[i], number_fractions[i]));
    }
    return old_number_fractions;
}

bool Composition::finalize(bool norm) {
    bool finalized = false;
    if (m_massFracMode) {
        finalized = finalizeMassFracMode(norm);
    } else {
        finalized = finalizeNumberFracMode(norm);
    }
    if (finalized) {
        m_finalized = true;
    }
    return finalized;
}

bool Composition::finalizeMassFracMode(bool norm) {
    std::vector<double> mass_fractions;
    mass_fractions.reserve(m_compositions.size());
    for (const auto& [_, entry] : m_compositions) {
        mass_fractions.push_back(entry.mass_fraction());
    }
    if (norm) {
        double sum = 0.0;
        for (const auto& mass_fraction : mass_fractions) {
            sum += mass_fraction;
        }
        for (int i = 0; i < mass_fractions.size(); ++i) {
            mass_fractions[i] /= sum;
        }
        for (auto& [symbol, entry] : m_compositions) {
            setMassFraction(symbol, entry.mass_fraction() / sum);
        }
    }
    try {
        validateComposition(mass_fractions);
    } catch (const std::runtime_error& e) {
        double massSum = 0.0;
        for (const auto& [_, entry] : m_compositions) {
            massSum += entry.mass_fraction();
        }
        LOG_ERROR(m_logger, "Composition is invalid (Total mass {}).", massSum);
        m_finalized = false;
        return false;
    }
    for (const auto& [_, entry] : m_compositions) {
        m_specificNumberDensity += entry.rel_abundance();
    }
    m_meanParticleMass = 1.0/m_specificNumberDensity;
    return true;
}

bool Composition::finalizeNumberFracMode(bool norm) {
    std::vector<double> number_fractions;
    number_fractions.reserve(m_compositions.size());
    for (const auto& [_, entry] : m_compositions) {
        number_fractions.push_back(entry.number_fraction());
    }
    if (norm) {
        double sum = 0.0;
        for (const auto& number_fraction : number_fractions) {
            sum += number_fraction;
        }
        for (auto& [symbol, entry] : m_compositions) {
            setNumberFraction(symbol, entry.number_fraction() / sum);
        }
    }
    try {
        validateComposition(number_fractions);
    } catch (const std::runtime_error& e) {
        double numberSum = 0.0;
        for (const auto& [_, entry] : m_compositions) {
            numberSum += entry.number_fraction();
        }
        LOG_ERROR(m_logger, "Composition is invalid (Total number {}).", numberSum);
        m_finalized = false;
        return false;
    }
    for (const auto& [_, entry] : m_compositions) {
        m_meanParticleMass += entry.rel_abundance();
    }
    m_specificNumberDensity = 1.0/m_meanParticleMass;
    return true;
}

double Composition::getMassFraction(const std::string& symbol) const {
    if (!m_finalized) {
        LOG_ERROR(m_logger, "Composition has not been finalized.");
        throw std::runtime_error("Composition has not been finalized (Consider running .finalize()).");
    }
    if (m_compositions.count(symbol) == 0) {
        LOG_ERROR(m_logger, "Symbol {} is not in the composition.", symbol);
        throw std::runtime_error("Symbol is not in the composition.");
    }
    if (m_massFracMode) {
        return m_compositions.at(symbol).mass_fraction();
    } else {
        return m_compositions.at(symbol).mass_fraction(m_meanParticleMass);
    }
}

std::unordered_map<std::string, double> Composition::getMassFraction() const {
    std::unordered_map<std::string, double> mass_fractions;
    for (const auto& [symbol, entry] : m_compositions) {
        mass_fractions[symbol] = getMassFraction(symbol);
    }
    return mass_fractions;
}


double Composition::getNumberFraction(const std::string& symbol) const {
    if (!m_finalized) {
        LOG_ERROR(m_logger, "Composition has not been finalized.");
        throw std::runtime_error("Composition has not been finalized (Consider running .finalize()).");
    }
    if (m_compositions.count(symbol) == 0) {
        LOG_ERROR(m_logger, "Symbol {} is not in the composition.", symbol);
        throw std::runtime_error("Symbol is not in the composition.");
    }
    if (!m_massFracMode) {
        return m_compositions.at(symbol).number_fraction();
    } else {
        return m_compositions.at(symbol).number_fraction(m_specificNumberDensity);
    }
}

std::unordered_map<std::string, double> Composition::getNumberFraction() const {
    std::unordered_map<std::string, double> number_fractions;
    for (const auto& [symbol, entry] : m_compositions) {
        number_fractions[symbol] = getNumberFraction(symbol);
    }
    return number_fractions;
}

std::pair<CompositionEntry, GlobalComposition> Composition::getComposition(const std::string& symbol) const {
    if (!m_finalized) {
        LOG_ERROR(m_logger, "Composition has not been finalized.");
        throw std::runtime_error("Composition has not been finalized (Consider running .finalize()).");
    }
    if (m_compositions.count(symbol) == 0) {
        LOG_ERROR(m_logger, "Symbol {} is not in the composition.", symbol);
        throw std::runtime_error("Symbol is not in the composition.");
    }
    return {m_compositions.at(symbol), {m_specificNumberDensity, m_meanParticleMass}};
}

std::pair<std::unordered_map<std::string, CompositionEntry>, GlobalComposition> Composition::getComposition() const {
    if (!m_finalized) {
        LOG_ERROR(m_logger, "Composition has not been finalized.");
        throw std::runtime_error("Composition has not been finalized (Consider running .finalize()).");
    }
    return {m_compositions, {m_specificNumberDensity, m_meanParticleMass}};
}

Composition Composition::subset(const std::vector<std::string>& symbols, std::string method) const {
    std::array<std::string, 2> methods = {"norm", "none"};

    if (std::find(methods.begin(), methods.end(), method) == methods.end()) {
        std::string errorMessage = "Invalid method: " + method + ". Valid methods are 'norm' and 'none'.";
        LOG_ERROR(m_logger, "Invalid method: {}. Valid methods are norm and none.", method);
        throw std::runtime_error(errorMessage);
    }

    Composition subsetComposition;
    for (const auto& symbol : symbols) {
        if (m_compositions.count(symbol) == 0) {
            LOG_ERROR(m_logger, "Symbol {} is not in the composition.", symbol);
            throw std::runtime_error("Symbol is not in the composition.");
        } else {
            subsetComposition.registerSymbol(symbol);
        }
        subsetComposition.setMassFraction(symbol, m_compositions.at(symbol).mass_fraction());
    }
    if (method == "norm") {
        bool isNorm = subsetComposition.finalize(true);
        if (!isNorm) {
            LOG_ERROR(m_logger, "Subset composition is invalid.");
            throw std::runtime_error("Subset composition is invalid.");
        }
    }
    return subsetComposition;
}

void Composition::setCompositionMode(bool massFracMode) {
    if (!m_finalized) {
        LOG_ERROR(m_logger, "Composition has not been finalized. Mode cannot be set unless composition is finalized.");
        throw std::runtime_error("Composition has not been finalized (Consider running .finalize()). The mode cannot be set unless the composition is finalized.");
    }

    bool okay = true;
    for (auto& [_, entry] : m_compositions) {
        if (massFracMode) {
            okay = entry.setMassFracMode(m_meanParticleMass);
        } else {
            okay = entry.setNumberFracMode(m_specificNumberDensity);
        }
        if (!okay) {
            LOG_ERROR(m_logger, "Composition mode could not be set.");
            throw std::runtime_error("Composition mode could not be set.");
        }
    }
    m_massFracMode = massFracMode;
}