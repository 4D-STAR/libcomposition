#ifndef COMPOSITION_H
#define COMPOSITION_H

#include <iomanip>
#include <iostream>
#include <string>
#include <unordered_map>
#include <set>

#include "quill/LogMacros.h"

#include "probe.h"
#include "config.h"

#include "atomicSpecies.h"

namespace composition{
    /**
     * @brief Represents an entry in the composition with a symbol and mass fraction.
     */
    struct CompositionEntry {
        std::string symbol; ///< The chemical symbol of the element.
        double mass_fraction; ///< The mass fraction of the element.

        /**
         * @brief Overloaded output stream operator for CompositionEntry.
         * @param os The output stream.
         * @param entry The CompositionEntry to output.
         * @return The output stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const CompositionEntry& entry) {
            os << "<" << entry.symbol << " : " << entry.mass_fraction << ">";
            return os;
        }
    };

    /**
     * @brief Manages the composition of elements.
     * @details The composition is a collection of elements with their respective mass fractions.
     * The general purpose of this class is to provide a standardized interface for managing the composition of
     * any part of 4DSSE. There are a few rules when using this class.
     * - Only species in the atomicSpecies.h database can be used. There are 1000s (All species from AME2020) in there so it should not be a problem.
     * - Before a mass fraction can be set with a particular instance of Composition, the symbol must be registered. (i.e. register He-3 before setting its mass fraction)
     * - Before any composition information can be retrived (e.g. getComposition), the composition must be finalized (call to .finalize()). This checks if the total mass fraction sums to approximatly 1 (within 1 part in 10^8)
     * - Any changes made to the composition after finalization will "unfinalize" the composition. This means that the composition must be finalized again before any information can be retrived.
     * - The mass fraction of any individual species must be no more than 1 and no less than 0.
     * - The only exception to the finalize rule is if the compositon was constructed with symbols and mass fractions at instantiation time. In this case, the composition is automatically finalized.  
     * however, this means that the composition passed to the constructor must be valid.
     * 
     * @example Constructing a finalized composition with symbols and mass fractions:
     * @code
     * std::vector<std::string> symbols = {"H", "He"};
     * std::vector<double> mass_fractions = {0.7, 0.3};
     * Composition comp(symbols, mass_fractions);
     * @endcode
     * @example Constructing a composition with symbols and finalizing it later:
     * @code
     * std::vector<std::string> symbols = {"H", "He"};
     * Composition comp(symbols);
     * comp.setComposition("H", 0.7);
     * comp.setComposition("He", 0.3);
     * comp.finalize();
     * @endcode
     */
    class Composition {
    private:
        Config& m_config = Config::getInstance();
        Probe::LogManager& m_logManager = Probe::LogManager::getInstance();
        quill::Logger* m_logger = m_logManager.getLogger("log");

        bool m_finalized = false;

        std::set<std::string> m_registeredSymbols;
        std::unordered_map<std::string, CompositionEntry> m_compositions;

        /**
         * @brief Checks if the given symbol is valid. 
         * @details A symbol is valid if it is in the atomic species database (species in atomicSpecies.h). These include all the isotopes from AME2020.
         * @param symbol The symbol to check.
         * @return True if the symbol is valid, false otherwise.
         */
        bool isValidSymbol(const std::string& symbol) const;

        /**
         * @brief Checks if the given mass fractions are valid.
         * @param mass_fractions The mass fractions to check.
         * @return True if the mass fractions are valid, false otherwise.
         */
        bool isValidComposition(const std::vector<double>& mass_fractions) const;

        /**
         * @brief Validates the given mass fractions.
         * @param mass_fractions The mass fractions to validate.
         * @throws std::invalid_argument if the mass fractions are invalid.
         */
        void validateComposition(const std::vector<double>& mass_fractions) const;

    public:
        /**
         * @brief Default constructor.
         */
        Composition() = default;

        /**
         * @brief Default destructor.
         */
        ~Composition() = default;

        /**
         * @brief Finalizes the composition.
         * @return True if the composition is successfully finalized, false otherwise.
         */
        bool finalize();

        /**
         * @brief Constructs a Composition with the given symbols.
         * @param symbols The symbols to initialize the composition with.
         * @example
         * @code
         * std::vector<std::string> symbols = {"H", "O"};
         * Composition comp(symbols);
         * @endcode
         */
        Composition(const std::vector<std::string>& symbols);

        /**
         * @brief Constructs a Composition with the given symbols and mass fractions.
         * @param symbols The symbols to initialize the composition with.
         * @param mass_fractions The mass fractions corresponding to the symbols.
         * @example
         * @code
         * std::vector<std::string> symbols = {"H", "O"};
         * std::vector<double> mass_fractions = {0.1, 0.9};
         * Composition comp(symbols, mass_fractions);
         * @endcode
         */
        Composition(const std::vector<std::string>& symbols, const std::vector<double>& mass_fractions);

        /**
         * @brief Registers a new symbol.
         * @param symbol The symbol to register.
         * @example
         * @code
         * Composition comp;
         * comp.registerSymbol("H");
         * @endcode
         */
        void registerSymbol(const std::string& symbol);

        /**
         * @brief Gets the registered symbols.
         * @return A set of registered symbols.
         */
        std::set<std::string> getRegisteredSymbols() const;

        /**
         * @brief Sets the composition for a given symbol.
         * @param symbol The symbol to set the composition for.
         * @param mass_fraction The mass fraction to set.
         * @return The mass fraction that was set.
         * @example
         * @code
         * Composition comp;
         * comp.setComposition("H", 0.1);
         * @endcode
         */
        double setComposition(const std::string& symbol, const double& mass_fraction);

        /**
         * @brief Sets the composition for multiple symbols.
         * @param symbols The symbols to set the composition for.
         * @param mass_fractions The mass fractions corresponding to the symbols.
         * @return A vector of mass fractions that were set.
         * @example
         * @code
         * std::vector<std::string> symbols = {"H", "O"};
         * std::vector<double> mass_fractions = {0.1, 0.9};
         * Composition comp;
         * comp.setComposition(symbols, mass_fractions);
         * @endcode
         */
        std::vector<double> setComposition(const std::vector<std::string>& symbols, const std::vector<double>& mass_fractions);

        /**
         * @brief Gets the compositions.
         * @return An unordered map of compositions.
         */
        std::unordered_map<std::string, CompositionEntry> getCompositions() const;

        /**
         * @brief Gets the composition for a given symbol.
         * @param symbol The symbol to get the composition for.
         * @return The CompositionEntry for the given symbol.
         */
        CompositionEntry getComposition(const std::string& symbol) const;

        /**
         * @brief Overloaded output stream operator for Composition.
         * @param os The output stream.
         * @param composition The Composition to output.
         * @return The output stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const Composition& composition) {
            os << "Composition: \n";
            for (const auto& [symbol, entry] : composition.m_compositions) {
                os << entry << "\n";
            }
            return os;
        }

    };
};

#endif // COMPOSITION_H