/* ***********************************************************************
//
//   Copyright (C) 2025 -- The 4D-STAR Collaboration
//   File Author: Emily Boudreaux
//   Last Modified: March 26, 2025
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
#pragma once

#include <string>
#include <unordered_map>
#include <set>

#include <optional>

#include "fourdst/config/config.h"
#include "fourdst/logging/logging.h"
#include "fourdst/composition/composition_abstract.h"
#include "fourdst/atomic/atomicSpecies.h"

namespace fourdst::composition {
    /**
     * @struct CanonicalComposition
     * @brief Represents the canonical (X, Y, Z) composition of stellar material.
     * @details This is a standard astrophysical representation where:
     * - X is the total mass fraction of all hydrogen isotopes.
     * - Y is the total mass fraction of all helium isotopes.
     * - Z is the total mass fraction of all other elements (metals).
     * By definition, X + Y + Z should sum to 1.0.
     */
    struct CanonicalComposition {
        double X = 0.0; ///< Mass fraction of Hydrogen.
        double Y = 0.0; ///< Mass fraction of Helium.
        double Z = 0.0; ///< Mass fraction of Metals.

        /**
         * @brief Overloads the stream insertion operator for easy printing.
         * @param os The output stream.
         * @param composition The CanonicalComposition object to print.
         * @return The output stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const CanonicalComposition& composition) {
            os << "<CanonicalComposition: "
               << "X = " << composition.X << ", "
               << "Y = " << composition.Y << ", "
               << "Z = " << composition.Z << ">";
            return os;
        }
    };

    /**
     * @class Composition
     * @brief Manages a collection of chemical species and their abundances.
     * @details This class is a primary interface for defining and manipulating material compositions.
     * It can operate in two modes: mass fraction or number fraction.
     *
     * **Key Rules and Workflow:**
     * 1.  **Registration:** Before setting an abundance for a species, its symbol (e.g., "He-4") must be registered using `registerSymbol()` or `registerSpecies()`. All registered species must conform to the same abundance mode (mass or number fraction).
     * 2.  **Setting Abundances:** Use `setMassFraction()` or `setNumberFraction()` to define the composition.
     * 3.  **Finalization:** Before querying any compositional data (e.g., `getMassFraction()`, `getMeanParticleMass()`), the composition must be **finalized** by calling `finalize()`. This step validates the composition (abundances sum to ~1.0) and computes global properties.
     * 4.  **Modification:** Any modification to abundances after finalization will un-finalize the composition, requiring another call to `finalize()` before data can be retrieved again.
     * 5.  **Construction:** A pre-finalized composition can be created by providing symbols and valid, normalized abundances to the constructor.
     *
     * @throws This class throws various exceptions from `fourdst::composition::exceptions` for invalid operations, such as using unregistered symbols, providing invalid abundances, or accessing data from a non-finalized composition.
     *
     * @par Mass Fraction Example:
     * @code
     * Composition comp;
     * comp.registerSymbol("H-1");
     * comp.registerSymbol("He-4");
     * comp.setMassFraction("H-1", 0.75);
     * comp.setMassFraction("He-4", 0.25);
     * if (comp.finalize()) {
     *     double he_mass_frac = comp.getMassFraction("He-4"); // Returns 0.25
     * }
     * @endcode
     *
     * @par Number Fraction Example:
     * @code
     * Composition comp;
     * comp.registerSymbol("H-1", false); // Register in number fraction mode
     * comp.registerSymbol("He-4", false);
     * comp.setNumberFraction("H-1", 0.9);
     * comp.setNumberFraction("He-4", 0.1);
     * if (comp.finalize()) {
     *     double he_num_frac = comp.getNumberFraction("He-4"); // Returns 0.1
     * }
     * @endcode
     */
    class Composition : public CompositionAbstract {
    private:
        struct CompositionCache {
            std::optional<CanonicalComposition> canonicalComp; ///< Cached canonical composition data.
            std::optional<std::vector<double>> massFractions; ///< Cached vector of mass fractions.
            std::optional<std::vector<double>> numberFractions; ///< Cached vector of number fractions.
            std::optional<std::vector<double>> molarAbundances; ///< Cached vector of molar abundances.
            std::optional<std::vector<atomic::Species>> sortedSpecies; ///< Cached vector of sorted species (by mass).
            std::optional<std::vector<std::string>> sortedSymbols; ///< Cached vector of sorted species (by mass).
            std::optional<double> Ye; ///< Cached electron abundance.

            void clear() {
                canonicalComp = std::nullopt;
                massFractions = std::nullopt;
                numberFractions = std::nullopt;
                molarAbundances = std::nullopt;
                sortedSymbols = std::nullopt;
                sortedSpecies = std::nullopt;
                Ye = std::nullopt;
            }

            [[nodiscard]] bool is_clear() const {
                return !canonicalComp.has_value() && !massFractions.has_value() &&
                       !numberFractions.has_value() && !molarAbundances.has_value() && !sortedSymbols.has_value() &&
                       !Ye.has_value() && !sortedSpecies.has_value();
            }
        };
    private:
        // logging::LogManager& m_logManager = logging::LogManager::getInstance();
        static quill::Logger* getLogger() {
            static quill::Logger* logger = logging::LogManager::getInstance().getLogger("log");
            return logger;
        }

        std::set<atomic::Species> m_registeredSpecies;
        std::map<atomic::Species, double> m_molarAbundances;

        mutable CompositionCache m_cache; ///< Cache for computed properties to avoid redundant calculations.

    public:
        /**
         * @brief Default constructor.
         */
        Composition() = default;

        /**
         * @brief Default destructor.
         */
        ~Composition() override = default;

        /**
         * @brief Constructs a Composition and registers the given symbols.
         * @param symbols The symbols to register. The composition will be in mass fraction mode by default.
         * @throws exceptions::InvalidSymbolError if any symbol is invalid.
         * @par Usage Example:
         * @code
         * std::vector<std::string> symbols = {"H-1", "O-16"};
         * Composition comp(symbols);
         * comp.setMassFraction("H-1", 0.11);
         * comp.setMassFraction("O-16", 0.89);
         * comp.finalize();
         * @endcode
         */
        explicit Composition(const std::vector<std::string>& symbols);

        explicit Composition(const std::vector<atomic::Species>& species);

        /**
         * @brief Constructs a Composition and registers the given symbols from a set.
         * @param symbols The symbols to register. The composition will be in mass fraction mode by default.
         * @throws exceptions::InvalidSymbolError if any symbol is invalid.
         * @par Usage Example:
         * @code
         * std::set<std::string> symbols = {"H-1", "O-16"};
         * Composition comp(symbols);
         * @endcode
         */
        explicit Composition(const std::set<std::string>& symbols);

        explicit Composition(const std::set<atomic::Species>& species);

        /**
         * @brief Constructs and finalizes a Composition with the given symbols and fractions.
         * @details This constructor provides a convenient way to create a fully-formed, finalized composition in one step.
         * The provided fractions must be valid and sum to 1.0.
         * @param symbols The symbols to initialize the composition with.
         * @param molarAbundances The corresponding molar abundances for each symbol.
         * @throws exceptions::InvalidCompositionError if the number of symbols and fractions do not match, or if the fractions do not sum to ~1.0.
         * @throws exceptions::InvalidSymbolError if any symbol is invalid.
         * @post The composition is immediately finalized.
         * @par Usage Example:
         * @code
         * std::vector<std::string> symbols = {"H-1", "O-16"};
         * std::vector<double> mass_fractions = {0.1119, 0.8881};
         * Composition comp(symbols, mass_fractions); // Finalized on construction
         * @endcode
         */
        Composition(const std::vector<std::string>& symbols, const std::vector<double>& molarAbundances);

        Composition(const std::vector<atomic::Species>& species, const std::vector<double>& molarAbundances);

        Composition(const std::set<std::string>& symbols, const std::vector<double>& molarAbundances);

        /**
         * @brief Constructs a Composition from another Composition.
         * @param composition The Composition to copy.
         */
        Composition(const Composition& composition);

        /**
         * @brief Assignment operator.
         * @param other The Composition to assign from.
         * @return A reference to this Composition.
         */
        Composition& operator=(Composition const& other);

        /**
         * @brief Registers a new symbol for inclusion in the composition.
         * @details A symbol must be registered before its abundance can be set. The first registration sets the mode (mass/number fraction) for the entire composition.
         * @param symbol The symbol to register (e.g., "Fe-56").
         * @throws exceptions::InvalidSymbolError if the symbol is not in the atomic species database.
         * @throws exceptions::CompositionModeError if attempting to register with a mode that conflicts with the existing mode.
         * @par Usage Example:
         * @code
         * Composition comp;
         * comp.registerSymbol("H-1"); // Now in mass fraction mode
         * comp.registerSymbol("He-4"); // Must also be mass fraction mode
         * @endcode
         */
        void registerSymbol(const std::string& symbol);

        /**
         * @brief Registers multiple new symbols.
         * @param symbols The symbols to register.
         * @throws exceptions::InvalidSymbolError if any symbol is invalid.
         * @throws exceptions::CompositionModeError if the mode conflicts with an already set mode.
         * @par Usage Example:
         * @code
         * std::vector<std::string> symbols = {"H-1", "O-16"};
         * Composition comp;
         * comp.registerSymbol(symbols);
         * @endcode
         */
        void registerSymbol(const std::vector<std::string>& symbols);

        /**
         * @brief Registers a new species by extracting its symbol.
         * @param species The species to register.
         * @throws exceptions::InvalidSymbolError if the species' symbol is invalid.
         * @throws exceptions::CompositionModeError if the mode conflicts.
         * @par Usage Example:
         * @code
         * #include "fourdst/composition/species.h" // Assuming species like H1 are defined here
         * Composition comp;
         * comp.registerSpecies(fourdst::atomic::species.at("H-1"));
         * @endcode
         */
        void registerSpecies(const atomic::Species& species);


        /**
         * @brief Registers a vector of new species.
         * @param species The vector of species to register.
         * @throws exceptions::InvalidSymbolError if any species' symbol is invalid.
         * @throws exceptions::CompositionModeError if the mode conflicts.
         * @par Usage Example:
         * @code
         * #include "fourdst/composition/species.h"
         * Composition comp;
         * std::vector<fourdst::atomic::Species> my_species = { ... };
         * comp.registerSpecies(my_species, false); // Number fraction mode
         * @endcode
         */
        void registerSpecies(const std::vector<atomic::Species>& species);

        /**
         * @brief Checks if a given isotope is present in the composition.
         * @pre The composition must be finalized.
         * @param species The isotope to check for.
         * @return True if the isotope is in the composition, false otherwise.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] bool contains(const atomic::Species& species) const override;

        [[nodiscard]] bool contains(const std::string& symbol) const override;

        [[nodiscard]] size_t size() const override;

        void setMolarAbundance(
            const std::string& symbol,
            const double& molar_abundance
        );

        void setMolarAbundance(
            const atomic::Species& species,
            const double& molar_abundance
        );

        void setMolarAbundance(
            const std::vector<std::string>& symbols,
            const std::vector<double>& molar_abundances
        );

        void setMolarAbundance(
            const std::vector<atomic::Species>& species,
            const std::vector<double>& molar_abundances
        );

        void setMolarAbundance(
            const std::set<std::string>& symbols,
            const std::vector<double>& molar_abundances
        );

        void setMolarAbundance(
            const std::set<atomic::Species>& species,
            const std::vector<double>& molar_abundances
        );

        /**
         * @brief Gets the registered symbols.
         * @return A set of registered symbols.
         */
        [[nodiscard]] std::set<std::string> getRegisteredSymbols() const override;

        /**
         * @brief Get a set of all species that are registered in the composition.
         * @return A set of `atomic::Species` objects registered in the composition.
         */
        [[nodiscard]] const std::set<atomic::Species> &getRegisteredSpecies() const override;

        /**
         * @brief Gets the mass fractions of all species in the composition.
         * @pre The composition must be finalized.
         * @return An unordered map of symbols to their mass fractions.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] std::unordered_map<atomic::Species, double> getMassFraction() const override;

        /**
         * @brief Gets the mass fraction for a given symbol.
         * @pre The composition must be finalized.
         * @param symbol The symbol to get the mass fraction for.
         * @return The mass fraction for the given symbol.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         * @throws exceptions::UnregisteredSymbolError if the symbol is not in the composition.
         */
        [[nodiscard]] double getMassFraction(const std::string& symbol) const override;

        /**
         * @brief Gets the mass fraction for a given isotope.
         * @pre The composition must be finalized.
         * @param species The isotope to get the mass fraction for.
         * @return The mass fraction for the given isotope.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         * @throws exceptions::UnregisteredSymbolError if the isotope is not registered in the composition.
         */
        [[nodiscard]] double getMassFraction(const fourdst::atomic::Species& species) const override;

        /**
         * @brief Gets the number fraction for a given symbol.
         * @pre The composition must be finalized.
         * @param symbol The symbol to get the number fraction for.
         * @return The number fraction for the given symbol.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         * @throws exceptions::UnregisteredSymbolError if the symbol is not in the composition.
         */
        [[nodiscard]] double getNumberFraction(const std::string& symbol) const override;

        /**
         * @brief Gets the number fraction for a given isotope.
         * @pre The composition must be finalized.
         * @param species The isotope to get the number fraction for.
         * @return The number fraction for the given isotope.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         * @throws exceptions::UnregisteredSymbolError if the isotope is not registered in the composition.
         */
        [[nodiscard]] double getNumberFraction(const fourdst::atomic::Species& species) const override;

        /**
         * @brief Gets the number fractions of all species in the composition.
         * @pre The composition must be finalized.
         * @return An unordered map of symbols to their number fractions.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] std::unordered_map<atomic::Species, double> getNumberFraction() const override;

        /**
        * @brief Gets the molar abundance (X_i / A_i) for a given symbol.
        * @pre The composition must be finalized.
        * @param symbol The symbol to get the molar abundance for.
        * @return The molar abundance for the given symbol.
        * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
        * @throws exceptions::UnregisteredSymbolError if the symbol is not in the composition.
        */
        [[nodiscard]] double getMolarAbundance(const std::string& symbol) const override;

        /**
         * @brief Gets the molar abundance for a given isotope.
         * @pre The composition must be finalized.
         * @param species The isotope to get the molar abundance for.
         * @return The molar abundance for the given isotope.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         * @throws exceptions::UnregisteredSymbolError if the isotope is not registered in the composition.
         */
        [[nodiscard]] double getMolarAbundance(const atomic::Species& species) const override;

        /**
         * @brief Compute the mean particle mass of the composition.
         * @pre The composition must be finalized.
         * @return Mean particle mass in atomic mass units (g/mol).
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] double getMeanParticleMass() const override;

        /**
         * @brief Compute the electron abundance of the composition.
         * @details Ye is defined as the sum over all species of (Z_i * X_i / A_i), where Z_i is the atomic number, X_i is the mass fraction, and A_i is the atomic mass of species i.
         * @return Ye (electron abundance).
         * @pre The composition must be finalized.
         */
        [[nodiscard]] double getElectronAbundance() const override;


        /**
         * @brief Gets the current canonical composition (X, Y, Z).
         * @details Calculates the total mass fractions for H, He, and metals.
         * @pre The composition must be finalized.
         * @param harsh If true, this will throw an error if `1 - (X + Y)` is not equal to the directly summed `Z` (within a tolerance). If false, it will only log a warning.
         * @return The `CanonicalComposition` struct.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         * @throws std::runtime_error if `harsh` is true and the canonical composition is not self-consistent.
         */
        [[nodiscard]] CanonicalComposition getCanonicalComposition(bool harsh=false) const;

        /**
         * @brief Get a uniform vector representation of the mass fraction stored in the composition object sorted such that the lightest species is at index 0 and the heaviest is at the last index.
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @return the vector of mass fractions sorted by species mass (lightest to heaviest).
         * @pre The composition must be finalized.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] std::vector<double> getMassFractionVector() const override;

        /**
         * @brief Get a uniform vector representation of the number fractions stored in the composition object sorted such that the lightest species is at index 0 and the heaviest is at the last index.
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @return the vector of number fractions sorted by species mass (lightest to heaviest).
         * @pre The composition must be finalized.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] std::vector<double> getNumberFractionVector() const override;

        /**
         * @brief Get a uniform vector representation of the molar abundances stored in the composition object sorted such that the lightest species is at index 0 and the heaviest is at the last index.
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @return the vector of molar abundances sorted by species mass (lightest to heaviest).
         * @pre The composition must be finalized.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] std::vector<double> getMolarAbundanceVector() const override;

        /**
         * @brief get the index in the sorted vector representation for a given symbol
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @pre The composition must be finalized.
         * @pre symbol must be registered in the composition
         * @param symbol the symbol to look up the index for. Note that this is the index species data will be at if you were to call getMolarAbundanceVector(), getMassFractionVector(), or getNumberFractionVector()
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         * @throws exceptions::UnregisteredSymbolError if the symbol is not registered in the composition
         * @return The index of the symbol in the sorted vector representation.
         */
        [[nodiscard]] size_t getSpeciesIndex(const std::string& symbol) const override;

        /**
         * @brief get the index in the sorted vector representation for a given symbol
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @pre The composition must be finalized.
         * @pre symbol must be registered in the composition
         * @param species the species to look up the index for. Note that this is the index species data will be at if you were to call getMolarAbundanceVector(), getMassFractionVector(), or getNumberFractionVector()
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         * @throws exceptions::UnregisteredSymbolError if the symbol is not registered in the composition
         * @return The index of the symbol in the sorted vector representation.
         */
        [[nodiscard]] size_t getSpeciesIndex(const atomic::Species& species) const override;

        /**
         * @brief Get the species at a given index in the sorted vector representation.
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @pre The composition must be finalized.
         * @param index The index in the sorted vector representation for which to return the species. Must be in [0, N-1] where N is the number of registered species.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         * @throws std::out_of_range if the index is out of range.
         * @return The species at the given index in the sorted vector representation.
         */
        [[nodiscard]] atomic::Species getSpeciesAtIndex(size_t index) const override;

        /**
         * @brief Overloaded output stream operator for Composition.
         * @param os The output stream.
         * @param composition The Composition to output.
         * @return The output stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const Composition& composition);

        /**
         * @brief Returns an iterator to the beginning of the composition map.
         * @return An iterator to the beginning.
         */
        auto begin() {
            return m_molarAbundances.begin();
        }

        /**
         * @brief Returns a const iterator to the beginning of the composition map.
         * @return A const iterator to the beginning.
         */
        [[nodiscard]] auto begin() const {
            return m_molarAbundances.cbegin();
        }

        /**
         * @brief Returns an iterator to the end of the composition map.
         * @return An iterator to the end.
         */
        auto end() {
            return m_molarAbundances.end();
        }

        /**
         * @brief Returns a const iterator to the end of the composition map.
         * @return A const iterator to the end.
         */
        [[nodiscard]] auto end() const {
            return m_molarAbundances.cend();
        }

    };
}; // namespace fourdst::composition
