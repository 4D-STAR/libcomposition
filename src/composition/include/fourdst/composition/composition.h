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

#include <utility>
#include <optional>

#include "fourdst/config/config.h"
#include "fourdst/logging/logging.h"
#include "fourdst/composition/composition_abstract.h"
#include "fourdst/composition/atomicSpecies.h"

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
     * @brief Represents global properties of a finalized composition.
     * @details This struct holds derived quantities that describe the entire composition,
     * such as mean particle mass. It is typically returned by `Composition` methods
     * after the composition has been finalized and is intended for internal or advanced use.
     */
    struct GlobalComposition {
        double specificNumberDensity; ///< The specific number density (moles per unit mass, sum of X_i/M_i), where X_i is mass fraction and M_i is molar mass. Units: mol/g.
        double meanParticleMass; ///< The mean mass per particle (inverse of specific number density). Units: g/mol.

        // Overload the output stream operator for GlobalComposition
        friend std::ostream& operator<<(std::ostream& os, const GlobalComposition& comp);
    };

    /**
     * @brief Represents a single entry (an isotope) within a composition.
     * @details This struct holds the properties of one component, including its symbol,
     * the corresponding `atomic::Species` object, and its abundance (either as a mass
     * fraction or number fraction). It manages the state and conversions for that single entry.
     */
    struct CompositionEntry {
        std::optional<std::string> m_symbol = std::nullopt; ///< The chemical symbol of the species (e.g., "H-1", "Fe-56").
        std::optional<atomic::Species> m_isotope = std::nullopt; ///< The `atomic::Species` object containing detailed isotope data.
        bool m_massFracMode = true; ///< The mode of the composition entry. True if mass fraction, false if number fraction.

        double m_massFraction = 0.0; ///< The mass fraction of the species. Valid only if `m_massFracMode` is true.
        double m_numberFraction = 0.0; ///< The number fraction (mole fraction) of the species. Valid only if `m_massFracMode` is false.
        double m_relAbundance = 0.0; ///< The relative abundance, used internally for conversions. For mass fraction mode, this is X_i / A_i; for number fraction mode, it's n_i * A_i.
        double m_molesPerMass = 0.0;
        double m_cachedNumberFraction = 0.0; ///< Cached number fraction for conversions when in mass fraction mode.


        bool m_initialized = false; ///< True if the composition entry has been initialized with a valid species.

        /**
         * @brief Default constructor. Initializes a default entry (H-1), but in an uninitialized state.
         */
        CompositionEntry();

        /**
         * @brief Constructs a CompositionEntry for a given symbol and abundance mode.
         * @param symbol The chemical symbol of the species (e.g., "He-4").
         * @param massFracMode True to operate in mass fraction mode, false for number fraction mode.
         * @throws exceptions::InvalidSpeciesSymbolError if the symbol does not exist in the atomic species database.
         * @throws exceptions::EntryAlreadyInitializedError if setSpecies is called on an already initialized entry.
         * @par Usage Example:
         * @code
         * CompositionEntry entry("H-1", true); // Entry for H-1 in mass fraction mode.
         * @endcode
         */
        explicit CompositionEntry(const std::string& symbol, bool massFracMode=true);

        /**
         * @brief Copy constructor.
         * @param entry The CompositionEntry to copy.
         */
        CompositionEntry(const CompositionEntry& entry);

        /**
         * @brief Sets the species for the composition entry. This can only be done once.
         * @param symbol The chemical symbol of the species.
         * @throws exceptions::EntryAlreadyInitializedError if the entry has already been initialized.
         * @throws exceptions::InvalidSpeciesSymbolError if the symbol is not found in the atomic species database.
         */
        void setSpecies(const std::string& symbol);

        /**
         * @brief Gets the chemical symbol of the species.
         * @return The chemical symbol.
         */
        [[nodiscard]] std::string symbol() const;

        /**
         * @brief Gets the mass fraction of the species.
         * @pre The entry must be in mass fraction mode.
         * @return The mass fraction of the species.
         * @throws exceptions::CompositionModeError if the entry is in number fraction mode.
         */
        [[nodiscard]] double mass_fraction() const;

        /**
         * @brief Gets the number fraction of the species.
         * @pre The entry must be in number fraction mode.
         * @return The number fraction of the species.
         * @throws exceptions::CompositionModeError if the entry is in mass fraction mode.
         */
        [[nodiscard]] double number_fraction() const;

        /**
         * @brief Gets the number fraction, converting from mass fraction if necessary.
         * @param totalMolesPerMass The total moles per unit mass (specific number density) of the entire composition.
         * @return The number fraction of the species.
         */
        [[nodiscard]] double number_fraction(double totalMolesPerMass) const;

        /**
         * @brief Gets the relative abundance of the species.
         * @return The relative abundance.
         */
        [[nodiscard]] double rel_abundance() const;

        /**
         * @brief Gets the isotope data for the species.
         * @return A const reference to the `atomic::Species` object.
         */
        [[nodiscard]] atomic::Species isotope() const;

        /**
         * @brief Gets the mode of the composition entry.
         * @return True if in mass fraction mode, false if in number fraction mode.
         */
        [[nodiscard]] bool getMassFracMode() const;

        /**
         * @brief Sets the mass fraction of the species.
         * @param mass_fraction The mass fraction to set. Must be in [0, 1].
         * @pre The entry must be in mass fraction mode.
         * @throws exceptions::CompositionModeError if the entry is in number fraction mode.
         */
        void setMassFraction(double mass_fraction);

        /**
         * @brief Sets the number fraction of the species.
         * @param number_fraction The number fraction to set. Must be in [0, 1].
         * @pre The entry must be in number fraction mode.
         * @throws exceptions::CompositionModeError if the entry is in mass fraction mode.
         */
        void setNumberFraction(double number_fraction);

        /**
         * @brief Switches the mode to mass fraction mode.
         * @param meanMolarMass The mean molar mass of the composition, required for conversion.
         * @return True if the mode was successfully set, false otherwise.
         */
        bool setMassFracMode(double meanMolarMass);

        /**
         * @brief Switches the mode to number fraction mode.
         * @param totalMolesPerMass The total moles per unit mass (specific number density) of the composition.
         * @return True if the mode was successfully set, false otherwise.
         */
        bool setNumberFracMode(double totalMolesPerMass);

        /**
         * @brief Overloaded output stream operator for CompositionEntry.
         * @param os The output stream.
         * @param entry The CompositionEntry to output.
         * @return The output stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const CompositionEntry& entry);
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
            std::optional<GlobalComposition> globalComp; ///< Cached global composition data.
            std::optional<CanonicalComposition> canonicalComp; ///< Cached canonical composition data.
            std::optional<std::vector<double>> massFractions; ///< Cached vector of mass fractions.
            std::optional<std::vector<double>> numberFractions; ///< Cached vector of number fractions.
            std::optional<std::vector<double>> molarAbundances; ///< Cached vector of molar abundances.
            std::optional<std::vector<atomic::Species>> sortedSpecies; ///< Cached vector of sorted species (by mass).
            std::optional<std::vector<std::string>> sortedSymbols; ///< Cached vector of sorted species (by mass).
            std::optional<double> Ye; ///< Cached electron abundance.

            void clear() {
                globalComp = std::nullopt;
                canonicalComp = std::nullopt;
                massFractions = std::nullopt;
                numberFractions = std::nullopt;
                molarAbundances = std::nullopt;
                sortedSymbols = std::nullopt;
                sortedSpecies = std::nullopt;
                Ye = std::nullopt;
            }

            [[nodiscard]] bool is_clear() const {
                return !globalComp.has_value() && !canonicalComp.has_value() && !massFractions.has_value() &&
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

        bool m_finalized = false; ///< True if the composition is finalized.
        double m_specificNumberDensity = 0.0; ///< The specific number density of the composition (\sum_{i} X_i m_i. Where X_i is the number fraction of the ith species and m_i is the mass of the ith species).
        double m_meanParticleMass = 0.0; ///< The mean particle mass of the composition (\sum_{i} \frac{n_i}{m_i}. where n_i is the number fraction of the ith species and m_i is the mass of the ith species).
        bool m_massFracMode = true; ///< True if mass fraction mode, false if number fraction mode.

        std::set<std::string> m_registeredSymbols; ///< The registered symbols.
        std::unordered_map<std::string, CompositionEntry> m_compositions; ///< The compositions.

        mutable CompositionCache m_cache; ///< Cache for computed properties to avoid redundant calculations.


        /**
         * @brief Checks if the given symbol is valid by checking against the global species database.
         * @param symbol The symbol to check.
         * @return True if the symbol is valid, false otherwise.
         */
        static bool isValidSymbol(const std::string& symbol);

        /**
         * @brief Checks if the given fractions are valid (sum to ~1.0).
         * @param fractions The fractions to check.
         * @return True if the fractions are valid, false otherwise.
         */
        [[nodiscard]] static bool isValidComposition(const std::vector<double>& fractions) ;

        /**
         * @brief Validates the given fractions, throwing an exception on failure.
         * @param fractions The fractions to validate.
         * @throws exceptions::InvalidCompositionError if the fractions are invalid.
         */
        static void validateComposition(const std::vector<double>& fractions) ;

        /**
         * @brief Finalizes the composition in mass fraction mode.
         * @param norm If true, the composition will be normalized to sum to 1.
         * @return True if the composition is successfully finalized, false otherwise.
         */
        bool finalizeMassFracMode(bool norm);

        /**
         * @brief Finalizes the composition in number fraction mode.
         * @param norm If true, the composition will be normalized to sum to 1.
         * @return True if the composition is successfully finalized, false otherwise.
         */
        bool finalizeNumberFracMode(bool norm);

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
         * @brief Finalizes the composition, making it ready for querying.
         * @details This method checks if the sum of all fractions (mass or number) is approximately 1.0.
         * It also computes global properties like mean particle mass. This **must** be called before
         * any `get...` method can be used.
         * @param norm If true, the composition will be normalized to sum to 1 before validation. [Default: false]
         * @return True if the composition is valid and successfully finalized, false otherwise.
         * @post If successful, `m_finalized` is true and global properties are computed.
         */
        [[nodiscard]] bool finalize(bool norm=false);

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

        /**
         * @brief Constructs and finalizes a Composition with the given symbols and fractions.
         * @details This constructor provides a convenient way to create a fully-formed, finalized composition in one step.
         * The provided fractions must be valid and sum to 1.0.
         * @param symbols The symbols to initialize the composition with.
         * @param fractions The fractions (mass or number) corresponding to the symbols.
         * @param massFracMode True if `fractions` are mass fractions, false if they are number fractions. [Default: true]
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
        Composition(const std::vector<std::string>& symbols, const std::vector<double>& fractions, bool massFracMode=true);

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
         * @param massFracMode True for mass fraction mode, false for number fraction mode. This is only effective for the first symbol registered.
         * @throws exceptions::InvalidSymbolError if the symbol is not in the atomic species database.
         * @throws exceptions::CompositionModeError if attempting to register with a mode that conflicts with the existing mode.
         * @par Usage Example:
         * @code
         * Composition comp;
         * comp.registerSymbol("H-1"); // Now in mass fraction mode
         * comp.registerSymbol("He-4"); // Must also be mass fraction mode
         * @endcode
         */
        void registerSymbol(const std::string& symbol, bool massFracMode=true);

        /**
         * @brief Registers multiple new symbols.
         * @param symbols The symbols to register.
         * @param massFracMode True for mass fraction mode, false for number fraction mode.
         * @throws exceptions::InvalidSymbolError if any symbol is invalid.
         * @throws exceptions::CompositionModeError if the mode conflicts with an already set mode.
         * @par Usage Example:
         * @code
         * std::vector<std::string> symbols = {"H-1", "O-16"};
         * Composition comp;
         * comp.registerSymbol(symbols);
         * @endcode
         */
        void registerSymbol(const std::vector<std::string>& symbols, bool massFracMode=true);

        /**
         * @brief Registers a new species by extracting its symbol.
         * @param species The species to register.
         * @param massFracMode True for mass fraction mode, false for number fraction mode.
         * @throws exceptions::InvalidSymbolError if the species' symbol is invalid.
         * @throws exceptions::CompositionModeError if the mode conflicts.
         * @par Usage Example:
         * @code
         * #include "fourdst/composition/species.h" // Assuming species like H1 are defined here
         * Composition comp;
         * comp.registerSpecies(fourdst::atomic::species.at("H-1"));
         * @endcode
         */
        void registerSpecies(const fourdst::atomic::Species& species, bool massFracMode=true);


        /**
         * @brief Registers a vector of new species.
         * @param species The vector of species to register.
         * @param massFracMode True for mass fraction mode, false for number fraction mode.
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
        void registerSpecies(const std::vector<fourdst::atomic::Species>& species, bool massFracMode=true);


        /**
         * @brief Gets the registered symbols.
         * @return A set of registered symbols.
         */
        [[nodiscard]] std::set<std::string> getRegisteredSymbols() const override;

        /**
         * @brief Get a set of all species that are registered in the composition.
         * @return A set of `atomic::Species` objects registered in the composition.
         */
        [[nodiscard]] std::set<fourdst::atomic::Species> getRegisteredSpecies() const override;

        /**
         * @brief Sets the mass fraction for a given symbol.
         * @param symbol The symbol to set the mass fraction for.
         * @param mass_fraction The mass fraction to set (must be in [0, 1]).
         * @return The previous mass fraction that was set for the symbol.
         * @throws exceptions::UnregisteredSymbolError if the symbol is not registered.
         * @throws exceptions::CompositionModeError if the composition is in number fraction mode.
         * @throws exceptions::InvalidCompositionError if the mass fraction is not between 0 and 1.
         * @post The composition is marked as not finalized.
         * @par Usage Example:
         * @code
         * Composition comp;
         * comp.registerSymbol("H-1");
         * comp.setMassFraction("H-1", 0.7);
         * @endcode
         */
        double setMassFraction(const std::string& symbol, const double& mass_fraction);

        /**
         * @brief Sets the mass fraction for multiple symbols.
         * @param symbols The symbols to set the mass fractions for.
         * @param mass_fractions The mass fractions corresponding to the symbols.
         * @return A vector of the previous mass fractions that were set.
         * @throws exceptions::InvalidCompositionError if symbol and fraction counts differ.
         * @throws See `setMassFraction(const std::string&, const double&)` for other exceptions.
         * @post The composition is marked as not finalized.
         */
        std::vector<double> setMassFraction(const std::vector<std::string>& symbols, const std::vector<double>& mass_fractions);

        /**
         * @brief Sets the mass fraction for a given species.
         * @param species The species to set the mass fraction for.
         * @param mass_fraction The mass fraction to set.
         * @return The previous mass fraction that was set for the species.
         * @throws exceptions::UnregisteredSymbolError if the species is not registered.
         * @throws exceptions::CompositionModeError if the composition is in number fraction mode.
         * @throws exceptions::InvalidCompositionError if the mass fraction is not between 0 and 1.
         */
        double setMassFraction(const fourdst::atomic::Species& species, const double& mass_fraction);

        /**
         * @brief Sets the mass fraction for multiple species.
         * @param species The vector of species to set the mass fractions for.
         * @param mass_fractions The vector of mass fractions corresponding to the species.
         * @return A vector of the previous mass fractions that were set.
         * @throws See `setMassFraction(const std::vector<std::string>&, const std::vector<double>&)` for exceptions.
         */
        std::vector<double> setMassFraction(const std::vector<fourdst::atomic::Species>& species, const std::vector<double>& mass_fractions);

        /**
         * @brief Sets the number fraction for a given symbol.
         * @param symbol The symbol to set the number fraction for.
         * @param number_fraction The number fraction to set (must be in [0, 1]).
         * @return The previous number fraction that was set.
         * @throws exceptions::UnregisteredSymbolError if the symbol is not registered.
         * @throws exceptions::CompositionModeError if the composition is in mass fraction mode.
         * @throws exceptions::InvalidCompositionError if the number fraction is not between 0 and 1.
         * @post The composition is marked as not finalized.
         */
        double setNumberFraction(const std::string& symbol, const double& number_fraction);

        /**
         * @brief Sets the number fraction for multiple symbols.
         * @param symbols The symbols to set the number fractions for.
         * @param number_fractions The number fractions corresponding to the symbols.
         * @return A vector of the previous number fractions that were set.
         * @throws exceptions::InvalidCompositionError if symbol and fraction counts differ.
         * @throws See `setNumberFraction(const std::string&, const double&)` for other exceptions.
         */
        std::vector<double> setNumberFraction(const std::vector<std::string>& symbols, const std::vector<double>& number_fractions);

        /**
         * @brief Sets the number fraction for a given species.
         * @param species The species to set the number fraction for.
         * @param number_fraction The number fraction to set for the species.
         * @return The previous number fraction that was set for the species.
         * @throws exceptions::UnregisteredSymbolError if the species is not registered.
         * @throws exceptions::CompositionModeError if the composition is in mass fraction mode.
         * @throws exceptions::InvalidCompositionError if the number fraction is not between 0 and 1.
         */
        double setNumberFraction(const fourdst::atomic::Species& species, const double& number_fraction);

        /**
         * @brief Sets the number fraction for multiple species.
         * @param species The vector of species to set the number fractions for.
         * @param number_fractions The vector of number fractions corresponding to the species.
         * @return The vector of the previous number fractions that were set.
         * @throws See `setNumberFraction(const std::vector<std::string>&, const std::vector<double>&)` for exceptions.
         */
        std::vector<double> setNumberFraction(const std::vector<fourdst::atomic::Species>& species, const std::vector<double>& number_fractions);

        /**
         * @brief Mixes this composition with another to produce a new composition.
         * @details The mixing is performed linearly on the mass fractions. The formula for each species is:
         * `new_X_i = fraction * this_X_i + (1 - fraction) * other_X_i`.
         * The resulting composition is automatically finalized.
         * @param other The other composition to mix with.
         * @param fraction The mixing fraction. A value of 1.0 means the new composition is 100% `this`, 0.0 means 100% `other`.
         * @return A new, finalized `Composition` object representing the mixture.
         * @pre Both `this` and `other` compositions must be finalized.
         * @throws exceptions::CompositionNotFinalizedError if either composition is not finalized.
         * @throws exceptions::InvalidCompositionError if the fraction is not between 0 and 1.
         */
        [[nodiscard]] Composition mix(const Composition& other, double fraction) const;

        /**
         * @brief Gets the mass fractions of all species in the composition.
         * @pre The composition must be finalized.
         * @return An unordered map of symbols to their mass fractions.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] std::unordered_map<std::string, double> getMassFraction() const override;

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
        [[nodiscard]] std::unordered_map<std::string, double> getNumberFraction() const override;

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
        [[nodiscard]] double getMolarAbundance(const fourdst::atomic::Species& species) const override;

        /**
         * @brief Gets the composition entry and global composition data for a given symbol.
         * @pre The composition must be finalized.
         * @param symbol The symbol to get the composition for.
         * @return A pair containing the CompositionEntry and GlobalComposition for the given symbol.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         * @throws exceptions::UnregisteredSymbolError if the symbol is not in the composition.
         */
        [[nodiscard]] std::pair<CompositionEntry, GlobalComposition> getComposition(const std::string& symbol) const;

        /**
         * @brief Gets the composition entry and global composition data for a given species.
         * @pre The composition must be finalized.
         * @param species The species to get the composition for.
         * @return A pair containing the CompositionEntry and GlobalComposition for the given species.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         * @throws exceptions::UnregisteredSymbolError if the species is not in the composition.
         */
        [[nodiscard]] std::pair<CompositionEntry, GlobalComposition> getComposition(const fourdst::atomic::Species& species) const;

        /**
         * @brief Gets all composition entries and the global composition data.
         * @pre The composition must be finalized.
         * @return A pair containing an unordered map of all CompositionEntries and the GlobalComposition.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] std::pair<std::unordered_map<std::string, CompositionEntry>, GlobalComposition> getComposition() const;

        /**
         * @brief Compute the mean particle mass of the composition.
         * @pre The composition must be finalized.
         * @return Mean particle mass in atomic mass units (g/mol).
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] double getMeanParticleMass() const override;

        /**
         * @brief Compute the mean atomic number of the composition.
         * @pre The composition must be finalized.
         * @return Mean atomic number <Z>.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] double getMeanAtomicNumber() const override;

        /**
         * @brief Compute the electron abundance of the composition.
         * @details Ye is defined as the sum over all species of (Z_i * X_i / A_i), where Z_i is the atomic number, X_i is the mass fraction, and A_i is the atomic mass of species i.
         * @return Ye (electron abundance).
         * @pre The composition must be finalized.
         */
        [[nodiscard]] double getElectronAbundance() const override;

        /**
         * @brief Creates a new Composition object containing a subset of species from this one.
         * @param symbols The symbols to include in the subset.
         * @param method The method for handling the abundances of the new subset. Can be "norm" (normalize abundances to sum to 1) or "none" (keep original abundances).
         * @return A new `Composition` object containing the subset.
         * @throws exceptions::UnregisteredSymbolError if any requested symbol is not in the original composition.
         * @throws exceptions::InvalidMixingMode if an invalid method is provided.
         * @throws exceptions::FailedToFinalizeCompositionError if normalization fails.
         */
        [[nodiscard]] Composition subset(const std::vector<std::string>& symbols, const std::string& method="norm") const;

        /**
         * @brief Checks if a symbol is registered in the composition.
         * @param symbol The symbol to check.
         * @return True if the symbol is registered, false otherwise.
         */
        [[nodiscard]] bool hasSymbol(const std::string& symbol) const override;

        /**
         * @brief Checks if a species is registered in the composition.
         * @param species The species to check.
         * @return True if the species is registered, false otherwise.
         */
        [[nodiscard]] bool hasSpecies(const fourdst::atomic::Species &species) const override;

        /**
         * @brief Checks if a given isotope is present in the composition.
         * @pre The composition must be finalized.
         * @param isotope The isotope to check for.
         * @return True if the isotope is in the composition, false otherwise.
         * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
         */
        [[nodiscard]] bool contains(const atomic::Species& isotope) const override;

        /**
        * @brief Sets the composition mode (mass fraction vs. number fraction).
        * @details This function converts all entries in the composition to the specified mode.
        * @pre The composition must be finalized before the mode can be switched.
        * @param massFracMode True to switch to mass fraction mode, false for number fraction mode.
        * @throws exceptions::CompositionNotFinalizedError if the composition is not finalized.
        * @throws std::runtime_error if the conversion fails for an unknown reason.
        */
        void setCompositionMode(bool massFracMode);

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
        * @brief Overloads the + operator to mix two compositions with a 50/50 fraction.
        * @details This is a convenience operator that calls `mix(other, 0.5)`.
        * @param other The other composition to mix with.
        * @return The new, mixed composition.
        * @pre Both compositions must be finalized.
        * @throws See `mix()` for exceptions.
        */
        Composition operator+(const Composition& other) const;

        /**
         * @brief Returns an iterator to the beginning of the composition map.
         * @return An iterator to the beginning.
         */
        auto begin() {
            return m_compositions.begin();
        }

        /**
         * @brief Returns a const iterator to the beginning of the composition map.
         * @return A const iterator to the beginning.
         */
        [[nodiscard]] auto begin() const {
            return m_compositions.cbegin();
        }

        /**
         * @brief Returns an iterator to the end of the composition map.
         * @return An iterator to the end.
         */
        auto end() {
            return m_compositions.end();
        }

        /**
         * @brief Returns a const iterator to the end of the composition map.
         * @return A const iterator to the end.
         */
        [[nodiscard]] auto end() const {
            return m_compositions.cend();
        }

    };
}; // namespace fourdst::composition
