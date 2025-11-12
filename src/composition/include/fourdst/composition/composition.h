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
     * In order to use the Composition class a user must first register symbols or species. Symbols are
     * the string representation of a species (i.e. deuterium would be "H-2" whereas Beryllium 7 would
     * be "Be-7") and then set the molar abundances. Species are the data structure
     * fourdst::atomic::Species version. Here Deuterium would be represented by the Species
     * fourdst::atomic::H_2 whereas Beryllium 7 would be fourdst::atomic::Be_7. Once the symbols/species have been registered
     * the user can then set molar abundances.
     *
     * Once the Composition object has been populated the user can query mass fractions, number fractions, electron
     * abundances, mean particle mass, molar abundance, and Canonical (X, Y, Z) composition.
     *
     * @note This class only accepts molar abundances as input. If you
     * wish to construct a Composition using a vector of mass fractions, those must first be converted
     * to molar abundances. There is a helper function `fourdst::composition::buildCompositionFromMassFractions` which
     * wll facilitate just that.
     *
     *
     * @throws This class throws various exceptions from `fourdst::composition::exceptions` for invalid operations, such as using unregistered symbols or providing invalid abundances.
     *
     * @par Basic Example:
     * @code
     * Composition comp;
     * comp.registerSymbol("H-1");
     * comp.registerSymbol("He-4");
     * comp.setMolarAbundance("H-1", 0.75);
     * comp.setMolarAbundance("He-4", 0.25); // Note Molar Abundances do not need to sum to 1
     * @endcode
     *
     */
    // ReSharper disable once CppClassCanBeFinal
    class Composition final : public CompositionAbstract {
    private:
        /**
         * @struct CompositionCache
         * @brief Caches computed properties of the composition to avoid redundant calculations.
         * @details This struct holds optional cached values for various computed properties of the composition,
         * such as canonical composition, mass fractions, number fractions, molar abundances, sorted species,
         * sorted symbols, and electron abundance. The cache can be cleared when the composition is modified.
         */
        struct CompositionCache {
            std::optional<CanonicalComposition> canonicalComp; ///< Cached canonical composition data.
            std::optional<std::vector<double>> massFractions; ///< Cached vector of mass fractions.
            std::optional<std::vector<double>> numberFractions; ///< Cached vector of number fractions.
            std::optional<std::vector<double>> molarAbundances; ///< Cached vector of molar abundances.
            std::optional<std::vector<atomic::Species>> sortedSpecies; ///< Cached vector of sorted species (by mass).
            std::optional<std::vector<std::string>> sortedSymbols; ///< Cached vector of sorted species (by mass).
            std::optional<double> Ye; ///< Cached electron abundance.

            /**
             * @brief Clears all cached values.
             */
            void clear() {
                canonicalComp = std::nullopt;
                massFractions = std::nullopt;
                numberFractions = std::nullopt;
                molarAbundances = std::nullopt;
                sortedSymbols = std::nullopt;
                sortedSpecies = std::nullopt;
                Ye = std::nullopt;
            }

            /**
             * @brief Checks if the cache is clear (i.e., all cached values are empty).
             * @return True if the cache is clear, false otherwise.
             */
            [[nodiscard]] bool is_clear() const {
                return !canonicalComp.has_value() && !massFractions.has_value() &&
                       !numberFractions.has_value() && !molarAbundances.has_value() && !sortedSymbols.has_value() &&
                       !Ye.has_value() && !sortedSpecies.has_value();
            }
        };
    private:
        /**
         * @brief Gets the logger instance for the Composition class. This is static to ensure that all composition
         * objects share the same logger instance.
         * @return pointer to the logger instance.
         */
        static quill::Logger* getLogger() {
            static quill::Logger* logger = logging::LogManager::getInstance().getLogger("log");
            return logger;
        }

        std::set<atomic::Species> m_registeredSpecies; ///< Set of registered species in the composition.
        std::map<atomic::Species, double> m_molarAbundances; ///< Map of species to their molar abundances.

        mutable CompositionCache m_cache; ///< Cache for computed properties to avoid redundant calculations.

    public:
        /**
         * @brief Default constructor.
         * @details Creates an empty Composition object. No symbols or species are registered initially; however,
         * the user can register symbols or species later using the provided methods.
         */
        Composition() = default;

        /**
         * @brief Default destructor.
         */
        ~Composition() override = default;

        /**
         * @brief Constructs a Composition and registers the given symbols from a vector.
         * @param symbols The symbols to register.
         * @throws exceptions::UnknownSymbolError if any symbol is invalid. Symbols are invalid if they are not registered at compile time in the atomic species database (`fourdst/atomic/species.h`).
         * @par Example:
         * @code
         * std::vector<std::string> symbols = {"H-1", "O-16"};
         * Composition comp(symbols);
         * @endcode
         */
        explicit Composition(const std::vector<std::string>& symbols);

        /**
         * @brief Constructs a Composition and registers the given species from a vector.
         * @param species The species to register.
         * @par Example:
         * @code
         * std::vector<fourdst::atomic::Species> species = {fourdst::atomic::H_1, fourdst::atomic::O_16};
         * Composition comp(species);
         * @endcode
         *
         * @note Because species are strongly typed, this constructor does not need to check if the species is valid.
         * that is to say that the compiler will only allow valid species to be passed in. Therefore, this constructor
         * is marked noexcept and may therefore be slightly more performant than the symbol-based constructor.
         */
        explicit Composition(const std::vector<atomic::Species>& species);

        /**
         * @brief Constructs a Composition and registers the given symbols from a set.
         * @param symbols The symbols to register.
         * @throws exceptions::UnknownSymbolError if any symbol is invalid. Symbols are invalid if they are not registered at compile time in the atomic species database (`fourdst/atomic/species.h`).
         * @par Example:
         * @code
         * std::set<std::string> symbols = {"H-1", "O-16"};
         * Composition comp(symbols);
         * @endcode
         */
        explicit Composition(const std::set<std::string>& symbols);

        /**
         * @brief Constructs a Composition and registers the given species from a set.
         * @param species The species to register.
         * @par Example:
         * @code
         * std::set<fourdst::atomic::Species> species = {fourdst::atomic::H_1, fourdst::atomic::O_16};
         * Composition comp(species);
         * @endcode
         *
         * @note Because species are strongly typed, this constructor does not need to check if the species is valid.
         * that is to say that the compiler will only allow valid species to be passed in. Therefore, this constructor
         * is marked noexcept and may therefore be slightly more performant than the symbol-based constructor.
         */
        explicit Composition(const std::set<atomic::Species>& species);

        /**
         * @brief Constructs a Composition from symbols and their corresponding molar abundances.
         * @param symbols The symbols to register.
         * @param molarAbundances The corresponding molar abundances for each symbol.
         * @throws exceptions::UnknownSymbolError if any symbol is invalid. Symbols are invalid if they are not registered at compile time in the atomic species database (`fourdst/atomic/species.h`).
         * @throws exceptions::InvalidCompositionError if the number of symbols does not match the number of molar abundances.
         * @par Example:
         * @code
         * std::vector<std::string> symbols = {"H-1", "O-16"};
         * std::vector<double> molarAbundances = {1.03, 0.6};
         * Composition comp(symbols, molarAbundances);
         * @endcode
         *
         * @note Molar abundances do not need to sum to 1.0, they are an absolute quantity.
         */
        Composition(const std::vector<std::string>& symbols, const std::vector<double>& molarAbundances);

        /**
         * @brief Constructs a Composition from species and their corresponding molar abundances.
         * @param species The species to register.
         * @param molarAbundances The corresponding molar abundances for each species.
         * @throws exceptions::InvalidCompositionError if the number of species does not match the number of molar abundances.
         * @par Example:
         * @code
         * std::vector<fourdst::atomic::Species> species = {fourdst::atomic::H_1, fourdst::atomic::O_16};
         * std::vector<double> molarAbundances = {1.03, 0.6};
         * Composition comp(species, molarAbundances);
         * @endcode
         *
         * @note Molar abundances do not need to sum to 1.0, they are an absolute quantity.
         */
        Composition(const std::vector<atomic::Species>& species, const std::vector<double>& molarAbundances);

        /**
         * @brief Constructs a Composition from symbols in a set and their corresponding molar abundances.
         * @param symbols The symbols to register.
         * @param molarAbundances The corresponding molar abundances for each symbol.
         * @throws exceptions::UnknownSymbolError if any symbol is invalid. Symbols are invalid if they are not registered at compile time in the atomic species database (`fourdst/atomic/species.h`).
         * @throws exceptions::InvalidCompositionError if the number of symbols does not match the number of molar abundances.
         * @par Example:
         * @code
         * std::set<std::string> symbols = {"H-1", "O-16"};
         * std::vector<double> molarAbundances = {1.03, 0.6};
         * Composition comp(symbols, molarAbundances);
         * @endcode
         *
         * @note Molar abundances do not need to sum to 1.0, they are an absolute quantity.
         */
        Composition(const std::set<std::string>& symbols, const std::vector<double>& molarAbundances);

        /**
         * @brief Constructs a Composition from another Composition.
         * @param composition The Composition to copy.
         */
        Composition(const Composition& composition);

        explicit Composition(const CompositionAbstract& composition);

        /**
         * @brief Assignment operator.
         * @param other The Composition to assign from.
         * @return A reference to this Composition.
         */
        Composition& operator=(Composition const& other);

        /**
         * @brief Registers a new symbol for inclusion in the composition.
         * @details A symbol must be registered before its abundance can be set.
         * @param symbol The symbol to register (e.g., "Fe-56").
         * @throws exceptions::UnknownSymbolError if the symbol is not in the atomic species database.
         * @par Example:
         * @code
         * Composition comp;
         * comp.registerSymbol("H-1");
         * comp.registerSymbol("He-4");
         * @endcode
         *
         * @note upon registering a symbol, its molar abundance is initialized to 0.0.
         */
        void registerSymbol(const std::string& symbol);

        /**
         * @brief Registers multiple new symbols.
         * @param symbols The symbols to register.
         * @throws exceptions::UnknownSymbolError if any symbol is invalid.
         * @par Example:
         * @code
         * std::vector<std::string> symbols = {"H-1", "O-16"};
         * Composition comp;
         * comp.registerSymbol(symbols);
         * @endcode
         *
         * @note upon registering a symbol, its molar abundance is initialized to 0.0. Therefore, registering a vector
         * of symbols will initialize all their molar abundances to 0.0.
         */
        void registerSymbol(const std::vector<std::string>& symbols);

        /**
         * @brief Registers a new species by extracting its symbol.
         * @param species The species to register.
         * @par Example:
         * @code
         * #include "fourdst/composition/species.h"
         *
         * Composition comp;
         * comp.registerSpecies(fourdst::atomic::C_12);
         * @endcode
         *
         * @note Because species are strongly typed, this method does not need to check if the species is valid.
         * that is to say that the compiler will only allow valid species to be passed in. Therefore, this method
         * is marked noexcept and may therefore be slightly more performant than the symbol-based method.
         *
         * @note upon registering a species, its molar abundance is initialized to 0.0.
         *
         * @note All species are in the `fourdst/atomic/species.h` header file. These can be accessed directly through
         * their fully qualified names (e.g., `fourdst::atomic::C_12` for Carbon-12). Alternatively, these can also
         * be accessed through a string-indexed map located in `fourdst/atomic/species.h` called `fourdst::atomic::species` (
         * e.g., `fourdst::atomic::species.at("C-12")` for Carbon-12).
         */
        void registerSpecies(const atomic::Species& species) noexcept;


        /**
         * @brief Registers a vector of new species.
         * @param species The vector of species to register.
         * @par Example:
         * @code
         * #include "fourdst/composition/species.h"
         * Composition comp;
         * std::vector<fourdst::atomic::Species> my_species = { ... };
         * comp.registerSpecies(my_species);
         * @endcode
         *
         * @note upon registering a species, its molar abundance is initialized to 0.0. Therefore, registering a vector
         * of species will initialize all their molar abundances to 0.0.
         *
         * @note All species are in the `fourdst/atomic/species.h` header file. These can be accessed directly through
         * their fully qualified names (e.g., `fourdst::atomic::C_12` for Carbon-12). Alternatively, these can also
         * be accessed through a string-indexed map located in `fourdst/atomic/species.h` called `fourdst::atomic::species` (
         * e.g., `fourdst::atomic::species.at("C-12")` for Carbon-12).
         */
        void registerSpecies(const std::vector<atomic::Species>& species) noexcept;

        /**
         * @brief Checks if a given species is present in the composition.
         * @param species The isotope to check for.
         * @return True if the species is in the composition, false otherwise.
         */
        [[nodiscard]] bool contains(const atomic::Species& species) const noexcept override;

        /**
         * @brief Checks if a given symbol is present in the composition.
         * @param symbol The symbol to check for.
         * @return True if the symbol is in the composition, false otherwise.
         * @throws exceptions::UnknownSymbolError if the symbol is not in the atomic species database.
         */
        [[nodiscard]] bool contains(const std::string& symbol) const override;

        /**
         * @brief Gets the number of registered species in the composition.
         * @return The number of registered species.
         */
        [[nodiscard]] size_t size() const noexcept override;

        /**
         * @brief Sets the molar abundance for a given symbol.
         * @param symbol The symbol to set the molar abundance for.
         * @param molar_abundance The molar abundance to set.
         *
         * @throws exceptions::UnknownSymbolError if the symbol is not in the atomic species database.
         * @throws exceptions::UnregisteredSymbolError if the symbol is not in the composition.
         * @throws exceptions::InvalidCompositionError if the molar abundance is negative.
         *
         * @par Example:
         * @code
         * Composition comp({"H-1", "He-4"});
         * comp.setMolarAbundance("H-1", 1.0);
         * comp.setMolarAbundance("He-4", 0.5);
         * @endcode
         */
        void setMolarAbundance(
            const std::string& symbol,
            const double& molar_abundance
        );

        /**
         * @brief Sets the molar abundance for a given isotope.
         * @param species The isotope to set the molar abundance for.
         * @param molar_abundance The molar abundance to set.
         *
         * @throws exceptions::UnregisteredSymbolError if the isotope is not registered in the composition.
         * @throws exceptions::InvalidCompositionError if the molar abundance is negative.
         *
         * @par Example:
         * @code
         * #include "fourdst/composition/species.h"
         * Composition comp({fourdst::atomic::H_1, fourdst::atomic::He_4});
         * comp.setMolarAbundance(fourdst::atomic::H_1, 1.0);
         * comp.setMolarAbundance(fourdst::atomic::He_4, 0.5);
         * @endcode
         *
         * @note Since this method does not need to validate the species exists in the database, it will generally be
         * slightly more performant than the symbol-based method.
         */
        void setMolarAbundance(
            const atomic::Species& species,
            const double& molar_abundance
        );

        /**
         * @brief Sets the molar abundances for a list of symbols.
         * @param symbols The symbols to set the molar abundances for.
         * @param molar_abundances The molar abundances to set.
         *
         * @throws exceptions::UnknownSymbolError if any symbol is not in the atomic species database.
         * @throws exceptions::UnregisteredSymbolError if any symbol is not in the composition.
         * @throws exceptions::InvalidCompositionError if any molar abundance is negative.
         *
         * @par Example:
         * @code
         * Composition comp({"H-1", "He-4"});
         * comp.setMolarAbundance({"H-1", "He-4"}, {1.0, 0.5});
         * @endcode
         */
        void setMolarAbundance(
            const std::vector<std::string>& symbols,
            const std::vector<double>& molar_abundances
        );

        /**
         * @brief Sets the molar abundances for a list of isotopes.
         * @param species The isotopes to set the molar abundances for.
         * @param molar_abundances The molar abundances to set.
         *
         * @throws exceptions::UnregisteredSymbolError if any isotope is not registered in the composition.
         * @throws exceptions::InvalidCompositionError if any molar abundance is negative.
         *
         * @par Example:
         * @code
         * #include "fourdst/composition/species.h"
         * Composition comp({fourdst::atomic::H_1, fourdst::atomic::He_4});
         * comp.setMolarAbundance({fourdst::atomic::H_1, fourdst::atomic::He_4}, {1.0, 0.5});
         * @endcode
         *
         * @note Since this method does not need to validate the species exists in the database, it will generally be
         * slightly more performant than the symbol-based method.
         */
        void setMolarAbundance(
            const std::vector<atomic::Species>& species,
            const std::vector<double>& molar_abundances
        );

        /**
         * @brief Sets the molar abundances for a set of symbols.
         * @param symbols The symbols to set the molar abundances for.
         * @param molar_abundances The molar abundances to set.
         *
         * @throws exceptions::UnknownSymbolError if any symbol is not in the atomic species database.
         * @throws exceptions::UnregisteredSymbolError if any symbol is not in the composition.
         * @throws exceptions::InvalidCompositionError if any molar abundance is negative.
         *
         * @par Example:
         * @code
         * std::set<std::string> symbols = {"H-1", "He-4"};
         * Composition comp(symbols);
         * comp.setMolarAbundance(symbols, {1.0, 0.5});
         * @endcode
         */
        void setMolarAbundance(
            const std::set<std::string>& symbols,
            const std::vector<double>& molar_abundances
        );

        /**
         * @brief Sets the molar abundances for a set of isotopes.
         * @param species The isotopes to set the molar abundances for.
         * @param molar_abundances The molar abundances to set.
         *
         * @throws exceptions::UnregisteredSymbolError if any isotope is not registered in the composition.
         * @throws exceptions::InvalidCompositionError if any molar abundance is negative.
         *
         * @par Example:
         * @code
         * #include "fourdst/composition/species.h"
         * std::set<fourdst::atomic::Species> species = {fourdst::atomic::H_1, fourdst::atomic::He_4};
         * Composition comp(species);
         * comp.setMolarAbundance(species, {1.0, 0.5});
         * @endcode
         *
         * @note Since this method does not need to validate the species exists in the database, it will generally be
         * slightly more performant than the symbol-based method.
         */
        void setMolarAbundance(
            const std::set<atomic::Species>& species,
            const std::vector<double>& molar_abundances
        );

        /**
         * @brief Gets the registered symbols.
         * @return A set of registered symbols.
         *
         * @note This method will construct a new set each time it is called. If you need just need access to the
         * registered species, consider using `getRegisteredSpecies()` instead which returns a constant reference
         * to the internal set.
         */
        [[nodiscard]] std::set<std::string> getRegisteredSymbols() const noexcept override;

        /**
         * @brief Get a set of all species that are registered in the composition.
         * @return A set of `atomic::Species` objects registered in the composition.
         *
         * @note This will return a constant reference to the internal m_registeredSpecies set, therefore the return
         * value of this method will only be valid as long as the Composition object is valid (i.e. it cannot
         * outlive the Composition object it was called on).
         */
        [[nodiscard]] const std::set<atomic::Species> &getRegisteredSpecies() const noexcept override;

        /**
         * @brief Gets the mass fractions of all species in the composition.
         * @return An unordered map of symbols to their mass fractions.
         *
         * @note This method will construct a new unordered map each time it is called.
         */
        [[nodiscard]] std::unordered_map<atomic::Species, double> getMassFraction() const noexcept override;

        /**
         * @brief Gets the mass fraction for a given symbol. See the overload for species-based lookup for more details
         * on how mass fractions are calculated.
         * @param symbol The symbol to get the mass fraction for.
         * @return The mass fraction for the given symbol.
         * @throws exceptions::UnknownSymbolError if the symbol is not in the atomic species database.
         * @throws exceptions::UnregisteredSymbolError if the symbol is not in the composition.
         */
        [[nodiscard]] double getMassFraction(const std::string& symbol) const override;

        /**
         * @brief Gets the mass fraction for a given species.
         * @details The mass fraction X_i for a species is calculated using the formula:
         * \f[
         * X_i = \frac{(Y_i \cdot A_i)}{\sum_j (Y_j \cdot A_j)}
         * \f]
         * where:
         * - \f$Y_i\f$ is the molar abundance of species i.
         * - \f$A_i\f$ is the atomic mass of species i.
         * - The denominator sums over all species j in the composition.
         *
         * This formula ensures that the mass fractions of all species sum to 1.0.
         *
         * @param species The species to get the mass fraction for.
         * @return The mass fraction for the given isotope.
         * @throws exceptions::UnregisteredSymbolError if the isotope is not registered in the composition.
         */
        [[nodiscard]] double getMassFraction(const atomic::Species& species) const override;

        /**
         * @brief Gets the number fraction for a given symbol. See the overload for species-based lookup for more details
         * on how number fractions are calculated.
         * @param symbol The symbol to get the number fraction for.
         * @return The number fraction for the given symbol.
         * @throws exceptions::UnknownSymbolError if the symbol is not in the atomic species database.
         * @throws exceptions::UnregisteredSymbolError if the symbol is not in the composition.
         */
        [[nodiscard]] double getNumberFraction(const std::string& symbol) const override;

        /**
         * @brief Gets the number fraction for a given species.
         * @details The number fraction Y_i for a species is calculated using the formula:
         * \f[
         * X_i = \frac{Y_i}{\sum_j Y_j}
         * \f]
         * where:
         * - \f$Y_i\f$ is the molar abundance of species i.
         * - The denominator sums over all species j in the composition.
         *
         * This formula ensures that the number fractions of all species sum to 1.0.
         *
         * @param species The species to get the number fraction for.
         * @return The number fraction for the given isotope.
         * @throws exceptions::UnregisteredSymbolError if the isotope is not registered in the composition.
         */
        [[nodiscard]] double getNumberFraction(const atomic::Species& species) const override;

        /**
         * @brief Gets the number fractions of all species in the composition.
         * @return An unordered map of symbols to their number fractions.
         *
         * @note This method will construct a new unordered map each time it is called.
         */
        [[nodiscard]] std::unordered_map<atomic::Species, double> getNumberFraction() const noexcept override;

        /**
         * @brief Gets the molar abundances of all species in the composition.
         * @throws exceptions::UnknownSymbolError if any symbol is not in the atomic species database.
         * @throws exceptions::UnregisteredSymbolError if any symbol is not in the composition.
         * @return The molar abundance of the symbol.
         *
         * @note These are the most performant quantities to retrieve since they are stored directly in the composition object and
         * require no computation. This overload is slightly less performant than the species-based overload since it
         * needs to validate the symbol exists in the atomic species database.
         */
        [[nodiscard]] double getMolarAbundance(const std::string& symbol) const override;

        /**
         * @brief Gets the molar abundance for a given species.
         * @param species The species to get the molar abundance for.
         * @return The molar abundance for the given isotope.
         * @throws exceptions::UnregisteredSymbolError if the isotope is not registered in the composition.
         *
         * @note These are the most performant quantities to retrieve since they are stored directly in the composition object and
         * require no computation.
         */
        [[nodiscard]] double getMolarAbundance(const atomic::Species& species) const override;

        /**
         * @brief Compute the mean particle mass of the composition.
         * @details The mean particle mass is calculated using the formula:
         * \f[
         * \bar{A} = \frac{\sum_i (Y_i \cdot A_i)}{\sum_j Y_j}
         * \f]
         * where:
         * - \f$Y_i\f$ is the molar abundance of species i.
         * - \f$A_i\f$ is the atomic mass of species i.
         * - The sums run over all species i in the composition.
         *
         * @return Mean particle mass in atomic mass units (g/mol).
         */
        [[nodiscard]] double getMeanParticleMass() const noexcept override;

        /**
         * @brief Compute the electron abundance of the composition.
         * @details The electron abundance is calculated using the formula:
         * \f[
         * Y_e = \sum_i (Y_i \cdot Z_i)
         * \f]
         * where:
         * - \f$Y_i\f$ is the molar abundance of species i.
         * - \f$Z_i\f$ is the atomic number (number of protons) of species i.
         *
         *
         * @return Ye (electron abundance).
         */
        [[nodiscard]] double getElectronAbundance() const noexcept override;


        /**
         * @brief Compute the canonical composition (X, Y, Z) of the composition.
         * @details The canonical composition is defined as:
         * - X: mass fraction of hydrogen (\f$\sum_{i=1}^{7}X_{^{i}H}\f$)
         * - Y: mass fraction of helium (\f$\sum_{i=3}^{10}X_{^{i}He}\f$)
         * - Z: mass fraction of all other elements (Everything else)
         *
         * The canonical composition is computed by summing the mass fractions of all registered species
         * in the composition according to their element type.
         *
         * @return A CanonicalComposition struct containing the X, Y, and Z values.
         *
         * @throws exceptions::InvalidCompositionError if, after constructing the canonical composition, the sum X + Y + Z is not approximately equal to 1.0 (within a tolerance of 1e-16)
         */
        [[nodiscard]] CanonicalComposition getCanonicalComposition() const;

        /**
         * @brief Get a uniform vector representation of the mass fraction stored in the composition object sorted such that the lightest species is at index 0 and the heaviest is at the last index.
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @return the vector of mass fractions sorted by species mass (lightest to heaviest).
         */
        [[nodiscard]] std::vector<double> getMassFractionVector() const noexcept override;

        /**
         * @brief Get a uniform vector representation of the number fractions stored in the composition object sorted such that the lightest species is at index 0 and the heaviest is at the last index.
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @return the vector of number fractions sorted by species mass (lightest to heaviest).
         */
        [[nodiscard]] std::vector<double> getNumberFractionVector() const noexcept override;

        /**
         * @brief Get a uniform vector representation of the molar abundances stored in the composition object sorted such that the lightest species is at index 0 and the heaviest is at the last index.
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @return the vector of molar abundances sorted by species mass (lightest to heaviest).
         */
        [[nodiscard]] std::vector<double> getMolarAbundanceVector() const noexcept override;

        /**
         * @brief get the index in the sorted vector representation for a given symbol
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @param symbol the symbol to look up the index for. Note that this is the index species data will be at if you were to call getMolarAbundanceVector(), getMassFractionVector(), or getNumberFractionVector()
         * @throws exceptions::UnknownSymbolError if the symbol is not in the atomic species database.
         * @throws exceptions::UnregisteredSymbolError if the symbol is not registered in the composition
         * @return The index of the symbol in the sorted vector representation.
         */
        [[nodiscard]] size_t getSpeciesIndex(const std::string& symbol) const override;

        /**
         * @brief get the index in the sorted vector representation for a given symbol
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @param species the species to look up the index for. Note that this is the index species data will be at if you were to call getMolarAbundanceVector(), getMassFractionVector(), or getNumberFractionVector()
         * @throws exceptions::UnregisteredSymbolError if the symbol is not registered in the composition
         * @return The index of the symbol in the sorted vector representation.
         */
        [[nodiscard]] size_t getSpeciesIndex(const atomic::Species& species) const override;

        /**
         * @brief Get the species at a given index in the sorted vector representation.
         * @details This is primarily useful for external libraries which need to ensure that vector representation uniformity is maintained
         * @param index The index in the sorted vector representation for which to return the species. Must be in [0, N-1] where N is the number of registered species.
         * @throws std::out_of_range if the index is out of range.
         * @return The species at the given index in the sorted vector representation.
         */
        [[nodiscard]] atomic::Species getSpeciesAtIndex(size_t index) const override;

        [[nodiscard]] std::unique_ptr<CompositionAbstract> clone() const override;

        /**
         * @brief Overloaded output stream operator for Composition.
         * @param os The output stream.
         * @param composition The Composition to output.
         * @return The output stream.
         */
        friend std::ostream& operator<<(std::ostream& os, const Composition& composition);

        /**
         * @brief Returns an iterator to the beginning of the molar abundance map
         * @return An iterator to the beginning.
         *
         * @par Example:
         * @code
         * Composition comp({"H-1", "He-4", "C-12"}, {1.02, 0.56, 0.02});
         *
         * for (const auto& [sp, y] : comp) {
         *      std::cout << "Species: " << sp << ", Molar Abundance: " << y << std::endl;
         * }
         * @endcode
         *
         * @note Because we store molar abundances as a sorted map, keyed by species. And Because the < and > operators
         * for species are defined based on their atomic mass. When iterating over the molar abundance map, species will be
         * seen in order from lightest to heaviest.
         */
        [[nodiscard]] std::map<atomic::Species, double>::iterator begin() override {
            return m_molarAbundances.begin();
        }

        /**
         * @brief Returns a const iterator to the beginning of the molar abundance map.
         * @return A const iterator to the beginning.
         *
         * @par Example:
         * @code
         * Composition comp({"H-1", "He-4", "C-12"}, {1.02, 0.56, 0.02});
         *
         * for (const auto& [sp, y] : comp) {
         *      std::cout << "Species: " << sp << ", Molar Abundance: " << y << std::endl;
         * }
         * @endcode
         *
         * @note Because we store molar abundances as a sorted map, keyed by species. And Because the < and > operators
         * for species are defined based on their atomic mass. When iterating over the molar abundance map, species will be
         * seen in order from lightest to heaviest.
         */
        [[nodiscard]] std::map<atomic::Species, double>::const_iterator begin() const override {
            return m_molarAbundances.cbegin();
        }

        /**
         * @brief Returns an iterator to the end of the molar abundance map.
         * @return An iterator to the end.
         *
         * @par Example:
         * @code
         * Composition comp({"H-1", "He-4", "C-12"}, {1.02, 0.56, 0.02});
         *
         * for (const auto& [sp, y] : comp) {
         *      std::cout << "Species: " << sp << ", Molar Abundance: " << y << std::endl;
         * }
         * @endcode
         *
         * @note Because we store molar abundances as a sorted map, keyed by species. And Because the < and > operators
         * for species are defined based on their atomic mass. When iterating over the molar abundance map, species will be
         * seen in order from lightest to heaviest.
         */
        [[nodiscard]] std::map<atomic::Species, double>::iterator end() override {
            return m_molarAbundances.end();
        }

        /**
         * @brief Returns a const iterator to the end of the molar abundance map.
         * @return A const iterator to the end.
         *
         * @par Example:
         * @code
         * Composition comp({"H-1", "He-4", "C-12"}, {1.02, 0.56, 0.02});
         *
         * for (const auto& [sp, y] : comp) {
         *      std::cout << "Species: " << sp << ", Molar Abundance: " << y << std::endl;
         * }
         * @endcode
         *
         * @note Because we store molar abundances as a sorted map, keyed by species. And Because the < and > operators
         * for species are defined based on their atomic mass. When iterating over the molar abundance map, species will be
         * seen in order from lightest to heaviest.
         */
        [[nodiscard]] std::map<atomic::Species, double>::const_iterator end() const override {
            return m_molarAbundances.cend();
        }

    };

    inline bool operator==(const Composition& a, const Composition& b) noexcept {
        if (a.size() != b.size()) return false;

        // Compare species sets quickly
        if (a.getRegisteredSpecies() != b.getRegisteredSpecies())
            return false;

        // Compare all abundances
        for (auto itA = a.begin(), itB = b.begin();
             itA != a.end() && itB != b.end(); ++itA, ++itB) {
            if (itA->first != itB->first)
                return false;
            if (itA->second != itB->second)
                return false;
             }
        return true;
    }
}; // namespace fourdst::composition
