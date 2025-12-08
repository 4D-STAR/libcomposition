#pragma once


#include <format>
#include <string_view>
#include <string>
#include <optional>
#include <limits>
#include <stdexcept>

/**
 * @namespace fourdst::atomic
 * @brief Contains classes and functions related to atomic data, such as properties of atomic species.
 */
namespace fourdst::atomic {
    /**
     * @brief Converts a spin-parity string (JPI string) to a double-precision floating-point number.
     * @param jpi_string The spin-parity string to convert (e.g., "1/2+", "5/2-", "0+").
     * @return The spin value as a double. Returns `NaN` for invalid or unparsable strings.
     */
    inline double convert_jpi_to_double(const std::string& jpi_string) noexcept;

    /**
     * @struct Species
     * @brief Represents an atomic species (isotope) with its fundamental physical properties.
     *
     * This struct holds data parsed from nuclear data libraries, such as atomic mass,
     * half-life, and spin. It is a fundamental data structure for representing the
     * components of a material composition.
     *
     * @note This struct is designed to be lightweight and is primarily a data container.
     *
     * @par Usage Example
     * @code
     * #include "fourdst/composition/atomicSpecies.h"
     * #include <iostream>
     *
     * int main() {
     *     // Create a species for Deuterium (H-2)
     *     fourdst::atomic::Species deuterium(
     *         "H2", "H", 1002, 1, 1, 2, 2224.52, "", 0.0, -1.0, "1+", "", 2.0141017781, 4.0e-11
     *     );
     *
     *     std::cout << "Species: " << deuterium.name() << std::endl;
     *     std::cout << "Atomic Mass: " << deuterium.mass() << " u" << std::endl;
     *     std::cout << "Spin: " << deuterium.spin() << std::endl;
     *
     *     return 0;
     * }
     * @endcode
     */
    struct Species {
        std::string m_name; ///< Name of the species (e.g., "Fe56").
        std::string m_el; ///< Element symbol (e.g., "Fe").
        int m_nz; ///< NZ identifier, typically 1000*Z + A.
        int m_n; ///< Number of neutrons.
        int m_z; ///< Atomic number (number of protons).
        int m_a; ///< Mass number (N + Z).
        double m_bindingEnergy; ///< Binding energy in keV.
        std::string m_betaCode; ///< Beta decay code.
        double m_betaDecayEnergy; ///< Beta decay energy in keV.
        double m_halfLife_s; ///< Half-life in seconds. A value of -1.0 typically indicates stability.
        std::string m_spinParity; ///< Spin and parity as a string (e.g., "1/2-").
        std::string m_decayModes; ///< Decay modes as a string.
        double m_atomicMass; ///< Atomic mass in atomic mass units (u).
        double m_atomicMassUnc; ///< Uncertainty in the atomic mass.
        mutable std::optional<double> m_spin = std::nullopt; ///< Nuclear spin as a double, derived from m_spinParity.

        /**
         * @brief Constructs a Species object with detailed properties.
         *
         * @param name Name of the species.
         * @param el Element symbol.
         * @param nz NZ identifier.
         * @param n Number of neutrons.
         * @param z Atomic number.
         * @param a Mass number.
         * @param bindingEnergy Binding energy.
         * @param betaCode Beta decay code.
         * @param betaDecayEnergy Beta decay energy.
         * @param halfLife_s Half-life in seconds.
         * @param spinParity Spin and parity string.
         * @param decayModes Decay modes string.
         * @param atomicMass Atomic mass.
         * @param atomicMassUnc Atomic mass uncertainty.
         *
         * @post The `m_spin` member is initialized by parsing `m_spinParity` using `convert_jpi_to_double`.
         */
        Species(
            const std::string_view name,
            const std::string_view el,
            const int nz,
            const int n,
            const int z,
            const int a,
            const double bindingEnergy,
            const std::string_view betaCode,
            const double betaDecayEnergy,
            const double halfLife_s,
            const std::string_view spinParity,
            const std::string_view decayModes,
            const double atomicMass,
            const double atomicMassUnc
        ) :
        m_name(name),
        m_el(el),
        m_nz(nz),
        m_n(n),
        m_z(z),
        m_a(a),
        m_bindingEnergy(bindingEnergy),
        m_betaCode(betaCode),
        m_betaDecayEnergy(betaDecayEnergy),
        m_halfLife_s(halfLife_s),
        m_spinParity(spinParity),
        m_decayModes(decayModes),
        m_atomicMass(atomicMass),
        m_atomicMassUnc(atomicMassUnc) {};

        /**
         * @brief Copy constructor for Species.
         * @param species The Species object to copy.
         * @post A new Species object is created as a deep copy of `species`. The `m_spin` member is re-calculated.
         */
        Species(const Species& species) {
            m_name = species.m_name;
            m_el = species.m_el;
            m_nz = species.m_nz;
            m_n = species.m_n;
            m_z = species.m_z;
            m_a = species.m_a;
            m_bindingEnergy = species.m_bindingEnergy;
            m_betaCode = species.m_betaCode;
            m_betaDecayEnergy = species.m_betaDecayEnergy;
            m_halfLife_s = species.m_halfLife_s;
            m_spinParity = species.m_spinParity;
            m_decayModes = species.m_decayModes;
            m_atomicMass = species.m_atomicMass;
            m_atomicMassUnc = species.m_atomicMassUnc;
        }


        /**
         * @brief Gets the atomic mass of the species.
         * @return The atomic mass in atomic mass units (u).
         */
        [[nodiscard]] double mass() const {
            return m_atomicMass;
        }

        /**
         * @brief Gets the uncertainty in the atomic mass.
         * @return The atomic mass uncertainty.
         */
        [[nodiscard]] double massUnc() const {
            return m_atomicMassUnc;
        }

        /**
         * @brief Gets the half-life of the species.
         * @return The half-life in seconds.
         */
        [[nodiscard]] double halfLife() const {
            return m_halfLife_s;
        }

        /**
         * @brief Gets the spin and parity as a string.
         * @return A string_view of the spin and parity (e.g., "1/2+").
         */
        [[nodiscard]] std::string_view spinParity() const {
            return m_spinParity;
        }

        /**
         * @brief Gets the decay modes as a string.
         * @return A string_view of the decay modes.
         */
        [[nodiscard]] std::string_view decayModes() const {
            return m_decayModes;
        }

        /**
         * @brief Gets the binding energy of the species.
         * @return The binding energy in keV.
         */
        [[nodiscard]] double bindingEnergy() const {
            return m_bindingEnergy;
        }

        /**
         * @brief Gets the beta decay energy of the species.
         * @return The beta decay energy in keV.
         */
        [[nodiscard]] double betaDecayEnergy() const {
            return m_betaDecayEnergy;
        }

        /**
         * @brief Gets the beta decay code.
         * @return A string_view of the beta decay code.
         */
        [[nodiscard]] std::string_view betaCode() const {
            return m_betaCode;
        }

        /**
         * @brief Gets the name of the species.
         * @return A string_view of the species name (e.g., "Fe56").
         */
        [[nodiscard]] std::string_view name() const {
            return m_name;
        }

        /**
         * @brief Gets the element symbol of the species.
         * @return A string_view of the element symbol (e.g., "Fe").
         */
        [[nodiscard]] std::string_view el() const {
            return m_el;
        }

        /**
         * @brief Gets the NZ identifier of the species.
         * @return The NZ identifier (1000*Z + A).
         */
        [[nodiscard]] int nz() const {
            return m_nz;
        }

        /**
         * @brief Gets the number of neutrons.
         * @return The number of neutrons (N).
         */
        [[nodiscard]] int n() const {
            return m_n;
        }

        /**
         * @brief Gets the atomic number (number of protons).
         * @return The atomic number (Z).
         */
        [[nodiscard]] int z() const {
            return m_z;
        }

        /**
         * @brief Gets the mass number.
         * @return The mass number (A = N + Z).
         */
        [[nodiscard]] int a() const {
            return m_a;
        }

        /**
         * @brief Gets the nuclear spin as a numeric value.
         * @return The spin as a double.
         */
        [[nodiscard]] double spin() const {
            if (!m_spin.has_value()) { // The spin calculation is very expensive, and we almost never need it so we only compute it the first time it is requested
                m_spin = convert_jpi_to_double(m_spinParity);
            }
            return m_spin.value();
        }

        /**
         * @brief Overloads the stream insertion operator for easy printing of a Species object.
         * @param os The output stream.
         * @param species The Species object to print.
         * @return The output stream with the species name.
         */
        friend std::ostream& operator<<(std::ostream& os, const Species& species) {
            os << species.m_name;
            return os;
        }

        friend bool operator==(const Species& lhs, const Species& rhs);
        friend bool operator!=(const Species& lhs, const Species& rhs);
        friend std::partial_ordering operator<=>(const Species &lhs, const Species &rhs);
    };
    /**
     * @brief Equality operator for Species. Compares based on name.
     * @param lhs The left-hand side Species.
     * @param rhs The right-hand side Species.
     * @return `true` if the names are identical, `false` otherwise.
     */
    inline bool operator==(const Species& lhs, const Species& rhs) {
        return (lhs.m_name == rhs.m_name);
    }
    /**
     * @brief Inequality operator for Species. Compares based on name.
     * @param lhs The left-hand side Species.
     * @param rhs The right-hand side Species.
     * @return `true` if the names are different, `false` otherwise.
     */
    inline bool operator!=(const Species& lhs, const Species& rhs) {
        return (lhs.m_name != rhs.m_name);
    }

    inline std::partial_ordering operator<=>(const Species &lhs, const Species &rhs) {
        if (const auto cmp = lhs.m_atomicMass <=> rhs.m_atomicMass; cmp != 0) {
            return cmp;
        }

        return lhs.m_name <=> rhs.m_name;
    }




    /**
     * @brief Converts a spin-parity string (JPI string) to a double-precision floating-point number representing the spin.
     *
     * @details
     * **Purpose and Usage:**
     * This function is a utility for converting the textual representation of nuclear spin and parity,
     * commonly found in nuclear data files (e.g., "5/2-"), into a numerical format (`double`)
     * that can be used in calculations. It is used internally by the `Species` constructor to
     * initialize the `m_spin` member.
     *
     * **Algorithm:**
     * The function robustly parses various formats of JPI strings:
     * 1.  **Sanitization:** It first removes common non-numeric characters like `(`, `)`, `*`, and `#`.
     * 2.  **Parity-Only:** Strings containing only `+` or `-` are treated as spin 0.
     * 3.  **Comma-Separated:** If multiple values are present (e.g., "3/2+,5/2+"), it considers only the first one.
     * 4.  **Parity Suffix:** It removes the trailing `+` or `-` parity indicator.
     * 5.  **Fractional vs. Integer:**
     *     - If the string contains a `/` (e.g., "5/2"), it parses the numerator and denominator to calculate the fraction.
     *     - Otherwise, it attempts to parse the string directly as a double (for integer spins like "1", "2").
     * 6.  **Error Handling:** If the string is empty, malformed, or results in division by zero,
     *     it returns `std::numeric_limits<double>::quiet_NaN()` to indicate a parsing failure.
     *     This function does not throw exceptions but relies on the return value for error signaling.
     *
     * @param jpi_string The spin-parity string to convert.
     *
     * @pre The input `jpi_string` is a standard C++ string.
     * @post The function returns a `double` representing the spin, or `NaN` if parsing fails. The input string is not modified.
     *
     * @return The spin value as a `double`. Returns `NaN` for invalid or unparsable strings.
     */
    inline double convert_jpi_to_double(const std::string& jpi_string) noexcept {
        std::string s = jpi_string;

        if (s.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        std::erase_if(s, [](const char c) {
            return c == '(' || c == ')' || c == '*' || c == '#';
        });

        if (s == "+" || s == "-") {
            return 0.0;
        }

        if (const size_t comma_pos = s.find(','); comma_pos != std::string::npos) {
            s = s.substr(0, comma_pos);
        }

        if (!s.empty() && (s.back() == '+' || s.back() == '-')) {
            s.pop_back();
        }

        if (s.empty()) {
            return std::numeric_limits<double>::quiet_NaN();
        }

        try {
            if (size_t slash_pos = s.find('/'); slash_pos != std::string::npos) {
                if (slash_pos == 0) {
                    s = "1" + s;
                    slash_pos = 1;
                }
                const std::string numerator_str = s.substr(0, slash_pos);
                const std::string denominator_str = s.substr(slash_pos + 1);
                if (denominator_str.empty()) {
                    return std::numeric_limits<double>::quiet_NaN();
                }
                const double numerator = std::stod(numerator_str);
                const double denominator = std::stod(denominator_str);
                if (denominator == 0.0) {
                    return std::numeric_limits<double>::quiet_NaN();
                }
                return numerator / denominator;
            } else {
                return std::stod(s);
            }
        } catch (const std::invalid_argument&) {
            return std::numeric_limits<double>::quiet_NaN();
        } catch (const std::out_of_range&) {
            return std::numeric_limits<double>::quiet_NaN();
        }
    }

}

/**
 * @brief Specialization of `std::hash` for `fourdst::atomic::Species`.
 *
 * @details
 * This allows `fourdst::atomic::Species` objects to be used as keys in unordered
 * associative containers like `std::unordered_map` and `std::unordered_set`.
 * The hash is computed based on the species' name (`m_name`), as it is expected
 * to be a unique identifier for each species.
 *
 * @par Usage Example
 * @code
 * #include "fourdst/composition/atomicSpecies.h"
 * #include <unordered_map>
 * #include <string>
 *
 * int main() {
 *     std::unordered_map<fourdst::atomic::Species, double> abundance;
 *     fourdst::atomic::Species h1("H1", ...);
 *     abundance[h1] = 0.999;
 *     return 0;
 * }
 * @endcode
 */

template<>
struct std::hash<fourdst::atomic::Species> {
    /**
         * @brief Computes the hash for a Species object.
         * @param s The Species object to hash.
         * @return The hash value of the species' name.
         */
    size_t operator()(const fourdst::atomic::Species& s) const noexcept {
        return std::hash<std::string>()(s.m_name);
    }
};
