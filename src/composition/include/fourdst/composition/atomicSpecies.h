#pragma once
#include <unordered_map>
#include <string_view>
#include <string>
#include <iostream>

namespace fourdst::atomic {
    inline double convert_jpi_to_double(const std::string& jpi_string);

    struct Species {
        std::string m_name; //< Name of the species
        std::string m_el; //< Element symbol
        int m_nz; //< NZ
        int m_n; //< N
        int m_z; //< Z
        int m_a; //< A
        double m_bindingEnergy; //< Binding energy
        std::string m_betaCode; //< Beta decay code
        double m_betaDecayEnergy; //< Beta decay energy
        double m_halfLife_s; //< Half-life in seconds
        std::string m_spinParity; //< Spin and parity
        std::string m_decayModes; //< Decay modes
        double m_atomicMass; //< Atomic mass
        double m_atomicMassUnc; //< Atomic mass uncertainty
        double m_spin = 0.0; //< Spin of the species, default is 0.0

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
        m_atomicMassUnc(atomicMassUnc) {
            m_spin = convert_jpi_to_double(m_spinParity);
        };

        //Copy constructor
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
            m_spin = convert_jpi_to_double(m_spinParity);
        }


        [[nodiscard]] double mass() const {
            return m_atomicMass;
        }

        [[nodiscard]] double massUnc() const {
            return m_atomicMassUnc;
        }

        [[nodiscard]] double halfLife() const {
            return m_halfLife_s;
        }

        [[nodiscard]] std::string_view spinParity() const {
            return m_spinParity;
        }

        [[nodiscard]] std::string_view decayModes() const {
            return m_decayModes;
        }

        [[nodiscard]] double bindingEnergy() const {
            return m_bindingEnergy;
        }

        [[nodiscard]] double betaDecayEnergy() const {
            return m_betaDecayEnergy;
        }

        [[nodiscard]] std::string_view betaCode() const {
            return m_betaCode;
        }

        [[nodiscard]] std::string_view name() const {
            return m_name;
        }

        [[nodiscard]] std::string_view el() const {
            return m_el;
        }

        [[nodiscard]] int nz() const {
            return m_nz;
        }

        [[nodiscard]] int n() const {
            return m_n;
        }

        [[nodiscard]] int z() const {
            return m_z;
        }

        [[nodiscard]] int a() const {
            return m_a;
        }

        [[nodiscard]] double spin() const {
            return m_spin;
        }

        friend std::ostream& operator<<(std::ostream& os, const Species& species) {
            os << species.m_name;
            return os;
        }

        friend bool operator==(const Species& lhs, const Species& rhs);
        friend bool operator!=(const Species& lhs, const Species& rhs);
        friend bool operator<(const Species& lhs, const Species& rhs);
        friend bool operator>(const Species& lhs, const Species& rhs);
    };
    inline bool operator==(const Species& lhs, const Species& rhs) {
        return (lhs.m_name == rhs.m_name);
    }
    inline bool operator!=(const Species& lhs, const Species& rhs) {
        return (lhs.m_name != rhs.m_name);
    }
    inline bool operator<(const Species& lhs, const Species& rhs) {
        return (lhs.m_atomicMass < rhs.m_atomicMass);
    }
    inline bool operator>(const Species& lhs, const Species& rhs) {
        return (lhs.m_atomicMass > rhs.m_atomicMass);
    }

    inline double convert_jpi_to_double(const std::string& jpi_string) {
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

template<>
struct std::hash<fourdst::atomic::Species> {
    size_t operator()(const fourdst::atomic::Species& s) const noexcept {
        return std::hash<std::string>()(s.m_name);
    }
}; // namespace std
