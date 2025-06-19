#pragma once
#include <unordered_map>
#include <string_view>
#include <string>
#include <iostream>
#include "atomicSpecies.h"

namespace serif::atomic {
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
        double m_atomicMass; //< Atomic mass
        double m_atomicMassUnc; //< Atomic mass uncertainty

        Species(const std::string_view name, const std::string_view el, const int nz, const int n, const int z, const int a, const double bindingEnergy, const std::string_view betaCode, const double betaDecayEnergy, const double atomicMass, const double atomicMassUnc)
            : m_name(name), m_el(el), m_nz(nz), m_n(n), m_z(z), m_a(a), m_bindingEnergy(bindingEnergy), m_betaCode(betaCode), m_betaDecayEnergy(betaDecayEnergy), m_atomicMass(atomicMass), m_atomicMassUnc(atomicMassUnc) {};

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
            m_atomicMass = species.m_atomicMass;
            m_atomicMassUnc = species.m_atomicMassUnc;
        }


        double mass() const {
            return m_atomicMass;
        }

        double massUnc() const {
            return m_atomicMassUnc;
        }

        double bindingEnergy() const {
            return m_bindingEnergy;
        }

        double betaDecayEnergy() const {
            return m_betaDecayEnergy;
        }

        std::string_view betaCode() const {
            return m_betaCode;
        }

        std::string_view name() const {
            return m_name;
        }

        std::string_view el() const {
            return m_el;
        }

        int nz() const {
            return m_nz;
        }

        int n() const {
            return m_n;
        }

        int z() const {
            return m_z;
        }

        int a() const {
            return m_a;
        }

        friend std::ostream& operator<<(std::ostream& os, const Species& species) {
            os << species.m_name << " (" << species.m_atomicMass << " u)";
            return os;
        }

        friend bool operator==(const Species& lhs, const Species& rhs);
        friend bool operator!=(const Species& lhs, const Species& rhs);
    };
    inline bool operator==(const Species& lhs, const Species& rhs) {
        return (lhs.m_name == rhs.m_name);
    }
    inline bool operator!=(const Species& lhs, const Species& rhs) {
        return (lhs.m_name != rhs.m_name);
    }

}

namespace std {
    template<>
    struct hash<serif::atomic::Species> {
        size_t operator()(const serif::atomic::Species& s) const noexcept {
            return std::hash<std::string>()(s.m_name);
        }
    };
} // namespace std