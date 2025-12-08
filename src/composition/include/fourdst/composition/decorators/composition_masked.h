#pragma once

#include <vector>
#include <set>
#include <unordered_map>
#include <memory>
#include <string>

#include "fourdst/composition/composition_abstract.h"
#include "fourdst/composition/decorators/composition_decorator_abstract.h"
#include "fourdst/composition/iterators/composition_abstract_iterator.h"
#include "fourdst/atomic/atomicSpecies.h"

namespace fourdst::composition {
    class MaskedComposition final : public CompositionDecorator {
    public:
        using iterator = detail::CompositionIterator<false>;
        using const_iterator = detail::CompositionIterator<true>;

        MaskedComposition(
            const CompositionAbstract& baseComposition,
            const std::vector<atomic::Species>& activeSpecies
        );

        [[nodiscard]] bool contains(const atomic::Species &species) const noexcept override;
        [[nodiscard]] bool contains(const std::string &symbol) const override;

        [[nodiscard]] const std::vector<atomic::Species>& getRegisteredSpecies() const noexcept override;
        [[nodiscard]] std::set<std::string> getRegisteredSymbols() const noexcept override;

        [[nodiscard]] size_t size() const noexcept override;

        [[nodiscard]] std::unordered_map<atomic::Species, double> getMassFraction() const noexcept override;
        [[nodiscard]] std::unordered_map<atomic::Species, double> getNumberFraction() const noexcept override;

        [[nodiscard]] double getMassFraction(const std::string &symbol) const override;
        [[nodiscard]] double getMassFraction(const atomic::Species &species) const override;
        [[nodiscard]] double getNumberFraction(const std::string &symbol) const override;
        [[nodiscard]] double getNumberFraction(const atomic::Species &species) const override;
        [[nodiscard]] double getMolarAbundance(const std::string &symbol) const override;
        [[nodiscard]] double getMolarAbundance(const atomic::Species &species) const override;
        [[nodiscard]] double getMeanParticleMass() const noexcept override;

        [[nodiscard]] double getElectronAbundance() const noexcept override;

        [[nodiscard]] std::vector<double> getMassFractionVector() const noexcept override;
        [[nodiscard]] std::vector<double> getNumberFractionVector() const noexcept override;
        [[nodiscard]] std::vector<double> getMolarAbundanceVector() const noexcept override;

        [[nodiscard]] size_t getSpeciesIndex(const std::string &symbol) const override;
        [[nodiscard]] size_t getSpeciesIndex(const atomic::Species &species) const override;
        [[nodiscard]] atomic::Species getSpeciesAtIndex(size_t index) const override;

        [[nodiscard]] std::unique_ptr<CompositionAbstract> clone() const override;

        [[nodiscard]] iterator begin() override;
        [[nodiscard]] iterator end() override;

        [[nodiscard]] const_iterator begin() const override;
        [[nodiscard]] const_iterator end() const override;

        [[nodiscard]] size_t hash() const override;

    private:
        std::vector<atomic::Species> m_activeSpecies;
        std::vector<double> m_molarAbundances;
    };

}