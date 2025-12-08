#pragma once

#include <vector>
#include <iterator>
#include <utility>
#include <compare>

#include "fourdst/atomic/atomicSpecies.h"

namespace fourdst::composition::detail {

    template <bool IsConst>
    class CompositionIterator {
    public:
        using iterator_category = std::random_access_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        // Returns a pair of references. Natively supports structured binding [sp, y]
        using value_type        = std::pair<const atomic::Species, double>;


        // Define reference types based on const-ness
        using SpeciesRef = const atomic::Species&;
        using AbundRef   = std::conditional_t<IsConst, const double&, double&>;
        using reference  = std::pair<SpeciesRef, AbundRef>;

        struct ArrowProxy {
            reference m_payload;
            const reference* operator->() const { return &m_payload; }
        };

        using pointer = ArrowProxy;

    private:
        using SpecIt = std::vector<atomic::Species>::const_iterator;
        using AbunIt = std::conditional_t<IsConst,
                                          std::vector<double>::const_iterator,
                                          std::vector<double>::iterator>;

        SpecIt m_sIt;
        AbunIt m_aIt;

    public:
        CompositionIterator() = default;
        CompositionIterator(SpecIt sIt, AbunIt aIt) : m_sIt(sIt), m_aIt(aIt) {}

        template <bool WasConst, typename = std::enable_if_t<IsConst && !WasConst>>
        CompositionIterator(const CompositionIterator<WasConst>& other)
            : m_sIt(other.getSpeciesIt()), m_aIt(other.getAbundanceIt()) {}

        [[nodiscard]] SpecIt getSpeciesIt() const { return m_sIt; }
        [[nodiscard]] AbunIt getAbundanceIt() const { return m_aIt; }

        reference operator*() const {
            return { *m_sIt, *m_aIt };
        }

        ArrowProxy operator->() const {
            return ArrowProxy{ **this };
        }

        reference operator[](difference_type n) const {
            return { *(m_sIt + n), *(m_aIt + n) };
        }

        // --- Movement ---
        CompositionIterator& operator++() { ++m_sIt; ++m_aIt; return *this; }
        CompositionIterator operator++(int) { auto tmp = *this; ++(*this); return tmp; }

        CompositionIterator& operator--() { --m_sIt; --m_aIt; return *this; } // FIXED
        CompositionIterator operator--(int) { auto tmp = *this; --(*this); return tmp; }

        CompositionIterator& operator+=(difference_type n) { m_sIt += n; m_aIt += n; return *this; }
        CompositionIterator& operator-=(difference_type n) { m_sIt -= n; m_aIt -= n; return *this; }

        // --- Arithmetic ---
        friend CompositionIterator operator+(CompositionIterator it, difference_type n) { return it += n; }

        // Commutative addition (n + it)
        friend CompositionIterator operator+(difference_type n, CompositionIterator it) { return it += n; }
        friend CompositionIterator operator-(CompositionIterator it, difference_type n) { return it -= n; }

        // Difference between iterators
        friend difference_type operator-(const CompositionIterator& lhs, const CompositionIterator& rhs) {
            return lhs.m_sIt - rhs.m_sIt;
        }

        template <bool R>
        bool operator==(const CompositionIterator<R>& other) const { return m_sIt == other.getSpeciesIt(); }

        template <bool R>
        bool operator!=(const CompositionIterator<R>& other) const { return m_sIt != other.getSpeciesIt(); }

        template <bool R>
        bool operator<(const CompositionIterator<R>& other) const { return m_sIt < other.getSpeciesIt(); }

        template <bool R>
        bool operator>(const CompositionIterator<R>& other) const { return m_sIt > other.getSpeciesIt(); }

        template <bool R>
        bool operator<=(const CompositionIterator<R>& other) const { return m_sIt <= other.getSpeciesIt(); }

        template <bool R>
        bool operator>=(const CompositionIterator<R>& other) const { return m_sIt >= other.getSpeciesIt(); }
    };

}