#pragma once

#include <cstring>
#include <cmath>
#include <vector>
#include <bit>

#include "xxhash64.h"
#include "fourdst/composition/composition.h"
#include "fourdst/composition/composition_abstract.h"

namespace fourdst::composition::utils {
    struct CompositionHash {
        template <typename CompositionT>
        static uint64_t hash_exact(const CompositionT& comp) {
            uint64_t h0 = kSeed;
            uint64_t h1 = kSeed ^ kPrime1;
            uint64_t h2 = kSeed ^ kPrime2;
            uint64_t h3 = kSeed ^ kPrime3;

            auto it = comp.begin();
            size_t remaining = comp.size();

            while (remaining >= 4) {
                const auto& p0 = *it;
                ++it;
                h0 ^= pack_species_id(p0.first);
                h0 = mum(h0, kPrime1);
                h0 ^= normalize_double_bits(p0.second);
                h0 = mum(h0, kPrime2);

                const auto& p1 = *it;
                ++it;
                h1 ^= pack_species_id(p1.first);
                h1 = mum(h1, kPrime1);
                h1 ^= normalize_double_bits(p1.second);
                h1 = mum(h1, kPrime2);

                const auto& p2 = *it;
                ++it;
                h2 ^= pack_species_id(p2.first);
                h2 = mum(h2, kPrime1);
                h2 ^= normalize_double_bits(p2.second);
                h2 = mum(h2, kPrime2);

                const auto& p3 = *it;
                ++it;
                h3 ^= pack_species_id(p3.first);
                h3 = mum(h3, kPrime1);
                h3 ^= normalize_double_bits(p3.second);
                h3 = mum(h3, kPrime2);

                remaining -= 4;
            }

            while (remaining > 0) {
                const auto& p = *it;
                ++it;
                h0 ^= pack_species_id(p.first);
                h0 = mum(h0, kPrime1);
                h0 ^= normalize_double_bits(p.second);
                h0 = mum(h0, kPrime2);
                --remaining;
            }

            return mum(h0 ^ h1 ^ h2 ^ h3, kPrime3);
        }

    private:
        static constexpr uint64_t kSeed = 0xC04D5EEDBEEFull;
        static constexpr uint64_t kPrime1 = 0xa0761d6478bd642fULL;
        static constexpr uint64_t kPrime2 = 0xe7037ed1a0b428dbULL;
        static constexpr uint64_t kPrime3 = 0x8ebc6af09c88c6e3ULL;

        // --- Helper: Fast integer mixing ---
        static inline uint64_t mum(const uint64_t a, const uint64_t b) noexcept {
            const unsigned __int128 r = static_cast<unsigned __int128>(a) * static_cast<unsigned __int128>(b);
            return static_cast<uint64_t>(r) ^ static_cast<uint64_t>(r >> 64);
        }

        static inline uint64_t mix(const uint64_t h) noexcept {
            return mum(h, kPrime1);
        }

        // --- Normalization Logic ---
        static inline uint64_t normalize_double_bits(double v) noexcept {
            if (v == 0.0) v = 0.0; // fold -0.0 -> +0.0
            if (std::isnan(v)) {
                return 0x7ff8000000000000ULL; // canonical quiet NaN
            }
            return std::bit_cast<uint64_t>(v);
        }

        static inline uint32_t pack_species_id(const auto& s) noexcept {
            const auto z = static_cast<uint16_t>(s.z());
            const auto a = static_cast<uint16_t>(s.a());
            return (static_cast<uint32_t>(z) << 16) | static_cast<uint32_t>(a);
        }
    };
}

namespace std {
    template<>
    struct hash<fourdst::composition::Composition> {
        std::size_t operator()(const fourdst::composition::Composition& c) const noexcept {
            return static_cast<std::size_t>(
                fourdst::composition::utils::CompositionHash::hash_exact(c)
            );
        }
    };
}