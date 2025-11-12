#pragma once

#include <cstring>
#include <cmath>
#include <vector>
#include <bit>

#include "xxhash64.h"

namespace fourdst::composition::utils {
    struct CompositionHash {
        static constexpr uint64_t kSeed = 0xC04D5EEDBEEFull;
        static constexpr char kTag[] = "4DSTAR:Composition";

        template <typename CompositionT>
        static uint64_t hash_exact(const CompositionT& comp) {
            std::vector<std::uint8_t> buf;
            reserve_bytes(comp, buf);
            write_header(comp, buf);

            for (auto it = comp.begin(); it != comp.end(); ++it) {
                const auto& species = it->first;
                const double abundance = it->second;

                const std::uint32_t spWord = pack_species(species);
                push_le32(buf, spWord);

                const std::uint64_t bits = normalize_double_bits(abundance);
                push_le64(buf, bits);
            }

            return XXHash64::hash(buf.data(), buf.size(), kSeed);
        }

        static inline bool is_finite(double v) noexcept {
            return std::isfinite(v);
        }

        static inline std::int64_t quantize_index(double v, double eps) noexcept {
            const auto ld_v = static_cast<long double>(v);
            const auto ld_eps = static_cast<long double>(eps);

            const long double scaled = ld_v / ld_eps;
            const long long idx = std::llroundl(scaled);
            return static_cast<std::int64_t>(idx);
        }

        template <typename CompositionT>
        static uint64_t hash_quantized(const CompositionT& comp, double eps) noexcept {
            std::vector<std::uint8_t> buf;
            reserve_bytes(comp, buf);
            write_header(comp, buf);
            push_bytes(buf, reinterpret_cast<const std::uint8_t*>("quantized"), 9);
            push_le64(buf, encode_fp64(eps));

            for (auto it = comp.begin(); it != comp.end(); ++it) {
                const auto& species = it->first;
                const double abundance = it->second;

                const std::uint32_t spWord = pack_species(species);
                push_le32(buf, spWord);

                if (!is_finite(abundance) || eps <= 0.0) {
                    const std::uint64_t bits = normalize_double_bits(abundance);
                    push_le64(buf, bits);
                } else {
                    const std::int64_t idx = quantize_index(abundance, eps);
                    push_le64(buf, static_cast<std::uint64_t>(idx));
                }
            }

            return XXHash64::hash(buf.data(), buf.size(), kSeed ^ 0x7319'BEEF'1234ull);
        }


    private:
        template <typename SpeciesT>
        static std::uint32_t pack_species(const SpeciesT& s) noexcept {
            // Adjust accessors if your Species API differs.
            const auto z = static_cast<std::uint16_t>(s.z());
            const auto a = static_cast<std::uint16_t>(s.a());
            return (static_cast<std::uint32_t>(z) << 16) | static_cast<std::uint32_t>(a);
        }

        static inline std::uint64_t normalize_double_bits(double v) noexcept {
            if (v == 0.0) v = 0.0; // fold -0.0 -> +0.0
            if (std::isnan(v)) {
                return 0x7ff8000000000000ULL; // canonical quiet NaN
            }
            return std::bit_cast<std::uint64_t>(v);
        }

        static inline double quantize(double v, double eps) noexcept {
            if (!std::isfinite(v) || eps <= 0.0) return v;
            const double q = std::nearbyint(v / eps) * eps;
            return (q == 0.0) ? 0.0 : q;
        }

        static inline std::uint64_t encode_fp64(double v) noexcept {
            return std::bit_cast<std::uint64_t>(v);
        }

        // ---------- byte helpers (explicit little-endian) ----------
        static inline void push_le32(std::vector<std::uint8_t>& b, std::uint32_t x) {
            b.push_back(static_cast<std::uint8_t>( x        & 0xFF));
            b.push_back(static_cast<std::uint8_t>((x >> 8 ) & 0xFF));
            b.push_back(static_cast<std::uint8_t>((x >> 16) & 0xFF));
            b.push_back(static_cast<std::uint8_t>((x >> 24) & 0xFF));
        }

        static inline void push_le64(std::vector<std::uint8_t>& b, std::uint64_t x) noexcept {
            b.push_back(static_cast<std::uint8_t>( x        & 0xFF));
            b.push_back(static_cast<std::uint8_t>((x >> 8 ) & 0xFF));
            b.push_back(static_cast<std::uint8_t>((x >> 16) & 0xFF));
            b.push_back(static_cast<std::uint8_t>((x >> 24) & 0xFF));
            b.push_back(static_cast<std::uint8_t>((x >> 32) & 0xFF));
            b.push_back(static_cast<std::uint8_t>((x >> 40) & 0xFF));
            b.push_back(static_cast<std::uint8_t>((x >> 48) & 0xFF));
            b.push_back(static_cast<std::uint8_t>((x >> 56) & 0xFF));
        }

        static inline void push_bytes(std::vector<std::uint8_t>& b, const std::uint8_t* p, std::size_t n) noexcept{
            b.insert(b.end(), p, p + n);
        }

        template <typename CompositionT>
        static void write_header(const CompositionT& comp, std::vector<std::uint8_t>& buf) noexcept {
            push_bytes(buf, reinterpret_cast<const std::uint8_t*>(kTag), sizeof(kTag) - 1);

            const std::size_t nRegistered = comp.getRegisteredSpecies().size();
            std::size_t nMolar = 0;
            for (auto it = comp.begin(); it != comp.end(); ++it) { ++nMolar; }

            push_le64(buf, static_cast<std::uint64_t>(nRegistered));
            push_le64(buf, static_cast<std::uint64_t>(nMolar));
        }

        template <typename CompositionT>
        static void reserve_bytes(const CompositionT& comp, std::vector<std::uint8_t>& buf) noexcept {
            std::size_t nMolar = 0;
            for (auto it = comp.begin(); it != comp.end(); ++it) { ++nMolar; }
            const std::size_t approx = (sizeof(kTag) - 1) + 16 + nMolar * (4 + 8 + 0 /*quantized flag optional*/);
            buf.reserve(approx);
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