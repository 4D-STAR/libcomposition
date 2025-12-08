#include "benchmark_utils.h"

#include "fourdst/composition/composition.h"
#include "fourdst/atomic/species.h"

#include <chrono>
#include <random>
#include <ranges>

std::chrono::duration<double, std::nano> benchmark_construction(const size_t iterations, const size_t nSpecies) {
    using namespace fourdst::composition;
    using namespace fourdst::atomic;

    // Setup random machine to get random double between 0 and 1 for molar abundances
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<Species> species_to_register;
    std::vector<double> molarAbundances;

    size_t count = 0;
    for (const auto& sp : species | std::views::values) {
        if (count >= nSpecies) {
            break;
        }
        species_to_register.push_back(sp);
        molarAbundances.push_back(dis(gen));
        count++;
    }

    const auto duration = fdst_benchmark_function([&]() {
        for (size_t i = 0; i < iterations; ++i) {
            fourdst::composition::Composition comp(species_to_register, molarAbundances);
        }
    });

    return duration / static_cast<double>(iterations);
}

std::chrono::duration<double, std::nano> benchmark_access(const size_t iterations, const size_t nSpecies) {
    using namespace fourdst::composition;
    using namespace fourdst::atomic;

    // Setup random machine to get random double between 0 and 1 for molar abundances
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0.0, 1.0);

    std::vector<Species> species_to_register;
    std::vector<double> molarAbundances;

    size_t count = 0;
    for (const auto& sp : species | std::views::values) {
        if (count >= nSpecies) {
            break;
        }
        species_to_register.push_back(sp);
        molarAbundances.push_back(dis(gen));
        count++;
    }

    const Composition comp(species_to_register, molarAbundances);

    std::uniform_int_distribution<>(0, nSpecies - 1);
    std::vector<Species> random_lookup_species;
    for (size_t i = 0; i < iterations; ++i) {
        random_lookup_species.push_back(species_to_register[static_cast<size_t>(dis(gen))]);
    }

    const auto duration = fdst_benchmark_function([&]() {
        for (size_t i = 0; i < iterations; ++i) {
            volatile double y = comp.getMolarAbundance(random_lookup_species[i]);
            do_not_optimize(y);
        }
    });


    return duration / static_cast<double>(iterations);
}

int main () {
    constexpr size_t nIterations = 1000;
    constexpr size_t nSpecies = 100;

    std::vector<double> durations;
    durations.resize(nIterations);

    for (size_t i = 0; i < nIterations; ++i) {
        std::print("Iteration {}/{}\r", i + 1, nIterations);
        auto duration = benchmark_construction(10, nSpecies);
        durations[i] = duration.count();
    }
    std::println("");

    std::println("Average time to construct composition over {} iterations: {} ns", nIterations,
                 std::accumulate(durations.begin(), durations.end(), 0.0) / nIterations);
    std::println("Max time to construct composition over {} iterations: {} ns", nIterations,
                 *std::ranges::max_element(durations));
    std::println("Min time to construct composition over {} iterations: {} ns", nIterations,
                 *std::ranges::min_element(durations));


    plot_ascii_histogram(durations, "Composition Construction Time Histogram");


    durations.clear();
    durations.resize(nIterations);
    for (size_t i = 0; i < nIterations; ++i) {
        std::print("Iteration {}/{}\r", i + 1, nIterations);
        auto duration = benchmark_access(1000, nSpecies);
        durations[i] = duration.count();
    }
    std::println("");
    std::println("Average time to access composition over {} iterations: {} ns", nIterations,
                 std::accumulate(durations.begin(), durations.end(), 0.0) / nIterations);
    std::println("Max time to access composition over {} iterations: {} ns", nIterations,
                 *std::ranges::max_element(durations));
    std::println("Min time to access composition over {} iterations: {} ns", nIterations,
                 *std::ranges::min_element(durations));

    plot_ascii_histogram(durations, "Composition Access Time Histogram");
}