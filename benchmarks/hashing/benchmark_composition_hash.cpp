#include "fourdst/composition/composition.h"
#include "fourdst/composition/utils/composition_hash.h"
#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/atomic/species.h"

#include <numeric>
#include <print>
#include <vector>
#include <ranges>
#include <chrono>

#include "benchmark_utils.h"


std::chrono::duration<double, std::nano> build_and_hash_compositions(const size_t iter, const size_t nSpecies = 8) {
    using namespace fourdst::composition;
    using namespace fourdst::atomic;

    Composition comp;
    size_t count = 0;
    for (const auto& sp :  species | std::views::values) {
        if (count >= nSpecies) {
            break;
        }
        comp.registerSpecies(sp);
        comp.setMolarAbundance(sp, 0.1);
        count++;
    }

    const auto duration = fdst_benchmark_function([&]() {
        for (size_t i = 0; i < iter; ++i) {
            uint64_t hashValue = utils::CompositionHash::hash_exact(comp);
            do_not_optimize(hashValue);
        }
    });

    return duration / static_cast<double>(iter);
}

int main() {
    using namespace fourdst::composition;
    using namespace fourdst::atomic;


    const size_t nIterations = 1000;
    std::vector<double> durations;
    durations.resize(nIterations);
    for (size_t i = 0; i < nIterations; ++i) {
        std::print("Iteration {}/{}\r", i + 1, nIterations);
        auto duration = build_and_hash_compositions(1000, 100);
        durations[i] = duration.count();
    }
    std::println("");

    std::println("Average time to build and hash composition over {} iterations: {} ns", nIterations,
                 std::accumulate(durations.begin(), durations.end(), 0.0) / nIterations);
    std::println("Max time to build and hash composition over {} iterations: {} ns", nIterations,
                 *std::ranges::max_element(durations));
    std::println("Min time to build and hash composition over {} iterations: {} ns", nIterations,
                 *std::ranges::min_element(durations));
    std::println("Standard deviation of time to build and hash composition over {} iterations: {} ns", nIterations,
                 [] (const std::vector<double>& data, const double mean) {
                     double sum = 0.0;
                     for (const auto& d : data) {
                         sum += (d - mean) * (d - mean);
                     }
                     return std::sqrt(sum / data.size());
                 } (durations, std::accumulate(durations.begin(), durations.end(), 0.0) / nIterations)
    );
    std::println("Index of max time: {}", std::distance(durations.begin(),
                 std::ranges::max_element(durations)));
    std::println("Index of min time: {}", std::distance(durations.begin(),
                 std::ranges::min_element(durations)));

    std::vector<double> log_duration = durations;
    std::ranges::transform(log_duration, log_duration.begin(), [](const double d) {
        return std::log10(d);
    });

    std::vector<double> filtered_durations;
    const double mean = std::accumulate(durations.begin(), durations.end(), 0.0) / durations.size();
    const double stddev = [] (const std::vector<double>& data, const double mean) {
        double sum = 0.0;
        for (const auto& d : data) {
            sum += (d - mean) * (d - mean);
        }
        return std::sqrt(sum / data.size());
    } (durations, mean);

    for (const auto& d : durations) {
        if (std::abs(d - mean) <= 3 * stddev) {
            filtered_durations.push_back(d);
        }
    }
    std::println("{}", plot_ascii_histogram(filtered_durations, "Build and Hash Composition Times (ns)"));
}