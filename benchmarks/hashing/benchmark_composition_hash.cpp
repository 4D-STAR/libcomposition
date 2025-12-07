#include "fourdst/composition/composition.h"
#include "fourdst/composition/utils/composition_hash.h"
#include "fourdst/composition/utils.h"
#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/atomic/species.h"

#include <chrono>
#include <numeric>
#include <print>
#include <string>
#include <vector>
#include <cstdint>
#include <ranges>

template <class T>
void do_not_optimize(T&& datum) {
    asm volatile("" : "+r" (datum));
}

uint32_t calc_num_bins(const std::vector<double>& data) {
    // Use Sturges' formula
    const size_t n = data.size();
    return static_cast<uint32_t>(std::ceil(std::log2(n) + 1));
}

std::string plot_ascii_histogram(std::vector<double> data, std::string title) {
    // Use std::format
    const uint32_t nBins = calc_num_bins(data);
    const double minVal = *std::ranges::min_element(data);
    const double maxVal = *std::ranges::max_element(data);

    std::string histogram;
    histogram += std::format("{:^60}\n", title);
    histogram += std::string(60, '=') + "\n";
    std::vector<uint32_t> bins(nBins, 0);
    const double binWidth = (maxVal - minVal) / nBins;
    for (const auto& value : data) {
        const uint32_t binIndex = static_cast<uint32_t>((value - minVal) / binWidth);
        if (binIndex < nBins) {
            bins[binIndex]++;
        } else {
            bins[nBins - 1]++;
        }
    }
    const uint32_t maxBinCount = *std::ranges::max_element(bins);
    for (uint32_t i = 0; i < nBins; ++i) {
        const double binStart = minVal + i * binWidth;
        const double binEnd = binStart + binWidth;
        const uint32_t barLength = static_cast<uint32_t>(std::round((static_cast<double>(bins[i]) / maxBinCount) * 50.0));
        histogram += std::format("[{:.2e}, {:.2e}): {:>15} | {:}\n",
                                 binStart, binEnd, bins[i], std::string(barLength, '*'));
    }
    return histogram;

}

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

    const auto start = std::chrono::high_resolution_clock::now();
    for (size_t i = 0; i < iter; ++i) {
        uint64_t hashValue = utils::CompositionHash::hash_exact(comp);
        do_not_optimize(hashValue);
    }
    const auto end = std::chrono::high_resolution_clock::now();

    const std::chrono::duration<double, std::nano> duration = (end - start)/iter;
    return duration;
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