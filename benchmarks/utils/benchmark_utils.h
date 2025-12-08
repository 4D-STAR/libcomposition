#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>
#include <string>
#include <cmath>
#include <format>
#include <chrono>


template <class T>
void do_not_optimize(T&& datum) {
    asm volatile("" : "+r" (datum));
}

inline uint32_t calc_num_bins(const std::vector<double>& data) {
    const size_t n = data.size();
    return static_cast<uint32_t>(std::ceil(std::log2(n) + 1));
}

inline std::string plot_ascii_histogram(std::vector<double> data, std::string title) {
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

template <typename Func>
auto fdst_benchmark_function(Func&& func_call) {
    auto start = std::chrono::high_resolution_clock::now();

    // Forward the callable
    std::forward<Func>(func_call)();

    auto end = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

    do_not_optimize(duration.count());

    return duration;
}