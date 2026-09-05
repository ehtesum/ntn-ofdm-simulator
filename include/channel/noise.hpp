#pragma once

#include "common/types.hpp"
#include <random>
#include <cstdint>

namespace ntn::channel {

class AWGNChannel {
public:
    // Construct AWGN generator with an optional deterministic seed (default: 42)
    explicit AWGNChannel(uint64_t seed = 42);

    // Set the RNG seed for deterministic reproducibility
    void set_seed(uint64_t seed);

    // Generates N complex Gaussian noise samples with total variance sigma^2 (N0)
    // Real and Imaginary components each have variance sigma^2 / 2
    [[nodiscard]] ComplexVector generate_noise(size_t num_samples, double noise_variance);

    // Generates a single complex Gaussian noise sample
    [[nodiscard]] Complex generate_sample(double noise_variance);

    // Adds AWGN noise to a signal in-place given an explicit noise variance
    void add_noise_variance(ComplexVector& signal, double noise_variance);

    // Adds AWGN noise to a signal in-place given a target SNR (in dB)
    // Automatically calculates signal power to determine the exact required noise variance
    void add_noise_snr(ComplexVector& signal, double snr_db);

private:
    std::mt19937_64 rng_;
    std::normal_distribution<double> std_normal_{0.0, 1.0};
};

} // namespace ntn::channel
