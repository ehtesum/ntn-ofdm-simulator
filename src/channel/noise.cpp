#include "channel/noise.hpp"
#include "signal/dsp_utils.hpp"
#include <cmath>

namespace ntn::channel {

AWGNChannel::AWGNChannel(uint64_t seed)
    : rng_(seed) {
}

void AWGNChannel::set_seed(uint64_t seed) {
    rng_.seed(seed);
}

Complex AWGNChannel::generate_sample(double noise_variance) {
    if (noise_variance <= 0.0) {
        return Complex(0.0, 0.0);
    }
    // For complex AWGN, total variance sigma^2 = var(I) + var(Q)
    // Therefore std_dev per dimension = sqrt(noise_variance / 2.0)
    const double std_dev = std::sqrt(noise_variance / 2.0);
    double real_part = std_normal_(rng_) * std_dev;
    double imag_part = std_normal_(rng_) * std_dev;
    return Complex(real_part, imag_part);
}

ComplexVector AWGNChannel::generate_noise(size_t num_samples, double noise_variance) {
    ComplexVector noise;
    noise.reserve(num_samples);
    if (noise_variance <= 0.0) {
        noise.assign(num_samples, Complex(0.0, 0.0));
        return noise;
    }
    const double std_dev = std::sqrt(noise_variance / 2.0);
    for (size_t i = 0; i < num_samples; ++i) {
        double real_part = std_normal_(rng_) * std_dev;
        double imag_part = std_normal_(rng_) * std_dev;
        noise.emplace_back(real_part, imag_part);
    }
    return noise;
}

void AWGNChannel::add_noise_variance(ComplexVector& signal, double noise_variance) {
    if (signal.empty() || noise_variance <= 0.0) {
        return;
    }
    const double std_dev = std::sqrt(noise_variance / 2.0);
    for (auto& sample : signal) {
        sample += Complex(std_normal_(rng_) * std_dev, std_normal_(rng_) * std_dev);
    }
}

void AWGNChannel::add_noise_snr(ComplexVector& signal, double snr_db) {
    if (signal.empty()) {
        return;
    }
    double signal_power = signal::calculate_power(signal);
    if (signal_power <= 1e-15) {
        return;
    }
    double noise_var = signal::calculate_noise_variance(signal_power, snr_db);
    add_noise_variance(signal, noise_var);
}

} // namespace ntn::channel
