#include "compensation/ici_mitigation.hpp"
#include "signal/fft.hpp"
#include "signal/dsp_utils.hpp"
#include <cmath>
#include <numbers>
#include <algorithm>
#include <stdexcept>

namespace ntn::compensation {

Complex IciKernel::leakage_coefficient(int k, int m, size_t N, double epsilon) {
    int d = m - k;
    double delta = static_cast<double>(d) + epsilon;

    if (std::abs(delta) < 1e-12) {
        return Complex(1.0, 0.0);
    }

    double pi = std::numbers::pi;
    double num = std::sin(pi * delta);
    double den = static_cast<double>(N) * std::sin(pi * delta / static_cast<double>(N));

    if (std::abs(den) < 1e-14) {
        return Complex(1.0, 0.0);
    }

    double magnitude = num / den;
    double phase = pi * delta * (static_cast<double>(N - 1)) / static_cast<double>(N);

    return std::polar(magnitude, phase);
}

IciMetrics IciKernel::measure_ici(const ofdm::OFDMConfig& config,
                                  int test_subcarrier,
                                  double epsilon) {
    size_t N = config.fft_size;
    ComplexVector tx_freq(N, Complex(0.0, 0.0));
    
    // Transmit unit impulse on test subcarrier
    int idx = (test_subcarrier % static_cast<int>(N) + static_cast<int>(N)) % static_cast<int>(N);
    tx_freq[idx] = Complex(1.0, 0.0);

    // Modulate to time-domain via IFFT
    ComplexVector tx_time = tx_freq;
    signal::FFT::inverse(tx_time);

    // Apply normalized frequency offset in time domain
    // phase = 2 * pi * epsilon * n / N
    ComplexVector rx_time(N);
    double phase_step = 2.0 * std::numbers::pi * epsilon / static_cast<double>(N);
    for (size_t n = 0; n < N; ++n) {
        double phase = phase_step * static_cast<double>(n);
        rx_time[n] = tx_time[n] * std::polar(1.0, phase);
    }

    // Demodulate back to frequency domain via FFT
    ComplexVector rx_freq = rx_time;
    signal::FFT::forward(rx_freq);

    // Compute powers
    double p_carrier = std::norm(rx_freq[idx]);
    double p_ici = 0.0;
    for (size_t k = 0; k < N; ++k) {
        if (static_cast<int>(k) != idx) {
            p_ici += std::norm(rx_freq[k]);
        }
    }

    IciMetrics metrics;
    metrics.carrier_power = p_carrier;
    metrics.ici_power = p_ici;
    metrics.cir_linear = (p_ici > 1e-18) ? (p_carrier / p_ici) : 1e12;
    metrics.cir_db = 10.0 * std::log10(metrics.cir_linear);

    // Spectrum around test subcarrier (+/- 8 subcarriers)
    int span = 8;
    metrics.leakage_spectrum.resize(2 * span + 1);
    for (int offset = -span; offset <= span; ++offset) {
        int target = (idx + offset + static_cast<int>(N)) % static_cast<int>(N);
        metrics.leakage_spectrum[offset + span] = std::norm(rx_freq[target]);
    }

    return metrics;
}

IciEqualizer::IciEqualizer(const ofdm::OFDMConfig& config)
    : config_(config) {}

ComplexVector IciEqualizer::equalize(const ComplexVector& rx_subcarriers, double estimated_epsilon) const {
    size_t N = rx_subcarriers.size();
    if (N != config_.fft_size) {
        throw std::invalid_argument("IciEqualizer: rx_subcarriers size mismatch with fft_size");
    }

    // If offset is negligible, return original subcarriers
    if (std::abs(estimated_epsilon) < 1e-6) {
        return rx_subcarriers;
    }

    ComplexVector equalized = rx_subcarriers;

    // Tridiagonal Banded Matrix Equalizer (Thomas Algorithm)
    // S_{k,m} with m = k-1, k, k+1
    Complex b_diag = IciKernel::leakage_coefficient(0, 0, N, estimated_epsilon);
    Complex a_sub  = IciKernel::leakage_coefficient(1, 0, N, estimated_epsilon); // d = -1
    Complex c_sup  = IciKernel::leakage_coefficient(0, 1, N, estimated_epsilon); // d = +1

    auto active_indices = config_.active_subcarrier_indices();
    if (active_indices.empty()) {
        return equalized;
    }

    // Group contiguous bands (handling negative and positive frequency clusters)
    std::vector<std::vector<size_t>> bands;
    std::vector<size_t> current_band;
    for (size_t idx : active_indices) {
        if (current_band.empty() || idx == current_band.back() + 1) {
            current_band.push_back(idx);
        } else {
            bands.push_back(current_band);
            current_band = {idx};
        }
    }
    if (!current_band.empty()) {
        bands.push_back(current_band);
    }

    // Equalize each contiguous band independently
    for (const auto& band : bands) {
        size_t M = band.size();
        if (M == 0) continue;
        if (M == 1) {
            equalized[band[0]] /= b_diag;
            continue;
        }

        std::vector<Complex> c_prime(M, Complex(0.0, 0.0));
        std::vector<Complex> d_prime(M, Complex(0.0, 0.0));

        // Forward elimination
        c_prime[0] = c_sup / b_diag;
        d_prime[0] = rx_subcarriers[band[0]] / b_diag;

        for (size_t i = 1; i < M; ++i) {
            Complex denom = b_diag - a_sub * c_prime[i - 1];
            if (i < M - 1) {
                c_prime[i] = c_sup / denom;
            }
            d_prime[i] = (rx_subcarriers[band[i]] - a_sub * d_prime[i - 1]) / denom;
        }

        // Back substitution
        equalized[band[M - 1]] = d_prime[M - 1];
        for (int i = static_cast<int>(M) - 2; i >= 0; --i) {
            equalized[band[i]] = d_prime[i] - c_prime[i] * equalized[band[i + 1]];
        }
    }

    return equalized;
}

ComplexVector IciEqualizer::compensate_cpe(const ComplexVector& rx_subcarriers,
                                           const std::vector<int>& pilot_indices,
                                           const ComplexVector& pilot_symbols) const {
    if (pilot_indices.empty() || pilot_indices.size() != pilot_symbols.size()) {
        return rx_subcarriers;
    }

    // Estimate common phase rotation: arg(sum(rx * conj(pilot)))
    Complex sum_cross(0.0, 0.0);
    for (size_t i = 0; i < pilot_indices.size(); ++i) {
        int idx = pilot_indices[i];
        if (idx >= 0 && idx < static_cast<int>(rx_subcarriers.size())) {
            sum_cross += rx_subcarriers[idx] * std::conj(pilot_symbols[i]);
        }
    }

    double cpe_angle = std::arg(sum_cross);
    Complex de_rotation = std::polar(1.0, -cpe_angle);

    ComplexVector compensated = rx_subcarriers;
    for (auto& s : compensated) {
        s *= de_rotation;
    }

    return compensated;
}

} // namespace ntn::compensation
