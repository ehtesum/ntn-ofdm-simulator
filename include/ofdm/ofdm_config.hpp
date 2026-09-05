#pragma once

#include "common/types.hpp"
#include <cstddef>
#include <vector>
#include <stdexcept>

namespace ntn::ofdm {

struct OFDMConfig {
    size_t fft_size{512};                   // IFFT/FFT size (must be power of 2, e.g. 64, 128, 256, 512, 1024)
    size_t num_active_subcarriers{300};     // Number of data subcarriers (e.g. 300 for 25 PRBs in 5G NR)
    size_t cp_length{36};                   // Cyclic prefix length in samples (~7% of FFT size)
    bool null_dc{true};                     // Whether to null the DC subcarrier (k=0) to prevent LO leakage
    double subcarrier_spacing_hz{15.0e3};   // Subcarrier spacing delta_f (15 kHz for 3GPP NR mu=0)

    // Total samples per OFDM symbol including Cyclic Prefix
    [[nodiscard]] constexpr size_t symbol_duration_samples() const noexcept {
        return fft_size + cp_length;
    }

    // Baseband sampling rate: Fs = N_fft * delta_f
    [[nodiscard]] constexpr double sampling_rate_hz() const noexcept {
        return static_cast<double>(fft_size) * subcarrier_spacing_hz;
    }

    // Time-domain duration of one OFDM symbol including CP in seconds
    [[nodiscard]] constexpr double symbol_duration_seconds() const noexcept {
        return static_cast<double>(symbol_duration_samples()) / sampling_rate_hz();
    }

    // Validates configuration parameters
    void validate() const {
        if (fft_size < 16 || (fft_size & (fft_size - 1)) != 0) {
            throw std::invalid_argument("FFT size must be a power of 2 and >= 16.");
        }
        size_t max_active = null_dc ? (fft_size - 1) : fft_size;
        if (num_active_subcarriers == 0 || num_active_subcarriers > max_active) {
            throw std::invalid_argument("Active subcarriers must be > 0 and <= " + std::to_string(max_active));
        }
        if (cp_length == 0 || cp_length >= fft_size) {
            throw std::invalid_argument("Cyclic prefix length must be > 0 and < fft_size.");
        }
    }

    // Computes the exact FFT bin indices [0, N_fft-1] allocated to active subcarriers.
    // Centered symmetrically around DC (k=0):
    // Negative frequencies: bins [N_fft - N_sc/2, N_fft - 1]
    // Positive frequencies: bins [1, N_sc/2] (if null_dc is true)
    [[nodiscard]] std::vector<size_t> active_subcarrier_indices() const {
        validate();
        std::vector<size_t> indices;
        indices.reserve(num_active_subcarriers);

        size_t half_sc = num_active_subcarriers / 2;
        size_t neg_sc = half_sc;
        size_t pos_sc = num_active_subcarriers - half_sc;

        // Negative subcarrier indices in standard DFT ordering
        for (size_t i = 0; i < neg_sc; ++i) {
            indices.push_back(fft_size - neg_sc + i);
        }

        // Positive subcarrier indices
        size_t start_pos = null_dc ? 1 : 0;
        for (size_t i = 0; i < pos_sc; ++i) {
            indices.push_back(start_pos + i);
        }

        return indices;
    }
};

} // namespace ntn::ofdm
