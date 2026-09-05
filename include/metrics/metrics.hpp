#pragma once

#include "common/types.hpp"
#include "modem/modulation.hpp"
#include <cstddef>

namespace ntn::metrics {

struct BERResult {
    size_t total_bits{0};
    size_t bit_errors{0};
    double ber{0.0};
};

struct SERResult {
    size_t total_symbols{0};
    size_t symbol_errors{0};
    double ser{0.0};
};

// Computes Bit Error Rate by comparing transmitted and received bitstreams
[[nodiscard]] BERResult calculate_ber(const ByteVector& tx_bits, const ByteVector& rx_bits);

// Computes Symbol Error Rate by comparing transmitted and sliced/received symbols
[[nodiscard]] SERResult calculate_ser(const ComplexVector& tx_symbols, const ComplexVector& rx_symbols, const modem::Modulator& mod);

// Mathematical Q-function: Q(x) = 0.5 * erfc(x / sqrt(2))
[[nodiscard]] double q_function(double x);

// Theoretical AWGN BER for Gray-coded BPSK: Q(sqrt(2 * Eb/N0))
[[nodiscard]] double theoretical_ber_bpsk(double eb_n0_db);

// Theoretical AWGN BER for Gray-coded QPSK: Q(sqrt(2 * Eb/N0)) (identical to BPSK)
[[nodiscard]] double theoretical_ber_qpsk(double eb_n0_db);

// Theoretical AWGN BER for Gray-coded 16-QAM approximation
[[nodiscard]] double theoretical_ber_16qam(double eb_n0_db);

} // namespace ntn::metrics
