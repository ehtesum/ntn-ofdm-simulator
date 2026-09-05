#include "metrics/metrics.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ntn::metrics {

BERResult calculate_ber(const ByteVector& tx_bits, const ByteVector& rx_bits) {
    if (tx_bits.size() != rx_bits.size()) {
        throw std::invalid_argument("Transmitted and received bit vectors must have identical length.");
    }
    BERResult result;
    result.total_bits = tx_bits.size();
    if (result.total_bits == 0) {
        return result;
    }

    size_t errors = 0;
    for (size_t i = 0; i < result.total_bits; ++i) {
        if (tx_bits[i] != rx_bits[i]) {
            errors++;
        }
    }
    result.bit_errors = errors;
    result.ber = static_cast<double>(errors) / static_cast<double>(result.total_bits);
    return result;
}

SERResult calculate_ser(const ComplexVector& tx_symbols, const ComplexVector& rx_symbols, const modem::Modulator& mod) {
    if (tx_symbols.size() != rx_symbols.size()) {
        throw std::invalid_argument("Transmitted and received symbol vectors must have identical length.");
    }
    SERResult result;
    result.total_symbols = tx_symbols.size();
    if (result.total_symbols == 0) {
        return result;
    }

    ByteVector tx_bits = mod.demodulate(tx_symbols);
    ByteVector rx_bits = mod.demodulate(rx_symbols);
    size_t bps = mod.bits_per_symbol();

    size_t symbol_errors = 0;
    for (size_t s = 0; s < result.total_symbols; ++s) {
        bool sym_err = false;
        for (size_t b = 0; b < bps; ++b) {
            if (tx_bits[s * bps + b] != rx_bits[s * bps + b]) {
                sym_err = true;
                break;
            }
        }
        if (sym_err) {
            symbol_errors++;
        }
    }

    result.symbol_errors = symbol_errors;
    result.ser = static_cast<double>(symbol_errors) / static_cast<double>(result.total_symbols);
    return result;
}

double q_function(double x) {
    return 0.5 * std::erfc(x / std::numbers::sqrt2);
}

double theoretical_ber_bpsk(double eb_n0_db) {
    double eb_n0_lin = db_to_linear(eb_n0_db);
    return q_function(std::sqrt(2.0 * eb_n0_lin));
}

double theoretical_ber_qpsk(double eb_n0_db) {
    // Gray-coded QPSK has the exact same bit error probability as BPSK per bit!
    return theoretical_ber_bpsk(eb_n0_db);
}

double theoretical_ber_16qam(double eb_n0_db) {
    double eb_n0_lin = db_to_linear(eb_n0_db);
    // Standard Gray-coded square 16-QAM BER approximation:
    // P_b ~ (3/8) * erfc(sqrt( (2/5) * Eb/N0 )) = (3/4) * Q(sqrt( (4/5) * Eb/N0 ))
    return 0.75 * q_function(std::sqrt(0.8 * eb_n0_lin));
}

} // namespace ntn::metrics
