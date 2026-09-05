#include "ofdm/ofdm_modem.hpp"
#include <stdexcept>
#include <cmath>

namespace ntn::ofdm {

// ----------------------------------------------------------------------
// OFDM Transmitter
// ----------------------------------------------------------------------

OFDMTransmitter::OFDMTransmitter(OFDMConfig config)
    : config_(config), active_indices_(config.active_subcarrier_indices()) {
}

ComplexVector OFDMTransmitter::modulate_symbol(const ComplexVector& freq_symbols) const {
    if (freq_symbols.size() != config_.num_active_subcarriers) {
        throw std::invalid_argument("Input freq_symbols size (" + std::to_string(freq_symbols.size()) +
                                    ") does not match active subcarriers (" +
                                    std::to_string(config_.num_active_subcarriers) + ").");
    }

    // Allocate frequency grid of size N_fft and initialize to zero (null subcarriers and DC)
    ComplexVector grid(config_.fft_size, Complex(0.0, 0.0));
    for (size_t i = 0; i < config_.num_active_subcarriers; ++i) {
        grid[active_indices_[i]] = freq_symbols[i];
    }

    // Transform from frequency to time domain via IFFT
    signal::FFT::inverse(grid);

    // Prepend Cyclic Prefix: copy the last cp_length samples to the beginning
    ComplexVector time_symbol;
    time_symbol.reserve(config_.symbol_duration_samples());

    // CP samples: grid[N_fft - cp_length ... N_fft - 1]
    time_symbol.insert(time_symbol.end(), grid.end() - config_.cp_length, grid.end());
    // Useful symbol: grid[0 ... N_fft - 1]
    time_symbol.insert(time_symbol.end(), grid.begin(), grid.end());

    return time_symbol;
}

ComplexVector OFDMTransmitter::modulate_symbols(const std::vector<ComplexVector>& symbols_per_ofdm) const {
    ComplexVector continuous_stream;
    continuous_stream.reserve(symbols_per_ofdm.size() * config_.symbol_duration_samples());

    for (const auto& sym : symbols_per_ofdm) {
        ComplexVector modulated = modulate_symbol(sym);
        continuous_stream.insert(continuous_stream.end(), modulated.begin(), modulated.end());
    }
    return continuous_stream;
}

ComplexVector OFDMTransmitter::modulate_stream(const ComplexVector& all_freq_symbols) const {
    if (all_freq_symbols.size() % config_.num_active_subcarriers != 0) {
        throw std::invalid_argument("Frequency symbols size must be a multiple of active subcarriers.");
    }
    size_t num_ofdm_syms = all_freq_symbols.size() / config_.num_active_subcarriers;
    std::vector<ComplexVector> per_symbol_vec(num_ofdm_syms);

    for (size_t s = 0; s < num_ofdm_syms; ++s) {
        per_symbol_vec[s].assign(
            all_freq_symbols.begin() + s * config_.num_active_subcarriers,
            all_freq_symbols.begin() + (s + 1) * config_.num_active_subcarriers
        );
    }
    return modulate_symbols(per_symbol_vec);
}

// ----------------------------------------------------------------------
// OFDM Receiver
// ----------------------------------------------------------------------

OFDMReceiver::OFDMReceiver(OFDMConfig config)
    : config_(config), active_indices_(config.active_subcarrier_indices()) {
}

ComplexVector OFDMReceiver::demodulate_symbol(const ComplexVector& time_samples) const {
    if (time_samples.size() != config_.symbol_duration_samples()) {
        throw std::invalid_argument("Input time_samples size (" + std::to_string(time_samples.size()) +
                                    ") does not match OFDM symbol duration (" +
                                    std::to_string(config_.symbol_duration_samples()) + ").");
    }

    // Strip Cyclic Prefix: discard first cp_length samples
    ComplexVector time_core(time_samples.begin() + config_.cp_length, time_samples.end());

    // Transform from time to frequency domain via FFT
    signal::FFT::forward(time_core);

    // Extract active subcarriers
    ComplexVector freq_symbols;
    freq_symbols.reserve(config_.num_active_subcarriers);
    for (size_t idx : active_indices_) {
        freq_symbols.push_back(time_core[idx]);
    }

    return freq_symbols;
}

std::vector<ComplexVector> OFDMReceiver::demodulate_symbols(const ComplexVector& time_stream, size_t num_ofdm_symbols) const {
    const size_t sym_len = config_.symbol_duration_samples();
    if (time_stream.size() < num_ofdm_symbols * sym_len) {
        throw std::invalid_argument("Time stream length is insufficient for requested OFDM symbols.");
    }

    std::vector<ComplexVector> result;
    result.reserve(num_ofdm_symbols);

    for (size_t s = 0; s < num_ofdm_symbols; ++s) {
        ComplexVector sym_time(
            time_stream.begin() + s * sym_len,
            time_stream.begin() + (s + 1) * sym_len
        );
        result.push_back(demodulate_symbol(sym_time));
    }
    return result;
}

ComplexVector OFDMReceiver::demodulate_stream(const ComplexVector& time_stream, size_t num_ofdm_symbols) const {
    auto per_symbol = demodulate_symbols(time_stream, num_ofdm_symbols);
    ComplexVector flat_stream;
    flat_stream.reserve(num_ofdm_symbols * config_.num_active_subcarriers);

    for (const auto& sym : per_symbol) {
        flat_stream.insert(flat_stream.end(), sym.begin(), sym.end());
    }
    return flat_stream;
}

// ----------------------------------------------------------------------
// High-Level OFDM Transceiver
// ----------------------------------------------------------------------

OFDMTransceiver::OFDMTransceiver(OFDMConfig config)
    : config_(config), tx_(config), rx_(config) {
}

ComplexVector OFDMTransceiver::transmit_bits(const ByteVector& bits, const modem::Modulator& mod) const {
    size_t bits_per_ofdm = bits_per_ofdm_symbol(mod);

    if (bits.size() % bits_per_ofdm != 0) {
        throw std::invalid_argument("Total bits must be an integer multiple of bits per OFDM symbol (" +
                                    std::to_string(bits_per_ofdm) + ").");
    }

    // 1. Modulate bits to QPSK / 16-QAM constellation symbols
    ComplexVector constell_symbols = mod.modulate(bits);

    // 2. Modulate constellation symbols onto OFDM subcarriers
    return tx_.modulate_stream(constell_symbols);
}

ByteVector OFDMTransceiver::receive_bits(const ComplexVector& time_signal, size_t num_ofdm_symbols, const modem::Modulator& mod) const {
    // 1. Demodulate OFDM time signal to frequency constellation symbols
    ComplexVector recovered_symbols = rx_.demodulate_stream(time_signal, num_ofdm_symbols);

    // 2. Slicing/demodulating constellation symbols to bits
    return mod.demodulate(recovered_symbols);
}

} // namespace ntn::ofdm
