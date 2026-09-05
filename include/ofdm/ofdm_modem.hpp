#pragma once

#include "ofdm/ofdm_config.hpp"
#include "modem/modulation.hpp"
#include "signal/fft.hpp"
#include <memory>

namespace ntn::ofdm {

// OFDM Transmitter: Maps frequency-domain symbols -> IFFT -> Cyclic Prefix insertion
class OFDMTransmitter {
public:
    explicit OFDMTransmitter(OFDMConfig config);

    // Modulates a single OFDM symbol from frequency-domain subcarrier symbols.
    // Input must contain exactly config.num_active_subcarriers complex symbols.
    // Output contains config.symbol_duration_samples() time-domain samples.
    [[nodiscard]] ComplexVector modulate_symbol(const ComplexVector& freq_symbols) const;

    // Modulates multiple OFDM symbols into a continuous time-domain waveform
    [[nodiscard]] ComplexVector modulate_symbols(const std::vector<ComplexVector>& symbols_per_ofdm) const;

    // Modulates a flat stream of frequency symbols into a continuous time-domain waveform
    [[nodiscard]] ComplexVector modulate_stream(const ComplexVector& all_freq_symbols) const;

    [[nodiscard]] const OFDMConfig& config() const noexcept { return config_; }

private:
    OFDMConfig config_;
    std::vector<size_t> active_indices_;
};

// OFDM Receiver: Cyclic Prefix removal -> FFT -> Subcarrier demapping
class OFDMReceiver {
public:
    explicit OFDMReceiver(OFDMConfig config);

    // Demodulates a single OFDM symbol from time-domain samples.
    // Input must contain config.symbol_duration_samples() time-domain samples.
    // Output contains config.num_active_subcarriers frequency-domain symbols.
    [[nodiscard]] ComplexVector demodulate_symbol(const ComplexVector& time_samples) const;

    // Demodulates a continuous time-domain stream into per-OFDM symbol subcarrier vectors
    [[nodiscard]] std::vector<ComplexVector> demodulate_symbols(const ComplexVector& time_stream, size_t num_ofdm_symbols) const;

    // Demodulates a continuous time-domain stream into a flat vector of frequency symbols
    [[nodiscard]] ComplexVector demodulate_stream(const ComplexVector& time_stream, size_t num_ofdm_symbols) const;

    [[nodiscard]] const OFDMConfig& config() const noexcept { return config_; }

private:
    OFDMConfig config_;
    std::vector<size_t> active_indices_;
};

// High-level OFDM Transceiver combining Transmitter, Receiver, and Digital Modulator
class OFDMTransceiver {
public:
    explicit OFDMTransceiver(OFDMConfig config);

    // Returns total payload bits accommodated per OFDM symbol for a given modulation
    [[nodiscard]] size_t bits_per_ofdm_symbol(const modem::Modulator& mod) const noexcept {
        return config_.num_active_subcarriers * mod.bits_per_symbol();
    }

    // Transmits binary data: Bits -> QAM/QPSK Symbols -> OFDM Modulation
    [[nodiscard]] ComplexVector transmit_bits(const ByteVector& bits, const modem::Modulator& mod) const;

    // Receives binary data: Time Samples -> OFDM Demodulation -> QAM/QPSK Slicing -> Bits
    [[nodiscard]] ByteVector receive_bits(const ComplexVector& time_signal, size_t num_ofdm_symbols, const modem::Modulator& mod) const;

    [[nodiscard]] const OFDMTransmitter& transmitter() const noexcept { return tx_; }
    [[nodiscard]] const OFDMReceiver& receiver() const noexcept { return rx_; }
    [[nodiscard]] const OFDMConfig& config() const noexcept { return config_; }

private:
    OFDMConfig config_;
    OFDMTransmitter tx_;
    OFDMReceiver rx_;
};

} // namespace ntn::ofdm
