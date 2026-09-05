#include "common/types.hpp"
#include "ofdm/ofdm_config.hpp"
#include "ofdm/ofdm_modem.hpp"
#include "modem/modulation.hpp"
#include "metrics/metrics.hpp"
#include "signal/dsp_utils.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include <random>

namespace {

int total_tests = 0;
int passed_tests = 0;

void run_test(const std::string& name, bool condition, const std::string& extra = "") {
    total_tests++;
    if (condition) {
        passed_tests++;
        std::cout << "  [PASS] " << name << "\n";
    } else {
        std::cout << "  [FAIL] " << name;
        if (!extra.empty()) {
            std::cout << " (" << extra << ")";
        }
        std::cout << "\n";
    }
}

} // namespace

int main() {
    std::cout << "\n=========================================\n";
    std::cout << " NTN-OFDM SIMULATOR - OFDM MODEM TESTS   \n";
    std::cout << "=========================================\n\n";

    // ----------------------------------------------------
    // Test Suite 1: OFDM Configuration & Subcarrier Grid
    // ----------------------------------------------------
    std::cout << "[Suite 1: OFDM Grid Configuration]\n";
    {
        ntn::ofdm::OFDMConfig cfg;
        cfg.fft_size = 512;
        cfg.num_active_subcarriers = 300;
        cfg.cp_length = 36;
        cfg.null_dc = true;

        auto indices = cfg.active_subcarrier_indices();
        run_test("Active subcarrier count matches configuration", indices.size() == 300);

        // Check that DC (index 0) is never allocated when null_dc is true
        bool dc_found = false;
        for (size_t idx : indices) {
            if (idx == 0) dc_found = true;
        }
        run_test("DC subcarrier (k=0) is strictly nulled", !dc_found);

        // Check sampling rate and symbol duration
        // 512 * 15 kHz = 7.68 MHz sampling rate (standard LTE/NR rate)
        run_test("Sampling rate for 512 FFT @ 15 kHz SCS is 7.68 MHz",
                 std::abs(cfg.sampling_rate_hz() - 7.68e6) < 1.0);

        // Check config validation error handling
        bool invalid_fft_caught = false;
        try {
            ntn::ofdm::OFDMConfig bad_fft;
            bad_fft.fft_size = 500; // not power of 2
            bad_fft.validate();
        } catch (const std::invalid_argument&) {
            invalid_fft_caught = true;
        }
        run_test("Invalid non-power-of-2 FFT size throws invalid_argument", invalid_fft_caught);
    }

    // ----------------------------------------------------
    // Test Suite 2: Cyclic Prefix Structure Verification
    // ----------------------------------------------------
    std::cout << "\n[Suite 2: Cyclic Prefix Structure]\n";
    {
        ntn::ofdm::OFDMConfig cfg;
        cfg.fft_size = 128;
        cfg.num_active_subcarriers = 72; // 6 PRBs
        cfg.cp_length = 16;
        cfg.null_dc = true;

        ntn::ofdm::OFDMTransmitter tx(cfg);

        // Deterministic pseudo-random symbols
        std::mt19937_64 rng(42);
        std::uniform_real_distribution<double> dist(-1.0, 1.0);
        ntn::ComplexVector freq_syms(cfg.num_active_subcarriers);
        for (auto& s : freq_syms) {
            s = ntn::Complex(dist(rng), dist(rng));
        }

        ntn::ComplexVector time_sym = tx.modulate_symbol(freq_syms);

        run_test("Total modulated samples = fft_size + cp_length",
                 time_sym.size() == cfg.symbol_duration_samples());

        // In an OFDM symbol, the first cp_length samples must be an exact copy
        // of the last cp_length samples of the useful FFT window:
        // time_sym[0 ... cp-1] == time_sym[fft_size ... fft_size + cp - 1]
        double max_cp_diff = 0.0;
        for (size_t i = 0; i < cfg.cp_length; ++i) {
            double diff = std::abs(time_sym[i] - time_sym[cfg.fft_size + i]);
            max_cp_diff = std::max(max_cp_diff, diff);
        }
        run_test("Cyclic Prefix perfectly matches end of FFT block (diff < 1e-14)",
                 max_cp_diff < 1e-14, "max_diff=" + std::to_string(max_cp_diff));
    }

    // ----------------------------------------------------
    // Test Suite 3: Single-Symbol Ideal Loopback
    // ----------------------------------------------------
    std::cout << "\n[Suite 3: Ideal Single-Symbol Loopback]\n";
    {
        ntn::ofdm::OFDMConfig cfg;
        cfg.fft_size = 256;
        cfg.num_active_subcarriers = 144;
        cfg.cp_length = 20;

        ntn::ofdm::OFDMTransmitter tx(cfg);
        ntn::ofdm::OFDMReceiver rx(cfg);

        ntn::ComplexVector tx_freq(cfg.num_active_subcarriers);
        for (size_t i = 0; i < cfg.num_active_subcarriers; ++i) {
            double angle = 2.0 * ntn::PI * static_cast<double>(i) / static_cast<double>(cfg.num_active_subcarriers);
            tx_freq[i] = ntn::Complex(std::cos(angle), std::sin(angle));
        }

        ntn::ComplexVector tx_time = tx.modulate_symbol(tx_freq);
        ntn::ComplexVector rx_freq = rx.demodulate_symbol(tx_time);

        run_test("Recovered frequency symbol vector size matches transmitted",
                 rx_freq.size() == tx_freq.size());

        double max_err = 0.0;
        for (size_t i = 0; i < tx_freq.size(); ++i) {
            max_err = std::max(max_err, std::abs(rx_freq[i] - tx_freq[i]));
        }
        run_test("Single-symbol loopback maximum error < 1e-12", max_err < 1e-12, "max_err=" + std::to_string(max_err));

        double evm = ntn::signal::calculate_evm_percent(tx_freq, rx_freq);
        run_test("Single-symbol loopback EVM is 0.0%", std::abs(evm) < 1e-10);
    }

    // ----------------------------------------------------
    // Test Suite 4: Multi-Symbol Slot Loopback (14 OFDM Symbols)
    // ----------------------------------------------------
    std::cout << "\n[Suite 4: 5G NR Slot (14 OFDM Symbols) QPSK Loopback]\n";
    {
        ntn::ofdm::OFDMConfig cfg;
        cfg.fft_size = 512;
        cfg.num_active_subcarriers = 300; // 25 PRBs
        cfg.cp_length = 36;
        cfg.null_dc = true;

        ntn::ofdm::OFDMTransceiver transceiver(cfg);
        ntn::modem::QPSKModulator qpsk;

        const size_t num_ofdm_symbols = 14; // Standard 5G NR slot duration
        const size_t total_payload_bits = transceiver.bits_per_ofdm_symbol(qpsk) * num_ofdm_symbols;
        // 300 subcarriers * 2 bits/sym * 14 symbols = 8,400 bits

        std::mt19937_64 rng(999);
        std::uniform_int_distribution<int> bit_dist(0, 1);
        ntn::ByteVector tx_bits(total_payload_bits);
        for (auto& b : tx_bits) {
            b = static_cast<uint8_t>(bit_dist(rng));
        }

        // Transmit -> Ideal Channel -> Receive
        ntn::ComplexVector time_signal = transceiver.transmit_bits(tx_bits, qpsk);

        run_test("Slot time-domain sample count equals 14 * symbol_duration",
                 time_signal.size() == num_ofdm_symbols * cfg.symbol_duration_samples());

        ntn::ByteVector rx_bits = transceiver.receive_bits(time_signal, num_ofdm_symbols, qpsk);

        run_test("Received bit count matches transmitted bit count (8,400 bits)",
                 rx_bits.size() == tx_bits.size());

        auto ber = ntn::metrics::calculate_ber(tx_bits, rx_bits);
        run_test("Ideal 14-symbol slot transmission has zero bit errors (BER = 0.0)",
                 ber.bit_errors == 0 && ber.ber == 0.0);
    }

    // ----------------------------------------------------
    // Test Suite 5: 16-QAM High-Spectral-Efficiency Loopback
    // ----------------------------------------------------
    std::cout << "\n[Suite 5: 16-QAM Multicarrier Loopback]\n";
    {
        ntn::ofdm::OFDMConfig cfg;
        cfg.fft_size = 256;
        cfg.num_active_subcarriers = 120; // 10 PRBs
        cfg.cp_length = 18;

        ntn::ofdm::OFDMTransceiver transceiver(cfg);
        ntn::modem::QAM16Modulator qam16;

        const size_t num_ofdm_symbols = 10;
        const size_t total_bits = transceiver.bits_per_ofdm_symbol(qam16) * num_ofdm_symbols;
        // 120 subcarriers * 4 bits/sym * 10 = 4,800 bits

        std::mt19937_64 rng(555);
        std::uniform_int_distribution<int> bit_dist(0, 1);
        ntn::ByteVector tx_bits(total_bits);
        for (auto& b : tx_bits) {
            b = static_cast<uint8_t>(bit_dist(rng));
        }

        ntn::ComplexVector time_signal = transceiver.transmit_bits(tx_bits, qam16);
        ntn::ByteVector rx_bits = transceiver.receive_bits(time_signal, num_ofdm_symbols, qam16);

        auto ber = ntn::metrics::calculate_ber(tx_bits, rx_bits);
        run_test("16-QAM multi-symbol loopback has zero bit errors (4,800 bits)",
                 ber.bit_errors == 0 && ber.ber == 0.0);
    }

    // ----------------------------------------------------
    // Test Suite 6: Numerology Scalability (64, 128, 512, 1024)
    // ----------------------------------------------------
    std::cout << "\n[Suite 6: Numerology Scalability across FFT Sizes]\n";
    {
        for (size_t n_fft : {64, 128, 512, 1024}) {
            ntn::ofdm::OFDMConfig cfg;
            cfg.fft_size = n_fft;
            cfg.num_active_subcarriers = n_fft / 2;
            cfg.cp_length = n_fft / 8;

            ntn::ofdm::OFDMTransmitter tx(cfg);
            ntn::ofdm::OFDMReceiver rx(cfg);

            ntn::ComplexVector input_syms(cfg.num_active_subcarriers, ntn::Complex(0.7071, -0.7071));
            auto time_out = tx.modulate_symbol(input_syms);
            auto freq_out = rx.demodulate_symbol(time_out);

            double err = 0.0;
            for (size_t k = 0; k < cfg.num_active_subcarriers; ++k) {
                err = std::max(err, std::abs(freq_out[k] - input_syms[k]));
            }

            std::string label = "Numerology N_fft=" + std::to_string(n_fft) + " loopback max error < 1e-12";
            run_test(label, err < 1e-12, "err=" + std::to_string(err));
        }
    }

    // ----------------------------------------------------
    // Summary
    // ----------------------------------------------------
    std::cout << "\n=========================================\n";
    std::cout << " TEST RESULTS: " << passed_tests << " / " << total_tests << " passed.\n";
    std::cout << "=========================================\n\n";

    return (passed_tests == total_tests) ? 0 : 1;
}
