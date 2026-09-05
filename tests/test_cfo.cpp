#include "common/types.hpp"
#include "ofdm/ofdm_config.hpp"
#include "ofdm/ofdm_modem.hpp"
#include "modem/modulation.hpp"
#include "metrics/metrics.hpp"
#include "signal/dsp_utils.hpp"
#include "synchronization/cfo.hpp"
#include "channel/noise.hpp"

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
    std::cout << " NTN-OFDM SIMULATOR - CFO UNIT TESTS     \n";
    std::cout << "=========================================\n\n";

    ntn::ofdm::OFDMConfig cfg;
    cfg.fft_size = 512;
    cfg.num_active_subcarriers = 300;
    cfg.cp_length = 36;
    cfg.subcarrier_spacing_hz = 15000.0; // 15 kHz
    const double fs = cfg.sampling_rate_hz(); // 7.68 MHz

    // ----------------------------------------------------
    // Test Suite 1: CFO Phase Rotation & Continuity
    // ----------------------------------------------------
    std::cout << "[Suite 1: CFO Channel Phase Physics]\n";
    {
        double f_cfo = 3000.0; // +3 kHz offset
        ntn::synchronization::CFOChannel cfo_chan(f_cfo, fs);

        // Constant amplitude DC signal
        ntn::ComplexVector sig1(100, ntn::Complex(1.0, 0.0));
        cfo_chan.apply_cfo(sig1);

        // Check phase of sample 50
        double expected_phase_50 = 2.0 * ntn::PI * f_cfo * 50.0 / fs;
        double actual_phase_50 = std::arg(sig1[50]);
        run_test("CFO phase rotation matches analytical formula at sample 50",
                 std::abs(actual_phase_50 - expected_phase_50) < 1e-9);

        // Check phase continuity across contiguous blocks
        ntn::ComplexVector sig2(100, ntn::Complex(1.0, 0.0));
        cfo_chan.apply_cfo(sig2);
        double expected_phase_100 = 2.0 * ntn::PI * f_cfo * 100.0 / fs;
        // sig2[0] corresponds to absolute sample index 100
        double actual_phase_100 = std::arg(sig2[0]);
        run_test("CFO phase accumulator preserves seamless continuity across blocks",
                 std::abs(actual_phase_100 - expected_phase_100) < 1e-9);
    }

    // ----------------------------------------------------
    // Test Suite 2: CP-Correlation Estimator Accuracy (Noiseless)
    // ----------------------------------------------------
    std::cout << "\n[Suite 2: CP-Correlation CFO Estimator Accuracy (Noiseless)]\n";
    {
        ntn::ofdm::OFDMTransmitter tx(cfg);
        ntn::synchronization::CFOEstimator estimator(cfg);

        // Random QPSK symbols
        std::mt19937_64 rng(42);
        std::uniform_int_distribution<int> bit_dist(0, 1);
        ntn::ByteVector bits(cfg.num_active_subcarriers * 2);
        for (auto& b : bits) b = static_cast<uint8_t>(bit_dist(rng));
        ntn::modem::QPSKModulator qpsk;
        auto freq_syms = qpsk.modulate(bits);
        auto time_sym = tx.modulate_symbol(freq_syms);

        // Test various normalized CFO values: epsilon = f_cfo / delta_f
        for (double eps : {0.0, 0.05, 0.15, 0.25, 0.40, -0.10, -0.30}) {
            double test_f_cfo = eps * cfg.subcarrier_spacing_hz;
            ntn::synchronization::CFOChannel cfo_impairment(test_f_cfo, fs);

            ntn::ComplexVector impaired = time_sym;
            cfo_impairment.apply_cfo(impaired);

            double estimated_eps = estimator.estimate_normalized_cfo(impaired);
            double estimated_cfo = estimator.estimate_cfo_hz(impaired);

            double err_eps = std::abs(estimated_eps - eps);
            double err_hz = std::abs(estimated_cfo - test_f_cfo);

            std::string label = "Estimator epsilon=" + std::to_string(eps) + " (Error < 1e-9)";
            run_test(label, err_eps < 1e-9 && err_hz < 1e-6, "err_hz=" + std::to_string(err_hz));
        }
    }

    // ----------------------------------------------------
    // Test Suite 3: Estimator Robustness Under AWGN Noise
    // ----------------------------------------------------
    std::cout << "\n[Suite 3: Multi-Symbol CFO Estimation Under AWGN]\n";
    {
        ntn::ofdm::OFDMTransmitter tx(cfg);
        ntn::synchronization::CFOEstimator estimator(cfg);
        ntn::channel::AWGNChannel awgn(777);

        const size_t num_symbols = 10;
        ntn::modem::QPSKModulator qpsk;
        std::vector<ntn::ComplexVector> frame_syms(num_symbols);
        std::mt19937_64 rng(101);
        std::uniform_int_distribution<int> bit_dist(0, 1);

        for (size_t s = 0; s < num_symbols; ++s) {
            ntn::ByteVector bits(cfg.num_active_subcarriers * 2);
            for (auto& b : bits) b = static_cast<uint8_t>(bit_dist(rng));
            frame_syms[s] = qpsk.modulate(bits);
        }

        auto clean_stream = tx.modulate_symbols(frame_syms);

        double target_cfo_hz = 2250.0; // epsilon = 0.15 on 15 kHz SCS
        ntn::synchronization::CFOChannel cfo_channel(target_cfo_hz, fs);
        auto impaired_stream = clean_stream;
        cfo_channel.apply_cfo(impaired_stream);

        // Add 12 dB AWGN
        awgn.add_noise_snr(impaired_stream, 12.0);

        // Multi-symbol estimation
        double estimated_cfo = estimator.estimate_cfo_hz_multisymbol(impaired_stream, num_symbols);
        double est_error = std::abs(estimated_cfo - target_cfo_hz);

        run_test("Multi-symbol CP estimator at 12 dB SNR within 35 Hz (< 0.25% error)",
                 est_error < 35.0, "error=" + std::to_string(est_error) + " Hz");
    }

    // ----------------------------------------------------
    // Test Suite 4: End-to-End Compensation & Constellation Restoration
    // ----------------------------------------------------
    std::cout << "\n[Suite 4: End-to-End CFO Compensation Performance]\n";
    {
        ntn::ofdm::OFDMTransceiver transceiver(cfg);
        ntn::modem::QPSKModulator qpsk;

        const size_t num_symbols = 14; // 1 5G NR slot
        const size_t total_bits = transceiver.bits_per_ofdm_symbol(qpsk) * num_symbols; // 8,400 bits

        std::mt19937_64 rng(888);
        std::uniform_int_distribution<int> bit_dist(0, 1);
        ntn::ByteVector tx_bits(total_bits);
        for (auto& b : tx_bits) b = static_cast<uint8_t>(bit_dist(rng));

        auto tx_time = transceiver.transmit_bits(tx_bits, qpsk);

        // Introduce moderate CFO: epsilon = 0.15 (2,250 Hz)
        const double cfo_hz = 2250.0;
        ntn::synchronization::CFOChannel cfo_impairment(cfo_hz, fs);
        auto rx_impaired = tx_time;
        cfo_impairment.apply_cfo(rx_impaired);

        // 1. Uncompensated reception
        auto rx_bits_uncomp = transceiver.receive_bits(rx_impaired, num_symbols, qpsk);
        auto ber_uncomp = ntn::metrics::calculate_ber(tx_bits, rx_bits_uncomp);

        run_test("Uncompensated CFO causes severe BER failure (BER > 0.20)",
                 ber_uncomp.ber > 0.20, "uncomp_ber=" + std::to_string(ber_uncomp.ber));

        // 2. Estimation and Compensation
        ntn::synchronization::CFOEstimator estimator(cfg);
        double est_cfo = estimator.estimate_cfo_hz_multisymbol(rx_impaired, num_symbols);

        ntn::synchronization::CFOCompensator compensator(est_cfo, fs);
        auto rx_corrected = rx_impaired;
        compensator.compensate(rx_corrected);

        auto rx_bits_comp = transceiver.receive_bits(rx_corrected, num_symbols, qpsk);
        auto ber_comp = ntn::metrics::calculate_ber(tx_bits, rx_bits_comp);

        run_test("CFO compensation completely restores reception (BER = 0.0, 0 bit errors)",
                 ber_comp.bit_errors == 0 && ber_comp.ber == 0.0);
    }

    // ----------------------------------------------------
    // Summary
    // ----------------------------------------------------
    std::cout << "\n=========================================\n";
    std::cout << " TEST RESULTS: " << passed_tests << " / " << total_tests << " passed.\n";
    std::cout << "=========================================\n\n";

    return (passed_tests == total_tests) ? 0 : 1;
}
