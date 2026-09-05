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
#include <fstream>
#include <vector>
#include <filesystem>
#include <random>

namespace fs = std::filesystem;

int main() {
    fs::create_directories("results/tables");
    fs::create_directories("results/figures");

    std::cout << "====================================================================\n";
    std::cout << " NTN-OFDM SIMULATOR — CARRIER FREQUENCY OFFSET (CFO) EXPERIMENT\n";
    std::cout << "====================================================================\n";

    ntn::ofdm::OFDMConfig cfg;
    cfg.fft_size = 512;
    cfg.num_active_subcarriers = 300;
    cfg.cp_length = 36;
    cfg.subcarrier_spacing_hz = 15000.0;
    const double fs_hz = cfg.sampling_rate_hz();

    ntn::ofdm::OFDMTransceiver transceiver(cfg);
    ntn::modem::QPSKModulator qpsk;
    ntn::synchronization::CFOEstimator estimator(cfg);

    const size_t num_symbols = 50; // 50 OFDM symbols * 300 sc * 2 bits = 30,000 bits
    const size_t total_bits = transceiver.bits_per_ofdm_symbol(qpsk) * num_symbols;

    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int> bit_dist(0, 1);
    ntn::ByteVector tx_bits(total_bits);
    for (auto& b : tx_bits) b = static_cast<uint8_t>(bit_dist(rng));

    auto tx_time = transceiver.transmit_bits(tx_bits, qpsk);

    std::vector<double> eps_values = {0.0, 0.02, 0.05, 0.08, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40};

    std::ofstream csv("results/tables/cfo_comparison.csv");
    csv << "Normalized_CFO,CFO_Hz,Est_CFO_Hz,CFO_Error_Hz,Uncomp_BER,Uncomp_EVM_Percent,Comp_BER,Comp_EVM_Percent\n";

    std::cout << "\n=====================================================================================================\n";
    std::cout << "  CFO IMPACT: UNCOMPENSATED vs. COMPENSATED PERFORMANCE (QPSK, Delta_f = 15 kHz)\n";
    std::cout << "=====================================================================================================\n";
    std::cout << std::setw(8) << "Epsilon"
              << std::setw(12) << "CFO (Hz)"
              << std::setw(14) << "Est CFO (Hz)"
              << std::setw(14) << "Uncomp BER"
              << std::setw(14) << "Uncomp EVM"
              << std::setw(14) << "Comp BER"
              << std::setw(14) << "Comp EVM"
              << "\n";
    std::cout << "-----------------------------------------------------------------------------------------------------\n";

    // Also export a representative constellation for epsilon = 0.10
    std::ofstream raw_uncomp("results/tables/cfo_constellation_uncomp.csv");
    std::ofstream raw_comp("results/tables/cfo_constellation_comp.csv");
    raw_uncomp << "I,Q\n";
    raw_comp << "I,Q\n";

    for (double eps : eps_values) {
        double cfo_hz = eps * cfg.subcarrier_spacing_hz;
        ntn::synchronization::CFOChannel channel(cfo_hz, fs_hz);

        auto rx_impaired = tx_time;
        channel.apply_cfo(rx_impaired);

        // Add moderate channel noise (SNR = 20 dB) to test realistic conditions
        ntn::channel::AWGNChannel awgn(1234);
        awgn.add_noise_snr(rx_impaired, 20.0);

        // --- 1. Uncompensated Reception ---
        auto rx_symbols_uncomp = transceiver.receiver().demodulate_stream(rx_impaired, num_symbols);
        auto rx_bits_uncomp = qpsk.demodulate(rx_symbols_uncomp);
        auto ber_uncomp = ntn::metrics::calculate_ber(tx_bits, rx_bits_uncomp);

        auto tx_symbols = qpsk.modulate(tx_bits);
        double evm_uncomp = ntn::signal::calculate_evm_percent(tx_symbols, rx_symbols_uncomp);

        // --- 2. Compensated Reception ---
        double est_cfo_hz = estimator.estimate_cfo_hz_multisymbol(rx_impaired, num_symbols);
        ntn::synchronization::CFOCompensator compensator(est_cfo_hz, fs_hz);

        auto rx_corrected = rx_impaired;
        compensator.compensate(rx_corrected);

        auto rx_symbols_comp = transceiver.receiver().demodulate_stream(rx_corrected, num_symbols);
        auto rx_bits_comp = qpsk.demodulate(rx_symbols_comp);
        auto ber_comp = ntn::metrics::calculate_ber(tx_bits, rx_bits_comp);
        double evm_comp = ntn::signal::calculate_evm_percent(tx_symbols, rx_symbols_comp);

        // Save representative constellation samples for eps = 0.10
        if (std::abs(eps - 0.10) < 1e-6) {
            for (size_t k = 0; k < 1200; ++k) {
                raw_uncomp << rx_symbols_uncomp[k].real() << "," << rx_symbols_uncomp[k].imag() << "\n";
                raw_comp << rx_symbols_comp[k].real() << "," << rx_symbols_comp[k].imag() << "\n";
            }
        }

        std::cout << std::fixed << std::setprecision(2) << std::setw(8) << eps
                  << std::setw(12) << cfo_hz
                  << std::setw(14) << est_cfo_hz
                  << std::scientific << std::setprecision(3)
                  << std::setw(14) << ber_uncomp.ber
                  << std::fixed << std::setprecision(1)
                  << std::setw(13) << evm_uncomp << "%"
                  << std::scientific << std::setprecision(3)
                  << std::setw(14) << ber_comp.ber
                  << std::fixed << std::setprecision(1)
                  << std::setw(13) << evm_comp << "%"
                  << "\n";

        csv << eps << "," << cfo_hz << "," << est_cfo_hz << "," << std::abs(est_cfo_hz - cfo_hz) << ","
            << ber_uncomp.ber << "," << evm_uncomp << "," << ber_comp.ber << "," << evm_comp << "\n";
    }

    std::cout << "=====================================================================================================\n";
    std::cout << "Saved comparison table to results/tables/cfo_comparison.csv\n";

    return 0;
}
