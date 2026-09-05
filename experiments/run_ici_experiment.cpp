#include "compensation/ici_mitigation.hpp"
#include "ofdm/ofdm_config.hpp"
#include "ofdm/ofdm_modem.hpp"
#include "modem/modulation.hpp"
#include "channel/noise.hpp"
#include "metrics/metrics.hpp"
#include "signal/dsp_utils.hpp"
#include "signal/fft.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <random>
#include <cmath>
#include <numbers>

using namespace ntn;

int main() {
    std::cout << "====================================================================\n";
    std::cout << " NTN-OFDM SIMULATOR — INTER-CARRIER INTERFERENCE (ICI) EXPERIMENT   \n";
    std::cout << "====================================================================\n\n";

    ofdm::OFDMConfig config;
    config.fft_size = 512;
    config.num_active_subcarriers = 300;
    config.cp_length = 36;
    config.subcarrier_spacing_hz = 15000.0;

    // --- EXPERIMENT 1: CIR vs NORMALIZED FREQUENCY OFFSET EPSILON ---
    std::cout << "1. Measuring Carrier-to-Interference Ratio (CIR) across normalized offset eps...\n";
    std::ofstream cir_csv("results/tables/ici_cir_curve.csv");
    cir_csv << "epsilon,frequency_offset_hz,empirical_cir_db,analytical_cir_db\n";

    double pi = std::numbers::pi;
    for (int step = 0; step <= 30; ++step) {
        double eps = step * 0.01;
        double f_offset_hz = eps * config.subcarrier_spacing_hz;

        auto metrics = compensation::IciKernel::measure_ici(config, 256, eps);
        double analytical_cir_db = 100.0;
        if (eps > 0.001) {
            double ici_approx = (pi * pi / 3.0) * eps * eps;
            analytical_cir_db = 10.0 * std::log10(1.0 / ici_approx);
        }

        cir_csv << std::fixed << std::setprecision(4)
                << eps << "," << f_offset_hz << ","
                << metrics.cir_db << "," << analytical_cir_db << "\n";
    }
    cir_csv.close();
    std::cout << "   Saved CIR measurements to results/tables/ici_cir_curve.csv\n";

    // --- EXPERIMENT 2: SPECTRAL LEAKAGE PROFILE ---
    std::cout << "2. Exporting spectral leakage profile around active subcarrier...\n";
    std::ofstream spec_csv("results/tables/ici_leakage_profile.csv");
    spec_csv << "subcarrier_offset,power_eps_000,power_eps_005,power_eps_015,power_eps_030\n";

    auto m0 = compensation::IciKernel::measure_ici(config, 256, 0.00);
    auto m1 = compensation::IciKernel::measure_ici(config, 256, 0.05);
    auto m2 = compensation::IciKernel::measure_ici(config, 256, 0.15);
    auto m3 = compensation::IciKernel::measure_ici(config, 256, 0.30);

    for (int offset = -8; offset <= 8; ++offset) {
        int idx = offset + 8;
        spec_csv << offset << ","
                 << m0.leakage_spectrum[idx] << ","
                 << m1.leakage_spectrum[idx] << ","
                 << m2.leakage_spectrum[idx] << ","
                 << m3.leakage_spectrum[idx] << "\n";
    }
    spec_csv.close();
    std::cout << "   Saved spectral leakage profiles to results/tables/ici_leakage_profile.csv\n";

    // --- EXPERIMENT 3: ICI EQUALIZATION PERFORMANCE ---
    std::cout << "3. Evaluating Banded Tridiagonal ICI Equalizer on 16-QAM OFDM slots...\n";
    std::ofstream eq_csv("results/tables/ici_equalization_results.csv");
    eq_csv << "residual_eps,raw_evm_pct,derotated_evm_pct,equalized_evm_pct,raw_ber,derotated_ber,equalized_ber\n";

    std::cout << "========================================================================================\n";
    std::cout << " Residual eps   Raw EVM (%)   De-rot EVM (%)   Equalized EVM (%)   Raw BER    De-rot BER   Eq BER\n";
    std::cout << "========================================================================================\n";

    modem::QAM16Modulator qam16;
    ofdm::OFDMTransmitter tx(config);
    compensation::IciEqualizer equalizer(config);
    channel::AWGNChannel awgn(42);

    std::mt19937 rng(1337);
    std::uniform_int_distribution<int> bit_dist(0, 1);
    auto active_indices = config.active_subcarrier_indices();

    for (int step = 1; step <= 10; ++step) {
        double eps = step * 0.01; // 0.01 to 0.10

        size_t num_slots = 5;
        double sum_raw_evm = 0.0;
        double sum_derot_evm = 0.0;
        double sum_eq_evm = 0.0;
        uint64_t total_bits = 0;
        uint64_t raw_errors = 0;
        uint64_t derot_errors = 0;
        uint64_t eq_errors = 0;

        for (size_t slot = 0; slot < num_slots; ++slot) {
            size_t bits_per_slot = 14 * config.num_active_subcarriers * 4;
            ByteVector tx_bits(bits_per_slot);
            for (auto& b : tx_bits) b = static_cast<uint8_t>(bit_dist(rng));

            // Transmit slot
            auto tx_symbols = qam16.modulate(tx_bits);
            auto tx_time = tx.modulate_stream(tx_symbols);

            // Channel: AWGN (high SNR 25 dB to isolate ICI) + residual CFO
            double snr_db = 25.0;
            double p_sig = signal::calculate_power(tx_time);
            double noise_var = signal::calculate_noise_variance(p_sig, snr_db);
            ComplexVector rx_time = tx_time;
            awgn.add_noise_variance(rx_time, noise_var);

            // Apply residual CFO across the time stream
            double phase_step = 2.0 * pi * eps / static_cast<double>(config.fft_size);
            for (size_t n = 0; n < rx_time.size(); ++n) {
                rx_time[n] *= std::polar(1.0, phase_step * static_cast<double>(n));
            }

            // Demodulate symbol-by-symbol
            size_t symbol_len = config.symbol_duration_samples();
            ComplexVector slot_raw_subcarriers;
            ComplexVector slot_derot_subcarriers;
            ComplexVector slot_eq_subcarriers;

            for (size_t s = 0; s < 14; ++s) {
                ComplexVector sym_time(rx_time.begin() + s * symbol_len, rx_time.begin() + (s + 1) * symbol_len);

                // --- 1. Raw extraction ---
                ComplexVector sym_no_cp(sym_time.begin() + config.cp_length, sym_time.end());
                ComplexVector sym_fft = sym_no_cp;
                signal::FFT::forward(sym_fft);
                for (size_t idx : active_indices) {
                    slot_raw_subcarriers.push_back(sym_fft[idx]);
                }

                // --- 2. De-rotated symbol (compensate symbol starting phase) ---
                double sym_start_phase = phase_step * static_cast<double>(s * symbol_len + config.cp_length);
                ComplexVector sym_derot = sym_no_cp;
                for (auto& sample : sym_derot) {
                    sample *= std::polar(1.0, -sym_start_phase);
                }
                ComplexVector derot_fft = sym_derot;
                signal::FFT::forward(derot_fft);
                for (size_t idx : active_indices) {
                    slot_derot_subcarriers.push_back(derot_fft[idx]);
                }

                // --- 3. Banded ICI Equalization on de-rotated FFT ---
                ComplexVector sym_eq_fft = equalizer.equalize(derot_fft, eps);
                for (size_t idx : active_indices) {
                    slot_eq_subcarriers.push_back(sym_eq_fft[idx]);
                }
            }

            // EVM
            sum_raw_evm += signal::calculate_evm_percent(tx_symbols, slot_raw_subcarriers);
            sum_derot_evm += signal::calculate_evm_percent(tx_symbols, slot_derot_subcarriers);
            sum_eq_evm += signal::calculate_evm_percent(tx_symbols, slot_eq_subcarriers);

            // Demodulate bits
            auto rx_raw_bits = qam16.demodulate(slot_raw_subcarriers);
            auto rx_derot_bits = qam16.demodulate(slot_derot_subcarriers);
            auto rx_eq_bits = qam16.demodulate(slot_eq_subcarriers);

            raw_errors += metrics::calculate_ber(tx_bits, rx_raw_bits).bit_errors;
            derot_errors += metrics::calculate_ber(tx_bits, rx_derot_bits).bit_errors;
            eq_errors += metrics::calculate_ber(tx_bits, rx_eq_bits).bit_errors;
            total_bits += tx_bits.size();
        }

        double avg_raw_evm = sum_raw_evm / num_slots;
        double avg_derot_evm = sum_derot_evm / num_slots;
        double avg_eq_evm = sum_eq_evm / num_slots;
        double raw_ber = static_cast<double>(raw_errors) / total_bits;
        double derot_ber = static_cast<double>(derot_errors) / total_bits;
        double eq_ber = static_cast<double>(eq_errors) / total_bits;

        std::cout << "      " << std::fixed << std::setprecision(2) << eps
                  << "        " << std::setprecision(1) << std::setw(5) << avg_raw_evm << "%"
                  << "        " << std::setprecision(1) << std::setw(5) << avg_derot_evm << "%"
                  << "            " << std::setprecision(1) << std::setw(5) << avg_eq_evm << "%"
                  << "      " << std::scientific << std::setprecision(1) << raw_ber
                  << "    " << std::scientific << std::setprecision(1) << derot_ber
                  << "    " << std::scientific << std::setprecision(1) << eq_ber << "\n";

        eq_csv << std::fixed << std::setprecision(4)
               << eps << "," << avg_raw_evm << "," << avg_derot_evm << "," << avg_eq_evm << ","
               << std::scientific << raw_ber << "," << derot_ber << "," << eq_ber << "\n";
    }
    std::cout << "========================================================================================\n";
    eq_csv.close();
    std::cout << "Saved ICI equalization results to results/tables/ici_equalization_results.csv\n\n";

    return 0;
}
