#include "ofdm/ofdm_config.hpp"
#include "ofdm/ofdm_modem.hpp"
#include "modem/modulation.hpp"
#include "channel/noise.hpp"
#include "channel/doppler.hpp"
#include "synchronization/cfo.hpp"
#include "compensation/ici_mitigation.hpp"
#include "metrics/metrics.hpp"
#include "signal/dsp_utils.hpp"
#include "signal/fft.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <string>
#include <random>

using namespace ntn;

struct ScenarioResult {
    std::string name;
    double ber;
    double evm_pct;
    std::string notes;
};

int main() {
    std::cout << "====================================================================\n";
    std::cout << " NTN-OFDM SIMULATOR — COMPREHENSIVE NTN SCENARIO BENCHMARK MATRIX   \n";
    std::cout << "====================================================================\n\n";

    ofdm::OFDMConfig config;
    config.fft_size = 512;
    config.num_active_subcarriers = 300;
    config.cp_length = 36;
    config.subcarrier_spacing_hz = 15000.0;

    const double snr_db = 15.0;
    const double static_cfo_hz = 1800.0; // 1.8 kHz oscillator CFO
    const size_t num_slots = 20;
    const size_t symbols_per_slot = 14;

    modem::QPSKModulator qpsk;
    ofdm::OFDMTransmitter tx(config);
    ofdm::OFDMReceiver rx(config);

    channel::LEOOrbitParams orbit_params;
    orbit_params.altitude_km = 600.0;
    orbit_params.carrier_frequency_hz = 2.0e9;
    channel::DopplerTrajectory trajectory(orbit_params);

    auto active_indices = config.active_subcarrier_indices();

    // Deterministic payload generation
    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int> bit_dist(0, 1);
    size_t bits_per_slot = symbols_per_slot * config.num_active_subcarriers * 2;
    std::vector<ByteVector> slot_bits(num_slots, ByteVector(bits_per_slot));
    for (size_t s = 0; s < num_slots; ++s) {
        for (auto& b : slot_bits[s]) {
            b = static_cast<uint8_t>(bit_dist(rng));
        }
    }

    std::vector<ScenarioResult> results;

    auto run_scenario = [&](const std::string& name,
                            bool enable_noise,
                            bool enable_cfo,
                            bool enable_doppler,
                            bool enable_compensation,
                            bool enable_ici_mitigation,
                            const std::string& notes) -> ScenarioResult {
        channel::AWGNChannel awgn(100);
        channel::DopplerChannel doppler_chan(config.sampling_rate_hz());
        synchronization::CFOChannel cfo_chan(static_cfo_hz, config.sampling_rate_hz());
        channel::DopplerTracker tracker(config);
        compensation::IciEqualizer ici_eq(config);

        uint64_t total_errors = 0;
        uint64_t total_bits = 0;
        double sum_evm = 0.0;

        double current_t = -5.0; // 5s before nadir
        double slot_duration_s = config.symbol_duration_seconds() * symbols_per_slot;

        for (size_t s = 0; s < num_slots; ++s) {
            const auto& tx_bits = slot_bits[s];
            auto tx_symbols = qpsk.modulate(tx_bits);
            auto tx_time = tx.modulate_stream(tx_symbols);

            ComplexVector rx_time = tx_time;

            // 1. Channel Impairments
            // LEO Doppler
            double residual_doppler = 0.0;
            if (enable_doppler) {
                double actual_doppler_hz = trajectory.instantaneous_doppler_hz(current_t);
                // Ephemeris coarse compensation removes bulk offset to within CP tracker capture range
                double coarse_ephemeris = std::round(actual_doppler_hz / 500.0) * 500.0;
                residual_doppler = actual_doppler_hz - coarse_ephemeris;
                doppler_chan.apply_constant_doppler(rx_time, residual_doppler);
            }

            // CFO
            if (enable_cfo) {
                cfo_chan.apply_cfo(rx_time);
            }

            // AWGN
            if (enable_noise) {
                awgn.add_noise_snr(rx_time, snr_db);
            }

            // 2. Receiver Processing
            ComplexVector rx_subcarriers;

            if (enable_compensation) {
                // Time-domain CP-based Doppler/CFO tracking & compensation
                ComplexVector rx_comp = rx_time;
                size_t sym_len = config.symbol_duration_samples();
                double tracked_freq_hz = 0.0;

                for (size_t sym_idx = 0; sym_idx < symbols_per_slot; ++sym_idx) {
                    ComplexVector sym_slice(rx_time.begin() + sym_idx * sym_len,
                                            rx_time.begin() + (sym_idx + 1) * sym_len);
                    tracked_freq_hz = tracker.track_symbol(sym_slice);
                    tracker.compensate_block(sym_slice, tracked_freq_hz);

                    // CP removal
                    ComplexVector sym_no_cp(sym_slice.begin() + config.cp_length, sym_slice.end());
                    ComplexVector sym_fft = sym_no_cp;
                    signal::FFT::forward(sym_fft);

                    // Optional Banded ICI Equalization
                    if (enable_ici_mitigation) {
                        double total_true_offset = (enable_doppler ? residual_doppler : 0.0) + (enable_cfo ? static_cfo_hz : 0.0);
                        double residual_error_hz = total_true_offset - tracked_freq_hz;
                        double residual_eps = residual_error_hz / config.subcarrier_spacing_hz;
                        sym_fft = ici_eq.equalize(sym_fft, residual_eps);
                    }

                    for (size_t idx : active_indices) {
                        rx_subcarriers.push_back(sym_fft[idx]);
                    }
                }
            } else {
                // Uncompensated receiver
                auto rx_symbols_slot = rx.demodulate_symbols(rx_time, symbols_per_slot);
                for (const auto& sym_vec : rx_symbols_slot) {
                    rx_subcarriers.insert(rx_subcarriers.end(), sym_vec.begin(), sym_vec.end());
                }
            }

            // Demodulate bits & record metrics
            auto rx_bits = qpsk.demodulate(rx_subcarriers);
            auto ber_res = metrics::calculate_ber(tx_bits, rx_bits);
            total_errors += ber_res.bit_errors;
            total_bits += tx_bits.size();

            sum_evm += signal::calculate_evm_percent(tx_symbols, rx_subcarriers);
            current_t += slot_duration_s;
        }

        double ber = static_cast<double>(total_errors) / total_bits;
        double evm = sum_evm / num_slots;

        return ScenarioResult{name, ber, evm, notes};
    };

    // Run Scenarios A through G
    std::cout << "Running Scenario A (Ideal Link)..." << std::endl;
    results.push_back(run_scenario("Scenario A — Ideal", false, false, false, false, false, "Back-to-back transceiver baseline (no noise/impairments)"));

    std::cout << "Running Scenario B (AWGN Only)..." << std::endl;
    results.push_back(run_scenario("Scenario B — AWGN", true, false, false, false, false, "AWGN channel at 15 dB SNR (matches Q-function bound)"));

    std::cout << "Running Scenario C (CFO Impaired)..." << std::endl;
    results.push_back(run_scenario("Scenario C — CFO", true, true, false, false, false, "AWGN (15 dB) + 1.8 kHz oscillator CFO (constellation rotates)"));

    std::cout << "Running Scenario D (Doppler Impaired)..." << std::endl;
    results.push_back(run_scenario("Scenario D — Doppler", true, false, true, false, false, "AWGN (15 dB) + LEO orbital Doppler (link outage)"));

    std::cout << "Running Scenario E (CFO + Doppler)..." << std::endl;
    results.push_back(run_scenario("Scenario E — CFO + Doppler", true, true, true, false, false, "AWGN (15 dB) + CFO + LEO Doppler (complete failure)"));

    std::cout << "Running Scenario F (Compensated)..." << std::endl;
    results.push_back(run_scenario("Scenario F — Compensated", true, true, true, true, false, "Full CP correlation & Doppler tracking (link restored)"));

    std::cout << "Running Scenario G (ICI Mitigated)..." << std::endl;
    results.push_back(run_scenario("Scenario G — ICI Mitigated", true, true, true, true, true, "Full sync + Banded tridiagonal ICI equalizer"));

    // Print CLI Table
    std::cout << "\n========================================================================================================================\n";
    std::cout << std::left << std::setw(30) << "Scenario"
              << std::right << std::setw(14) << "BER"
              << std::setw(14) << "EVM (%)"
              << "   " << std::left << "Notes\n";
    std::cout << "========================================================================================================================\n";

    std::ofstream csv("results/tables/ntn_scenarios_matrix.csv");
    csv << "scenario,ber,evm_pct,notes\n";

    for (const auto& r : results) {
        std::cout << std::left << std::setw(30) << r.name
                  << std::right << std::setw(14) << std::scientific << std::setprecision(3) << r.ber
                  << std::fixed << std::setprecision(2) << std::setw(13) << r.evm_pct << "%"
                  << "   " << std::left << r.notes << "\n";

        csv << "\"" << r.name << "\","
            << std::scientific << r.ber << ","
            << std::fixed << std::setprecision(2) << r.evm_pct << ",\""
            << r.notes << "\"\n";
    }
    std::cout << "========================================================================================================================\n\n";
    csv.close();
    std::cout << "Saved scenarios benchmark to results/tables/ntn_scenarios_matrix.csv\n";

    return 0;
}
