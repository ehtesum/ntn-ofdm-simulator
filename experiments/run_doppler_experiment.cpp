#include "common/types.hpp"
#include "ofdm/ofdm_config.hpp"
#include "ofdm/ofdm_modem.hpp"
#include "modem/modulation.hpp"
#include "metrics/metrics.hpp"
#include "signal/dsp_utils.hpp"
#include "channel/doppler.hpp"
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
    std::cout << " NTN-OFDM SIMULATOR — LEO SATELLITE DOPPLER EXPERIMENT\n";
    std::cout << "====================================================================\n";

    ntn::channel::LEOOrbitParams orbit;
    orbit.altitude_km = 600.0;
    orbit.carrier_frequency_hz = 2.0e9; // 2.0 GHz S-band
    orbit.min_slant_range_km = 600.0;
    orbit.pass_duration_s = 600.0;

    ntn::channel::DopplerTrajectory trajectory(orbit);

    std::cout << "\n--- LEO ORBIT PARAMETERS (3GPP NTN S-Band) ---\n"
              << "Orbital Altitude:       " << orbit.altitude_km << " km\n"
              << "Orbital Velocity:       " << (orbit.orbital_velocity_mps() / 1e3) << " km/s\n"
              << "Carrier Frequency:      " << (orbit.carrier_frequency_hz / 1e9) << " GHz\n"
              << "Max Doppler Shift:      " << (orbit.max_doppler_shift_hz() / 1e3) << " kHz\n"
              << "Max Doppler Rate:       " << orbit.max_doppler_rate_hz_per_sec() << " Hz/s\n"
              << "-----------------------------------------------\n\n";

    // -------------------------------------------------------------
    // 1. Export Full Pass S-Curve Data (-300s to +300s)
    // -------------------------------------------------------------
    std::ofstream s_curve_csv("results/tables/doppler_trajectory.csv");
    s_curve_csv << "Time_sec,Doppler_Hz,DopplerRate_Hz_per_sec\n";

    for (double t = -300.0; t <= 300.0; t += 2.0) {
        double fd = trajectory.instantaneous_doppler_hz(t);
        double dfd = trajectory.instantaneous_doppler_rate_hz_per_sec(t);
        s_curve_csv << t << "," << fd << "," << dfd << "\n";
    }
    std::cout << "Saved 600-second LEO pass trajectory to results/tables/doppler_trajectory.csv\n";

    // -------------------------------------------------------------
    // 2. Simulate Active Slot-by-Slot Communications Session
    // -------------------------------------------------------------
    ntn::ofdm::OFDMConfig cfg;
    cfg.fft_size = 512;
    cfg.num_active_subcarriers = 300;
    cfg.cp_length = 36;
    cfg.subcarrier_spacing_hz = 15000.0;
    const double fs = cfg.sampling_rate_hz();

    ntn::ofdm::OFDMTransceiver transceiver(cfg);
    ntn::modem::QPSKModulator qpsk;
    ntn::channel::DopplerChannel channel(fs);
    ntn::channel::DopplerTracker tracker(cfg);

    // Run over 50 consecutive OFDM slots (each slot = 14 OFDM symbols)
    const size_t symbols_per_slot = 14;
    const size_t num_slots = 30;
    const double slot_duration_s = cfg.symbol_duration_seconds() * symbols_per_slot; // ~1 ms

    std::cout << "\nSimulating " << num_slots << " consecutive 5G NR slots under time-varying Doppler...\n";

    std::ofstream track_csv("results/tables/doppler_tracking.csv");
    track_csv << "Slot,Time_sec,Actual_Doppler_Hz,Tracked_Doppler_Hz,Estimation_Error_Hz,Uncomp_BER,Comp_BER,Uncomp_EVM,Comp_EVM\n";

    std::cout << "========================================================================================================\n";
    std::cout << " Slot     Time(s)    Actual fd (Hz)   Tracked fd (Hz)   Error (Hz)   Uncomp BER   Comp BER   Comp EVM\n";
    std::cout << "========================================================================================================\n";

    std::mt19937_64 rng(12345);
    std::uniform_int_distribution<int> bit_dist(0, 1);

    // Scenario: Satellite is approaching nadir (highest Doppler rate region)
    // Let starting time t = -5.0 seconds
    double current_t = -5.0;

    for (size_t slot = 0; slot < num_slots; ++slot) {
        // Transmit 1 slot of payload bits
        const size_t slot_bits = transceiver.bits_per_ofdm_symbol(qpsk) * symbols_per_slot;
        ntn::ByteVector tx_bits(slot_bits);
        for (auto& b : tx_bits) b = static_cast<uint8_t>(bit_dist(rng));

        auto tx_time = transceiver.transmit_bits(tx_bits, qpsk);

        // Pre-compensation: gNB/UE pre-compensates known orbital coarse offset from ephemeris
        // Leaving residual Doppler within the tracking range
        double actual_doppler_hz = trajectory.instantaneous_doppler_hz(current_t);
        double coarse_ephemeris_compensation = std::round(actual_doppler_hz / 500.0) * 500.0;
        double residual_doppler_hz = actual_doppler_hz - coarse_ephemeris_compensation;

        auto rx_impaired = tx_time;
        channel.apply_constant_doppler(rx_impaired, residual_doppler_hz);

        // Add 18 dB AWGN
        ntn::channel::AWGNChannel awgn(999 + slot);
        awgn.add_noise_snr(rx_impaired, 18.0);

        // --- Uncompensated Reception ---
        auto rx_bits_uncomp = transceiver.receive_bits(rx_impaired, symbols_per_slot, qpsk);
        auto ber_uncomp = ntn::metrics::calculate_ber(tx_bits, rx_bits_uncomp);
        auto rx_syms_uncomp = transceiver.receiver().demodulate_stream(rx_impaired, symbols_per_slot);
        auto tx_syms = qpsk.modulate(tx_bits);
        double evm_uncomp = ntn::signal::calculate_evm_percent(tx_syms, rx_syms_uncomp);

        // --- Doppler Tracking & Compensation ---
        auto rx_corrected = rx_impaired;
        const size_t sym_len = cfg.symbol_duration_samples();
        double slot_tracked_cfo = 0.0;

        for (size_t s = 0; s < symbols_per_slot; ++s) {
            ntn::ComplexVector sym_slice(rx_impaired.begin() + s * sym_len,
                                         rx_impaired.begin() + (s + 1) * sym_len);
            slot_tracked_cfo = tracker.track_symbol(sym_slice);
            tracker.compensate_block(sym_slice, slot_tracked_cfo);
            std::copy(sym_slice.begin(), sym_slice.end(), rx_corrected.begin() + s * sym_len);
        }

        auto rx_bits_comp = transceiver.receive_bits(rx_corrected, symbols_per_slot, qpsk);
        auto ber_comp = ntn::metrics::calculate_ber(tx_bits, rx_bits_comp);
        auto rx_syms_comp = transceiver.receiver().demodulate_stream(rx_corrected, symbols_per_slot);
        double evm_comp = ntn::signal::calculate_evm_percent(tx_syms, rx_syms_comp);

        double est_error = std::abs(slot_tracked_cfo - residual_doppler_hz);

        if (slot % 3 == 0 || slot == num_slots - 1) {
            std::cout << std::setw(5) << slot
                      << std::fixed << std::setprecision(3) << std::setw(12) << current_t
                      << std::setprecision(1) << std::setw(18) << actual_doppler_hz
                      << std::setw(18) << (coarse_ephemeris_compensation + slot_tracked_cfo)
                      << std::setw(13) << est_error
                      << std::scientific << std::setprecision(2) << std::setw(13) << ber_uncomp.ber
                      << std::setw(11) << ber_comp.ber
                      << std::fixed << std::setprecision(1) << std::setw(10) << evm_comp << "%\n";
        }

        track_csv << slot << "," << current_t << "," << actual_doppler_hz << ","
                  << (coarse_ephemeris_compensation + slot_tracked_cfo) << ","
                  << est_error << ","
                  << ber_uncomp.ber << "," << ber_comp.ber << ","
                  << evm_uncomp << "," << evm_comp << "\n";

        current_t += slot_duration_s;
    }

    std::cout << "========================================================================================================\n";
    std::cout << "Saved Doppler tracking results to results/tables/doppler_tracking.csv\n";

    return 0;
}
