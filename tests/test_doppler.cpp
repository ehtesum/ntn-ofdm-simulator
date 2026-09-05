#include "common/types.hpp"
#include "ofdm/ofdm_config.hpp"
#include "ofdm/ofdm_modem.hpp"
#include "modem/modulation.hpp"
#include "metrics/metrics.hpp"
#include "channel/doppler.hpp"

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
    std::cout << " NTN-OFDM SIMULATOR - DOPPLER UNIT TESTS \n";
    std::cout << "=========================================\n\n";

    // ----------------------------------------------------
    // Test Suite 1: LEO Orbital Mechanics & Physics
    // ----------------------------------------------------
    std::cout << "[Suite 1: LEO Orbital Mechanics & Doppler Equations]\n";
    {
        ntn::channel::LEOOrbitParams params;
        params.altitude_km = 600.0;
        params.carrier_frequency_hz = 2.0e9; // 2.0 GHz
        params.min_slant_range_km = 600.0;

        double v_orb = params.orbital_velocity_mps();
        run_test("LEO orbital velocity at 600 km is ~7.56 km/s",
                 v_orb > 7500.0 && v_orb < 7600.0, "v=" + std::to_string(v_orb));

        double fd_max = params.max_doppler_shift_hz();
        run_test("Max Doppler shift at 2.0 GHz is ~50.45 kHz",
                 std::abs(fd_max - 50448.0) < 200.0, "fd_max=" + std::to_string(fd_max));

        double dfd_dt_max = params.max_doppler_rate_hz_per_sec();
        run_test("Max Doppler rate at closest approach is ~ -636 Hz/s",
                 std::abs(dfd_dt_max - (-636.0)) < 15.0, "dfd_dt=" + std::to_string(dfd_dt_max));
    }

    // ----------------------------------------------------
    // Test Suite 2: Doppler Trajectory (S-Curve) Properties
    // ----------------------------------------------------
    std::cout << "\n[Suite 2: Doppler Trajectory S-Curve Mathematical Symmetry]\n";
    {
        ntn::channel::LEOOrbitParams params;
        params.altitude_km = 600.0;
        params.carrier_frequency_hz = 2.0e9;
        params.min_slant_range_km = 600.0;

        ntn::channel::DopplerTrajectory traj(params);

        // t = 0: Doppler must be precisely 0 Hz at closest approach
        double fd_0 = traj.instantaneous_doppler_hz(0.0);
        run_test("Doppler at closest approach (t=0) is exactly 0 Hz", std::abs(fd_0) < 1e-9);

        // Approaching (t < 0): positive Doppler
        double fd_neg = traj.instantaneous_doppler_hz(-50.0);
        run_test("Approaching satellite (t=-50s) exhibits positive Doppler shift", fd_neg > 0.0,
                 "fd(-50)=" + std::to_string(fd_neg));

        // Receding (t > 0): negative Doppler
        double fd_pos = traj.instantaneous_doppler_hz(50.0);
        run_test("Receding satellite (t=+50s) exhibits negative Doppler shift", fd_pos < 0.0,
                 "fd(+50)=" + std::to_string(fd_pos));

        // Odd function symmetry: f_d(-t) == -f_d(t)
        run_test("Doppler trajectory is an exact odd function: fd(-t) == -fd(t)",
                 std::abs(fd_neg - (-fd_pos)) < 1e-9);
    }

    // ----------------------------------------------------
    // Test Suite 3: Doppler Channel Phase Integration
    // ----------------------------------------------------
    std::cout << "\n[Suite 3: Constant & Time-Varying Doppler Channel]\n";
    {
        const double fs = 1.92e6; // 1.92 MHz sampling rate
        ntn::channel::DopplerChannel channel(fs);

        // Apply constant 10 kHz Doppler to DC signal
        ntn::ComplexVector sig(1920, ntn::Complex(1.0, 0.0)); // 1 ms duration
        channel.apply_constant_doppler(sig, 10000.0);

        // Over 1 ms at 10 kHz, signal should complete exactly 10 full cycles (20*pi phase rotation)
        double total_phase = 2.0 * ntn::PI * 10000.0 * (1919.0 / fs);
        double sample_phase = std::arg(sig[1919]);
        double expected_mod_phase = std::fmod(total_phase, 2.0 * ntn::PI);
        if (expected_mod_phase > ntn::PI) expected_mod_phase -= 2.0 * ntn::PI;
        if (expected_mod_phase < -ntn::PI) expected_mod_phase += 2.0 * ntn::PI;

        run_test("Constant Doppler phase matches 10 kHz analytical rotation",
                 std::abs(sample_phase - expected_mod_phase) < 1e-6);
    }

    // ----------------------------------------------------
    // Test Suite 4: End-to-End Doppler Tracking & BER Restoration
    // ----------------------------------------------------
    std::cout << "\n[Suite 4: End-to-End Doppler Tracking & Recovery]\n";
    {
        ntn::ofdm::OFDMConfig cfg;
        cfg.fft_size = 512;
        cfg.num_active_subcarriers = 300;
        cfg.cp_length = 36;
        cfg.subcarrier_spacing_hz = 15000.0;

        ntn::ofdm::OFDMTransceiver transceiver(cfg);
        ntn::modem::QPSKModulator qpsk;
        ntn::channel::DopplerChannel channel(cfg.sampling_rate_hz());
        ntn::channel::DopplerTracker tracker(cfg);

        const size_t num_symbols = 20;
        const size_t total_bits = transceiver.bits_per_ofdm_symbol(qpsk) * num_symbols; // 12,000 bits

        std::mt19937_64 rng(54321);
        std::uniform_int_distribution<int> bit_dist(0, 1);
        ntn::ByteVector tx_bits(total_bits);
        for (auto& b : tx_bits) b = static_cast<uint8_t>(bit_dist(rng));

        auto tx_time = transceiver.transmit_bits(tx_bits, qpsk);

        // Residual uncompensated Doppler: 2,500 Hz
        const double residual_doppler_hz = 2500.0;
        auto rx_impaired = tx_time;
        channel.apply_constant_doppler(rx_impaired, residual_doppler_hz);

        // Uncompensated reception
        auto rx_bits_uncomp = transceiver.receive_bits(rx_impaired, num_symbols, qpsk);
        auto ber_uncomp = ntn::metrics::calculate_ber(tx_bits, rx_bits_uncomp);
        run_test("Uncompensated Doppler causes complete link failure (BER > 0.25)",
                 ber_uncomp.ber > 0.25, "ber=" + std::to_string(ber_uncomp.ber));

        // Doppler Tracking: symbol-by-symbol estimation and compensation
        auto rx_corrected = rx_impaired;
        const size_t sym_len = cfg.symbol_duration_samples();

        for (size_t s = 0; s < num_symbols; ++s) {
            ntn::ComplexVector sym_slice(rx_impaired.begin() + s * sym_len,
                                         rx_impaired.begin() + (s + 1) * sym_len);
            double tracked_cfo = tracker.track_symbol(sym_slice);
            tracker.compensate_block(sym_slice, tracked_cfo);
            std::copy(sym_slice.begin(), sym_slice.end(), rx_corrected.begin() + s * sym_len);
        }

        auto rx_bits_comp = transceiver.receive_bits(rx_corrected, num_symbols, qpsk);
        auto ber_comp = ntn::metrics::calculate_ber(tx_bits, rx_bits_comp);

        run_test("Doppler Tracker restores reception with 0 errors (BER = 0.0)",
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
