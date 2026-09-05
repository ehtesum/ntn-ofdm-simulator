#include "common/types.hpp"
#include "signal/dsp_utils.hpp"
#include "signal/fft.hpp"
#include "channel/noise.hpp"
#include "channel/link_budget.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>
#include <cassert>

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
    std::cout << " NTN-OFDM SIMULATOR - DSP UNIT TESTS     \n";
    std::cout << "=========================================\n\n";

    // ----------------------------------------------------
    // Test 1: Signal Power, RMS, Normalization
    // ----------------------------------------------------
    std::cout << "[Suite 1: Baseband DSP Utilities]\n";
    {
        // Pure DC signal of magnitude 2.0
        ntn::ComplexVector dc_signal(100, ntn::Complex(2.0, 0.0));
        double pwr = ntn::signal::calculate_power(dc_signal);
        double rms = ntn::signal::calculate_rms(dc_signal);
        run_test("DC Signal Power is 4.0", std::abs(pwr - 4.0) < 1e-9);
        run_test("DC Signal RMS is 2.0", std::abs(rms - 2.0) < 1e-9);

        // Power normalization to 1.0
        ntn::signal::normalize_power(dc_signal, 1.0);
        double norm_pwr = ntn::signal::calculate_power(dc_signal);
        run_test("Normalized Power is 1.0", std::abs(norm_pwr - 1.0) < 1e-9);

        // Constant envelope complex exponential: |s[n]| = 1
        ntn::ComplexVector exp_signal(64);
        for (size_t n = 0; n < 64; ++n) {
            double angle = 2.0 * ntn::PI * 3.0 * static_cast<double>(n) / 64.0;
            exp_signal[n] = ntn::Complex(std::cos(angle), std::sin(angle));
        }
        double exp_pwr = ntn::signal::calculate_power(exp_signal);
        double papr_lin = ntn::signal::calculate_papr_linear(exp_signal);
        double papr_db = ntn::signal::calculate_papr_db(exp_signal);
        run_test("Complex exponential power is 1.0", std::abs(exp_pwr - 1.0) < 1e-9);
        run_test("Complex exponential PAPR is 1.0 (0 dB)", std::abs(papr_lin - 1.0) < 1e-9 && std::abs(papr_db) < 1e-9);
    }

    // ----------------------------------------------------
    // Test 2: FFT / IFFT Round-Trip Numerical Precision
    // ----------------------------------------------------
    std::cout << "\n[Suite 2: Radix-2 Cooley-Tukey FFT / IFFT]\n";
    {
        for (size_t n : {64, 128, 512, 1024}) {
            ntn::ComplexVector orig(n);
            for (size_t i = 0; i < n; ++i) {
                double re = std::sin(static_cast<double>(i) * 0.3) + std::cos(static_cast<double>(i) * 0.7);
                double im = std::cos(static_cast<double>(i) * 0.5) - std::sin(static_cast<double>(i) * 0.2);
                orig[i] = ntn::Complex(re, im);
            }

            ntn::ComplexVector spectrum = ntn::signal::FFT::forward_copy(orig);
            ntn::ComplexVector recovered = ntn::signal::FFT::inverse_copy(spectrum);

            double max_err = 0.0;
            for (size_t i = 0; i < n; ++i) {
                max_err = std::max(max_err, std::abs(recovered[i] - orig[i]));
            }
            std::string label = "FFT/IFFT Round-Trip N=" + std::to_string(n) + " (Max error < 1e-12)";
            run_test(label, max_err < 1e-12, "max_err=" + std::to_string(max_err));
        }
    }

    // ----------------------------------------------------
    // Test 3: Parseval's Theorem Verification
    // ----------------------------------------------------
    std::cout << "\n[Suite 3: Parseval's Energy Conservation]\n";
    {
        const size_t n = 256;
        ntn::ComplexVector time_signal(n);
        for (size_t i = 0; i < n; ++i) {
            time_signal[i] = ntn::Complex(std::sin(static_cast<double>(i) * 0.1), std::cos(static_cast<double>(i) * 0.2));
        }

        double time_energy = 0.0;
        for (const auto& s : time_signal) {
            time_energy += std::norm(s);
        }

        ntn::ComplexVector freq_signal = ntn::signal::FFT::forward_copy(time_signal);
        double freq_energy = 0.0;
        for (const auto& s : freq_signal) {
            freq_energy += std::norm(s);
        }
        freq_energy /= static_cast<double>(n); // Unnormalized FFT scaling factor 1/N

        double energy_diff = std::abs(time_energy - freq_energy);
        run_test("Parseval's Theorem holds (energy conservation)", energy_diff < 1e-10);
    }

    // ----------------------------------------------------
    // Test 4: Single-Tone Frequency Pinpointing
    // ----------------------------------------------------
    std::cout << "\n[Suite 4: FFT Spectral Peak Accuracy]\n";
    {
        const size_t n = 128;
        const size_t target_bin = 17;
        ntn::ComplexVector tone(n);
        for (size_t i = 0; i < n; ++i) {
            double angle = 2.0 * ntn::PI * static_cast<double>(target_bin * i) / static_cast<double>(n);
            tone[i] = ntn::Complex(std::cos(angle), std::sin(angle));
        }

        ntn::ComplexVector spectrum = ntn::signal::FFT::forward_copy(tone);

        size_t peak_bin = 0;
        double peak_mag = 0.0;
        for (size_t i = 0; i < n; ++i) {
            double mag = std::abs(spectrum[i]);
            if (mag > peak_mag) {
                peak_mag = mag;
                peak_bin = i;
            }
        }
        run_test("Peak frequency bin precisely matches target k=17", peak_bin == target_bin);
        run_test("Peak magnitude equals N (128.0)", std::abs(peak_mag - 128.0) < 1e-9);
    }

    // ----------------------------------------------------
    // Test 5: FFTShift / IFFTShift Symmetry
    // ----------------------------------------------------
    std::cout << "\n[Suite 5: FFTShift / IFFTShift]\n";
    {
        const size_t n = 64;
        ntn::ComplexVector data(n);
        for (size_t i = 0; i < n; ++i) {
            data[i] = ntn::Complex(static_cast<double>(i), -static_cast<double>(i));
        }
        ntn::ComplexVector shifted = data;
        ntn::signal::FFT::fftshift(shifted);
        run_test("DC element moved to index N/2=32", shifted[32] == data[0]);

        ntn::signal::FFT::ifftshift(shifted);
        bool match = true;
        for (size_t i = 0; i < n; ++i) {
            if (shifted[i] != data[i]) match = false;
        }
        run_test("ifftshift perfectly inverts fftshift", match);
    }

    // ----------------------------------------------------
    // Test 6: EVM (Error Vector Magnitude) Metrics
    // ----------------------------------------------------
    std::cout << "\n[Suite 6: EVM Calculations]\n";
    {
        ntn::ComplexVector ref = {ntn::Complex(1, 1), ntn::Complex(-1, 1), ntn::Complex(-1, -1), ntn::Complex(1, -1)};
        // Identical received signal
        double evm_zero = ntn::signal::calculate_evm_percent(ref, ref);
        run_test("EVM with zero error is 0.0%", std::abs(evm_zero) < 1e-9);

        // Perturbed received signal (10% error vector magnitude)
        ntn::ComplexVector rx = ref;
        for (auto& s : rx) {
            s += ntn::Complex(0.1 * s.real(), 0.1 * s.imag());
        }
        double evm_pert = ntn::signal::calculate_evm_percent(ref, rx);
        run_test("EVM with 10% offset is 10.0%", std::abs(evm_pert - 10.0) < 1e-6);
    }

    // ----------------------------------------------------
    // Test 7: AWGN Channel Statistics
    // ----------------------------------------------------
    std::cout << "\n[Suite 7: AWGN Statistical Properties]\n";
    {
        ntn::channel::AWGNChannel awgn(12345);
        const size_t num_samples = 200000;
        const double target_variance = 2.0;

        ntn::ComplexVector noise = awgn.generate_noise(num_samples, target_variance);

        double mean_re = 0.0;
        double mean_im = 0.0;
        double empirical_var = 0.0;
        for (const auto& s : noise) {
            mean_re += s.real();
            mean_im += s.imag();
            empirical_var += std::norm(s);
        }
        mean_re /= static_cast<double>(num_samples);
        mean_im /= static_cast<double>(num_samples);
        empirical_var /= static_cast<double>(num_samples);

        run_test("AWGN Mean(Re) ~ 0 (|mean| < 0.015)", std::abs(mean_re) < 0.015, "mean=" + std::to_string(mean_re));
        run_test("AWGN Mean(Im) ~ 0 (|mean| < 0.015)", std::abs(mean_im) < 0.015, "mean=" + std::to_string(mean_im));
        run_test("AWGN Empirical Variance ~ target 2.0 (tolerance 2.5%)", std::abs(empirical_var - target_variance) < 0.05,
                 "var=" + std::to_string(empirical_var));
    }

    // ----------------------------------------------------
    // Test 8: Satellite Link Budget Calculation
    // ----------------------------------------------------
    std::cout << "\n[Suite 8: Satellite Free-Space Path Loss & Link Budget]\n";
    {
        // Analytical test: Carrier = 2.0 GHz, Distance = 600 km
        // FSPL = 20*log10(6e5) + 20*log10(2e9) + 20*log10(4*pi/299792458)
        // FSPL = 115.563 + 186.021 - 147.552 = ~154.03 dB
        double fspl = ntn::channel::LinkBudget::calculate_fspl(600.0e3, 2.0e9);
        run_test("FSPL at 2 GHz, 600 km is ~154.03 dB", std::abs(fspl - 154.03) < 0.1, "fspl=" + std::to_string(fspl));

        ntn::channel::LinkBudgetParams params;
        params.carrier_frequency_hz = 2.0e9;
        params.distance_m = 600.0e3;
        params.tx_power_dbm = 23.0; // Handheld / NTN UE (200 mW)
        params.tx_antenna_gain_dbi = 0.0;
        params.rx_antenna_gain_dbi = 30.0; // LEO satellite antenna
        params.bandwidth_hz = 180.0e3;    // 1 PRB (180 kHz)
        params.system_noise_temp_k = 290.0;
        params.misc_losses_db = 2.0;

        auto result = ntn::channel::LinkBudget::evaluate(params);
        run_test("Link Budget calculation completes and returns valid SNR", !std::isnan(result.snr_db));
        run_test("Thermal noise density N0 ~ -174 dBm/Hz at 290 K",
                 std::abs(result.noise_power_density_dbm_per_hz - (-174.0)) < 0.2);
    }

    // ----------------------------------------------------
    // Summary
    // ----------------------------------------------------
    std::cout << "\n=========================================\n";
    std::cout << " TEST RESULTS: " << passed_tests << " / " << total_tests << " passed.\n";
    std::cout << "=========================================\n\n";

    return (passed_tests == total_tests) ? 0 : 1;
}
