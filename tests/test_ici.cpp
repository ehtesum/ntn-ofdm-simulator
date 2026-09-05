#include "compensation/ici_mitigation.hpp"
#include "ofdm/ofdm_config.hpp"
#include "ofdm/ofdm_modem.hpp"
#include "modem/modulation.hpp"
#include "metrics/metrics.hpp"
#include "signal/dsp_utils.hpp"
#include "signal/fft.hpp"
#include <iostream>
#include <cassert>
#include <cmath>
#include <numbers>
#include <vector>

using namespace ntn;

int main() {
    std::cout << "\n======================================\n";
    std::cout << " NTN-OFDM SIMULATOR - ICI UNIT TESTS \n";
    std::cout << "======================================\n\n";

    int tests_passed = 0;
    int total_tests = 0;

    auto test_assert = [&](bool condition, const std::string& name) {
        total_tests++;
        if (condition) {
            std::cout << "  [PASS] " << name << "\n";
            tests_passed++;
        } else {
            std::cout << "  [FAIL] " << name << "\n";
            assert(false);
        }
    };

    ofdm::OFDMConfig config;
    config.fft_size = 64;
    config.num_active_subcarriers = 48;
    config.cp_length = 16;
    config.subcarrier_spacing_hz = 15000.0;

    // Suite 1: ICI Kernel Orthogonality
    std::cout << "[Suite 1: OFDM Orthogonality & Kernel Properties]\n";
    {
        // When epsilon = 0, diagonal is 1.0, off-diagonal is 0.0
        Complex s_diag = compensation::IciKernel::leakage_coefficient(10, 10, 64, 0.0);
        test_assert(std::abs(s_diag - Complex(1.0, 0.0)) < 1e-10, "S_{k,k}(0) == 1.0");

        Complex s_off = compensation::IciKernel::leakage_coefficient(10, 11, 64, 0.0);
        test_assert(std::abs(s_off) < 1e-10, "S_{k, m!=k}(0) == 0.0 (Orthogonality preserved)");

        // Symmetry: S_{k, k+1}(eps) and S_{k, k-1}(-eps)
        Complex s_plus = compensation::IciKernel::leakage_coefficient(10, 11, 64, 0.05);
        Complex s_minus = compensation::IciKernel::leakage_coefficient(10, 9, 64, -0.05);
        test_assert(std::abs(std::abs(s_plus) - std::abs(s_minus)) < 1e-6, "Leakage magnitude symmetric under signed offset");
    }

    // Suite 2: Empirical CIR vs Analytical Theory
    std::cout << "\n[Suite 2: Empirical CIR & Spectral Leakage]\n";
    {
        // Ideal case epsilon = 0.0
        auto m_zero = compensation::IciKernel::measure_ici(config, 32, 0.0);
        test_assert(m_zero.cir_db > 100.0, "Zero offset yields CIR > 100 dB");

        // Small offset epsilon = 0.05 (5% subcarrier spacing)
        // Analytical CIR ~ 10 * log10( 1 / (pi^2 / 3 * eps^2) ) ~ 20.85 dB
        auto m_small = compensation::IciKernel::measure_ici(config, 32, 0.05);
        double pi = std::numbers::pi;
        double expected_cir_db = 10.0 * std::log10(1.0 / ((pi * pi / 3.0) * (0.05 * 0.05)));
        test_assert(std::abs(m_small.cir_db - expected_cir_db) < 2.0, "Empirical CIR at eps=0.05 matches analytical ~20.8 dB");

        // Larger offset epsilon = 0.20 -> CIR drops to ~9 dB
        auto m_large = compensation::IciKernel::measure_ici(config, 32, 0.20);
        test_assert(m_large.cir_db < 12.0 && m_large.cir_db > 7.0, "Empirical CIR at eps=0.20 drops to ~9 dB");
    }

    // Suite 3: Common Phase Error (CPE) Mitigation
    std::cout << "\n[Suite 3: Common Phase Error (CPE) Mitigation]\n";
    {
        compensation::IciEqualizer eq(config);

        ComplexVector subcarriers(config.fft_size, Complex(1.0, 0.0));
        // Rotate all subcarriers by 35 degrees
        double true_rot = 35.0 * std::numbers::pi / 180.0;
        for (auto& s : subcarriers) {
            s *= std::polar(1.0, true_rot);
        }

        // Known pilot subcarriers at indices {12, 24, 36, 48}
        std::vector<int> pilot_idx = {12, 24, 36, 48};
        ComplexVector pilots = {Complex(1.0, 0.0), Complex(1.0, 0.0), Complex(1.0, 0.0), Complex(1.0, 0.0)};

        auto corrected = eq.compensate_cpe(subcarriers, pilot_idx, pilots);
        // Error after CPE correction should be nearly 0
        double residual_err = std::abs(corrected[20] - Complex(1.0, 0.0));
        test_assert(residual_err < 1e-6, "CPE compensation perfectly corrects uniform constellation phase rotation");
    }

    // Suite 4: Banded Tridiagonal Equalizer
    std::cout << "\n[Suite 4: Banded Frequency-Domain ICI Equalizer]\n";
    {
        modem::QPSKModulator qpsk;
        ofdm::OFDMTransmitter tx(config);
        ofdm::OFDMReceiver rx(config);
        compensation::IciEqualizer eq(config);

        // Generate known symbols
        std::vector<uint8_t> bits;
        for (size_t i = 0; i < config.num_active_subcarriers * 2; ++i) {
            bits.push_back(static_cast<uint8_t>(i % 2));
        }
        auto symbols = qpsk.modulate(bits);
        auto tx_time = tx.modulate_symbol(symbols);

        // Apply small residual offset epsilon = 0.04
        double eps = 0.04;
        double phase_step = 2.0 * std::numbers::pi * eps / static_cast<double>(config.fft_size);
        ComplexVector rx_time = tx_time;
        for (size_t n = 0; n < rx_time.size(); ++n) {
            rx_time[n] *= std::polar(1.0, phase_step * static_cast<double>(n));
        }

        // Demodulate raw without ICI equalizer
        auto rx_raw_subcarriers = rx.demodulate_symbol(rx_time);
        double evm_before = signal::calculate_evm_percent(symbols, rx_raw_subcarriers);

        // Perform FFT manually to equalize before demapping
        ComplexVector sym_no_cp(rx_time.begin() + config.cp_length, rx_time.end());
        ComplexVector rx_fft = sym_no_cp;
        signal::FFT::forward(rx_fft);

        // Equalize using known residual epsilon
        ComplexVector equalized_fft = eq.equalize(rx_fft, eps);
        // Demap active subcarriers
        auto active_indices = config.active_subcarrier_indices();
        ComplexVector rx_eq_subcarriers;
        rx_eq_subcarriers.reserve(active_indices.size());
        for (size_t idx : active_indices) {
            rx_eq_subcarriers.push_back(equalized_fft[idx]);
        }
        double evm_after = signal::calculate_evm_percent(symbols, rx_eq_subcarriers);

        test_assert(evm_after < evm_before, "Banded ICI equalizer reduces EVM caused by residual frequency offset");
    }

    std::cout << "\n======================================\n";
    std::cout << " TEST RESULTS: " << tests_passed << " / " << total_tests << " passed.\n";
    std::cout << "======================================\n\n";

    return (tests_passed == total_tests) ? 0 : 1;
}
