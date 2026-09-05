#include "common/types.hpp"
#include "modem/modulation.hpp"
#include "metrics/metrics.hpp"
#include "signal/dsp_utils.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>

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
    std::cout << " NTN-OFDM SIMULATOR - MODULATION TESTS   \n";
    std::cout << "=========================================\n\n";

    // ----------------------------------------------------
    // Test Suite 1: BPSK Modulator / Demodulator
    // ----------------------------------------------------
    std::cout << "[Suite 1: BPSK Modem]\n";
    {
        ntn::modem::BPSKModulator bpsk;
        ntn::ByteVector bits = {0, 1, 0, 0, 1, 1, 0, 1};
        auto syms = bpsk.modulate(bits);

        run_test("BPSK symbol count equals bit count", syms.size() == bits.size());
        run_test("Bit 0 maps to +1.0", syms[0].real() == 1.0 && syms[0].imag() == 0.0);
        run_test("Bit 1 maps to -1.0", syms[1].real() == -1.0 && syms[1].imag() == 0.0);

        double pwr = ntn::signal::calculate_power(syms);
        run_test("BPSK constellation average power is 1.0", std::abs(pwr - 1.0) < 1e-12);

        auto recovered = bpsk.demodulate(syms);
        run_test("BPSK ideal loopback produces identical bits", recovered == bits);

        // Perturbation within decision region
        ntn::ComplexVector noisy_syms = {ntn::Complex(0.05, 0.8), ntn::Complex(-0.1, -0.4)};
        auto noisy_rec = bpsk.demodulate(noisy_syms);
        run_test("BPSK hard decision slices positive Re to 0, negative Re to 1",
                 noisy_rec[0] == 0 && noisy_rec[1] == 1);
    }

    // ----------------------------------------------------
    // Test Suite 2: QPSK Modulator / Demodulator & Gray Coding
    // ----------------------------------------------------
    std::cout << "\n[Suite 2: QPSK Modem & Gray Mapping]\n";
    {
        ntn::modem::QPSKModulator qpsk;
        // All 4 2-bit combinations: 00, 01, 10, 11
        ntn::ByteVector bits = {0, 0,  0, 1,  1, 0,  1, 1};
        auto syms = qpsk.modulate(bits);

        run_test("QPSK produces N_bits / 2 symbols", syms.size() == 4);

        double pwr = ntn::signal::calculate_power(syms);
        run_test("QPSK constellation average power is 1.0", std::abs(pwr - 1.0) < 1e-12);

        auto recovered = qpsk.demodulate(syms);
        run_test("QPSK ideal loopback recovers all 4 quadrants with 0 errors", recovered == bits);

        // Gray Coding Property:
        // Quadrant 1 (00) -> Quadrant 4 (01): only 1 bit flips (second bit)
        // Quadrant 1 (00) -> Quadrant 2 (10): only 1 bit flips (first bit)
        // Diagonal: Quadrant 1 (00) -> Quadrant 3 (11): 2 bits flip
        auto sym_00 = qpsk.modulate({0, 0})[0];
        auto sym_01 = qpsk.modulate({0, 1})[0];
        auto sym_10 = qpsk.modulate({1, 0})[0];
        run_test("Adjacent quadrant 00 and 01 share Re axis", sym_00.real() == sym_01.real());
        run_test("Adjacent quadrant 00 and 10 share Im axis", sym_00.imag() == sym_10.imag());

        // Error handling: odd number of bits
        bool threw = false;
        try {
            [[maybe_unused]] auto unused = qpsk.modulate({0, 1, 0});
        } catch (const std::invalid_argument&) {
            threw = true;
        }
        run_test("QPSK throws invalid_argument on odd bit length", threw);
    }

    // ----------------------------------------------------
    // Test Suite 3: 16-QAM Modulator / Demodulator
    // ----------------------------------------------------
    std::cout << "\n[Suite 3: 16-QAM Modem]\n";
    {
        ntn::modem::QAM16Modulator qam16;

        // Generate all 16 possible 4-bit symbols: 0 to 15 (64 bits)
        ntn::ByteVector all_bits;
        all_bits.reserve(64);
        for (uint8_t i = 0; i < 16; ++i) {
            all_bits.push_back((i >> 3) & 1);
            all_bits.push_back((i >> 2) & 1);
            all_bits.push_back((i >> 1) & 1);
            all_bits.push_back(i & 1);
        }

        auto syms = qam16.modulate(all_bits);
        run_test("16-QAM produces 16 constellation symbols from 64 bits", syms.size() == 16);

        double pwr = ntn::signal::calculate_power(syms);
        run_test("16-QAM constellation average power is exactly 1.0", std::abs(pwr - 1.0) < 1e-12,
                 "pwr=" + std::to_string(pwr));

        auto recovered = qam16.demodulate(syms);
        run_test("16-QAM ideal loopback recovers all 16 symbols without error", recovered == all_bits);

        // Verification of decision boundaries under noise within margin
        // Add small perturbation (< 1.0 * inv_norm)
        ntn::ComplexVector perturbed = syms;
        for (auto& s : perturbed) {
            s += ntn::Complex(0.05, -0.05);
        }
        auto recovered_noisy = qam16.demodulate(perturbed);
        run_test("16-QAM hard decision handles small perturbations correctly", recovered_noisy == all_bits);
    }

    // ----------------------------------------------------
    // Test Suite 4: BER, SER, and Theoretical Curves
    // ----------------------------------------------------
    std::cout << "\n[Suite 4: Error Rate Metrics & Analytical Bounds]\n";
    {
        ntn::ByteVector tx = {0, 1, 0, 1, 1, 1, 0, 0, 1, 0};
        ntn::ByteVector rx = {0, 1, 0, 0, 1, 1, 0, 0, 0, 0}; // 2 errors (indices 3 and 8)

        auto ber = ntn::metrics::calculate_ber(tx, rx);
        run_test("BER calculation: 2 errors out of 10 bits is 0.2",
                 ber.bit_errors == 2 && std::abs(ber.ber - 0.2) < 1e-12);

        // Q-function checks
        // Q(0) = 0.5
        run_test("Q(0) = 0.5", std::abs(ntn::metrics::q_function(0.0) - 0.5) < 1e-9);
        // Q(3.0) ~ 0.0013499
        run_test("Q(3.0) ~ 0.00135", std::abs(ntn::metrics::q_function(3.0) - 0.001349898) < 1e-6);

        // Theoretical BER monotonic decrease
        double ber_4db = ntn::metrics::theoretical_ber_qpsk(4.0);
        double ber_8db = ntn::metrics::theoretical_ber_qpsk(8.0);
        double ber_12db = ntn::metrics::theoretical_ber_qpsk(12.0);
        run_test("Theoretical QPSK BER strictly decreases with higher SNR",
                 ber_4db > ber_8db && ber_8db > ber_12db);

        // 16-QAM vs QPSK at same Eb/N0: QPSK has lower BER than 16-QAM
        double ber_qam_8db = ntn::metrics::theoretical_ber_16qam(8.0);
        run_test("QPSK has lower BER than 16-QAM at same Eb/N0", ber_8db < ber_qam_8db);
    }

    // ----------------------------------------------------
    // Summary
    // ----------------------------------------------------
    std::cout << "\n=========================================\n";
    std::cout << " TEST RESULTS: " << passed_tests << " / " << total_tests << " passed.\n";
    std::cout << "=========================================\n\n";

    return (passed_tests == total_tests) ? 0 : 1;
}
