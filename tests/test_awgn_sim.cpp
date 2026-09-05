#include "common/types.hpp"
#include "simulation/awgn_sim.hpp"

#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>

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
    std::cout << " NTN-OFDM SIMULATOR - AWGN SIM TESTS     \n";
    std::cout << "=========================================\n\n";

    ntn::ofdm::OFDMConfig ofdm_cfg;
    ofdm_cfg.fft_size = 256;
    ofdm_cfg.num_active_subcarriers = 128;
    ofdm_cfg.cp_length = 16;
    ofdm_cfg.null_dc = true;

    ntn::simulation::AWGNSimConfig sim_cfg;
    sim_cfg.ofdm_config = ofdm_cfg;
    sim_cfg.modulation = ntn::modem::ModulationType::QPSK;
    sim_cfg.random_seed = 12345;

    ntn::simulation::AWGNSimulator sim(sim_cfg);

    // ----------------------------------------------------
    // Test 1: Low SNR vs High SNR Error Monotonicity
    // ----------------------------------------------------
    std::cout << "[Suite 1: SNR vs Error Monotonicity]\n";
    {
        auto pt_low  = sim.run_point(2.0, 50);
        auto pt_mid  = sim.run_point(6.0, 50);
        auto pt_high = sim.run_point(10.0, 50);

        run_test("BER monotonically decreases with SNR (2 dB > 6 dB > 10 dB)",
                 pt_low.ber > pt_mid.ber && pt_mid.ber > pt_high.ber,
                 "ber_low=" + std::to_string(pt_low.ber) + ", ber_high=" + std::to_string(pt_high.ber));

        run_test("EVM (%) strictly decreases as SNR increases",
                 pt_low.evm_percent > pt_mid.evm_percent && pt_mid.evm_percent > pt_high.evm_percent,
                 "evm_low=" + std::to_string(pt_low.evm_percent) + "%, evm_high=" + std::to_string(pt_high.evm_percent) + "%");
    }

    // ----------------------------------------------------
    // Test 2: Convergence to Theoretical Bounds at Moderate SNR
    // ----------------------------------------------------
    std::cout << "\n[Suite 2: Theoretical Agreement]\n";
    {
        // At SNR = 4.0 dB, QPSK Eb/N0 = 4.0 - 10*log10(2) = 0.99 dB
        // Theoretical BER ~ 0.056
        auto pt = sim.run_point(4.0, 200); // 200 * 128 * 2 = 51,200 bits
        double relative_error = std::abs(pt.ber - pt.theoretical_ber) / pt.theoretical_ber;
        run_test("Simulated QPSK BER at 4 dB matches theory within 15% relative error",
                 relative_error < 0.15,
                 "sim_ber=" + std::to_string(pt.ber) + ", theory=" + std::to_string(pt.theoretical_ber));
    }

    // ----------------------------------------------------
    // Test 3: Simulation Determinism & Reproducibility
    // ----------------------------------------------------
    std::cout << "\n[Suite 3: Deterministic Reproducibility]\n";
    {
        ntn::simulation::AWGNSimulator sim1(sim_cfg);
        ntn::simulation::AWGNSimulator sim2(sim_cfg);

        auto pt1 = sim1.run_point(6.0, 30);
        auto pt2 = sim2.run_point(6.0, 30);

        run_test("Same random seed produces identical bit error count", pt1.bit_errors == pt2.bit_errors);
        run_test("Same random seed produces identical empirical EVM", pt1.evm_percent == pt2.evm_percent);
    }

    // ----------------------------------------------------
    // Summary
    // ----------------------------------------------------
    std::cout << "\n=========================================\n";
    std::cout << " TEST RESULTS: " << passed_tests << " / " << total_tests << " passed.\n";
    std::cout << "=========================================\n\n";

    return (passed_tests == total_tests) ? 0 : 1;
}
