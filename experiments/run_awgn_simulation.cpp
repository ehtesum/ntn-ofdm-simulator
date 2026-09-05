#include "simulation/awgn_sim.hpp"
#include <iostream>
#include <iomanip>
#include <filesystem>

namespace fs = std::filesystem;

void print_table(const std::string& title, const std::vector<ntn::simulation::AWGNSimPoint>& results) {
    std::cout << "\n========================================================================================\n";
    std::cout << " " << title << "\n";
    std::cout << "========================================================================================\n";
    std::cout << std::setw(8) << "SNR(dB)"
              << std::setw(10) << "Eb/N0(dB)"
              << std::setw(12) << "Bits"
              << std::setw(10) << "Errors"
              << std::setw(14) << "Sim BER"
              << std::setw(14) << "Theory BER"
              << std::setw(12) << "SER"
              << std::setw(12) << "EVM(%)"
              << "\n";
    std::cout << "----------------------------------------------------------------------------------------\n";

    for (const auto& pt : results) {
        std::cout << std::fixed << std::setprecision(1) << std::setw(8) << pt.snr_db
                  << std::setw(10) << pt.eb_n0_db
                  << std::setw(12) << pt.total_bits
                  << std::setw(10) << pt.bit_errors
                  << std::scientific << std::setprecision(3)
                  << std::setw(14) << pt.ber
                  << std::setw(14) << pt.theoretical_ber
                  << std::setw(12) << pt.ser
                  << std::fixed << std::setprecision(2)
                  << std::setw(12) << pt.evm_percent
                  << "\n";
    }
    std::cout << "========================================================================================\n";
}

int main() {
    fs::create_directories("results/tables");

    std::cout << "====================================================================\n";
    std::cout << " NTN-OFDM SIMULATOR — MONTE CARLO AWGN PERFORMANCE SWEEP\n";
    std::cout << "====================================================================\n";

    ntn::ofdm::OFDMConfig ofdm_cfg;
    ofdm_cfg.fft_size = 512;
    ofdm_cfg.num_active_subcarriers = 300;
    ofdm_cfg.cp_length = 36;
    ofdm_cfg.null_dc = true;

    // -------------------------------------------------------------
    // 1. QPSK Simulation (0 dB to 16 dB)
    // -------------------------------------------------------------
    {
        ntn::simulation::AWGNSimConfig qpsk_sim_cfg;
        qpsk_sim_cfg.ofdm_config = ofdm_cfg;
        qpsk_sim_cfg.modulation = ntn::modem::ModulationType::QPSK;
        qpsk_sim_cfg.snr_range_db = {0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0};
        qpsk_sim_cfg.num_ofdm_symbols_per_snr = 300; // 300 * 300 * 2 = 180,000 bits per SNR
        qpsk_sim_cfg.random_seed = 42;

        std::cout << "\nRunning QPSK Simulation (" << qpsk_sim_cfg.num_ofdm_symbols_per_snr
                  << " OFDM symbols / point, 180,000 bits/SNR)..." << std::flush;

        ntn::simulation::AWGNSimulator sim(qpsk_sim_cfg);
        auto qpsk_results = sim.run_sweep();
        std::cout << " Done.\n";

        print_table("QPSK over AWGN: Empirical vs. Theoretical BER", qpsk_results);
        ntn::simulation::AWGNSimulator::export_csv("results/tables/awgn_results_qpsk.csv", qpsk_results, "QPSK");
        std::cout << "Saved QPSK results to results/tables/awgn_results_qpsk.csv\n";
    }

    // -------------------------------------------------------------
    // 2. 16-QAM Simulation (4 dB to 20 dB)
    // -------------------------------------------------------------
    {
        ntn::simulation::AWGNSimConfig qam_sim_cfg;
        qam_sim_cfg.ofdm_config = ofdm_cfg;
        qam_sim_cfg.modulation = ntn::modem::ModulationType::QAM16;
        qam_sim_cfg.snr_range_db = {4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
        qam_sim_cfg.num_ofdm_symbols_per_snr = 300; // 300 * 300 * 4 = 360,000 bits per SNR
        qam_sim_cfg.random_seed = 123;

        std::cout << "\nRunning 16-QAM Simulation (" << qam_sim_cfg.num_ofdm_symbols_per_snr
                  << " OFDM symbols / point, 360,000 bits/SNR)..." << std::flush;

        ntn::simulation::AWGNSimulator sim(qam_sim_cfg);
        auto qam_results = sim.run_sweep();
        std::cout << " Done.\n";

        print_table("16-QAM over AWGN: Empirical vs. Theoretical BER", qam_results);
        ntn::simulation::AWGNSimulator::export_csv("results/tables/awgn_results_16qam.csv", qam_results, "16-QAM");
        std::cout << "Saved 16-QAM results to results/tables/awgn_results_16qam.csv\n";
    }

    return 0;
}
