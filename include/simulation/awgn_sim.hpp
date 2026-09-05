#pragma once

#include "common/types.hpp"
#include "ofdm/ofdm_config.hpp"
#include "ofdm/ofdm_modem.hpp"
#include "modem/modulation.hpp"
#include "metrics/metrics.hpp"
#include "channel/noise.hpp"

#include <vector>
#include <string>

namespace ntn::simulation {

struct AWGNSimPoint {
    double snr_db{0.0};
    double eb_n0_db{0.0};
    double ber{0.0};
    double ser{0.0};
    double evm_percent{0.0};
    double evm_db{0.0};
    double theoretical_ber{0.0};
    size_t total_bits{0};
    size_t bit_errors{0};
    size_t total_symbols{0};
    size_t symbol_errors{0};
};

struct AWGNSimConfig {
    ofdm::OFDMConfig ofdm_config{};
    modem::ModulationType modulation{modem::ModulationType::QPSK};
    std::vector<double> snr_range_db{0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
    size_t num_ofdm_symbols_per_snr{200}; // e.g., 200 symbols * 300 sc * 2 bits = 120,000 bits
    uint64_t random_seed{42};
};

class AWGNSimulator {
public:
    explicit AWGNSimulator(AWGNSimConfig config);

    // Executes a single SNR point and measures empirical BER, SER, and EVM
    [[nodiscard]] AWGNSimPoint run_point(double snr_db, size_t num_ofdm_symbols);

    // Executes the full SNR sweep defined in config_
    [[nodiscard]] std::vector<AWGNSimPoint> run_sweep();

    // Exports results to a CSV file
    static void export_csv(const std::string& filepath, const std::vector<AWGNSimPoint>& results, const std::string& mod_name);

    [[nodiscard]] const AWGNSimConfig& config() const noexcept { return config_; }

private:
    AWGNSimConfig config_;
    std::unique_ptr<modem::Modulator> modulator_;
    ofdm::OFDMTransmitter tx_;
    ofdm::OFDMReceiver rx_;
    channel::AWGNChannel channel_;
};

} // namespace ntn::simulation
