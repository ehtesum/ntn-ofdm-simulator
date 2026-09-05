#include "simulation/awgn_sim.hpp"
#include "signal/dsp_utils.hpp"
#include <fstream>
#include <iomanip>
#include <random>
#include <cmath>
#include <iostream>

namespace ntn::simulation {

AWGNSimulator::AWGNSimulator(AWGNSimConfig config)
    : config_(config),
      modulator_(modem::create_modulator(config.modulation)),
      tx_(config.ofdm_config),
      rx_(config.ofdm_config),
      channel_(config.random_seed) {
}

AWGNSimPoint AWGNSimulator::run_point(double snr_db, size_t num_ofdm_symbols) {
    AWGNSimPoint pt;
    pt.snr_db = snr_db;

    const size_t num_active = config_.ofdm_config.num_active_subcarriers;
    const size_t bps = modulator_->bits_per_symbol();
    const size_t total_bits = num_ofdm_symbols * num_active * bps;

    pt.total_bits = total_bits;
    pt.total_symbols = num_ofdm_symbols * num_active;

    // Generate random bits
    std::mt19937_64 bit_rng(config_.random_seed + static_cast<uint64_t>(std::abs(snr_db) * 100.0));
    std::uniform_int_distribution<int> bit_dist(0, 1);
    ByteVector tx_bits(total_bits);
    for (auto& b : tx_bits) {
        b = static_cast<uint8_t>(bit_dist(bit_rng));
    }

    // 1. Modulate bits to QPSK/QAM constellation symbols
    ComplexVector tx_freq_symbols = modulator_->modulate(tx_bits);

    // 2. Modulate frequency symbols onto OFDM time-domain signal
    ComplexVector time_tx = tx_.modulate_stream(tx_freq_symbols);

    // 3. Pass through AWGN channel
    ComplexVector time_rx = time_tx;
    channel_.add_noise_snr(time_rx, snr_db);

    // 4. OFDM Demodulation: CP removal -> FFT -> subcarrier extraction
    ComplexVector rx_freq_symbols = rx_.demodulate_stream(time_rx, num_ofdm_symbols);

    // 5. Hard-decision demodulation to bits
    ByteVector rx_bits = modulator_->demodulate(rx_freq_symbols);

    // 6. Compute metrics
    auto ber_res = metrics::calculate_ber(tx_bits, rx_bits);
    auto ser_res = metrics::calculate_ser(tx_freq_symbols, rx_freq_symbols, *modulator_);

    pt.bit_errors = ber_res.bit_errors;
    pt.ber = ber_res.ber;
    pt.symbol_errors = ser_res.symbol_errors;
    pt.ser = ser_res.ser;
    pt.evm_percent = signal::calculate_evm_percent(tx_freq_symbols, rx_freq_symbols);
    pt.evm_db = signal::calculate_evm_db(tx_freq_symbols, rx_freq_symbols);

    // In OFDM, receiver FFT acts as a subcarrier matched filter bank.
    // Noise falling outside active subcarriers (guard bands) is rejected.
    // Hence, subcarrier SNR = SNR_channel + 10 * log10(N_fft / N_active).
    // Eb/N0 = SNR_subcarrier - 10 * log10(bits_per_symbol).
    const double guard_gain_db = 10.0 * std::log10(static_cast<double>(config_.ofdm_config.fft_size) /
                                                   static_cast<double>(num_active));
    const double subcarrier_snr_db = snr_db + guard_gain_db;
    pt.eb_n0_db = subcarrier_snr_db - 10.0 * std::log10(static_cast<double>(bps));

    // Theoretical closed-form comparison evaluated at exact subcarrier Eb/N0
    if (config_.modulation == modem::ModulationType::BPSK) {
        pt.theoretical_ber = metrics::theoretical_ber_bpsk(pt.eb_n0_db);
    } else if (config_.modulation == modem::ModulationType::QPSK) {
        pt.theoretical_ber = metrics::theoretical_ber_qpsk(pt.eb_n0_db);
    } else {
        pt.theoretical_ber = metrics::theoretical_ber_16qam(pt.eb_n0_db);
    }

    return pt;
}

std::vector<AWGNSimPoint> AWGNSimulator::run_sweep() {
    std::vector<AWGNSimPoint> results;
    results.reserve(config_.snr_range_db.size());

    for (double snr : config_.snr_range_db) {
        results.push_back(run_point(snr, config_.num_ofdm_symbols_per_snr));
    }
    return results;
}

void AWGNSimulator::export_csv(const std::string& filepath, const std::vector<AWGNSimPoint>& results, const std::string& mod_name) {
    std::ofstream ofs(filepath);
    if (!ofs.is_open()) {
        throw std::runtime_error("Could not open CSV file for writing: " + filepath);
    }

    ofs << "Modulation,SNR_dB,EbN0_dB,TotalBits,BitErrors,BER,Theoretical_BER,TotalSymbols,SymbolErrors,SER,EVM_Percent,EVM_dB\n";
    ofs << std::scientific << std::setprecision(6);

    for (const auto& pt : results) {
        ofs << mod_name << ","
            << std::fixed << std::setprecision(2) << pt.snr_db << ","
            << pt.eb_n0_db << ","
            << std::scientific
            << pt.total_bits << ","
            << pt.bit_errors << ","
            << pt.ber << ","
            << pt.theoretical_ber << ","
            << pt.total_symbols << ","
            << pt.symbol_errors << ","
            << pt.ser << ","
            << std::fixed << std::setprecision(3)
            << pt.evm_percent << ","
            << pt.evm_db << "\n";
    }
}

} // namespace ntn::simulation
