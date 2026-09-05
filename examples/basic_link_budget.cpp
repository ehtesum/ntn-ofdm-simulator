#include "channel/link_budget.hpp"
#include "signal/dsp_utils.hpp"
#include <iostream>

int main() {
    std::cout << "=========================================================\n";
    std::cout << " 3GPP NTN S-Band (2 GHz) Satellite Link Budget Demo\n";
    std::cout << "=========================================================\n\n";

    ntn::channel::LinkBudgetParams params;
    params.carrier_frequency_hz = 2.0e9;         // 2.0 GHz S-band
    params.distance_m = 600.0e3;                 // 600 km LEO orbit (nadir)
    params.tx_power_dbm = 23.0;                  // Handheld UE (23 dBm = 200 mW)
    params.tx_antenna_gain_dbi = 0.0;            // Omnidirectional UE antenna
    params.rx_antenna_gain_dbi = 28.0;           // High-gain satellite spot beam
    params.bandwidth_hz = 180.0e3;               // 1 5G PRB (12 subcarriers * 15 kHz)
    params.system_noise_temp_k = 290.0;          // Standard reference temperature
    params.misc_losses_db = 2.5;                 // Atmospheric + polarization loss

    double bit_rate = 150.0e3;                   // 150 kbps QPSK payload
    auto result = ntn::channel::LinkBudget::evaluate(params, bit_rate);

    std::cout << result.summary << "\n";

    std::cout << "Key Takeaway for Skylo:\n";
    std::cout << "At 600 km slant range in S-band, FSPL is ~154 dB.\n";
    std::cout << "With a 23 dBm handheld terminal and a 28 dBi satellite beam,\n";
    std::cout << "the uplink SNR over 1 PRB (180 kHz) is " << result.snr_db << " dB,\n";
    std::cout << "sufficient for robust QPSK modulation.\n";

    return 0;
}
