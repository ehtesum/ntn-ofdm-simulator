#pragma once

#include "common/types.hpp"
#include <string>

namespace ntn::channel {

// Input parameters for a satellite link budget
struct LinkBudgetParams {
    double carrier_frequency_hz{2.0e9};      // e.g., 2.0 GHz (S-band 3GPP NTN)
    double distance_m{600.0e3};               // Distance / Slant range in meters (e.g. 600 km LEO)
    double tx_power_dbm{23.0};                // User Equipment / Satellite Tx Power (dBm) (e.g. 23 dBm = 200 mW)
    double tx_antenna_gain_dbi{0.0};          // Tx antenna gain (dBi)
    double rx_antenna_gain_dbi{25.0};         // Rx antenna gain (dBi) (e.g. satellite reflector or phased array)
    double bandwidth_hz{15.0e3 * 12 * 4};     // Channel bandwidth (e.g. 4 PRBs @ 15 kHz SCS = 720 kHz)
    double system_noise_temp_k{290.0};        // System noise temperature in Kelvin
    double misc_losses_db{3.0};               // Atmospheric, polarization, and implementation losses (dB)
};

// Computed link budget output
struct LinkBudgetResult {
    double wavelength_m{0.0};
    double fspl_db{0.0};                      // Free-Space Path Loss (dB)
    double eirp_dbm{0.0};                     // Equivalent Isotropically Radiated Power (dBm)
    double rx_power_dbm{0.0};                 // Received power at antenna terminals (dBm)
    double noise_power_density_dbm_per_hz{0.0}; // Thermal noise density N0 (dBm/Hz)
    double total_noise_power_dbm{0.0};        // Total noise power across bandwidth (dBm)
    double snr_db{0.0};                       // Received SNR / C/N (dB)
    double eb_n0_db{0.0};                     // Eb/N0 (dB) for a given bit rate
    std::string summary;
};

class LinkBudget {
public:
    // Calculates slant range from satellite altitude (m) and ground elevation angle (degrees)
    [[nodiscard]] static double calculate_slant_range(double altitude_m, double elevation_deg);

    // Calculates Free-Space Path Loss (FSPL) in dB
    [[nodiscard]] static double calculate_fspl(double distance_m, double frequency_hz);

    // Evaluates full end-to-end link budget
    [[nodiscard]] static LinkBudgetResult evaluate(const LinkBudgetParams& params, double bit_rate_bps = 0.0);
};

} // namespace ntn::channel
