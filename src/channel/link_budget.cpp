#include "channel/link_budget.hpp"
#include <sstream>
#include <iomanip>
#include <cmath>

namespace ntn::channel {

double LinkBudget::calculate_slant_range(double altitude_m, double elevation_deg) {
    if (elevation_deg >= 90.0) {
        return altitude_m;
    }
    const double r_earth = 6371.0e3; // Earth radius in meters
    const double elev_rad = elevation_deg * (PI / 180.0);
    const double cos_el = std::cos(elev_rad);
    const double sin_el = std::sin(elev_rad);

    // Slant range equation based on spherical Earth geometry
    double term = (r_earth + altitude_m) / r_earth;
    double d = r_earth * (std::sqrt(term * term - cos_el * cos_el) - sin_el);
    return std::max(altitude_m, d);
}

double LinkBudget::calculate_fspl(double distance_m, double frequency_hz) {
    if (distance_m <= 0.0 || frequency_hz <= 0.0) {
        return 0.0;
    }
    // FSPL (dB) = 20 * log10(d) + 20 * log10(f) + 20 * log10(4 * pi / c)
    // 20 * log10(4 * pi / c) = 20 * log10(4 * pi / 299792458) = -147.5522
    const double c = SPEED_OF_LIGHT;
    return 20.0 * std::log10(distance_m) + 20.0 * std::log10(frequency_hz) + 20.0 * std::log10(4.0 * PI / c);
}

LinkBudgetResult LinkBudget::evaluate(const LinkBudgetParams& params, double bit_rate_bps) {
    LinkBudgetResult result;

    const double c = SPEED_OF_LIGHT;
    result.wavelength_m = c / params.carrier_frequency_hz;
    result.fspl_db = calculate_fspl(params.distance_m, params.carrier_frequency_hz);
    result.eirp_dbm = params.tx_power_dbm + params.tx_antenna_gain_dbi;

    // Prx = EIRP - FSPL + Grx - Misc_losses
    result.rx_power_dbm = result.eirp_dbm - result.fspl_db + params.rx_antenna_gain_dbi - params.misc_losses_db;

    // Noise spectral density: N0 = k_B * T_sys (Watts/Hz)
    // In dBm/Hz: 10 * log10(k_B * T_sys) + 30
    const double n0_watts_per_hz = BOLTZMANN_CONSTANT * params.system_noise_temp_k;
    result.noise_power_density_dbm_per_hz = 10.0 * std::log10(n0_watts_per_hz) + 30.0;

    // Total thermal noise power: N = N0 * B (in dBm)
    const double total_noise_watts = n0_watts_per_hz * params.bandwidth_hz;
    result.total_noise_power_dbm = 10.0 * std::log10(total_noise_watts) + 30.0;

    // SNR = Prx - N_total
    result.snr_db = result.rx_power_dbm - result.total_noise_power_dbm;

    // Eb/N0 = SNR + 10 * log10(Bandwidth / BitRate)
    if (bit_rate_bps > 0.0) {
        result.eb_n0_db = result.snr_db + 10.0 * std::log10(params.bandwidth_hz / bit_rate_bps);
    } else {
        result.eb_n0_db = result.snr_db;
    }

    std::ostringstream oss;
    oss << std::fixed << std::setprecision(2);
    oss << "=== SATELLITE LINK BUDGET SUMMARY ===\n"
        << "Carrier Frequency:     " << (params.carrier_frequency_hz / 1e9) << " GHz\n"
        << "Wavelength:            " << (result.wavelength_m * 100.0) << " cm\n"
        << "Slant Range (d):       " << (params.distance_m / 1e3) << " km\n"
        << "Tx Power:              " << params.tx_power_dbm << " dBm (" << dbm_to_watts(params.tx_power_dbm) << " W)\n"
        << "Tx Antenna Gain:       " << params.tx_antenna_gain_dbi << " dBi\n"
        << "Tx EIRP:               " << result.eirp_dbm << " dBm\n"
        << "Free-Space Path Loss:  " << result.fspl_db << " dB\n"
        << "Rx Antenna Gain:       " << params.rx_antenna_gain_dbi << " dBi\n"
        << "Misc Losses:           " << params.misc_losses_db << " dB\n"
        << "Received Power (Prx):  " << result.rx_power_dbm << " dBm\n"
        << "Noise Temp (Tsys):     " << params.system_noise_temp_k << " K\n"
        << "Noise Density (N0):    " << result.noise_power_density_dbm_per_hz << " dBm/Hz\n"
        << "Bandwidth (B):         " << (params.bandwidth_hz / 1e3) << " kHz\n"
        << "Total Noise Power:     " << result.total_noise_power_dbm << " dBm\n"
        << "Received SNR (C/N):    " << result.snr_db << " dB\n";
    if (bit_rate_bps > 0.0) {
        oss << "Bit Rate (Rb):         " << (bit_rate_bps / 1e3) << " kbps\n"
            << "Eb/N0:                 " << result.eb_n0_db << " dB\n";
    }
    oss << "=====================================\n";
    result.summary = oss.str();

    return result;
}

} // namespace ntn::channel
