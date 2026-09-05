#include "signal/dsp_utils.hpp"
#include <algorithm>
#include <numeric>
#include <cmath>

namespace ntn::signal {

double calculate_power(const ComplexVector& signal) {
    if (signal.empty()) {
        return 0.0;
    }
    double total_energy = 0.0;
    for (const auto& sample : signal) {
        total_energy += std::norm(sample); // |sample|^2 = real^2 + imag^2
    }
    return total_energy / static_cast<double>(signal.size());
}

double calculate_rms(const ComplexVector& signal) {
    return std::sqrt(calculate_power(signal));
}

double calculate_peak_power(const ComplexVector& signal) {
    if (signal.empty()) {
        return 0.0;
    }
    double peak = 0.0;
    for (const auto& sample : signal) {
        peak = std::max(peak, std::norm(sample));
    }
    return peak;
}

double calculate_papr_linear(const ComplexVector& signal) {
    double avg_pwr = calculate_power(signal);
    if (avg_pwr <= 1e-15) {
        return 1.0;
    }
    return calculate_peak_power(signal) / avg_pwr;
}

double calculate_papr_db(const ComplexVector& signal) {
    return linear_to_db(calculate_papr_linear(signal));
}

void normalize_power(ComplexVector& signal, double target_power) {
    if (signal.empty() || target_power <= 0.0) {
        return;
    }
    double current_power = calculate_power(signal);
    if (current_power <= 1e-15) {
        return;
    }
    double scale = std::sqrt(target_power / current_power);
    scale_signal(signal, scale);
}

void scale_signal(ComplexVector& signal, Complex scalar) {
    for (auto& sample : signal) {
        sample *= scalar;
    }
}

void scale_signal(ComplexVector& signal, double scalar) {
    for (auto& sample : signal) {
        sample *= scalar;
    }
}

double calculate_evm_percent(const ComplexVector& reference, const ComplexVector& received) {
    if (reference.empty() || reference.size() != received.size()) {
        throw std::invalid_argument("Reference and received vectors must have equal, non-zero length.");
    }
    double error_energy = 0.0;
    double reference_energy = 0.0;

    for (size_t i = 0; i < reference.size(); ++i) {
        Complex error = received[i] - reference[i];
        error_energy += std::norm(error);
        reference_energy += std::norm(reference[i]);
    }

    if (reference_energy <= 1e-15) {
        return 0.0;
    }
    return std::sqrt(error_energy / reference_energy) * 100.0;
}

double calculate_evm_db(const ComplexVector& reference, const ComplexVector& received) {
    double evm_linear = calculate_evm_percent(reference, received) / 100.0;
    return (evm_linear > 1e-15) ? 20.0 * std::log10(evm_linear) : -300.0;
}

double calculate_noise_variance(double signal_power, double snr_db) {
    double snr_linear = db_to_linear(snr_db);
    if (snr_linear <= 1e-15) {
        throw std::invalid_argument("SNR too low or negative infinity");
    }
    return signal_power / snr_linear;
}

} // namespace ntn::signal
