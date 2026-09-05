#pragma once

#include "common/types.hpp"
#include <cstddef>
#include <stdexcept>

namespace ntn::signal {

// Calculates the average power of a complex baseband signal: P = (1/N) * sum(|s[n]|^2)
[[nodiscard]] double calculate_power(const ComplexVector& signal);

// Calculates the Root Mean Square (RMS) value of a signal: RMS = sqrt(P)
[[nodiscard]] double calculate_rms(const ComplexVector& signal);

// Calculates peak power: max(|s[n]|^2)
[[nodiscard]] double calculate_peak_power(const ComplexVector& signal);

// Calculates Peak-to-Average Power Ratio (PAPR) in linear scale
[[nodiscard]] double calculate_papr_linear(const ComplexVector& signal);

// Calculates Peak-to-Average Power Ratio (PAPR) in decibels (dB)
[[nodiscard]] double calculate_papr_db(const ComplexVector& signal);

// Normalizes a signal in-place so its average power equals target_power (default 1.0)
void normalize_power(ComplexVector& signal, double target_power = 1.0);

// Scales a signal by a complex scalar factor
void scale_signal(ComplexVector& signal, Complex scalar);

// Scales a signal by a real scalar factor
void scale_signal(ComplexVector& signal, double scalar);

// Calculates Error Vector Magnitude (EVM) RMS in percent:
// EVM_rms = sqrt( sum(|rx - tx|^2) / sum(|tx|^2) ) * 100%
[[nodiscard]] double calculate_evm_percent(const ComplexVector& reference, const ComplexVector& received);

// Calculates EVM in dB: 20 * log10(EVM_linear)
[[nodiscard]] double calculate_evm_db(const ComplexVector& reference, const ComplexVector& received);

// Computes required complex baseband noise variance (N0) for a given signal power and SNR (dB)
// Total noise variance sigma^2 = N0 = P_signal / 10^(SNR/10)
// For complex AWGN, real and imaginary components each have variance sigma^2 / 2
[[nodiscard]] double calculate_noise_variance(double signal_power, double snr_db);

} // namespace ntn::signal
