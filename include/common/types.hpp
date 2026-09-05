#pragma once

#include <complex>
#include <vector>
#include <cstdint>
#include <numbers>
#include <cmath>

namespace ntn {

// Primary complex sample precision for simulation fidelity
using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;
using RealVector = std::vector<double>;
using ByteVector = std::vector<uint8_t>;

// Physical constants
inline constexpr double PI = 3.141592653589793238462643383279502884;
inline constexpr double TWO_PI = 2.0 * PI;
inline constexpr double SPEED_OF_LIGHT = 299792458.0; // m/s
inline constexpr double BOLTZMANN_CONSTANT = 1.380649e-23; // J/K

// Mathematical utilities
[[nodiscard]] constexpr double db_to_linear(double db) noexcept {
    return std::pow(10.0, db / 10.0);
}

[[nodiscard]] inline double linear_to_db(double linear) noexcept {
    return (linear > 1e-30) ? 10.0 * std::log10(linear) : -300.0;
}

[[nodiscard]] constexpr double dbm_to_watts(double dbm) noexcept {
    return std::pow(10.0, (dbm - 30.0) / 10.0);
}

[[nodiscard]] inline double watts_to_dbm(double watts) noexcept {
    return (watts > 1e-30) ? 10.0 * std::log10(watts) + 30.0 : -270.0;
}

} // namespace ntn
