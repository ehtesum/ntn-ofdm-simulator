#include "channel/doppler.hpp"
#include <cmath>
#include <stdexcept>

namespace ntn::channel {

// Standard Earth gravitational parameter and equatorial radius
inline constexpr double GM_EARTH = 3.986004418e14; // m^3 / s^2
inline constexpr double R_EARTH_M = 6371.0e3;        // meters

// ----------------------------------------------------------------------
// LEO Orbit Parameters & Physics
// ----------------------------------------------------------------------

double LEOOrbitParams::orbital_velocity_mps() const noexcept {
    const double r = R_EARTH_M + (altitude_km * 1.0e3);
    return std::sqrt(GM_EARTH / r);
}

double LEOOrbitParams::max_doppler_shift_hz() const noexcept {
    const double v = orbital_velocity_mps();
    return (v / SPEED_OF_LIGHT) * carrier_frequency_hz;
}

double LEOOrbitParams::max_doppler_rate_hz_per_sec() const noexcept {
    const double v = orbital_velocity_mps();
    const double d_min = min_slant_range_km * 1.0e3;
    return -(v * v * carrier_frequency_hz) / (SPEED_OF_LIGHT * d_min);
}

// ----------------------------------------------------------------------
// Doppler Trajectory (S-Curve)
// ----------------------------------------------------------------------

DopplerTrajectory::DopplerTrajectory(LEOOrbitParams params)
    : params_(params) {
}

double DopplerTrajectory::instantaneous_doppler_hz(double t_sec) const noexcept {
    const double v = params_.orbital_velocity_mps();
    const double d_min = params_.min_slant_range_km * 1.0e3;
    const double c = SPEED_OF_LIGHT;
    const double fc = params_.carrier_frequency_hz;

    // Slant range rho(t) = sqrt(d_min^2 + (v * t)^2)
    // Radial velocity v_r(t) = (v^2 * t) / rho(t)
    // Doppler shift f_d(t) = - (v_r(t) / c) * fc
    double vt = v * t_sec;
    double rho = std::sqrt(d_min * d_min + vt * vt);
    double v_radial = (v * vt) / rho;

    return -(v_radial / c) * fc;
}

double DopplerTrajectory::instantaneous_doppler_rate_hz_per_sec(double t_sec) const noexcept {
    const double v = params_.orbital_velocity_mps();
    const double d_min = params_.min_slant_range_km * 1.0e3;
    const double c = SPEED_OF_LIGHT;
    const double fc = params_.carrier_frequency_hz;

    double vt = v * t_sec;
    double rho2 = d_min * d_min + vt * vt;
    double rho3 = rho2 * std::sqrt(rho2);

    return -(v * v * d_min * d_min * fc) / (c * rho3);
}

// ----------------------------------------------------------------------
// Doppler Channel
// ----------------------------------------------------------------------

DopplerChannel::DopplerChannel(double sampling_rate_hz)
    : sampling_rate_hz_(sampling_rate_hz) {
    if (sampling_rate_hz_ <= 0.0) {
        throw std::invalid_argument("Sampling rate must be positive.");
    }
}

void DopplerChannel::apply_constant_doppler(ComplexVector& signal, double doppler_hz) {
    if (signal.empty() || std::abs(doppler_hz) < 1e-12) {
        return;
    }
    const double phase_step = TWO_PI * (doppler_hz / sampling_rate_hz_);
    double phase = current_phase_;

    for (auto& sample : signal) {
        sample *= Complex(std::cos(phase), std::sin(phase));
        phase += phase_step;
    }

    current_phase_ = std::fmod(phase, TWO_PI);
    if (current_phase_ > PI) current_phase_ -= TWO_PI;
    if (current_phase_ < -PI) current_phase_ += TWO_PI;
}

void DopplerChannel::apply_time_varying_doppler(ComplexVector& signal, const DopplerTrajectory& trajectory, double t_start_sec) {
    if (signal.empty()) {
        return;
    }
    const double dt = 1.0 / sampling_rate_hz_;
    double t = t_start_sec;
    double phase = current_phase_;

    for (auto& sample : signal) {
        double f_d = trajectory.instantaneous_doppler_hz(t);
        double dphi = TWO_PI * f_d * dt;
        phase += dphi;
        sample *= Complex(std::cos(phase), std::sin(phase));
        t += dt;
    }

    current_phase_ = std::fmod(phase, TWO_PI);
    if (current_phase_ > PI) current_phase_ -= TWO_PI;
    if (current_phase_ < -PI) current_phase_ += TWO_PI;
}

// ----------------------------------------------------------------------
// Doppler Tracker
// ----------------------------------------------------------------------

DopplerTracker::DopplerTracker(ofdm::OFDMConfig ofdm_config)
    : ofdm_config_(ofdm_config) {
}

double DopplerTracker::track_symbol(const ComplexVector& time_samples) {
    const size_t n_fft = ofdm_config_.fft_size;
    const size_t n_cp = ofdm_config_.cp_length;
    const size_t sym_len = ofdm_config_.symbol_duration_samples();

    if (time_samples.size() < sym_len) {
        throw std::invalid_argument("Symbol length too short for Doppler tracking.");
    }

    // Cyclic Prefix correlation for instantaneous frequency measurement
    Complex correlation(0.0, 0.0);
    for (size_t m = 0; m < n_cp; ++m) {
        correlation += time_samples[m + n_fft] * std::conj(time_samples[m]);
    }

    double angle = std::arg(correlation);
    double instant_cfo_hz = (angle / TWO_PI) * ofdm_config_.subcarrier_spacing_hz;

    // First measurement immediately initializes the filter without transient lag
    if (!initialized_) {
        current_freq_est_hz_ = instant_cfo_hz;
        initialized_ = true;
    } else {
        // Subsequent measurements filtered using exponential moving average
        current_freq_est_hz_ = alpha_ * current_freq_est_hz_ + (1.0 - alpha_) * instant_cfo_hz;
    }

    return current_freq_est_hz_;
}

void DopplerTracker::compensate_block(ComplexVector& time_samples, double estimated_cfo_hz) {
    if (time_samples.empty() || std::abs(estimated_cfo_hz) < 1e-12) {
        return;
    }
    const double fs = ofdm_config_.sampling_rate_hz();
    const double phase_step = -TWO_PI * (estimated_cfo_hz / fs);
    double phase = current_phase_rad_;

    for (auto& sample : time_samples) {
        sample *= Complex(std::cos(phase), std::sin(phase));
        phase += phase_step;
    }

    current_phase_rad_ = std::fmod(phase, TWO_PI);
    if (current_phase_rad_ > PI) current_phase_rad_ -= TWO_PI;
    if (current_phase_rad_ < -PI) current_phase_rad_ += TWO_PI;
}

void DopplerTracker::reset() noexcept {
    current_freq_est_hz_ = 0.0;
    current_phase_rad_ = 0.0;
    initialized_ = false;
}

} // namespace ntn::channel
