#pragma once

#include "common/types.hpp"
#include "ofdm/ofdm_config.hpp"
#include <cstddef>
#include <vector>

namespace ntn::channel {

// Physical orbit and carrier parameters for LEO NTN satellite scenarios
struct LEOOrbitParams {
    double altitude_km{600.0};              // Satellite orbital altitude (e.g. 600 km)
    double carrier_frequency_hz{2.0e9};     // RF Carrier frequency (2.0 GHz S-band)
    double min_slant_range_km{600.0};       // Slant range at closest approach / nadir (km)
    double pass_duration_s{600.0};          // Full horizon-to-horizon pass duration (~10 minutes)

    // Computes circular Keplerian orbital velocity: v = sqrt(G*M_E / (R_E + h))
    [[nodiscard]] double orbital_velocity_mps() const noexcept;

    // Computes maximum Doppler shift at horizon: f_d_max = (v / c) * f_c
    [[nodiscard]] double max_doppler_shift_hz() const noexcept;

    // Computes maximum Doppler rate (frequency slope) at nadir: df_d/dt = -v^2 * f_c / (c * d_min)
    [[nodiscard]] double max_doppler_rate_hz_per_sec() const noexcept;
};

// Generates time-varying Doppler trajectory for a LEO satellite pass
class DopplerTrajectory {
public:
    explicit DopplerTrajectory(LEOOrbitParams params = LEOOrbitParams{});

    // Evaluates instantaneous Doppler shift f_d(t) in Hz at time t (seconds) relative to closest approach (t=0)
    [[nodiscard]] double instantaneous_doppler_hz(double t_sec) const noexcept;

    // Evaluates instantaneous Doppler rate df_d/dt in Hz/s at time t (seconds)
    [[nodiscard]] double instantaneous_doppler_rate_hz_per_sec(double t_sec) const noexcept;

    [[nodiscard]] const LEOOrbitParams& params() const noexcept { return params_; }

private:
    LEOOrbitParams params_;
};

// Channel impairment model applying time-varying or constant Doppler shift to complex baseband signal
class DopplerChannel {
public:
    // sampling_rate_hz: Baseband sampling rate Fs
    explicit DopplerChannel(double sampling_rate_hz);

    // Applies a constant Doppler shift to a signal vector in-place: r[n] = s[n] * exp(j * (2*pi*f_d*n/Fs + phase))
    void apply_constant_doppler(ComplexVector& signal, double doppler_hz);

    // Applies time-varying Doppler according to a trajectory starting at t_start_sec
    void apply_time_varying_doppler(ComplexVector& signal, const DopplerTrajectory& trajectory, double t_start_sec);

    // Resets phase accumulator
    void reset_phase(double initial_phase_rad = 0.0) noexcept { current_phase_ = initial_phase_rad; }

    [[nodiscard]] double current_phase() const noexcept { return current_phase_; }
    [[nodiscard]] double sampling_rate_hz() const noexcept { return sampling_rate_hz_; }

private:
    double sampling_rate_hz_;
    double current_phase_{0.0};
};

// Two-stage Doppler Tracker for NTN OFDM links:
// Stage 1: Coarse frequency offset correction (e.g. from orbital ephemeris or sync raster)
// Stage 2: Slot-by-slot fine phase and frequency tracking
class DopplerTracker {
public:
    explicit DopplerTracker(ofdm::OFDMConfig ofdm_config);

    // Updates tracker estimate from a received OFDM symbol
    // Returns estimated instantaneous frequency offset in Hz
    [[nodiscard]] double track_symbol(const ComplexVector& time_samples);

    // Compensates a block of time-domain samples using the currently tracked Doppler frequency and phase
    void compensate_block(ComplexVector& time_samples, double estimated_cfo_hz);

    // Resets tracker state
    void reset() noexcept;

    [[nodiscard]] double current_frequency_estimate_hz() const noexcept { return current_freq_est_hz_; }
    [[nodiscard]] double current_phase_rad() const noexcept { return current_phase_rad_; }

private:
    ofdm::OFDMConfig ofdm_config_;
    double current_freq_est_hz_{0.0};
    double current_phase_rad_{0.0};
    double alpha_{0.7}; // Exponential smoothing factor for tracker noise filtering
    bool initialized_{false};
};

} // namespace ntn::channel
