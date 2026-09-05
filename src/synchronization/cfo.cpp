#include "synchronization/cfo.hpp"
#include <cmath>
#include <stdexcept>

namespace ntn::synchronization {

// ----------------------------------------------------------------------
// CFO Channel
// ----------------------------------------------------------------------

CFOChannel::CFOChannel(double f_cfo_hz, double sampling_rate_hz, double initial_phase_rad)
    : f_cfo_hz_(f_cfo_hz),
      sampling_rate_hz_(sampling_rate_hz),
      current_phase_(initial_phase_rad) {
    if (sampling_rate_hz_ <= 0.0) {
        throw std::invalid_argument("Sampling rate must be positive.");
    }
}

void CFOChannel::apply_cfo(ComplexVector& signal) {
    if (signal.empty() || std::abs(f_cfo_hz_) < 1e-12) {
        return;
    }

    const double phase_step = TWO_PI * (f_cfo_hz_ / sampling_rate_hz_);
    double phase = current_phase_;

    for (auto& sample : signal) {
        // Complex phase rotator: exp(j * phase) = cos(phase) + j * sin(phase)
        Complex rotator(std::cos(phase), std::sin(phase));
        sample *= rotator;
        phase += phase_step;
    }

    // Keep phase bounded in [-pi, pi) to prevent numerical precision loss over time
    current_phase_ = std::fmod(phase, TWO_PI);
    if (current_phase_ > PI) current_phase_ -= TWO_PI;
    if (current_phase_ < -PI) current_phase_ += TWO_PI;
}

void CFOChannel::reset_phase(double initial_phase_rad) noexcept {
    current_phase_ = initial_phase_rad;
}

// ----------------------------------------------------------------------
// CFO Estimator (Cyclic Prefix Correlation)
// ----------------------------------------------------------------------

CFOEstimator::CFOEstimator(ofdm::OFDMConfig ofdm_config)
    : ofdm_config_(ofdm_config) {
}

double CFOEstimator::estimate_normalized_cfo(const ComplexVector& time_samples) const {
    const size_t sym_len = ofdm_config_.symbol_duration_samples();
    const size_t n_fft = ofdm_config_.fft_size;
    const size_t n_cp = ofdm_config_.cp_length;

    if (time_samples.size() < sym_len) {
        throw std::invalid_argument("Input time_samples length (" + std::to_string(time_samples.size()) +
                                    ") is smaller than symbol duration (" + std::to_string(sym_len) + ").");
    }

    // Cross-correlate the Cyclic Prefix (indices 0 ... cp-1) with the useful symbol tail (indices fft ... fft + cp - 1)
    Complex correlation(0.0, 0.0);
    for (size_t m = 0; m < n_cp; ++m) {
        Complex cp_sample = time_samples[m];
        Complex tail_sample = time_samples[m + n_fft];
        // tail * conj(cp)
        correlation += tail_sample * std::conj(cp_sample);
    }

    // Phase angle of the correlation metric
    double angle = std::arg(correlation); // in (-pi, pi]

    // epsilon = angle / (2 * pi)
    return angle / TWO_PI;
}

double CFOEstimator::estimate_cfo_hz(const ComplexVector& time_samples) const {
    double norm_cfo = estimate_normalized_cfo(time_samples);
    // f_cfo = epsilon * delta_f
    return norm_cfo * ofdm_config_.subcarrier_spacing_hz;
}

double CFOEstimator::estimate_cfo_hz_multisymbol(const ComplexVector& time_stream, size_t num_ofdm_symbols) const {
    const size_t sym_len = ofdm_config_.symbol_duration_samples();
    const size_t n_fft = ofdm_config_.fft_size;
    const size_t n_cp = ofdm_config_.cp_length;

    if (time_stream.size() < num_ofdm_symbols * sym_len) {
        throw std::invalid_argument("Time stream length is insufficient for requested OFDM symbols.");
    }

    Complex total_correlation(0.0, 0.0);
    for (size_t s = 0; s < num_ofdm_symbols; ++s) {
        const size_t offset = s * sym_len;
        for (size_t m = 0; m < n_cp; ++m) {
            Complex cp_sample = time_stream[offset + m];
            Complex tail_sample = time_stream[offset + m + n_fft];
            total_correlation += tail_sample * std::conj(cp_sample);
        }
    }

    double angle = std::arg(total_correlation);
    double norm_cfo = angle / TWO_PI;
    return norm_cfo * ofdm_config_.subcarrier_spacing_hz;
}

// ----------------------------------------------------------------------
// CFO Compensator
// ----------------------------------------------------------------------

CFOCompensator::CFOCompensator(double estimated_cfo_hz, double sampling_rate_hz)
    : estimated_cfo_hz_(estimated_cfo_hz),
      sampling_rate_hz_(sampling_rate_hz) {
    if (sampling_rate_hz_ <= 0.0) {
        throw std::invalid_argument("Sampling rate must be positive.");
    }
}

void CFOCompensator::compensate(ComplexVector& signal) {
    if (signal.empty() || std::abs(estimated_cfo_hz_) < 1e-12) {
        return;
    }

    // Negative phase step to counteract frequency offset: exp(-j * 2*pi*f_cfo*n/Fs)
    const double phase_step = -TWO_PI * (estimated_cfo_hz_ / sampling_rate_hz_);
    double phase = current_phase_;

    for (auto& sample : signal) {
        Complex de_rotator(std::cos(phase), std::sin(phase));
        sample *= de_rotator;
        phase += phase_step;
    }

    current_phase_ = std::fmod(phase, TWO_PI);
    if (current_phase_ > PI) current_phase_ -= TWO_PI;
    if (current_phase_ < -PI) current_phase_ += TWO_PI;
}

void CFOCompensator::reset_phase(double initial_phase_rad) noexcept {
    current_phase_ = initial_phase_rad;
}

} // namespace ntn::synchronization
