#pragma once

#include "common/types.hpp"
#include "ofdm/ofdm_config.hpp"
#include <cstddef>
#include <vector>

namespace ntn::synchronization {

// Models Carrier Frequency Offset (CFO) impairment in complex baseband
class CFOChannel {
public:
    // f_cfo_hz: frequency offset in Hertz
    // sampling_rate_hz: baseband sampling rate (Fs)
    // initial_phase_rad: initial phase offset at sample n=0
    explicit CFOChannel(double f_cfo_hz, double sampling_rate_hz, double initial_phase_rad = 0.0);

    // Applies CFO phase rotation in-place: r[n] = s[n] * exp(j * (2*pi*f_cfo*n/Fs + phase))
    // Updates internal phase accumulator to maintain phase continuity across multiple symbols/frames
    void apply_cfo(ComplexVector& signal);

    // Resets internal phase accumulator to initial phase
    void reset_phase(double initial_phase_rad = 0.0) noexcept;

    [[nodiscard]] double f_cfo_hz() const noexcept { return f_cfo_hz_; }
    [[nodiscard]] double normalized_cfo(double subcarrier_spacing_hz) const noexcept {
        return f_cfo_hz_ / subcarrier_spacing_hz;
    }
    [[nodiscard]] double current_phase() const noexcept { return current_phase_; }

    void set_f_cfo_hz(double f_cfo_hz) noexcept { f_cfo_hz_ = f_cfo_hz; }

private:
    double f_cfo_hz_;
    double sampling_rate_hz_;
    double current_phase_;
};

// Estimates fractional Carrier Frequency Offset using Cyclic Prefix correlation (Moose / van de Beek method)
class CFOEstimator {
public:
    explicit CFOEstimator(ofdm::OFDMConfig ofdm_config);

    // Estimates fractional CFO epsilon = f_cfo / delta_f from a single OFDM symbol.
    // Operating range: |epsilon| <= 0.5 (unambiguous range of CP correlation).
    // time_samples must contain ofdm_config.symbol_duration_samples() samples.
    [[nodiscard]] double estimate_normalized_cfo(const ComplexVector& time_samples) const;

    // Estimates CFO in Hertz from a single OFDM symbol
    [[nodiscard]] double estimate_cfo_hz(const ComplexVector& time_samples) const;

    // Estimates CFO averaged across multiple OFDM symbols for enhanced noise resilience
    [[nodiscard]] double estimate_cfo_hz_multisymbol(const ComplexVector& time_stream, size_t num_ofdm_symbols) const;

    [[nodiscard]] const ofdm::OFDMConfig& config() const noexcept { return ofdm_config_; }

private:
    ofdm::OFDMConfig ofdm_config_;
};

// Compensates Carrier Frequency Offset by de-rotating baseband samples
class CFOCompensator {
public:
    explicit CFOCompensator(double estimated_cfo_hz, double sampling_rate_hz);

    // Applies de-rotation in-place: r_comp[n] = r[n] * exp(-j * (2*pi*f_cfo*n/Fs + phase))
    // Maintains phase continuity across successive calls
    void compensate(ComplexVector& signal);

    // Resets phase accumulator
    void reset_phase(double initial_phase_rad = 0.0) noexcept;

    void set_estimated_cfo_hz(double cfo_hz) noexcept { estimated_cfo_hz_ = cfo_hz; }
    [[nodiscard]] double estimated_cfo_hz() const noexcept { return estimated_cfo_hz_; }

private:
    double estimated_cfo_hz_;
    double sampling_rate_hz_;
    double current_phase_{0.0};
};

} // namespace ntn::synchronization
