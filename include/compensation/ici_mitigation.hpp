#pragma once

#include "common/types.hpp"
#include "ofdm/ofdm_config.hpp"
#include <vector>
#include <complex>

namespace ntn::compensation {

/**
 * @brief Frequency-domain ICI metrics structure.
 */
struct IciMetrics {
    double carrier_power{0.0};       ///< Desired subcarrier signal power
    double ici_power{0.0};           ///< Total inter-carrier interference power
    double cir_linear{0.0};          ///< Carrier-to-Interference Ratio (linear)
    double cir_db{0.0};              ///< Carrier-to-Interference Ratio (dB)
    std::vector<double> leakage_spectrum; ///< Power leakage profile around target subcarrier
};

/**
 * @brief Computes analytical and empirical ICI kernel values.
 *
 * S(d, epsilon) = (1/N) * (1 - e^{j 2pi (d + epsilon)}) / (1 - e^{j 2pi (d + epsilon)/N})
 */
class IciKernel {
public:
    /**
     * @brief Computes the leakage coefficient from subcarrier m to subcarrier k
     *        under normalized frequency offset epsilon.
     * @param k Target subcarrier index
     * @param m Interfering subcarrier index
     * @param N FFT size
     * @param epsilon Normalized frequency offset (f_offset / subcarrier_spacing)
     * @return Complex leakage coefficient S_{k,m}(epsilon)
     */
    static Complex leakage_coefficient(int k, int m, size_t N, double epsilon);

    /**
     * @brief Measures empirical ICI and CIR by transmitting an isolated active subcarrier
     *        through a frequency-offset channel and observing adjacent subcarrier leakage.
     * @param config OFDM configuration
     * @param test_subcarrier Index of the active subcarrier (e.g., center subcarrier)
     * @param epsilon Normalized frequency offset
     * @return IciMetrics structure with CIR and leakage spectrum
     */
    static IciMetrics measure_ici(const ofdm::OFDMConfig& config,
                                  int test_subcarrier,
                                  double epsilon);
};

/**
 * @brief Frequency-Domain Banded ICI Equalizer.
 *
 * Mitigates inter-carrier interference caused by residual frequency offset
 * using tridiagonal banded matrix inversion (Thomas algorithm, O(N) complexity).
 */
class IciEqualizer {
public:
    explicit IciEqualizer(const ofdm::OFDMConfig& config);

    /**
     * @brief Equalizes received frequency-domain subcarriers to remove ICI.
     * @param rx_subcarriers Received frequency-domain subcarriers (size N_fft)
     * @param estimated_epsilon Residual normalized frequency offset
     * @return Equalized frequency-domain subcarriers
     */
    ComplexVector equalize(const ComplexVector& rx_subcarriers, double estimated_epsilon) const;

    /**
     * @brief Mitigates Common Phase Error (CPE) across OFDM subcarriers.
     * @param rx_subcarriers Received subcarriers
     * @param pilot_indices Indices of known pilot subcarriers
     * @param pilot_symbols Known pilot symbol values
     * @return CPE-compensated subcarriers
     */
    ComplexVector compensate_cpe(const ComplexVector& rx_subcarriers,
                                 const std::vector<int>& pilot_indices,
                                 const ComplexVector& pilot_symbols) const;

private:
    ofdm::OFDMConfig config_;
};

} // namespace ntn::compensation
