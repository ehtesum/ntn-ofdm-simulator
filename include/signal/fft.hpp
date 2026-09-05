#pragma once

#include "common/types.hpp"
#include <cstddef>
#include <stdexcept>

namespace ntn::signal {

class FFT {
public:
    // Computes forward discrete Fourier transform in-place (unnormalized):
    // X[k] = sum_{n=0}^{N-1} x[n] * exp(-j * 2 * pi * k * n / N)
    static void forward(ComplexVector& data);

    // Computes inverse discrete Fourier transform in-place (normalized by 1/N):
    // x[n] = (1/N) * sum_{k=0}^{N-1} X[k] * exp(+j * 2 * pi * k * n / N)
    static void inverse(ComplexVector& data);

    // Out-of-place forward FFT
    [[nodiscard]] static ComplexVector forward_copy(const ComplexVector& data) {
        ComplexVector copy = data;
        forward(copy);
        return copy;
    }

    // Out-of-place inverse FFT
    [[nodiscard]] static ComplexVector inverse_copy(const ComplexVector& data) {
        ComplexVector copy = data;
        inverse(copy);
        return copy;
    }

    // Shifts the zero-frequency component (DC) to the center of the spectrum
    static void fftshift(ComplexVector& data);

    // Inverse of fftshift
    static void ifftshift(ComplexVector& data);

    // Helper to check if a number is a positive power of 2
    [[nodiscard]] static constexpr bool is_power_of_two(size_t n) noexcept {
        return (n > 0) && ((n & (n - 1)) == 0);
    }

private:
    // Core Radix-2 Cooley-Tukey implementation with bit reversal
    static void transform_internal(ComplexVector& data, bool inverse);
    static void bit_reverse_reorder(ComplexVector& data);
};

} // namespace ntn::signal
