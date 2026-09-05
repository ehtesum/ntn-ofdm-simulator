#include "signal/fft.hpp"
#include <algorithm>
#include <cmath>

namespace ntn::signal {

void FFT::forward(ComplexVector& data) {
    if (data.empty()) {
        return;
    }
    if (!is_power_of_two(data.size())) {
        throw std::invalid_argument("FFT size must be a power of 2.");
    }
    transform_internal(data, false);
}

void FFT::inverse(ComplexVector& data) {
    if (data.empty()) {
        return;
    }
    size_t n = data.size();
    if (!is_power_of_two(n)) {
        throw std::invalid_argument("IFFT size must be a power of 2.");
    }
    transform_internal(data, true);

    // Normalize inverse transform by 1/N
    const double inv_n = 1.0 / static_cast<double>(n);
    for (auto& sample : data) {
        sample *= inv_n;
    }
}

void FFT::bit_reverse_reorder(ComplexVector& data) {
    const size_t n = data.size();
    size_t j = 0;
    for (size_t i = 0; i < n - 1; ++i) {
        if (i < j) {
            std::swap(data[i], data[j]);
        }
        size_t k = n >> 1;
        while (k <= j) {
            j -= k;
            k >>= 1;
        }
        j += k;
    }
}

void FFT::transform_internal(ComplexVector& data, bool inverse) {
    const size_t n = data.size();
    if (n <= 1) {
        return;
    }

    bit_reverse_reorder(data);

    // Cooley-Tukey butterfly stages
    for (size_t len = 2; len <= n; len <<= 1) {
        const double angle = (inverse ? TWO_PI : -TWO_PI) / static_cast<double>(len);
        const Complex wlen(std::cos(angle), std::sin(angle));

        for (size_t i = 0; i < n; i += len) {
            Complex w(1.0, 0.0);
            const size_t half_len = len >> 1;
            for (size_t j = 0; j < half_len; ++j) {
                Complex u = data[i + j];
                Complex v = data[i + j + half_len] * w;
                data[i + j] = u + v;
                data[i + j + half_len] = u - v;
                w *= wlen;
            }
        }
    }
}

void FFT::fftshift(ComplexVector& data) {
    if (data.empty()) {
        return;
    }
    const size_t n = data.size();
    const size_t pivot = (n + 1) / 2;
    std::rotate(data.begin(), data.begin() + pivot, data.end());
}

void FFT::ifftshift(ComplexVector& data) {
    if (data.empty()) {
        return;
    }
    const size_t n = data.size();
    const size_t pivot = n / 2;
    std::rotate(data.begin(), data.begin() + pivot, data.end());
}

} // namespace ntn::signal
