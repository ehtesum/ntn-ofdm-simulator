#pragma once

#include "common/types.hpp"
#include <memory>
#include <string>
#include <string_view>

namespace ntn::modem {

enum class ModulationType {
    BPSK,
    QPSK,
    QAM16
};

[[nodiscard]] std::string to_string(ModulationType type);
[[nodiscard]] ModulationType from_string(std::string_view name);

// Abstract base class for digital modulators/demodulators
class Modulator {
public:
    virtual ~Modulator() = default;

    // Number of bits represented per complex symbol
    [[nodiscard]] virtual size_t bits_per_symbol() const noexcept = 0;

    // Type of modulation
    [[nodiscard]] virtual ModulationType type() const noexcept = 0;

    // Maps binary input stream (0 and 1 values) to complex baseband symbols
    // If bits.size() is not a multiple of bits_per_symbol(), an exception is thrown or padded
    [[nodiscard]] virtual ComplexVector modulate(const ByteVector& bits) const = 0;

    // Performs hard-decision demodulation (minimum Euclidean distance slicing)
    [[nodiscard]] virtual ByteVector demodulate(const ComplexVector& symbols) const = 0;

    // Returns all ideal reference constellation points (normalized to average power = 1.0)
    [[nodiscard]] virtual ComplexVector reference_constellation() const = 0;

    // Name identifier
    [[nodiscard]] virtual std::string name() const = 0;
};

// Concrete implementations:

// Binary Phase Shift Keying (1 bit/symbol: '0' -> +1, '1' -> -1)
class BPSKModulator final : public Modulator {
public:
    [[nodiscard]] size_t bits_per_symbol() const noexcept override { return 1; }
    [[nodiscard]] ModulationType type() const noexcept override { return ModulationType::BPSK; }
    [[nodiscard]] ComplexVector modulate(const ByteVector& bits) const override;
    [[nodiscard]] ByteVector demodulate(const ComplexVector& symbols) const override;
    [[nodiscard]] ComplexVector reference_constellation() const override;
    [[nodiscard]] std::string name() const override { return "BPSK"; }
};

// Quadrature Phase Shift Keying with Gray Coding (2 bits/symbol)
// Mapping: 00 -> (1+j)/sqrt(2), 01 -> (-1+j)/sqrt(2), 11 -> (-1-j)/sqrt(2), 10 -> (1-j)/sqrt(2)
class QPSKModulator final : public Modulator {
public:
    [[nodiscard]] size_t bits_per_symbol() const noexcept override { return 2; }
    [[nodiscard]] ModulationType type() const noexcept override { return ModulationType::QPSK; }
    [[nodiscard]] ComplexVector modulate(const ByteVector& bits) const override;
    [[nodiscard]] ByteVector demodulate(const ComplexVector& symbols) const override;
    [[nodiscard]] ComplexVector reference_constellation() const override;
    [[nodiscard]] std::string name() const override { return "QPSK"; }

    static constexpr double NORMALIZATION_FACTOR = 0.70710678118654752440; // 1 / sqrt(2)
};

// 16-Quadrature Amplitude Modulation with Gray Coding (4 bits/symbol)
// Normalized by 1 / sqrt(10) so average power is 1.0
class QAM16Modulator final : public Modulator {
public:
    [[nodiscard]] size_t bits_per_symbol() const noexcept override { return 4; }
    [[nodiscard]] ModulationType type() const noexcept override { return ModulationType::QAM16; }
    [[nodiscard]] ComplexVector modulate(const ByteVector& bits) const override;
    [[nodiscard]] ByteVector demodulate(const ComplexVector& symbols) const override;
    [[nodiscard]] ComplexVector reference_constellation() const override;
    [[nodiscard]] std::string name() const override { return "16-QAM"; }

    static constexpr double NORMALIZATION_FACTOR = 0.31622776601683793320; // 1 / sqrt(10)
};

// Factory function
[[nodiscard]] std::unique_ptr<Modulator> create_modulator(ModulationType type);

} // namespace ntn::modem
