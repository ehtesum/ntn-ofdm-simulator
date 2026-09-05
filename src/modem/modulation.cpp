#include "modem/modulation.hpp"
#include <algorithm>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace ntn::modem {

std::string to_string(ModulationType type) {
    switch (type) {
        case ModulationType::BPSK:  return "BPSK";
        case ModulationType::QPSK:  return "QPSK";
        case ModulationType::QAM16: return "16-QAM";
    }
    return "UNKNOWN";
}

ModulationType from_string(std::string_view name) {
    std::string lower;
    lower.reserve(name.size());
    for (char c : name) {
        lower.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
    }
    if (lower == "bpsk") return ModulationType::BPSK;
    if (lower == "qpsk") return ModulationType::QPSK;
    if (lower == "16qam" || lower == "16-qam" || lower == "qam16") return ModulationType::QAM16;
    throw std::invalid_argument("Unknown modulation type: " + std::string(name));
}

// ----------------------------------------------------------------------
// BPSK Implementation
// ----------------------------------------------------------------------

ComplexVector BPSKModulator::modulate(const ByteVector& bits) const {
    ComplexVector symbols;
    symbols.reserve(bits.size());
    for (uint8_t bit : bits) {
        // '0' -> +1.0, '1' -> -1.0
        double re = (bit == 0) ? 1.0 : -1.0;
        symbols.emplace_back(re, 0.0);
    }
    return symbols;
}

ByteVector BPSKModulator::demodulate(const ComplexVector& symbols) const {
    ByteVector bits;
    bits.reserve(symbols.size());
    for (const auto& sym : symbols) {
        // Decision boundary at Re(s) = 0
        bits.push_back((sym.real() >= 0.0) ? 0 : 1);
    }
    return bits;
}

ComplexVector BPSKModulator::reference_constellation() const {
    return { Complex(1.0, 0.0), Complex(-1.0, 0.0) };
}

// ----------------------------------------------------------------------
// QPSK Implementation (Gray Coded)
// ----------------------------------------------------------------------

ComplexVector QPSKModulator::modulate(const ByteVector& bits) const {
    if (bits.size() % 2 != 0) {
        throw std::invalid_argument("Bit vector size must be a multiple of 2 for QPSK.");
    }
    ComplexVector symbols;
    symbols.reserve(bits.size() / 2);

    for (size_t i = 0; i < bits.size(); i += 2) {
        // Gray mapping:
        // b0=0 -> +1, b0=1 -> -1
        // b1=0 -> +1, b1=1 -> -1
        double i_val = (bits[i] == 0) ? 1.0 : -1.0;
        double q_val = (bits[i + 1] == 0) ? 1.0 : -1.0;
        symbols.emplace_back(i_val * NORMALIZATION_FACTOR, q_val * NORMALIZATION_FACTOR);
    }
    return symbols;
}

ByteVector QPSKModulator::demodulate(const ComplexVector& symbols) const {
    ByteVector bits;
    bits.reserve(symbols.size() * 2);

    for (const auto& sym : symbols) {
        // Hard-decision boundaries along I and Q axes:
        // Re >= 0 -> b0 = 0, Re < 0 -> b0 = 1
        // Im >= 0 -> b1 = 0, Im < 0 -> b1 = 1
        bits.push_back((sym.real() >= 0.0) ? 0 : 1);
        bits.push_back((sym.imag() >= 0.0) ? 0 : 1);
    }
    return bits;
}

ComplexVector QPSKModulator::reference_constellation() const {
    const double k = NORMALIZATION_FACTOR;
    return {
        Complex( k,  k), // 00
        Complex(-k,  k), // 10
        Complex(-k, -k), // 11
        Complex( k, -k)  // 01
    };
}

// ----------------------------------------------------------------------
// 16-QAM Implementation (Gray Coded)
// ----------------------------------------------------------------------

namespace {

// 2-bit Gray mapping per dimension:
// 00 -> +3, 01 -> +1, 11 -> -1, 10 -> -3
inline double map_2bits_to_amplitude(uint8_t b0, uint8_t b1) {
    if (b0 == 0 && b1 == 0) return  3.0;
    if (b0 == 0 && b1 == 1) return  1.0;
    if (b0 == 1 && b1 == 1) return -1.0;
    return -3.0; // 10
}

// Hard slicing along 1D dimension with decision boundaries at -2, 0, +2
inline void slice_amplitude_to_2bits(double val, uint8_t& b0, uint8_t& b1) {
    if (val > 2.0) {
        b0 = 0; b1 = 0; // +3
    } else if (val > 0.0) {
        b0 = 0; b1 = 1; // +1
    } else if (val > -2.0) {
        b0 = 1; b1 = 1; // -1
    } else {
        b0 = 1; b1 = 0; // -3
    }
}

} // namespace

ComplexVector QAM16Modulator::modulate(const ByteVector& bits) const {
    if (bits.size() % 4 != 0) {
        throw std::invalid_argument("Bit vector size must be a multiple of 4 for 16-QAM.");
    }
    ComplexVector symbols;
    symbols.reserve(bits.size() / 4);

    for (size_t i = 0; i < bits.size(); i += 4) {
        double i_amp = map_2bits_to_amplitude(bits[i], bits[i + 1]);
        double q_amp = map_2bits_to_amplitude(bits[i + 2], bits[i + 3]);
        symbols.emplace_back(i_amp * NORMALIZATION_FACTOR, q_amp * NORMALIZATION_FACTOR);
    }
    return symbols;
}

ByteVector QAM16Modulator::demodulate(const ComplexVector& symbols) const {
    ByteVector bits;
    bits.reserve(symbols.size() * 4);

    // Invert normalization: unscale by 1 / NORMALIZATION_FACTOR = sqrt(10)
    const double inv_norm = 1.0 / NORMALIZATION_FACTOR; // sqrt(10) ~ 3.16227766

    for (const auto& sym : symbols) {
        double i_val = sym.real() * inv_norm;
        double q_val = sym.imag() * inv_norm;

        uint8_t b0, b1, b2, b3;
        slice_amplitude_to_2bits(i_val, b0, b1);
        slice_amplitude_to_2bits(q_val, b2, b3);

        bits.push_back(b0);
        bits.push_back(b1);
        bits.push_back(b2);
        bits.push_back(b3);
    }
    return bits;
}

ComplexVector QAM16Modulator::reference_constellation() const {
    ComplexVector constell;
    constell.reserve(16);
    const double k = NORMALIZATION_FACTOR;
    const double levels[4] = {-3.0, -1.0, 1.0, 3.0};
    for (double i_lvl : levels) {
        for (double q_lvl : levels) {
            constell.emplace_back(i_lvl * k, q_lvl * k);
        }
    }
    return constell;
}

// ----------------------------------------------------------------------
// Modulator Factory
// ----------------------------------------------------------------------

std::unique_ptr<Modulator> create_modulator(ModulationType type) {
    switch (type) {
        case ModulationType::BPSK:
            return std::make_unique<BPSKModulator>();
        case ModulationType::QPSK:
            return std::make_unique<QPSKModulator>();
        case ModulationType::QAM16:
            return std::make_unique<QAM16Modulator>();
    }
    throw std::invalid_argument("Unsupported modulation type");
}

} // namespace ntn::modem
