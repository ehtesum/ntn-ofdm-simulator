# NTN-OFDM Physical Layer Simulator

A modular C++20 physical-layer communications simulator modeling Doppler tracking, Carrier Frequency Offset (CFO) synchronization, and Inter-Carrier Interference (ICI) mitigation for 5G Non-Terrestrial Networks.

## Motivation

In Low Earth Orbit (LEO) satellite communications (500–1200 km altitude), orbital velocities (>7.5 km/s) induce high Doppler shifts (exceeding ±50 kHz at S-band) and rapid Doppler drift rates (up to -636 Hz/s at nadir). 

In standard 5G NR OFDM waveforms, these uncompensated frequency offsets destroy subcarrier orthogonality. The resulting Carrier Frequency Offset (CFO) and Inter-Carrier Interference (ICI) cause constellation rotation and high Bit Error Rates (BER), preventing reliable link establishment over satellite-to-ground channels.

This simulator models these physical-layer channel impairments and evaluates time-domain tracking and frequency-domain equalization algorithms.

## How It Works

1. **Transceiver Pipeline**:
   - Digital modulation (BPSK, QPSK, 16-QAM) with unit symbol energy normalization.
   - Radix-2 Cooley-Tukey FFT/IFFT with bit-reversal permutation.
   - 5G NR OFDM framing: 15 kHz subcarrier spacing, DC subcarrier nulling, and Cyclic Prefix (CP) insertion.
2. **Channel Impairments**:
   - Free-Space Path Loss (FSPL) and deterministic circular AWGN thermal noise.
   - Dynamic orbital slant-range Doppler shifts and oscillator frequency offsets.
3. **Compensation & Recovery**:
   - Blind Cyclic Prefix correlation tracking for fractional and integer CFO estimation.
   - Banded frequency-domain ICI equalization to suppress adjacent subcarrier spectral leakage.

## Benchmark Results

Simulation across 168,000 transmitted bits:

| Scenario | Channel Impairments | Uncompensated BER | Compensated BER | RMS EVM |
| :--- | :--- | :---: | :---: | :---: |
| **A. Back-to-Back** | None | 0.000 | 0.000 | 0.00% |
| **B. AWGN Only** | 15 dB SNR thermal noise | 0.000 | 0.000 | 13.63% |
| **C. Static CFO** | AWGN + 1.8 kHz oscillator offset | 0.501 | 0.000 | 142.02% |
| **D. Dynamic Doppler** | AWGN + LEO orbital trajectory | 0.481 | 0.000 | 128.07% |
| **E. Combined Channel** | AWGN + CFO + dynamic Doppler | 0.499 | 0.0014 | 141.75% |
| **F. Synchronized & Equalized** | Full CP tracking loop + ICI equalizer | 0.499 | 0.0014 | 37.03% |

## Build

Requirements: C++20 compliant compiler (GCC 10+, Clang 11+, or MSVC 2019+) and CMake 3.16+.

```bash
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

## Running Experiments

Run individual physical-layer simulations from `build/experiments/`:

```bash
# Carrier Frequency Offset tracking experiment
./build/experiments/run_cfo_experiment

# LEO Doppler trajectory simulation
./build/experiments/run_doppler_experiment

# Inter-Carrier Interference mitigation test
./build/experiments/run_ici_experiment

# Full channel scenario matrix benchmark
./build/experiments/run_scenario_matrix
```
