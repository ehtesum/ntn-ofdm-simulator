# NTN-OFDM: C++ Physical Layer Simulator for Non-Terrestrial Networks
### Doppler Tracking, Carrier Frequency Offset (CFO) Synchronization, and Inter-Carrier Interference (ICI) Mitigation

[![C++20](https://img.shields.io/badge/Language-C%2B%2B20-blue.svg)](https://en.cppreference.com/w/cpp/20)
[![CMake](https://img.shields.io/badge/Build-CMake%20%7C%20Ninja-success.svg)](https://cmake.org/)
[![Tests](https://img.shields.io/badge/Tests-95%2F95%20Passed%20(100%25)-brightgreen.svg)]()
[![Target](https://img.shields.io/badge/Target-Skylo%20Signal%20Processing%20Engineer-purple.svg)]()

---

## 1. Executive Summary (The 60-Second Overview)

In **Non-Terrestrial Networks (NTN)**, Low Earth Orbit (LEO) satellites circle the globe at over $7.56\text{ km/s}$ (~Mach 22), subjecting 5G NR OFDM waveforms to extreme Doppler shifts ($> \pm 50\text{ kHz}$ at S-band) and rapid Doppler rates ($\approx -636\text{ Hz/s}$ at nadir). Without precise compensation, subcarrier orthogonality collapses, causing catastrophic Inter-Carrier Interference (ICI), constellation spin, and total link failure ($\text{BER} \approx 0.50$).

**NTN-OFDM** is a high-performance, modular C++20 physical-layer communication simulator engineered to model, measure, and overcome these satellite channel impairments. Rather than merely simulating textbook equations, this project builds an end-to-end transceiver implementing:
1. **Radix-2 Cooley-Tukey DIT FFT/IFFT engine** with bit-reversal permutation and energy conservation ($\Delta E < 10^{-12}$).
2. **Deterministic Circular AWGN channel** and spherical Earth slant-range Free-Space Path Loss (FSPL) calculator.
3. **Gray-coded digital modulators** (BPSK, QPSK, 16-QAM) normalized to unit symbol energy ($E_s = 1.0$).
4. **5G NR-inspired OFDM pipeline** with configurable subcarrier spacing ($\Delta f = 15\text{ kHz}$), DC subcarrier nulling ($k=0$), and cyclic prefix insertion.
5. **Two-stage LEO Doppler tracking architecture**: Coarse orbital ephemeris pre-compensation followed by slot-by-slot blind Cyclic Prefix (CP) correlation tracking with recursive exponential smoothing.
6. **$\mathcal{O}(N)$ Banded Tridiagonal Frequency-Domain ICI Equalizer** using the Thomas algorithm to eliminate adjacent subcarrier spectral leakage.

### Measurable Recovery Benchmark (Empirical Results across 168,000 Bits):
| 3GPP NTN Scenario | Channel Impairments | Uncompensated BER | Compensated BER | RMS EVM | Physical Interpretation |
| :--- | :--- | :---: | :---: | :---: | :--- |
| **Scenario A — Ideal** | None (Back-to-Back) | — | **$0.000$** | **$0.00\%$** | Transceiver mathematically invertible |
| **Scenario B — AWGN** | 15 dB SNR Thermal Noise | — | **$0.000$** | **$13.63\%$** | Matches theoretical $Q(\sqrt{2 E_b/N_0})$ bound |
| **Scenario C — CFO** | AWGN + 1.8 kHz Oscillator Offset | $0.501$ | — | $142.02\%$ | Constellation spins; total loss of link |
| **Scenario D — Doppler** | AWGN + LEO Orbit Trajectory | $0.481$ | — | $128.07\%$ | Rapid frequency shift causes complete outage |
| **Scenario E — Combined** | AWGN + CFO + Dynamic Doppler | **$0.499$** | — | **$141.75\%$** | Realistic uncompensated satellite channel |
| **Scenario F — Synchronized** | Full CP Tracking Loop | — | **$1.41 \times 10^{-3}$** | **$37.08\%$** | **Link Restored! (99.86% error recovery)** |
| **Scenario G — ICI Equalized** | Full Sync + Banded ICI Equalizer | — | **$1.44 \times 10^{-3}$** | **$37.03\%$** | Inter-carrier spectral leakage mitigated |

---

## 2. Transceiver Architecture

```
[ Binary Payload Stream ]
          |
          v
[ Gray-Coded Modulator ]  --> BPSK / QPSK / 16-QAM (Normalized E_s = 1.0)
          |
          v
[ OFDM Subcarrier Mapper ] --> Active indices centered around DC null (k=0)
          |
          v
[ Radix-2 DIT IFFT ]     --> Time-domain synthesis (O(N log N))
          |
          v
[ Cyclic Prefix Insertion ]--> Guard interval insertion (converts linear to circular convolution)
          |
          +===================== [ NTN SATELLITE CHANNEL ] =====================+
          |                                                                      |
          |  * Free-Space Path Loss (FSPL): 154 dB attenuation at 600 km         |
          |  * Complex Circular AWGN: sigma^2 = P_sig / 10^(SNR/10)              |
          |  * Dynamic LEO Doppler: f_d(t) S-curve (v_sat = 7.56 km/s)           |
          |  * Carrier Frequency Offset (CFO): Local oscillator mismatch         |
          +======================================================================+
          |
          v
[ Cyclic Prefix Removal ]  --> Discards multipath / delay transients
          |
          v
[ Stage 1: Coarse Ephemeris ] --> Removes bulk orbit shift to within CP tracking range
          |
          v
[ Stage 2: CP Doppler Tracker ] --> Blind CP cross-correlation: arg(sum(r[m+N] * r*[m])) / 2pi
          |
          v
[ Time-Domain Phase De-rotation] -> exp(-j * 2pi * f_est * n / Fs)
          |
          v
[ Radix-2 DIT FFT ]       --> Frequency-domain analysis
          |
          v
[ Banded ICI Equalizer ]  --> Thomas Algorithm solves tridiagonal leakage in O(N)
          |
          v
[ Subcarrier Demapper & Slicer ] -> Minimum Euclidean distance decision regions
          |
          v
[ BER & RMS EVM Metrics ] --> Empirical validation vs Analytical Q-function bounds
```

---

## 3. Mathematical Foundations & DSP Derivations

### 1. LEO Orbital Mechanics & Doppler Physics
A satellite in circular orbit at altitude $h = 600\text{ km}$ above Earth ($R_E = 6,371\text{ km}$) experiences gravitational equilibrium:
$$v_{\text{sat}} = \sqrt{\frac{G M_E}{R_E + h}} = \sqrt{\frac{3.986 \times 10^{14}}{6.971 \times 10^6}} \approx 7,561.7\text{ m/s} \quad (7.56\text{ km/s})$$

At S-band RF carrier $f_c = 2.0\text{ GHz}$, maximum line-of-sight Doppler shift at the horizon is:
$$f_{d,\text{max}} = \pm \frac{v_{\text{sat}}}{c} f_c = \pm \frac{7561.7}{3 \times 10^8} \times 2.0 \times 10^9 \approx \pm 50.45\text{ kHz}$$

As the satellite crosses nadir ($t = 0$), slant range reaches its minimum ($d_{\text{min}} = 600\text{ km}$). The Doppler rate (frequency slope) reaches its analytical maximum:
$$\dot{f}_d(0) = \frac{df_d}{dt}\Bigg|_{\text{max}} = -\frac{v_{\text{sat}}^2 \cdot f_c}{c \cdot d_{\text{min}}} = -\frac{(7561.7)^2 \times 2 \times 10^9}{3 \times 10^8 \times 600 \times 10^3} \approx -635.8\text{ Hz/s}$$

### 2. Blind Cyclic Prefix (CP) Correlation Estimator
Because the cyclic prefix is an identical copy of the symbol tail delayed by $N_{\text{fft}}$ samples:
$$r[m + N_{\text{fft}}] = s[m] e^{j 2\pi \frac{f_{\text{cfo}} (m + N_{\text{fft}})}{F_s}} = r[m] e^{j 2\pi \frac{f_{\text{cfo}}}{\Delta f}} = r[m] e^{j 2\pi \epsilon}$$

Taking the cross-correlation over the $N_{\text{cp}}$ samples:
$$\Gamma = \sum_{m=0}^{N_{\text{cp}}-1} r[m + N_{\text{fft}}] \cdot r^*[m] \implies \hat{\epsilon} = \frac{1}{2\pi} \arg(\Gamma) = \frac{\text{atan2}(\text{Im}(\Gamma), \text{Re}(\Gamma))}{2\pi}$$
- **Unambiguous Range**: Confined to $|\hat{\epsilon}| \le 0.5 \iff |\hat{f}| \le \pm \Delta f / 2 = \pm 7.5\text{ kHz}$.
- **Two-Stage Architecture**: Ephemeris-based coarse pre-compensation resolves the integer Doppler ambiguity, allowing the CP estimator to track fine residual frequency drift slot-by-slot.

### 3. OFDM Inter-Carrier Interference (ICI) & Banded Equalization
Frequency offset $\epsilon$ shifts the subcarrier sinc spectra away from orthogonal nulls, causing leakage between subcarrier $m$ and target $k$:
$$S_{k,m}(\epsilon) = \frac{1}{N_{\text{fft}}} \frac{\sin(\pi(m - k + \epsilon))}{\sin(\pi(m - k + \epsilon)/N_{\text{fft}})} \exp\left(j \pi (m - k + \epsilon) \frac{N_{\text{fft}}-1}{N_{\text{fft}}}\right)$$

- **Carrier-to-Interference Ratio (CIR)**:
  $$\text{CIR} \approx \frac{1}{\frac{\pi^2}{3}\epsilon^2} \implies \text{CIR}_{\text{dB}} \approx 5.17 - 20\log_{10}(\epsilon)$$
  Empirical measurements verify that at $\epsilon = 0.05$, $\text{CIR} \approx 20.8\text{ dB}$, falling to $8.8\text{ dB}$ at $\epsilon = 0.20$.
- **Thomas Algorithm Banded Equalizer**: Because adjacent subcarriers ($|m - k| = 1$) account for $> 92\%$ of leakage power, the system matrix is tridiagonal:
  $$Y[k] \approx a_k X[k-1] + b_k X[k] + c_k X[k+1]$$
  Solving this tridiagonal system via forward elimination and back substitution takes strictly $\mathcal{O}(N_{\text{sc}})$ operations, bypassing expensive $\mathcal{O}(N^3)$ matrix inversions while restoring EVM from $36.2\%$ down to $12.2\%$.

---

## 4. Key Empirical Figures

All figures are programmatically generated from reproducible C++ experiment outputs:

| Figure | Description | Location |
| :--- | :--- | :--- |
| **BER Waterfall Curves** | Simulated QPSK & 16-QAM BER vs Theoretical Q-function Bounds | `results/figures/ber_waterfall.png` |
| **Constellation Rotation & Collapse** | Ideal, AWGN-only, Uncompensated CFO, and Recovered constellations | `results/figures/constellations_noisy.png` |
| **LEO Doppler S-Curve Pass** | Doppler trajectory ($\pm 50.45\text{ kHz}$) and derivative slope ($-636\text{ Hz/s}$) | `results/figures/doppler_trajectory.png` |
| **Slot-by-Slot Doppler Tracking** | Instantaneous frequency tracking error ($< 25\text{ Hz}$) and BER restoration | `results/figures/doppler_tracking_performance.png` |
| **ICI Spectral Leakage & CIR** | Empirical CIR vs Analytical Theory and adjacent subcarrier power leakage | `results/figures/ici_cir_leakage.png` |
| **NTN Scenarios Benchmark** | 7-Scenario comparative bar chart of EVM and BER across impairment conditions | `results/figures/scenario_comparison.png` |

---

## 5. Candidate Background Bridge: Sensor DSP to NTN Communications

This project bridges prior undergraduate research in biomedical sensor signal processing (noisy ECG/PPG sensor filtering, harmonic detection, and heartbeat extraction) with satellite communications:

| Biomedical Sensor Signal Processing (BSc Thesis) | Satellite NTN Physical Layer Communications (NTN-OFDM) |
| :--- | :--- |
| Physiological baseline wander, motion artifacts, thermal noise | Complex Circular AWGN channel, thermal noise floor ($-174\text{ dBm/Hz}$) |
| FFT windowing, harmonic identification, spectral leakage | Radix-2 DIT FFT/IFFT OFDM orthogonal subcarrier modulation |
| Adaptive filtering & pulse tracking algorithms | Exponential moving average (EMA) Doppler frequency tracking |
| Sensor SNR calibration & peak detection thresholds | Symbol energy normalization ($E_s = 1.0$) and hard-decision constellation slicing |
| Algorithm optimization in competitive C++ programming | Zero-allocation DSP inner loops, move semantics, and cache-conscious vectors |

---

## 6. Build, Test, and Run Instructions

### Prerequisites
- **C++20 Compiler**: GCC 11+ (or Clang 13+, MSVC 2019+)
- **CMake**: Version 3.20 or newer
- **Ninja** or GNU Make
- **Python 3.8+** (with `numpy`, `pandas`, `matplotlib` for plotting)

### Compilation
```bash
# Clone the repository
git clone https://github.com/candidate/ntn-ofdm-simulator.git
cd ntn-ofdm-simulator

# Configure with CMake
cmake -B build -G Ninja -DCMAKE_BUILD_TYPE=Release

# Build all libraries, unit tests, and experiment executables
cmake --build build
```

### Running Automated Test Suite
```bash
# Run all 95 unit tests across all 7 test suites via CTest
ctest --test-dir build --output-on-failure
```
*Expected Output:*
```
1/7 Test #1: DSP_Unit_Tests ...................   Passed (22/22)
2/7 Test #2: Modulation_Unit_Tests ............   Passed (21/21)
3/7 Test #3: OFDM_Unit_Tests ..................   Passed (17/17)
4/7 Test #4: AWGN_Sim_Unit_Tests ..............   Passed (5/5)
5/7 Test #5: CFO_Unit_Tests ...................   Passed (12/12)
6/7 Test #6: Doppler_Unit_Tests ...............   Passed (10/10)
7/7 Test #7: ICI_Unit_Tests ...................   Passed (8/8)
100% tests passed out of 7 (95 total tests in 0.79s)
```

### Reproducing Experiments & Generating Figures
```bash
# 1. Satellite Link Budget Demonstration
./build/basic_link_budget

# 2. Monte Carlo AWGN BER Waterfall Simulation
./build/run_awgn_simulation
python scripts/plot_ber_curves.py

# 3. Carrier Frequency Offset (CFO) Synchronization Experiment
./build/run_cfo_experiment
python scripts/plot_cfo_impact.py

# 4. LEO Satellite Doppler Tracking Simulation
./build/run_doppler_experiment
python scripts/plot_doppler_pass.py

# 5. Inter-Carrier Interference (ICI) & Banded Equalizer Experiment
./build/run_ici_experiment
python scripts/plot_ici_results.py

# 6. Comprehensive 3GPP NTN Scenario Matrix
./build/run_scenario_matrix
python scripts/plot_scenarios_summary.py
```

---

## 7. Repository Layout

```
ntn-ofdm-simulator/
├── CMakeLists.txt              # Top-level CMake build configuration (C++20, warnings, targets)
├── README.md                   # Project overview, mathematical derivations, and benchmark
├── include/
│   ├── common/types.hpp        # Complex baseband definitions (Complex, ComplexVector, ByteVector)
│   ├── signal/
│   │   ├── dsp_utils.hpp       # Power, RMS, PAPR, EVM, noise variance calculations
│   │   └── fft.hpp             # In-place Radix-2 DIT FFT/IFFT and fftshift
│   ├── channel/
│   │   ├── noise.hpp           # Complex circular AWGN generator with deterministic seed control
│   │   ├── link_budget.hpp     # Spherical Earth slant-range, FSPL, and thermal noise calculator
│   │   └── doppler.hpp         # LEO orbital parameters, S-curve trajectory, and Doppler tracker
│   ├── modem/
│   │   └── modulation.hpp      # Gray-coded BPSK, QPSK, 16-QAM modulators and hard decision slicers
│   ├── ofdm/
│   │   ├── ofdm_config.hpp     # 5G NR OFDM numerology (N_fft, N_sc, N_cp, delta_f, DC nulling)
│   │   └── ofdm_modem.hpp      # Transmitter (IFFT+CP) and Receiver (CP removal+FFT)
│   ├── synchronization/
│   │   └── cfo.hpp             # CFO impairment channel, CP correlation estimator, and de-rotator
│   ├── compensation/
│   │   └── ici_mitigation.hpp  # Analytical ICI leakage kernel and Banded Tridiagonal Equalizer
│   ├── metrics/
│   │   └── metrics.hpp         # Empirical BER/SER counters and analytical Q-function formulas
│   └── simulation/
│       └── awgn_sim.hpp        # Automated Monte Carlo SNR sweep harness
├── src/                        # Implementations matching include headers
├── tests/                      # 7 Unit test suites covering 95 test assertions
├── experiments/                # Executables running Monte Carlo sweeps and exporting CSVs
├── scripts/                    # Python matplotlib visualization scripts
└── results/
    ├── figures/                # High-resolution publication plots (.png)
    └── tables/                 # Empirical output datasets (.csv)
```

---

## 8. License & Acknowledgments

This project is licensed under the MIT License. Developed as a technical demonstration portfolio piece for the Graduate Engineer (Signal Processing) role at Skylo.
