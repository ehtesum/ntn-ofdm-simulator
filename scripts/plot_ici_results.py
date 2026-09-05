#!/usr/bin/env python3
"""
NTN-OFDM Simulator — ICI Visualization Script
Generates publication-quality figures for:
  1. results/figures/ici_cir_leakage.png (Empirical vs Analytical CIR and Spectral Leakage)
  2. results/figures/ici_mitigation_evm.png (Equalization EVM and BER improvement)
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11

os.makedirs("results/figures", exist_ok=True)

# -------------------------------------------------------------
# Figure 1: CIR Curve and Spectral Leakage
# -------------------------------------------------------------
cir_file = "results/tables/ici_cir_curve.csv"
leakage_file = "results/tables/ici_leakage_profile.csv"

if os.path.exists(cir_file) and os.path.exists(leakage_file):
    df_cir = pd.read_csv(cir_file)
    df_leak = pd.read_csv(leakage_file)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Subplot 1: CIR vs Normalized Frequency Offset
    valid_cir = df_cir[df_cir['epsilon'] > 0.005]
    ax1.plot(valid_cir['epsilon'], valid_cir['empirical_cir_db'], 'bo-', label='Empirical CIR (Simulator)', linewidth=2, markersize=5)
    ax1.plot(valid_cir['epsilon'], valid_cir['analytical_cir_db'], 'r--', label=r'Analytical Bound: $10\log_{10}\left(\frac{1}{\frac{\pi^2}{3}\epsilon^2}\right)$', linewidth=2.2)

    ax1.set_xlabel(r'Normalized Frequency Offset $\epsilon = \Delta f / \Delta f_{\mathrm{sc}}$', fontweight='bold')
    ax1.set_ylabel('Carrier-to-Interference Ratio (CIR) [dB]', fontweight='bold')
    ax1.set_title('OFDM Orthogonality Degradation under Frequency Offset', fontweight='bold', pad=12)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(frameon=True, facecolor='white', framealpha=0.9)
    ax1.set_ylim(0, 50)

    # Subplot 2: Spectral Leakage Profile across Adjacent Subcarriers
    offsets = df_leak['subcarrier_offset']
    ax2.plot(offsets, 10 * np.log10(np.maximum(df_leak['power_eps_000'], 1e-12)), 'k-o', label=r'$\epsilon = 0.00$ (Ideal Orthogonality)', linewidth=1.8, markersize=5)
    ax2.plot(offsets, 10 * np.log10(np.maximum(df_leak['power_eps_005'], 1e-12)), 'g-s', label=r'$\epsilon = 0.05$ (5% subcarrier offset)', linewidth=1.8, markersize=5)
    ax2.plot(offsets, 10 * np.log10(np.maximum(df_leak['power_eps_015'], 1e-12)), 'm-^', label=r'$\epsilon = 0.15$ (15% subcarrier offset)', linewidth=1.8, markersize=5)
    ax2.plot(offsets, 10 * np.log10(np.maximum(df_leak['power_eps_030'], 1e-12)), 'r-d', label=r'$\epsilon = 0.30$ (30% subcarrier offset)', linewidth=1.8, markersize=5)

    ax2.set_xlabel('Subcarrier Offset $(k - k_0)$', fontweight='bold')
    ax2.set_ylabel('Normalized Subcarrier Power [dB]', fontweight='bold')
    ax2.set_title('Inter-Carrier Interference (ICI) Spectral Leakage', fontweight='bold', pad=12)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(frameon=True, facecolor='white', framealpha=0.9)
    ax2.set_ylim(-60, 5)

    plt.tight_layout()
    plt.savefig("results/figures/ici_cir_leakage.png", dpi=300)
    plt.close()
    print("Saved ICI CIR and leakage plot to results/figures/ici_cir_leakage.png")

# -------------------------------------------------------------
# Figure 2: Equalization Performance (EVM and BER)
# -------------------------------------------------------------
eq_file = "results/tables/ici_equalization_results.csv"
if os.path.exists(eq_file):
    df_eq = pd.read_csv(eq_file)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

    # Subplot 1: EVM (%) vs Residual Epsilon
    ax1.plot(df_eq['residual_eps'], df_eq['raw_evm_pct'], 'r-s', label='Raw (No Compensation)', linewidth=2, markersize=6)
    ax1.plot(df_eq['residual_eps'], df_eq['derotated_evm_pct'], 'orange', marker='^', linestyle='--', label='Symbol De-rotation Only', linewidth=2, markersize=6)
    ax1.plot(df_eq['residual_eps'], df_eq['equalized_evm_pct'], 'g-o', label='Banded Tridiagonal Equalizer (O(N))', linewidth=2.2, markersize=6)

    ax1.set_xlabel(r'Residual Normalized Offset $\epsilon$', fontweight='bold')
    ax1.set_ylabel('RMS EVM [%]', fontweight='bold')
    ax1.set_title('16-QAM OFDM EVM vs Residual Frequency Offset', fontweight='bold', pad=12)
    ax1.grid(True, linestyle='--', alpha=0.6)
    ax1.legend(frameon=True, facecolor='white', framealpha=0.9)

    # Subplot 2: BER vs Residual Epsilon
    ax2.semilogy(df_eq['residual_eps'], np.maximum(df_eq['raw_ber'], 1e-6), 'r-s', label='Raw (No Compensation)', linewidth=2, markersize=6)
    ax2.semilogy(df_eq['residual_eps'], np.maximum(df_eq['derotated_ber'], 1e-6), 'orange', marker='^', linestyle='--', label='Symbol De-rotation Only', linewidth=2, markersize=6)
    ax2.semilogy(df_eq['residual_eps'], np.maximum(df_eq['equalized_ber'], 1e-6), 'g-o', label='Banded Tridiagonal Equalizer', linewidth=2.2, markersize=6)

    ax2.set_xlabel(r'Residual Normalized Offset $\epsilon$', fontweight='bold')
    ax2.set_ylabel('Bit Error Rate (BER)', fontweight='bold')
    ax2.set_title('16-QAM OFDM BER vs Residual Frequency Offset', fontweight='bold', pad=12)
    ax2.grid(True, linestyle='--', alpha=0.6)
    ax2.legend(frameon=True, facecolor='white', framealpha=0.9)
    ax2.set_ylim(5e-7, 1.0)

    plt.tight_layout()
    plt.savefig("results/figures/ici_mitigation_evm.png", dpi=300)
    plt.close()
    print("Saved ICI mitigation performance plot to results/figures/ici_mitigation_evm.png")
