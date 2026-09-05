#!/usr/bin/env python3
"""
Plots CFO impact on OFDM:
1. Constellation comparison (Uncompensated smeared ring vs. Compensated clean clusters).
2. BER & EVM vs. Normalized Carrier Frequency Offset (epsilon = f_cfo / delta_f).
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

os.makedirs("results/figures", exist_ok=True)

bg_color = "#0f172a"      # Slate 900
card_color = "#1e293b"    # Slate 800
text_color = "#f8fafc"    # Slate 50
grid_color = "#334155"    # Slate 700
accent_red = "#f87171"    # Red 400
accent_emerald = "#10b981"# Emerald 500
accent_cyan = "#38bdf8"   # Sky 400
accent_amber = "#f59e0b"  # Amber 500

# -------------------------------------------------------------
# 1. Constellation Before vs After CFO Compensation
# -------------------------------------------------------------
uncomp_path = "results/tables/cfo_constellation_uncomp.csv"
comp_path = "results/tables/cfo_constellation_comp.csv"

if os.path.exists(uncomp_path) and os.path.exists(comp_path):
    df_uncomp = pd.read_csv(uncomp_path)
    df_comp = pd.read_csv(comp_path)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5.5))
    fig.patch.set_facecolor(bg_color)

    # Before Compensation
    ax1 = axes[0]
    ax1.set_facecolor(card_color)
    ax1.scatter(df_uncomp["I"], df_uncomp["Q"], color=accent_red, alpha=0.5, s=20, edgecolors='none')
    ax1.set_title("Before Compensation ($\epsilon = 0.10, f_{\\mathrm{cfo}} = 1.5\\mathrm{\\ kHz}$)\nSevere Phase Rotation & Smearing",
                  color=text_color, fontsize=11, pad=12, fontweight='bold')
    ax1.axhline(0, color=grid_color, linestyle='--', linewidth=0.9)
    ax1.axvline(0, color=grid_color, linestyle='--', linewidth=0.9)
    ax1.set_xlim(-1.8, 1.8)
    ax1.set_ylim(-1.8, 1.8)
    ax1.set_xlabel("In-Phase (I)", color=text_color)
    ax1.set_ylabel("Quadrature (Q)", color=text_color)
    ax1.tick_params(colors=text_color)
    ax1.grid(True, color=grid_color, alpha=0.4, linestyle=":")
    for spine in ax1.spines.values():
        spine.set_color(grid_color)

    # After Compensation
    ax2 = axes[1]
    ax2.set_facecolor(card_color)
    ax2.scatter(df_comp["I"], df_comp["Q"], color=accent_emerald, alpha=0.5, s=20, edgecolors='none')
    # Ideal centroids
    k = 1.0 / np.sqrt(2.0)
    ideal_pts = np.array([k + 1j*k, -k + 1j*k, -k - 1j*k, k - 1j*k])
    ax2.scatter(ideal_pts.real, ideal_pts.imag, color="#ffffff", marker='x', s=90, linewidths=2.5, zorder=6, label="Ideal Centroids")
    ax2.set_title("After CP-Correlation Compensation\nClean Constellation Restored (BER $\\to 0.0$)",
                  color=text_color, fontsize=11, pad=12, fontweight='bold')
    ax2.axhline(0, color=grid_color, linestyle='--', linewidth=0.9)
    ax2.axvline(0, color=grid_color, linestyle='--', linewidth=0.9)
    ax2.set_xlim(-1.8, 1.8)
    ax2.set_ylim(-1.8, 1.8)
    ax2.set_xlabel("In-Phase (I)", color=text_color)
    ax2.set_ylabel("Quadrature (Q)", color=text_color)
    ax2.tick_params(colors=text_color)
    ax2.grid(True, color=grid_color, alpha=0.4, linestyle=":")
    for spine in ax2.spines.values():
        spine.set_color(grid_color)
    ax2.legend(facecolor=card_color, edgecolor=grid_color, labelcolor=text_color, loc="upper right")

    plt.tight_layout()
    fig_constell_path = "results/figures/cfo_constellation_before_after.png"
    plt.savefig(fig_constell_path, dpi=200, facecolor=fig.get_facecolor(), bbox_inches='tight')
    plt.close()
    print(f"Saved constellation comparison figure to {fig_constell_path}")

# -------------------------------------------------------------
# 2. BER and EVM vs Normalized CFO Sweep
# -------------------------------------------------------------
cfo_csv_path = "results/tables/cfo_comparison.csv"
if os.path.exists(cfo_csv_path):
    df_cfo = pd.read_csv(cfo_csv_path)

    fig, (ax_ber, ax_evm) = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.patch.set_facecolor(bg_color)
    ax_ber.set_facecolor(card_color)
    ax_evm.set_facecolor(card_color)

    # BER vs Epsilon
    uncomp_ber = np.clip(df_cfo["Uncomp_BER"], 1e-5, 1.0)
    comp_ber = np.clip(df_cfo["Comp_BER"], 1e-5, 1.0)

    ax_ber.semilogy(df_cfo["Normalized_CFO"], uncomp_ber, 'o-', color=accent_red,
                    linewidth=2.2, markersize=7, label="Uncompensated (Catastrophic Failure)")
    ax_ber.semilogy(df_cfo["Normalized_CFO"], comp_ber, 's-', color=accent_emerald,
                    linewidth=2.2, markersize=7, label="With CP-Correlation Compensation")

    ax_ber.set_title("Bit Error Rate (BER) vs. Carrier Frequency Offset\n$\\Delta f = 15\\mathrm{\\ kHz}$, SNR = $20\\mathrm{\\ dB}$",
                     color=text_color, fontsize=11, pad=12, fontweight='bold')
    ax_ber.set_xlabel("Normalized CFO $\\epsilon = f_{\\mathrm{cfo}} / \\Delta f$", color=text_color, fontsize=10)
    ax_ber.set_ylabel("Bit Error Rate (BER)", color=text_color, fontsize=10)
    ax_ber.set_ylim(1e-5, 0.6)
    ax_ber.grid(True, which="both", color=grid_color, alpha=0.4, linestyle=":")
    ax_ber.tick_params(colors=text_color, which="both")
    for spine in ax_ber.spines.values():
        spine.set_color(grid_color)
    ax_ber.legend(facecolor=card_color, edgecolor=grid_color, labelcolor=text_color, loc="center left")

    # EVM vs Epsilon
    ax_evm.plot(df_cfo["Normalized_CFO"], df_cfo["Uncomp_EVM_Percent"], 'o-', color=accent_red,
                linewidth=2.2, markersize=7, label="Uncompensated EVM (%)")
    ax_evm.plot(df_cfo["Normalized_CFO"], df_cfo["Comp_EVM_Percent"], 's-', color=accent_emerald,
                linewidth=2.2, markersize=7, label="Compensated EVM (%)")

    ax_evm.set_title("EVM (%) vs. Carrier Frequency Offset",
                     color=text_color, fontsize=11, pad=12, fontweight='bold')
    ax_evm.set_xlabel("Normalized CFO $\\epsilon = f_{\\mathrm{cfo}} / \\Delta f$", color=text_color, fontsize=10)
    ax_evm.set_ylabel("EVM RMS (%)", color=text_color, fontsize=10)
    ax_evm.grid(True, color=grid_color, alpha=0.4, linestyle=":")
    ax_evm.tick_params(colors=text_color)
    for spine in ax_evm.spines.values():
        spine.set_color(grid_color)
    ax_evm.legend(facecolor=card_color, edgecolor=grid_color, labelcolor=text_color, loc="center left")

    plt.tight_layout()
    fig_sweep_path = "results/figures/cfo_ber_impact.png"
    plt.savefig(fig_sweep_path, dpi=200, facecolor=fig.get_facecolor(), bbox_inches='tight')
    plt.close()
    print(f"Saved CFO performance curves to {fig_sweep_path}")
