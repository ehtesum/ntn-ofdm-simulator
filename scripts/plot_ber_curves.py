#!/usr/bin/env python3
"""
Generates publication-quality BER waterfall curves and EVM vs SNR plots
from simulated Monte Carlo AWGN data.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

os.makedirs("results/figures", exist_ok=True)

# Dark modern styling
bg_color = "#0f172a"      # Slate 900
card_color = "#1e293b"    # Slate 800
text_color = "#f8fafc"    # Slate 50
grid_color = "#334155"    # Slate 700
accent_cyan = "#38bdf8"   # Sky 400
accent_amber = "#f59e0b"  # Amber 500
accent_emerald = "#10b981"# Emerald 500
accent_purple = "#c084fc" # Purple 400

qpsk_path = "results/tables/awgn_results_qpsk.csv"
qam_path = "results/tables/awgn_results_16qam.csv"

if not os.path.exists(qpsk_path) or not os.path.exists(qam_path):
    print("CSV data files not found. Please run the simulation executable first.")
    exit(0)

df_qpsk = pd.read_csv(qpsk_path)
df_qam = pd.read_csv(qam_path)

# Filter out points with zero errors for log-scale plotting, or set floor
ber_floor = 1e-6

# -------------------------------------------------------------
# 1. BER vs Eb/N0 Waterfall Curves
# -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 6))
fig.patch.set_facecolor(bg_color)
ax.set_facecolor(card_color)

# Theoretical curves
ebn0_fine = np.linspace(-2, 14, 200)
# Q-function: 0.5 * erfc(x / sqrt(2))
from scipy.special import erfc
theory_qpsk = 0.5 * erfc(np.sqrt(10.0 ** (ebn0_fine / 10.0)))
theory_qam = 0.75 * 0.5 * erfc(np.sqrt(0.8 * 10.0 ** (ebn0_fine / 10.0)))

ax.plot(ebn0_fine, theory_qpsk, color=accent_cyan, linestyle="--", linewidth=1.8, label="Theory QPSK (AWGN)")
ax.plot(ebn0_fine, theory_qam, color=accent_purple, linestyle="--", linewidth=1.8, label="Theory 16-QAM (AWGN)")

# Simulated QPSK
qpsk_ber_plot = np.clip(df_qpsk["BER"], ber_floor, 1.0)
ax.semilogy(df_qpsk["EbN0_dB"], qpsk_ber_plot, 'o-', color=accent_cyan,
            markersize=7, linewidth=2.0, label="Simulated QPSK (OFDM)", zorder=5)

# Simulated 16-QAM
qam_ber_plot = np.clip(df_qam["BER"], ber_floor, 1.0)
ax.semilogy(df_qam["EbN0_dB"], qam_ber_plot, 's-', color=accent_purple,
            markersize=7, linewidth=2.0, label="Simulated 16-QAM (OFDM)", zorder=5)

ax.set_yscale("log")
ax.set_ylim(1e-5, 0.5)
ax.set_xlim(-2, 14)
ax.set_title("OFDM Bit Error Rate (BER) Waterfall in AWGN Channel\nSimulated Monte Carlo vs. Theoretical Bounds",
             color=text_color, fontsize=13, pad=14, fontweight='bold')
ax.set_xlabel("Energy per Bit to Noise Ratio: $E_b / N_0$ (dB)", color=text_color, fontsize=11)
ax.set_ylabel("Bit Error Rate (BER)", color=text_color, fontsize=11)
ax.grid(True, which="both", color=grid_color, alpha=0.4, linestyle=":")
ax.tick_params(colors=text_color, which="both")
for spine in ax.spines.values():
    spine.set_color(grid_color)
ax.legend(facecolor=card_color, edgecolor=grid_color, labelcolor=text_color, fontsize=10, loc="lower left")

plt.tight_layout()
ber_fig_path = "results/figures/ber_waterfall.png"
plt.savefig(ber_fig_path, dpi=200, facecolor=fig.get_facecolor(), bbox_inches='tight')
plt.close()
print(f"Saved BER waterfall plot to {ber_fig_path}")

# -------------------------------------------------------------
# 2. EVM (%) vs SNR (dB) Curve
# -------------------------------------------------------------
fig, ax = plt.subplots(figsize=(9, 5.5))
fig.patch.set_facecolor(bg_color)
ax.set_facecolor(card_color)

ax.plot(df_qpsk["SNR_dB"], df_qpsk["EVM_Percent"], 'o-', color=accent_emerald,
        linewidth=2.2, markersize=7, label="Measured EVM (%)", zorder=5)

# Theoretical EVM line: EVM_rms ~ 1 / sqrt(SNR_lin) * 100%
snr_fine = np.linspace(0, 20, 100)
theory_evm = 100.0 / np.sqrt(10.0 ** (snr_fine / 10.0))
ax.plot(snr_fine, theory_evm, color=accent_amber, linestyle="--", linewidth=1.8, label="Theoretical Bound $1/\\sqrt{\\mathrm{SNR}}$")

ax.set_title("Error Vector Magnitude (EVM) vs. Channel SNR",
             color=text_color, fontsize=13, pad=14, fontweight='bold')
ax.set_xlabel("Channel SNR (dB)", color=text_color, fontsize=11)
ax.set_ylabel("EVM RMS (%)", color=text_color, fontsize=11)
ax.grid(True, color=grid_color, alpha=0.4, linestyle=":")
ax.tick_params(colors=text_color)
for spine in ax.spines.values():
    spine.set_color(grid_color)
ax.legend(facecolor=card_color, edgecolor=grid_color, labelcolor=text_color, fontsize=10)

plt.tight_layout()
evm_fig_path = "results/figures/evm_vs_snr.png"
plt.savefig(evm_fig_path, dpi=200, facecolor=fig.get_facecolor(), bbox_inches='tight')
plt.close()
print(f"Saved EVM plot to {evm_fig_path}")
