#!/usr/bin/env python3
"""
Generates clean, publication-quality constellation diagrams for
BPSK, QPSK, and 16-QAM under ideal conditions and with calibrated AWGN.
Outputs are saved into results/figures/.
"""

import os
import numpy as np
import matplotlib.pyplot as plt

os.makedirs("results/figures", exist_ok=True)

# -------------------------------------------------------------
# 1. Generate Constellation Diagrams (Ideal with Bit Labels)
# -------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Dark modern styling for rich aesthetics
bg_color = "#0f172a"      # Slate 900
card_color = "#1e293b"    # Slate 800
text_color = "#f8fafc"    # Slate 50
grid_color = "#334155"    # Slate 700
accent_cyan = "#38bdf8"   # Sky 400
accent_amber = "#f59e0b"  # Amber 500
accent_emerald = "#10b981"# Emerald 500

fig.patch.set_facecolor(bg_color)

# --- BPSK ---
ax = axes[0]
ax.set_facecolor(card_color)
bpsk_pts = np.array([1.0 + 0.0j, -1.0 + 0.0j])
bpsk_labels = ["0", "1"]
ax.scatter(bpsk_pts.real, bpsk_pts.imag, color=accent_cyan, s=120, zorder=5, edgecolors="#ffffff", linewidths=1.5)
for pt, lbl in zip(bpsk_pts, bpsk_labels):
    ax.annotate(lbl, (pt.real, pt.imag), textcoords="offset points", xytext=(0, 14),
                ha='center', fontsize=12, fontweight='bold', color=accent_amber)
ax.set_title("BPSK Constellation (1 bit/sym)\n$E_s = 1.0$", color=text_color, fontsize=12, pad=12)

# --- QPSK ---
ax = axes[1]
ax.set_facecolor(card_color)
k_qpsk = 1.0 / np.sqrt(2.0)
qpsk_pts = np.array([
    k_qpsk + 1j * k_qpsk,   # 00
    -k_qpsk + 1j * k_qpsk,  # 10
    -k_qpsk - 1j * k_qpsk,  # 11
    k_qpsk - 1j * k_qpsk    # 01
])
qpsk_labels = ["00", "10", "11", "01"]
ax.scatter(qpsk_pts.real, qpsk_pts.imag, color=accent_emerald, s=120, zorder=5, edgecolors="#ffffff", linewidths=1.5)
for pt, lbl in zip(qpsk_pts, qpsk_labels):
    ax.annotate(lbl, (pt.real, pt.imag), textcoords="offset points", xytext=(0, 12),
                ha='center', fontsize=11, fontweight='bold', color=accent_amber)
ax.set_title("QPSK Gray Constellation (2 bits/sym)\n$E_s = 1.0$, $d_{min} = \\sqrt{2}$", color=text_color, fontsize=12, pad=12)

# --- 16-QAM ---
ax = axes[2]
ax.set_facecolor(card_color)
k_qam = 1.0 / np.sqrt(10.0)
levels = [-3.0, -1.0, 1.0, 3.0]
bit_map_1d = {-3.0: "10", -1.0: "11", 1.0: "01", 3.0: "00"}
qam_pts = []
qam_labels = []
for i_lvl in levels:
    for q_lvl in levels:
        qam_pts.append((i_lvl + 1j * q_lvl) * k_qam)
        qam_labels.append(bit_map_1d[i_lvl] + bit_map_1d[q_lvl])
qam_pts = np.array(qam_pts)
ax.scatter(qam_pts.real, qam_pts.imag, color="#a855f7", s=90, zorder=5, edgecolors="#ffffff", linewidths=1.2)
for pt, lbl in zip(qam_pts, qam_labels):
    ax.annotate(lbl, (pt.real, pt.imag), textcoords="offset points", xytext=(0, 9),
                ha='center', fontsize=8, fontweight='bold', color=accent_amber)
ax.set_title("16-QAM Gray Constellation (4 bits/sym)\n$E_s = 1.0$, $d_{min} = 2/\\sqrt{10}$", color=text_color, fontsize=12, pad=12)

for ax in axes:
    ax.axhline(0, color=grid_color, linestyle='--', linewidth=0.9)
    ax.axvline(0, color=grid_color, linestyle='--', linewidth=0.9)
    ax.grid(True, color=grid_color, alpha=0.4, linestyle=':')
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_xlabel("In-Phase (I)", color=text_color)
    ax.set_ylabel("Quadrature (Q)", color=text_color)
    ax.tick_params(colors=text_color)
    for spine in ax.spines.values():
        spine.set_color(grid_color)

plt.tight_layout()
ideal_path = "results/figures/constellations_ideal.png"
plt.savefig(ideal_path, dpi=200, facecolor=fig.get_facecolor(), bbox_inches='tight')
plt.close()
print(f"Saved ideal constellation diagram to {ideal_path}")

# -------------------------------------------------------------
# 2. Generate Noisy Constellation Diagram (SNR = 16 dB)
# -------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
fig.patch.set_facecolor(bg_color)
np.random.seed(42)
num_syms = 1500
snr_db = 16.0
snr_lin = 10.0 ** (snr_db / 10.0)
noise_sigma = np.sqrt(1.0 / snr_lin / 2.0)

for idx, (pts, title, col) in enumerate([
    (bpsk_pts, "BPSK @ 16 dB SNR", accent_cyan),
    (qpsk_pts, "QPSK @ 16 dB SNR", accent_emerald),
    (qam_pts, "16-QAM @ 16 dB SNR", "#a855f7")
]):
    ax = axes[idx]
    ax.set_facecolor(card_color)
    
    # Sample random reference points and add AWGN
    rand_indices = np.random.choice(len(pts), size=num_syms)
    clean_syms = pts[rand_indices]
    noise = np.random.normal(0, noise_sigma, num_syms) + 1j * np.random.normal(0, noise_sigma, num_syms)
    rx_syms = clean_syms + noise
    
    # Plot received noisy cloud
    ax.scatter(rx_syms.real, rx_syms.imag, color=col, alpha=0.4, s=18, edgecolors='none', label='Received')
    # Plot ideal centroids
    ax.scatter(pts.real, pts.imag, color="#ef4444", s=80, marker='x', linewidths=2.5, zorder=6, label='Ideal Centroid')
    
    ax.axhline(0, color=grid_color, linestyle='--', linewidth=0.9)
    ax.axvline(0, color=grid_color, linestyle='--', linewidth=0.9)
    ax.grid(True, color=grid_color, alpha=0.4, linestyle=':')
    ax.set_xlim(-1.6, 1.6)
    ax.set_ylim(-1.6, 1.6)
    ax.set_title(title, color=text_color, fontsize=12, pad=12)
    ax.set_xlabel("In-Phase (I)", color=text_color)
    ax.set_ylabel("Quadrature (Q)", color=text_color)
    ax.tick_params(colors=text_color)
    for spine in ax.spines.values():
        spine.set_color(grid_color)
    ax.legend(facecolor=card_color, edgecolor=grid_color, labelcolor=text_color, loc='upper right')

plt.tight_layout()
noisy_path = "results/figures/constellations_noisy.png"
plt.savefig(noisy_path, dpi=200, facecolor=fig.get_facecolor(), bbox_inches='tight')
plt.close()
print(f"Saved noisy constellation diagram to {noisy_path}")
