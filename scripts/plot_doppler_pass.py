#!/usr/bin/env python3
"""
Plots LEO Satellite Doppler Dynamics:
1. Doppler S-Curve and Doppler Rate over a 10-minute overhead satellite pass.
2. Slot-by-slot dynamic Doppler tracking performance and BER recovery.
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
accent_cyan = "#38bdf8"   # Sky 400
accent_amber = "#f59e0b"  # Amber 500
accent_emerald = "#10b981"# Emerald 500
accent_red = "#f87171"    # Red 400

# -------------------------------------------------------------
# 1. Full Satellite Pass Doppler Trajectory S-Curve
# -------------------------------------------------------------
traj_path = "results/tables/doppler_trajectory.csv"
if os.path.exists(traj_path):
    df_traj = pd.read_csv(traj_path)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    fig.patch.set_facecolor(bg_color)
    ax1.set_facecolor(card_color)
    ax2.set_facecolor(card_color)

    # Upper: Doppler Shift (kHz)
    ax1.plot(df_traj["Time_sec"], df_traj["Doppler_Hz"] / 1e3, color=accent_cyan, linewidth=2.4)
    ax1.axhline(0, color=grid_color, linestyle="--", linewidth=1.0)
    ax1.axvline(0, color=accent_amber, linestyle=":", linewidth=1.2, label="Closest Approach (Nadir)")
    ax1.set_title("LEO Satellite Doppler Shift Over Full Pass (600 km Altitude, 2.0 GHz Carrier)",
                  color=text_color, fontsize=12, pad=12, fontweight='bold')
    ax1.set_ylabel("Doppler Shift (kHz)", color=text_color, fontsize=11)
    ax1.grid(True, color=grid_color, alpha=0.4, linestyle=":")
    ax1.tick_params(colors=text_color)
    for s in ax1.spines.values(): s.set_color(grid_color)
    ax1.legend(facecolor=card_color, edgecolor=grid_color, labelcolor=text_color, loc="upper right")

    # Lower: Doppler Rate (Hz/s)
    ax2.plot(df_traj["Time_sec"], df_traj["DopplerRate_Hz_per_sec"], color=accent_amber, linewidth=2.4)
    ax2.axvline(0, color=accent_amber, linestyle=":", linewidth=1.2)
    ax2.set_title("Doppler Rate of Change ($df_d / dt$)", color=text_color, fontsize=11, pad=10)
    ax2.set_xlabel("Time Relative to Closest Approach (seconds)", color=text_color, fontsize=11)
    ax2.set_ylabel("Doppler Rate (Hz/s)", color=text_color, fontsize=11)
    ax2.grid(True, color=grid_color, alpha=0.4, linestyle=":")
    ax2.tick_params(colors=text_color)
    for s in ax2.spines.values(): s.set_color(grid_color)

    plt.tight_layout()
    traj_fig = "results/figures/doppler_trajectory.png"
    plt.savefig(traj_fig, dpi=200, facecolor=fig.get_facecolor(), bbox_inches='tight')
    plt.close()
    print(f"Saved Doppler trajectory plot to {traj_fig}")

# -------------------------------------------------------------
# 2. Dynamic Slot-by-Slot Tracking Performance
# -------------------------------------------------------------
track_path = "results/tables/doppler_tracking.csv"
if os.path.exists(track_path):
    df_track = pd.read_csv(track_path)

    fig, (ax_track, ax_ber) = plt.subplots(1, 2, figsize=(13, 5.5))
    fig.patch.set_facecolor(bg_color)
    ax_track.set_facecolor(card_color)
    ax_ber.set_facecolor(card_color)

    # Left: Actual vs Tracked Doppler
    ax_track.plot(df_track["Slot"], df_track["Actual_Doppler_Hz"], 'o-', color=accent_cyan,
                  linewidth=2.0, markersize=5, label="Actual Orbital Doppler")
    ax_track.plot(df_track["Slot"], df_track["Tracked_Doppler_Hz"], 's--', color=accent_emerald,
                  linewidth=1.8, markersize=5, label="Receiver Tracked Doppler")
    ax_track.set_title("Slot-by-Slot Doppler Tracking", color=text_color, fontsize=12, pad=12, fontweight='bold')
    ax_track.set_xlabel("5G NR Slot Number", color=text_color, fontsize=11)
    ax_track.set_ylabel("Carrier Doppler (Hz)", color=text_color, fontsize=11)
    ax_track.grid(True, color=grid_color, alpha=0.4, linestyle=":")
    ax_track.tick_params(colors=text_color)
    for s in ax_track.spines.values(): s.set_color(grid_color)
    ax_track.legend(facecolor=card_color, edgecolor=grid_color, labelcolor=text_color, loc="best")

    # Right: BER Uncompensated vs Compensated
    uncomp_ber = np.clip(df_track["Uncomp_BER"], 1e-5, 1.0)
    comp_ber = np.clip(df_track["Comp_BER"], 1e-5, 1.0)

    ax_ber.semilogy(df_track["Slot"], uncomp_ber, 'o-', color=accent_red,
                    linewidth=2.0, markersize=5, label="Uncompensated (Failure)")
    ax_ber.semilogy(df_track["Slot"], comp_ber, 's-', color=accent_emerald,
                    linewidth=2.0, markersize=5, label="With Doppler Tracker (BER=0.0)")
    ax_ber.set_title("Bit Error Rate (BER) Before vs. After Tracking", color=text_color, fontsize=12, pad=12, fontweight='bold')
    ax_ber.set_xlabel("5G NR Slot Number", color=text_color, fontsize=11)
    ax_ber.set_ylabel("Bit Error Rate (BER)", color=text_color, fontsize=11)
    ax_ber.set_ylim(1e-5, 0.6)
    ax_ber.grid(True, which="both", color=grid_color, alpha=0.4, linestyle=":")
    ax_ber.tick_params(colors=text_color, which="both")
    for s in ax_ber.spines.values(): s.set_color(grid_color)
    ax_ber.legend(facecolor=card_color, edgecolor=grid_color, labelcolor=text_color, loc="center left")

    plt.tight_layout()
    track_fig = "results/figures/doppler_tracking_performance.png"
    plt.savefig(track_fig, dpi=200, facecolor=fig.get_facecolor(), bbox_inches='tight')
    plt.close()
    print(f"Saved Doppler tracking performance plot to {track_fig}")
