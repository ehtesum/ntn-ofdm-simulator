#!/usr/bin/env python3
"""
NTN-OFDM Simulator — Scenario Benchmark Matrix Visualization
Generates results/figures/scenario_comparison.png
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.size'] = 11

csv_file = "results/tables/ntn_scenarios_matrix.csv"
if os.path.exists(csv_file):
    df = pd.read_csv(csv_file)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Shorten scenario labels for plot readability
    labels = [
        "A: Ideal",
        "B: AWGN (15 dB)",
        "C: +CFO",
        "D: +Doppler",
        "E: +CFO + Doppler",
        "F: Compensated",
        "G: ICI Mitigated"
    ]

    colors = ['#2ca02c', '#1f77b4', '#d62728', '#e377c2', '#8c564b', '#ff7f0e', '#2ca02c']

    # Subplot 1: EVM Comparison (Bar chart)
    bars1 = ax1.bar(labels, df['evm_pct'], color=colors, alpha=0.85, edgecolor='black', linewidth=1.2)
    ax1.set_ylabel('RMS EVM [%]', fontweight='bold')
    ax1.set_title('EVM across 3GPP NTN Transmission Scenarios', fontweight='bold', pad=12)
    ax1.grid(True, linestyle='--', alpha=0.5, axis='y')
    ax1.set_xticklabels(labels, rotation=30, ha='right', fontweight='bold')
    
    # Add values on top of bars
    for bar in bars1:
        yval = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2.0, yval + 2.0, f'{yval:.1f}%', ha='center', va='bottom', fontsize=10, fontweight='bold')
    ax1.set_ylim(0, 165)

    # Subplot 2: BER Comparison (Log scale bar chart)
    # Clamp zero BER to 1e-5 for log-scale visualization
    ber_values = np.maximum(df['ber'], 1e-5)
    bars2 = ax2.bar(labels, ber_values, color=colors, alpha=0.85, edgecolor='black', linewidth=1.2)
    ax2.set_yscale('log')
    ax2.set_ylabel('Bit Error Rate (BER) [Log Scale]', fontweight='bold')
    ax2.set_title('Bit Error Rate (BER) across Scenarios', fontweight='bold', pad=12)
    ax2.grid(True, linestyle='--', alpha=0.5, which='both', axis='y')
    ax2.set_xticklabels(labels, rotation=30, ha='right', fontweight='bold')
    ax2.set_ylim(1e-5, 1.2)

    for i, bar in enumerate(bars2):
        raw_ber = df['ber'][i]
        text_str = f"{raw_ber:.1e}" if raw_ber > 0 else "0.0 (error-free)"
        ypos = ber_values[i] * 1.3
        ax2.text(bar.get_x() + bar.get_width()/2.0, ypos, text_str, ha='center', va='bottom', fontsize=9, fontweight='bold')

    plt.tight_layout()
    os.makedirs("results/figures", exist_ok=True)
    plt.savefig("results/figures/scenario_comparison.png", dpi=300)
    plt.close()
    print("Saved scenario benchmark figure to results/figures/scenario_comparison.png")
