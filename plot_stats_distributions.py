#!/usr/bin/env python3
"""
Plot distributions of all statistics across windows.

Usage:
    python plot_stats_distributions.py results.tsv
    python plot_stats_distributions.py results.tsv --output distributions.png
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def plot_distributions(input_file, output_file=None):
    """Plot histogram distributions for all statistics."""
    
    # Read data
    print(f"Reading {input_file}...")
    # Auto-detect separator by reading first line
    with open(input_file, 'r') as f:
        first_line = f.readline()
    
    # Check which delimiter is used
    if '\t' in first_line:
        sep = '\t'
    else:
        sep = ','
    
    df = pd.read_csv(input_file, sep=sep)
    
    # Identify statistic columns (exclude metadata columns)
    metadata_cols = ['region_contig', 'region_start', 'region_end', 
                     'window_start', 'window_stop', 'window_contig', 'n_variants']
    stat_cols = [col for col in df.columns if col not in metadata_cols]
    
    print(f"Found {len(stat_cols)} statistics to plot")
    print(f"Total windows: {len(df):,}")
    
    # Calculate grid dimensions (aim for roughly square layout)
    n_stats = len(stat_cols)
    n_cols = int(np.ceil(np.sqrt(n_stats)))
    n_rows = int(np.ceil(n_stats / n_cols))
    
    # Create figure
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5*n_cols, 4*n_rows))
    fig.suptitle('Distribution of Statistics Across All Windows', 
                 fontsize=24, fontweight='bold', y=0.995)
    
    # Flatten axes array for easier indexing
    if n_stats == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    # Bold colors for histograms
    colors = ['#e41a1c', '#377eb8', '#4daf4a', '#984ea3', '#ff7f00', 
              '#ffff33', '#a65628', '#f781bf', '#999999']
    
    # Plot each statistic
    for idx, stat in enumerate(stat_cols):
        ax = axes[idx]
        
        # Get non-NaN values
        values = df[stat].dropna()
        n_valid = len(values)
        n_nan = len(df[stat]) - n_valid
        
        if n_valid == 0:
            ax.text(0.5, 0.5, f'{stat}\n(No valid values)', 
                   ha='center', va='center', fontsize=14, fontweight='bold')
            ax.set_xticks([])
            ax.set_yticks([])
            continue
        
        # Plot histogram with bold color
        color = colors[idx % len(colors)]
        ax.hist(values, bins=50, color=color, alpha=0.8, edgecolor='black', linewidth=0.5)
        
        # Title with summary stats
        mean_val = values.mean()
        median_val = values.median()
        title = f'{stat}\n'
        title += f'mean={mean_val:.3g}, median={median_val:.3g}\n'
        title += f'n={n_valid:,}'
        if n_nan > 0:
            title += f' ({n_nan:,} NaN)'
        
        ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
        ax.set_xlabel('Value', fontsize=11, fontweight='bold')
        ax.set_ylabel('Count', fontsize=11, fontweight='bold')
        
        # Bold tick labels
        ax.tick_params(labelsize=10)
        for label in ax.get_xticklabels() + ax.get_yticklabels():
            label.set_fontweight('bold')
        
        # Add grid for easier reading
        ax.grid(True, alpha=0.3, linestyle='--', linewidth=0.5)
        ax.set_axisbelow(True)
    
    # Hide unused subplots
    for idx in range(n_stats, len(axes)):
        axes[idx].axis('off')
    
    # Adjust layout
    plt.tight_layout()
    
    # Save or show
    if output_file:
        print(f"Saving plot to {output_file}...")
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Done!")
    else:
        print("Displaying plot...")
        plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = None
    
    # Handle --output flag
    if len(sys.argv) >= 4 and sys.argv[2] == '--output':
        output_file = sys.argv[3]
    elif len(sys.argv) == 3:
        # Assume second arg is output file if no --output flag
        output_file = sys.argv[2]
    
    if not Path(input_file).exists():
        print(f"Error: File not found: {input_file}")
        sys.exit(1)
    
    plot_distributions(input_file, output_file)

