#!/usr/bin/env python3
"""
Analyze RMSD matrix and produce per-structure statistics.

Outputs a table with mean, max, std RMSD for each structure, ranked by mean.
"""

import pandas as pd
import numpy as np
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent


def main():
    # Load RMSD matrix
    rmsd_file = PROJECT_ROOT / "opsin_output" / "rmsd_matrix.csv"

    if not rmsd_file.exists():
        print(f"[ERROR] RMSD matrix not found: {rmsd_file}")
        return

    print(f"Loading RMSD matrix from {rmsd_file}")
    df = pd.read_csv(rmsd_file, index_col=0)

    # Ensure index is string (avoid 1e12 issue)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    print(f"Loaded {len(df)} x {len(df.columns)} matrix")

    # Calculate statistics for each structure
    stats = []
    for struct_id in df.index:
        # Get all RMSD values for this structure (excluding self-comparison)
        row = df.loc[struct_id].drop(struct_id, errors='ignore')
        row = row.astype(float)
        row = row[~np.isnan(row)]  # Remove NaN values

        if len(row) == 0:
            continue

        stats.append({
            'structure': struct_id,
            'mean': row.mean(),
            'std': row.std(),
            'max': row.max(),
            'min': row.min(),
            'n_comparisons': len(row)
        })

    # Create DataFrame and sort by mean
    stats_df = pd.DataFrame(stats)
    stats_df = stats_df.sort_values('mean', ascending=True)

    # Print table
    print("\n" + "=" * 80)
    print("RMSD Statistics per Structure (ranked by mean)")
    print("=" * 80)
    print(f"{'Structure':<35} {'Mean':>8} {'Std':>8} {'Max':>8} {'Min':>8} {'N':>5}")
    print("-" * 80)

    for _, row in stats_df.iterrows():
        print(f"{row['structure']:<35} {row['mean']:>8.2f} {row['std']:>8.2f} {row['max']:>8.2f} {row['min']:>8.2f} {row['n_comparisons']:>5.0f}")

    print("-" * 80)
    print(f"{'Overall':<35} {stats_df['mean'].mean():>8.2f} {stats_df['std'].mean():>8.2f} {stats_df['max'].mean():>8.2f} {stats_df['min'].mean():>8.2f}")

    # Save to CSV
    output_file = PROJECT_ROOT / "opsin_output" / "rmsd_stats.csv"
    stats_df.to_csv(output_file, index=False)
    print(f"\nSaved statistics to {output_file}")


if __name__ == "__main__":
    main()
