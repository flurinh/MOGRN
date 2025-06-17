#!/usr/bin/env python3
"""Check available GRN positions in the MSA table."""

import pandas as pd
import os
from pathlib import Path

def check_grn_positions():
    """Load GRN table and display available positions."""
    # Check for MSA table
    msa_path = Path("opsin_output/msa_table_grn.csv")
    
    if not msa_path.exists():
        print(f"ERROR: MSA table not found at {msa_path}")
        print("Please run the analysis workflow first.")
        return
    
    # Load the MSA table
    print(f"Loading MSA table from: {msa_path}")
    msa_df = pd.read_csv(msa_path)
    
    # Display basic info
    print(f"\nTable shape: {msa_df.shape}")
    print(f"Number of structures: {len(msa_df)}")
    
    # Get GRN columns (excluding non-GRN columns)
    non_grn_cols = ['structure_id', 'structure_type', 'family', 'group']
    grn_columns = [col for col in msa_df.columns if col not in non_grn_cols]
    
    print(f"\nTotal GRN positions: {len(grn_columns)}")
    print("\nFirst 20 GRN positions:")
    for i, col in enumerate(grn_columns[:20]):
        print(f"  {col}")
    
    # Check for specific positions
    test_positions = ['5.5', '5.50', '6.50', '7.50', '1.50', '2.50', '3.50', '4.50']
    print("\nChecking common GRN positions:")
    for pos in test_positions:
        if pos in grn_columns:
            # Count non-empty values
            non_empty = msa_df[pos].notna().sum()
            print(f"  {pos}: FOUND ({non_empty} non-empty values)")
        else:
            print(f"  {pos}: NOT FOUND")
    
    # Find positions with most data
    print("\nTop 10 positions with most data:")
    position_counts = {}
    for col in grn_columns:
        count = msa_df[col].notna().sum()
        if count > 0:
            position_counts[col] = count
    
    sorted_positions = sorted(position_counts.items(), key=lambda x: x[1], reverse=True)
    for pos, count in sorted_positions[:10]:
        print(f"  {pos}: {count} structures")
    
    # Sample some actual residue values
    print("\nSample residue values from top positions:")
    for pos, _ in sorted_positions[:5]:
        values = msa_df[pos].dropna().unique()[:5]
        print(f"  {pos}: {', '.join(values)}")

if __name__ == "__main__":
    check_grn_positions()