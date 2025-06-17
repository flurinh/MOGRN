#!/usr/bin/env python3
"""
Debug script to understand why carbon mapping is failing.
"""

import os
import sys
import pandas as pd

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from src.data_processing import load_opsin_structures
from src.retinal_carbon_mapping_correct import assign_retinal_carbons_correct
from src.retinal_carbon_mapping_fixed import assign_retinal_carbons_robust, test_ring_finding


def debug_retinal_structure(struct_id, struct_data):
    """Debug a single retinal structure."""
    print(f"\n{'='*80}")
    print(f"Debugging structure: {struct_id}")
    print('='*80)
    
    df = struct_data['df']
    
    # Filter for retinal atoms
    retinal_mask = df['res_name3l'].isin(['RET', 'LYR'])
    retinal_df = df[retinal_mask].copy()
    
    # If we have LYR, exclude lysine atoms
    if 'LYR' in retinal_df['res_name3l'].values:
        lys_atoms = {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'}
        retinal_df = retinal_df[~retinal_df['res_atom_name'].isin(lys_atoms)]
    
    print(f"\nRetinal DataFrame shape: {retinal_df.shape}")
    print(f"Columns: {list(retinal_df.columns)}")
    
    # Check atom_name column
    print(f"\nUnique atom_name values: {sorted(retinal_df['atom_name'].unique())}")
    print(f"Count by atom_name:")
    print(retinal_df['atom_name'].value_counts().sort_index())
    
    # Check res_atom_name column
    print(f"\nUnique res_atom_name values: {sorted(retinal_df['res_atom_name'].unique())}")
    
    # Look for carbon atoms
    carbon_mask = retinal_df['atom_name'] == 'C'
    carbon_count = carbon_mask.sum()
    print(f"\nCarbon atoms found (atom_name == 'C'): {carbon_count}")
    
    # Try different ways to identify carbons
    if carbon_count == 0:
        # Maybe carbons are labeled differently
        print("\nNo carbons found with atom_name == 'C'. Checking other possibilities...")
        
        # Check if carbon atoms have specific names
        carbon_pattern_mask = retinal_df['res_atom_name'].str.match(r'^C\d+$')
        if carbon_pattern_mask.any():
            print(f"Found {carbon_pattern_mask.sum()} carbons with pattern C1, C2, etc. in res_atom_name")
            print("Carbon names:", sorted(retinal_df[carbon_pattern_mask]['res_atom_name'].unique()))
        
        # Check element column if it exists
        if 'element' in retinal_df.columns:
            element_carbon_mask = retinal_df['element'] == 'C'
            print(f"Carbons found via element column: {element_carbon_mask.sum()}")
        
        # Check type_symbol if it exists
        if 'type_symbol' in retinal_df.columns:
            type_carbon_mask = retinal_df['type_symbol'] == 'C'
            print(f"Carbons found via type_symbol column: {type_carbon_mask.sum()}")
    
    # Show first few rows of retinal atoms
    print(f"\nFirst 5 retinal atoms:")
    display_cols = ['res_atom_name', 'atom_name', 'x', 'y', 'z']
    if 'element' in retinal_df.columns:
        display_cols.append('element')
    if 'type_symbol' in retinal_df.columns:
        display_cols.append('type_symbol')
    print(retinal_df[display_cols].head())
    
    return retinal_df


def main():
    """Debug carbon mapping issues."""
    
    print("Loading opsin structures...")
    
    # Set up paths
    output_dir = os.path.join(os.path.dirname(__file__), 'opsin_output')
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    
    # Load structures
    data = load_opsin_structures(
        data_dir=data_dir,
        output_dir=output_dir,
        chain_id='A',
        use_cache=True,
        visualize=False
    )
    
    processed_structures = data.get('processed_structures', {})
    
    if not processed_structures:
        print("No structures loaded!")
        return
    
    # Find experimental structures
    experimental_structures = {
        sid: sdata for sid, sdata in processed_structures.items()
        if '_pred' not in sid and '_model' not in sid
    }
    
    print(f"\nFound {len(experimental_structures)} experimental structures")
    
    # Debug the first few structures
    for i, (struct_id, struct_data) in enumerate(experimental_structures.items()):
        if i >= 3:  # Just debug first 3
            break
        
        retinal_df = debug_retinal_structure(struct_id, struct_data)
        
        # Try running the carbon assignment
        if not retinal_df.empty:
            print(f"\nTrying original carbon assignment for {struct_id}...")
            try:
                assignments = assign_retinal_carbons_correct(retinal_df)
                print(f"Carbon assignments: {list(assignments.keys())}")
            except Exception as e:
                print(f"Error during carbon assignment: {e}")
            
            print(f"\nTrying robust carbon assignment for {struct_id}...")
            try:
                assignments_robust = assign_retinal_carbons_robust(retinal_df)
                print(f"Robust carbon assignments: {list(assignments_robust.keys())}")
            except Exception as e:
                print(f"Error during robust assignment: {e}")
                import traceback
                traceback.print_exc()
            
            print(f"\nTesting ring finding for {struct_id}...")
            try:
                ring = test_ring_finding(retinal_df)
            except Exception as e:
                print(f"Error during ring test: {e}")
                import traceback
                traceback.print_exc()


if __name__ == "__main__":
    main()