#!/usr/bin/env python3
"""
Debug script to check if cached structures have been centered/normalized
"""

import pickle
import numpy as np
import pandas as pd
import os

def load_processed_structures():
    """Load processed structures from cache"""
    cache_files = [
        "opsin_output/cache/structure_comparison_A.pkl",
        "opsin_output/cache/grn_assignment_A.pkl", 
        "opsin_output/cache/processed_structures_A.pkl"
    ]
    
    for cache_file in cache_files:
        if os.path.exists(cache_file):
            print(f"Trying to load from: {cache_file}")
            try:
                with open(cache_file, 'rb') as f:
                    data = pickle.load(f)
                
                if isinstance(data, dict) and 'processed_structures' in data:
                    return data['processed_structures']
                elif isinstance(data, dict) and len(data) > 100:
                    return data
            except Exception as e:
                print(f"Error loading {cache_file}: {e}")
                continue
    
    raise FileNotFoundError("No suitable processed structures cache found")

def get_grn_assigned_structures():
    """Get list of structures that have GRN assignments"""
    grn_file = "opsin_output/opsin_grn_tables/residue_table_grn.csv"
    grn_df = pd.read_csv(grn_file, index_col=0)
    return list(grn_df.index)

def main():
    print("=== Loading Structures ===")
    processed_structures = load_processed_structures()
    grn_structures = get_grn_assigned_structures()
    
    print(f"Found {len(processed_structures)} processed structures")
    print(f"Found {len(grn_structures)} GRN structures")
    
    # Check MerMAID1_model_0 and a few other structures
    test_structures = ['MerMAID1_model_0', 'OtHKR_model_0', '6CSN']
    
    for struct_id in test_structures:
        if struct_id in processed_structures:
            print(f"\n=== {struct_id} ===")
            struct_data = processed_structures[struct_id]
            
            # Check different dataframes
            for df_name in ['df', 'df_norm', 'df_ca_norm']:
                if df_name in struct_data:
                    df = struct_data[df_name]
                    if 'res_atom_name' in df.columns:
                        ca_df = df[df['res_atom_name'] == 'CA']
                    else:
                        ca_df = df
                    
                    if not ca_df.empty and all(col in ca_df.columns for col in ['x', 'y', 'z']):
                        coords = ca_df[['x', 'y', 'z']].astype(float).values
                        center = np.mean(coords, axis=0)
                        
                        print(f"  {df_name}: {len(coords)} atoms")
                        print(f"    Center: [{center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}]")
                        print(f"    Range X: [{np.min(coords[:, 0]):.1f}, {np.max(coords[:, 0]):.1f}]")
                        print(f"    Range Y: [{np.min(coords[:, 1]):.1f}, {np.max(coords[:, 1]):.1f}]")
                        print(f"    Range Z: [{np.min(coords[:, 2]):.1f}, {np.max(coords[:, 2]):.1f}]")
                    else:
                        print(f"  {df_name}: No CA atoms or missing coordinates")
                else:
                    print(f"  {df_name}: Not found")
        else:
            print(f"\n{struct_id}: Not found in processed structures")

if __name__ == "__main__":
    main()