#!/usr/bin/env python3
"""
Inspect the clustered MO structures to understand their format.
"""

import os
import pandas as pd
from pathlib import Path
import sys

# Add project root to path
project_root = Path(__file__).resolve().parent
sys.path.append(str(project_root))

# Add protos to path if needed
protos_path = project_root / "protos" / "src"
if protos_path.exists():
    sys.path.insert(0, str(protos_path))

# Import PROTOS modules
from protos.processing.structure.struct_base_processor import CifBaseProcessor


def inspect_structures():
    """Inspect the clustered structures to see what's in them."""
    
    user_data_root = project_root / "data_clustered_mo"
    
    # Initialize processor
    cp = CifBaseProcessor(
        name="clustered_mo_processor",
        data_root=str(user_data_root),
        processor_data_dir="structure",
        preload=False
    )
    
    # Load the dataset
    cp.load_dataset('clustered_mo')
    print(f"Loaded dataset with {len(cp.pdb_ids)} structures")
    
    # Inspect first structure in detail
    if len(cp.pdb_ids) > 0:
        first_id = cp.pdb_ids[0]
        print(f"\nInspecting structure: {first_id}")
        
        # Get data for this structure
        df = cp.data[cp.data['pdb_id'] == first_id].copy()
        print(f"Total atoms: {len(df)}")
        
        # Check columns
        print(f"\nColumns: {df.columns.tolist()}")
        
        # Check chain IDs
        if 'auth_chain_id' in df.columns:
            chains = df['auth_chain_id'].unique()
            print(f"\nChain IDs (auth_chain_id): {chains}")
            for chain in chains:
                n_atoms = len(df[df['auth_chain_id'] == chain])
                n_ca = len(df[(df['auth_chain_id'] == chain) & (df['atom_name'] == 'CA')])
                print(f"  Chain {chain}: {n_atoms} atoms, {n_ca} CA atoms")
        
        # Check alternative chain column names
        for col in ['chain_id', 'label_asym_id', 'struct_asym_id', 'pdbx_strand_id']:
            if col in df.columns:
                chains = df[col].unique()
                print(f"\nChain IDs ({col}): {chains}")
                for chain in chains[:3]:  # Show first 3
                    n_atoms = len(df[df[col] == chain])
                    print(f"  Chain {chain}: {n_atoms} atoms")
        
        # Check atom names
        atom_names = df['atom_name'].value_counts().head(10)
        print(f"\nTop 10 atom names:")
        print(atom_names)
        
        # Show sample of CA atoms
        ca_atoms = df[df['atom_name'] == 'CA']
        if len(ca_atoms) > 0:
            print(f"\nFound {len(ca_atoms)} CA atoms")
            print("Sample CA atoms:")
            print(ca_atoms[['atom_name', 'auth_chain_id', 'auth_seq_id', 'res_name3l', 'x', 'y', 'z']].head())
        else:
            print("\nNo CA atoms found!")
            # Check for alternative naming
            for alt_name in ['C-alpha', 'CA ', ' CA', 'C_alpha']:
                alt_ca = df[df['atom_name'] == alt_name]
                if len(alt_ca) > 0:
                    print(f"Found {len(alt_ca)} atoms with name '{alt_name}'")
        
        # Check res_atom_name column
        if 'res_atom_name' in df.columns:
            print("\nChecking res_atom_name column:")
            res_atom_names = df['res_atom_name'].value_counts().head(10)
            print(res_atom_names)
            
            # Check for CA in res_atom_name
            ca_in_res = df[df['res_atom_name'].str.contains('CA', na=False)]
            if len(ca_in_res) > 0:
                print(f"\nFound {len(ca_in_res)} atoms with 'CA' in res_atom_name")
                print("Sample:")
                print(ca_in_res[['res_atom_name', 'atom_name', 'auth_chain_id', 'res_name3l']].head())
        
        # Show sample of all atoms in chain A
        print("\nSample of atoms in chain A:")
        chain_a = df[df['auth_chain_id'] == 'A']
        print(chain_a[['atom_name', 'res_atom_name', 'res_name3l', 'auth_seq_id']].head(20))
    
    # Summary for all structures
    print("\n" + "="*60)
    print("Summary for all structures:")
    
    for i, pdb_id in enumerate(cp.pdb_ids[:5]):  # First 5
        df = cp.data[cp.data['pdb_id'] == pdb_id].copy()
        
        # Get chain info
        if 'auth_chain_id' in df.columns:
            chains = df['auth_chain_id'].unique()
            ca_counts = {}
            for chain in chains:
                n_ca = len(df[(df['auth_chain_id'] == chain) & (df['atom_name'] == 'CA')])
                if n_ca > 0:
                    ca_counts[chain] = n_ca
            
            print(f"\n{pdb_id}:")
            print(f"  Chains: {list(chains)}")
            print(f"  CA atoms per chain: {ca_counts}")


if __name__ == "__main__":
    inspect_structures()