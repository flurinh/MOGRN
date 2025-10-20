#!/usr/bin/env python3
"""
Step-by-step alignment of clustered MO structures to processed structures.
This version carefully validates data at each step.
"""

import os
import numpy as np
import pandas as pd
import pickle
import json
from pathlib import Path
import sys
from typing import Dict, List, Tuple
from tqdm import tqdm

# Add project root to path
project_root = Path(__file__).resolve().parent
sys.path.append(str(project_root))

# Add protos to path if needed
protos_path = project_root / "protos" / "src"
if protos_path.exists():
    sys.path.insert(0, str(protos_path))

# Import PROTOS modules
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.processing.structure.struct_alignment import get_structure_alignment
from Bio.PDB.qcprot import QCPSuperimposer

from src.data_processing import (
    load_opsin_structures,
    filter_structures_by_chain_and_retinal
)


def validate_structure_data(df: pd.DataFrame, pdb_id: str) -> bool:
    """Validate that structure data has required columns and content."""
    required_cols = ['auth_chain_id', 'x', 'y', 'z', 'auth_seq_id']
    
    # Check columns exist
    for col in required_cols:
        if col not in df.columns:
            print(f"  WARNING: {pdb_id} missing column: {col}")
            return False
    
    # Check which column has CA atoms
    ca_atoms = None
    if 'res_atom_name' in df.columns:
        ca_atoms = df[(df['res_atom_name'] == 'CA') & (df['auth_chain_id'] == 'A')]
    elif 'atom_name' in df.columns:
        ca_atoms = df[(df['atom_name'] == 'CA') & (df['auth_chain_id'] == 'A')]
    
    if ca_atoms is None or len(ca_atoms) == 0:
        print(f"  WARNING: {pdb_id} has no CA atoms in chain A")
        return False
    
    # Check coordinates are numeric
    try:
        coords = ca_atoms[['x', 'y', 'z']].values.astype(float)
        if np.any(np.isnan(coords)):
            print(f"  WARNING: {pdb_id} has NaN coordinates")
            return False
    except:
        print(f"  WARNING: {pdb_id} has non-numeric coordinates")
        return False
    
    return True


def load_and_validate_clustered_structures(user_data_root: Path):
    """Load clustered MO structures and validate them."""
    print("Loading clustered MO structures...")
    
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
    
    # Extract and validate each structure
    structures = {}
    for pdb_id in tqdm(cp.pdb_ids, desc="Validating clustered structures"):
        df = cp.data[cp.data['pdb_id'] == pdb_id].copy()
        
        if validate_structure_data(df, pdb_id):
            # Filter to chain A only
            df_chain_a = df[df['auth_chain_id'] == 'A'].copy()
            structures[pdb_id] = df_chain_a
    
    print(f"Successfully validated {len(structures)}/{len(cp.pdb_ids)} structures")
    return structures


def load_and_validate_processed_structures(cache_dir: Path, chain_id: str = 'A'):
    """Load and validate processed structures from cache."""
    print("\nLoading processed structures from cache...")
    
    processed_cache_path = cache_dir / f"processed_structures_{chain_id}.pkl"
    
    if not processed_cache_path.exists():
        print("ERROR: No processed structures cache found!")
        print(f"Expected at: {processed_cache_path}")
        return None
    
    with open(processed_cache_path, 'rb') as f:
        data = pickle.load(f)
    
    processed_structures = data.get('processed_structures', {})
    print(f"Loaded {len(processed_structures)} processed structures from cache")
    
    # Validate each structure
    validated = {}
    for pdb_id, struct_data in tqdm(processed_structures.items(), desc="Validating processed structures"):
        # Extract dataframe
        if isinstance(struct_data, dict) and 'df' in struct_data:
            df = struct_data['df']
        else:
            df = struct_data
        
        if validate_structure_data(df, pdb_id):
            validated[pdb_id] = df
    
    print(f"Successfully validated {len(validated)}/{len(processed_structures)} structures")
    return validated


def calculate_rmsd_simple(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Calculate RMSD between two coordinate sets."""
    if len(coords1) != len(coords2) or len(coords1) < 3:
        return float('inf')
    
    try:
        sup = QCPSuperimposer()
        sup.set(coords2, coords1)
        sup.run()
        return sup.get_rms()
    except Exception as e:
        print(f"    RMSD calculation error: {e}")
        return float('inf')


def get_ca_atoms(df: pd.DataFrame) -> pd.DataFrame:
    """Get CA atoms from structure dataframe, handling different column names."""
    if 'res_atom_name' in df.columns:
        return df[df['res_atom_name'] == 'CA'].copy()
    elif 'atom_name' in df.columns:
        return df[df['atom_name'] == 'CA'].copy()
    else:
        return pd.DataFrame()  # Empty dataframe


def align_single_pair(struct1_df: pd.DataFrame, struct2_df: pd.DataFrame, 
                     struct1_id: str, struct2_id: str, debug: bool = False) -> float:
    """Align a single pair of structures and return RMSD."""
    
    # Get CA atoms using the helper function
    ca1 = get_ca_atoms(struct1_df)
    ca2 = get_ca_atoms(struct2_df)
    
    if len(ca1) == 0 or len(ca2) == 0:
        return float('inf')
    
    # Sort by residue number
    ca1 = ca1.sort_values('auth_seq_id')
    ca2 = ca2.sort_values('auth_seq_id')
    
    # Get coordinates
    coords1 = ca1[['x', 'y', 'z']].values.astype(float)
    coords2 = ca2[['x', 'y', 'z']].values.astype(float)
    
    if debug:
        print(f"  {struct1_id}: {len(coords1)} CA atoms")
        print(f"  {struct2_id}: {len(coords2)} CA atoms")
    
    # Method 1: Direct alignment (if same length)
    if len(coords1) == len(coords2):
        rmsd = calculate_rmsd_simple(coords1, coords2)
        if debug:
            print(f"  Direct alignment RMSD: {rmsd:.2f}")
        if rmsd < 50.0:  # Reasonable RMSD
            return rmsd
    
    # Method 2: Truncate to shorter length
    min_len = min(len(coords1), len(coords2))
    if min_len >= 50:  # At least 50 residues
        rmsd = calculate_rmsd_simple(coords1[:min_len], coords2[:min_len])
        if debug:
            print(f"  Truncated alignment RMSD ({min_len} residues): {rmsd:.2f}")
        if rmsd < 50.0:
            return rmsd
    
    # Method 3: Try sequence-based alignment
    try:
        alignment = get_structure_alignment(ca1, ca2, method='sequence')
        if alignment is not None and len(alignment) >= 50:
            # Extract aligned coordinates
            aligned_coords1 = []
            aligned_coords2 = []
            
            for idx1, idx2 in alignment:
                if idx1 < len(coords1) and idx2 < len(coords2):
                    aligned_coords1.append(coords1[idx1])
                    aligned_coords2.append(coords2[idx2])
            
            if len(aligned_coords1) >= 50:
                aligned_coords1 = np.array(aligned_coords1)
                aligned_coords2 = np.array(aligned_coords2)
                rmsd = calculate_rmsd_simple(aligned_coords1, aligned_coords2)
                if debug:
                    print(f"  Sequence alignment RMSD ({len(aligned_coords1)} residues): {rmsd:.2f}")
                return rmsd
    except Exception as e:
        if debug:
            print(f"  Sequence alignment failed: {e}")
    
    return float('inf')


def perform_all_vs_all_alignment(clustered_structures: Dict, processed_structures: Dict):
    """Perform all vs all alignment with progress tracking."""
    
    clustered_ids = list(clustered_structures.keys())
    processed_ids = list(processed_structures.keys())
    
    n_clustered = len(clustered_ids)
    n_processed = len(processed_ids)
    
    print(f"\nPerforming {n_clustered} x {n_processed} = {n_clustered * n_processed} alignments")
    
    # Initialize RMSD matrix
    rmsd_matrix = np.full((n_clustered, n_processed), np.inf)
    
    # Results dictionary
    results = {}
    
    # Create progress bar for total alignments
    with tqdm(total=n_clustered * n_processed, desc="Aligning structures") as pbar:
        for i, clustered_id in enumerate(clustered_ids):
            clustered_df = clustered_structures[clustered_id]
            
            best_match_id = None
            best_rmsd = float('inf')
            valid_alignments = []
            
            # Test first alignment in detail
            debug = (i == 0)
            
            for j, processed_id in enumerate(processed_ids):
                processed_df = processed_structures[processed_id]
                
                # Perform alignment
                rmsd = align_single_pair(
                    clustered_df, processed_df,
                    clustered_id, processed_id,
                    debug=(debug and j == 0)
                )
                
                if rmsd != float('inf'):
                    rmsd_matrix[i, j] = rmsd
                    valid_alignments.append((processed_id, rmsd))
                    
                    if rmsd < best_rmsd:
                        best_rmsd = rmsd
                        best_match_id = processed_id
                
                pbar.update(1)
            
            # Store results for this clustered structure
            results[clustered_id] = {
                'best_match_id': best_match_id,
                'best_rmsd': best_rmsd,
                'n_valid_alignments': len(valid_alignments),
                'all_alignments': dict(valid_alignments)
            }
            
            # Print summary for this structure
            if best_match_id:
                print(f"{clustered_id}: Best match = {best_match_id} (RMSD = {best_rmsd:.2f} Å)")
                print(f"  Valid alignments: {len(valid_alignments)}/{n_processed}")
    
    return results, rmsd_matrix


def main():
    """Main function."""
    
    # Set up paths
    cache_dir = project_root / "opsin_output" / "cache"
    user_data_root = project_root / "data_clustered_mo"
    
    # Step 1: Load and validate clustered structures
    print("STEP 1: Load clustered structures")
    print("=" * 60)
    clustered_structures = load_and_validate_clustered_structures(user_data_root)
    
    if not clustered_structures:
        print("ERROR: No valid clustered structures found!")
        return
    
    # Step 2: Load and validate processed structures
    print("\nSTEP 2: Load processed structures")
    print("=" * 60)
    processed_structures = load_and_validate_processed_structures(cache_dir)
    
    if not processed_structures:
        print("ERROR: No valid processed structures found!")
        return
    
    # Step 3: Perform alignments
    print("\nSTEP 3: Perform alignments")
    print("=" * 60)
    results, rmsd_matrix = perform_all_vs_all_alignment(clustered_structures, processed_structures)
    
    # Step 4: Save results
    print("\nSTEP 4: Save results")
    print("=" * 60)
    
    # Save detailed results
    output_file = project_root / "clustered_mo_alignment_results_v2.json"
    output_data = {
        'summary': {
            'n_clustered': len(clustered_structures),
            'n_processed': len(processed_structures),
            'clustered_ids': list(clustered_structures.keys()),
            'processed_ids': list(processed_structures.keys())
        },
        'alignments': results
    }
    
    with open(output_file, 'w') as f:
        json.dump(output_data, f, indent=2)
    print(f"Results saved to: {output_file}")
    
    # Save RMSD matrix
    rmsd_csv_file = project_root / "clustered_mo_rmsd_matrix_v2.csv"
    rmsd_df = pd.DataFrame(
        rmsd_matrix,
        index=list(clustered_structures.keys()),
        columns=list(processed_structures.keys())
    )
    rmsd_df.to_csv(rmsd_csv_file)
    print(f"RMSD matrix saved to: {rmsd_csv_file}")
    
    # Print summary
    print("\nSUMMARY")
    print("=" * 60)
    n_successful = sum(1 for r in results.values() if r['best_match_id'] is not None)
    print(f"Successful alignments: {n_successful}/{len(results)}")
    
    # Show distribution of valid alignments per structure
    valid_counts = [r['n_valid_alignments'] for r in results.values()]
    if valid_counts:
        print(f"Valid alignments per structure: min={min(valid_counts)}, "
              f"max={max(valid_counts)}, avg={np.mean(valid_counts):.1f}")
    
    # Show best matches
    print("\nBest matches:")
    for clustered_id, result in results.items():
        if result['best_match_id']:
            print(f"  {clustered_id} -> {result['best_match_id']} (RMSD: {result['best_rmsd']:.2f} Å)")


if __name__ == "__main__":
    main()