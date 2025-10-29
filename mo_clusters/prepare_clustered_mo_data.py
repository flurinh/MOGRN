#!/usr/bin/env python3
"""
Process clustered MO structures:
Step 1: Load structures using protos
Step 2: Parse mapping and extract sequences
"""

import os
import pandas as pd
import numpy as np
import json
import shutil
from pathlib import Path
import sys
from typing import Dict, List, Tuple

# Add project root to path
project_root = Path(__file__).resolve().parent
sys.path.append(str(project_root))

# Import PROTOS modules (same as prepare_data.py)
from protos.io.paths.path_config import ProtosPaths, DataSource
from protos.processing.structure.struct_base_processor import CifBaseProcessor


def setup_processor(user_data_root):
    """
    Set up the CifBaseProcessor (simplified from prepare_data.py)
    """
    print("=== Setting up CifBaseProcessor ===")
    
    # Ensure user_data_root is a Path object and is absolute
    if isinstance(user_data_root, str):
        user_data_root = Path(user_data_root)
    user_data_root = user_data_root.resolve()
    
    print(f"Using user_data_root: {user_data_root}")
    
    # Initialize the processor
    cp = CifBaseProcessor(
        name="clustered_mo_processor",
        data_root=str(user_data_root),
        processor_data_dir="structure",
        preload=False
    )
    
    return cp


def process_clustered_mo_dataset(cp, user_data_root):
    """Process clustered MO structures from local files"""
    
    dataset_name = 'clustered_mo'
    source_path = project_root / "structures" / "clustered_mo"
    
    if not source_path.exists():
        print(f"Source path not found: {source_path}")
        return
    
    # Get all CIF files
    cif_files = list(source_path.glob("*.cif"))
    if not cif_files:
        print(f"No CIF files found in {source_path}")
        return
    
    print(f"Found {len(cif_files)} CIF files for {dataset_name}")
    
    # Create the destination directory in the PROTOS structure
    mmcif_dir = user_data_root / "structure" / "mmcif"
    mmcif_dir.mkdir(exist_ok=True, parents=True)
    
    # Copy files to the PROTOS structure directory
    pdb_ids = []
    successful_copies = 0
    
    for cif_file in cif_files:
        # Use the full filename as PDB ID to avoid issues
        pdb_id = cif_file.stem
        
        # Copy the file to the PROTOS directory
        dest_file = mmcif_dir / cif_file.name
        try:
            shutil.copy2(cif_file, dest_file)
            print(f"Copied {cif_file.name} to {dest_file}")
            successful_copies += 1
            pdb_ids.append(pdb_id)
        except Exception as e:
            print(f"Error copying {cif_file.name}: {e}")
    
    print(f"Copied {successful_copies} files")
    
    if not pdb_ids:
        print(f"No structures were successfully copied for {dataset_name}")
        return
    
    # Load structures into the processor
    try:
        # Reset processor to ensure a clean state
        cp.reset_data()
        
        # Load the structures
        cp.load_structures(pdb_ids)
        
        # Check if any structures were loaded
        print(f"Loaded {len(cp.pdb_ids)} structures")
        
        if len(cp.pdb_ids) > 0:
            # Create dataset with metadata
            metadata = {
                "description": "Clustered microbial opsin structures",
                "source": str(source_path),
                "creation_date": pd.Timestamp.now().strftime('%Y-%m-%d'),
                "structures_count": len(cp.pdb_ids)
            }
            
            # Create the dataset
            cp.create_dataset(
                dataset_id=dataset_name,
                name=dataset_name,
                description=metadata['description'],
                content=cp.pdb_ids,
                metadata=metadata
            )
            
            print(f"Successfully created dataset {dataset_name} with {len(cp.pdb_ids)} structures")
            
            # List the PDB IDs that were loaded
            print("\nLoaded PDB IDs:")
            for pdb_id in cp.pdb_ids[:5]:  # Show first 5
                print(f"  - {pdb_id}")
            if len(cp.pdb_ids) > 5:
                print(f"  ... and {len(cp.pdb_ids) - 5} more")
                
        else:
            print(f"No structures were successfully loaded for {dataset_name}")
            
    except Exception as e:
        print(f"Error processing dataset {dataset_name}: {e}")
        import traceback
        traceback.print_exc()


def parse_mapping_file(mapping_file: Path) -> Dict[str, str]:
    """Parse the mapping file to get MO_XXX to full name mapping."""
    mapping = {}
    
    with open(mapping_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            # Split by tab to get MO_ID and full name
            parts = line.split('\t')
            if len(parts) >= 2:
                mo_id = parts[0]
                full_name = parts[1]
                mapping[mo_id] = full_name
    
    print(f"Parsed {len(mapping)} mappings from {mapping_file}")
    return mapping


def extract_sequence_from_structure(cp, pdb_id: str, chain_id: str = 'A') -> str:
    """Extract amino acid sequence from a structure for a specific chain."""
    try:
        # Get data for this PDB ID
        df = cp.data[cp.data['pdb_id'] == pdb_id].copy()
        
        # Filter by chain
        df_chain = df[df['auth_chain_id'] == chain_id].copy()
        
        if df_chain.empty:
            print(f"Warning: Chain {chain_id} not found in {pdb_id}")
            return ""
        
        # Three-letter to one-letter amino acid code mapping
        aa_3to1 = {
            'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D',
            'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 'GLY': 'G',
            'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K',
            'MET': 'M', 'PHE': 'F', 'PRO': 'P', 'SER': 'S',
            'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
            'UNK': 'X'  # Unknown amino acid
        }
        
        # Get unique residues sorted by residue number
        # Group by residue number and name to get unique residues
        residues = df_chain.groupby(['auth_seq_id', 'res_name3l']).first().reset_index()
        residues = residues.sort_values('auth_seq_id')
        
        sequence = []
        for _, res in residues.iterrows():
            res_name = res['res_name3l'].upper()
            # Only include standard amino acids
            if res_name in aa_3to1:
                sequence.append(aa_3to1[res_name])
            elif len(res_name) == 3 and res_name not in ['HOH', 'WAT']:  # Exclude water
                sequence.append('X')
        
        return ''.join(sequence)
        
    except Exception as e:
        print(f"Error extracting sequence from {pdb_id}: {e}")
        import traceback
        traceback.print_exc()
        return ""


def process_structures_step2(cp, mapping_file: Path):
    """Step 2: Extract sequences and prepare data for property table."""
    
    # Parse the mapping file
    name_mapping = parse_mapping_file(mapping_file)
    
    # Process each structure
    results = []
    
    print("\nExtracting sequences from structures...")
    for i, pdb_id in enumerate(cp.pdb_ids):
        if i % 5 == 0:
            print(f"Processing {i+1}/{len(cp.pdb_ids)}...")
        
        # Extract MO_ID from pdb_id (e.g., MO_001_model_0 -> MO_001)
        mo_id = pdb_id.split('_model_')[0]
        
        # Get full name from mapping
        full_name = name_mapping.get(mo_id, "")
        
        # Extract sequence
        sequence = extract_sequence_from_structure(cp, pdb_id, 'A')
        
        # Store results
        result = {
            'pdb_id': pdb_id,
            'mo_id': mo_id,
            'full_name': full_name,
            'sequence': sequence,
            'seq_length': len(sequence)
        }
        results.append(result)
    
    return results


def main():
    """Main function - Process clustered MO structures"""
    
    # Set up paths
    user_data_root = project_root / "data_clustered_mo"
    mapping_file = project_root / "property" / "mo_small_name_mapping.txt"
    
    # Step 1: Setup and load structures
    print("Step 1: Setting up and loading clustered MO structures")
    print("=" * 60)
    
    # Ensure user_data_root exists
    user_data_root.mkdir(parents=True, exist_ok=True)
    
    # Setup processor
    try:
        cp = setup_processor(user_data_root)
        print("Processor setup complete")
    except Exception as e:
        print(f"Failed to initialize processor: {e}")
        return
    
    # Process the clustered_mo dataset
    process_clustered_mo_dataset(cp, user_data_root)
    
    # Test loading the dataset back
    print("\n" + "=" * 60)
    print("Testing: Loading the dataset back")
    try:
        cp.reset_data()
        cp.load_dataset('clustered_mo')
        print(f"Successfully loaded dataset 'clustered_mo' with {len(cp.pdb_ids)} structures")
    except Exception as e:
        print(f"Error loading dataset: {e}")
        return
    
    print("\nStep 1 complete!")
    
    # Step 2: Parse mapping and extract sequences
    print("\n" + "=" * 60)
    print("Step 2: Parse mapping and extract sequences")
    print("=" * 60)
    
    results = process_structures_step2(cp, mapping_file)
    
    # Display some results
    print(f"\nProcessed {len(results)} structures")
    print("\nSample results:")
    for r in results[:3]:
        print(f"  {r['mo_id']}: {r['full_name'][:50]}... (seq_length: {r['seq_length']})")
    
    print("\nStep 2 complete!")
    
    # Save intermediate results for inspection
    intermediate_file = project_root / "clustered_mo_intermediate.json"
    with open(intermediate_file, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nIntermediate results saved to: {intermediate_file}")


if __name__ == "__main__":
    main()