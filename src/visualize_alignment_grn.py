#!/usr/bin/env python3
"""
Interactive GRN-based visualization with slider for exploring opsin structures.
Shows all aligned structures with GRN position highlighting via slider.
"""

import pickle
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly import colors
import os
from pathlib import Path
from sklearn.decomposition import PCA
from collections import Counter, defaultdict

# Import protos for structure loading (using new StructureProcessor API)
from protos.processing.structure import StructureProcessor
from src.data_processing import load_experimental_dataset

# Import helix color scheme
from src.opsin_color_scheme import HELIX_NUMBER_COLORS, get_categorical_colors


# Compatibility class for pickle loading
class DatasetCompat:
    """Compatibility class for loading pickled data created by opsin_analysis_workflow."""
    def __init__(self, pdb_ids, data=None):
        self.pdb_ids = list(pdb_ids)
        self.data = data if data is not None else pd.DataFrame()

    def format_data_types(self):
        if self.data.empty:
            return
        for col in ['x', 'y', 'z']:
            if col in self.data.columns:
                self.data[col] = pd.to_numeric(self.data[col], errors='coerce')


# Basic 3-letter to 1-letter amino acid mapping for hover text enrichment
THREE_TO_ONE_AA = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "RET": "RET", "HOH": "HOH"
}


def load_rmsd_cache(cache_dir="opsin_output/cache"):
    """Load alignment paths and rotation matrices from RMSD cache"""
    cache_path = Path(cache_dir)
    if not cache_path.exists():
        raise FileNotFoundError(f"Cache directory not found: {cache_path}")

    cache_files = sorted(cache_path.glob("rmsd_cache_*.pkl"))
    if not cache_files:
        raise FileNotFoundError(f"No RMSD cache files found in {cache_path}")

    cache_file = cache_files[0]
    print(f"Loading cache from: {cache_file}")

    with open(cache_file, 'rb') as f:
        cache_data = pickle.load(f)

    return cache_data


def load_processed_structures(cache_dir="opsin_output/cache"):
    """Load processed structures from cache"""
    cache_path = Path(cache_dir)
    candidate_files = [
        cache_path / "structure_comparison_A.pkl",
        cache_path / "grn_assignment_A.pkl",
        cache_path / "processed_structures_A.pkl"
    ]

    for cache_file in candidate_files:
        if cache_file.exists():
            print(f"Trying to load from: {cache_file}")
            try:
                with open(cache_file, 'rb') as f:
                    data = pickle.load(f)

                if isinstance(data, dict) and 'processed_structures' in data:
                    processed_structures = data['processed_structures']
                    print(f"Found {len(processed_structures)} processed structures in {cache_file}")
                    return processed_structures
                elif isinstance(data, dict) and len(data) > 100:
                    print(f"Found {len(data)} structures directly in {cache_file}")
                    return data
            except Exception as e:
                print(f"Error loading {cache_file}: {e}")
                continue

    raise FileNotFoundError(f"No suitable processed structures cache found in {cache_path}")


def load_alignment_transformations(cache_dir="opsin_output/cache", global_ref="7bmh"):
    """Load pairwise alignment transformations from RMSD cache"""
    cache_path = Path(cache_dir)
    rmsd_cache_files = list(cache_path.glob("rmsd_cache*.pkl"))

    if not rmsd_cache_files:
        print("[WARN] No RMSD cache found for alignment transformations")
        return {}

    with open(rmsd_cache_files[0], 'rb') as f:
        rmsd_cache = pickle.load(f)

    alignment_paths = rmsd_cache.get('alignment_paths', {})

    # Extract transformations to global reference
    transformations = {}
    for (src, tgt), data in alignment_paths.items():
        if tgt == global_ref:
            transformations[src] = {
                'rotation': np.array(data['rotation']),
                'translation': np.array(data['translation']),
                'rmsd': data.get('rmsd', np.nan)
            }

    # Add identity for global ref itself
    transformations[global_ref] = {
        'rotation': np.eye(3),
        'translation': np.zeros(3),
        'rmsd': 0.0
    }

    print(f"Loaded {len(transformations)} alignment transformations to {global_ref}")
    return transformations


def load_grn_table(grn_file="opsin_output/curated_grn_postprocessed.csv"):
    """Load and parse the GRN table"""
    grn_path = Path(grn_file)

    print(f"Reading GRN table from: {grn_path}")
    # Use dtype={0: str} to prevent scientific notation parsing (e.g., "1e12" -> 1000000000000)
    grn_df = pd.read_csv(grn_path, index_col=0, dtype={0: str})
    grn_df.index = grn_df.index.astype(str)

    print(f"GRN table: {len(grn_df)} structures, {len(grn_df.columns)} GRN positions")
    return grn_df


def parse_grn_residue(cell_value):
    """
    Parse GRN table cell value to extract amino acid and sequence position

    Args:
        cell_value: Cell value like "M1", "V2", "P3", or "-"

    Returns:
        tuple: (amino_acid, seq_pos) or (None, None) for gaps
    """
    if pd.isna(cell_value) or cell_value == '-':
        return None, None

    cell_str = str(cell_value).strip()
    if not cell_str or cell_str == '-':
        return None, None

    # Extract amino acid (first character) and sequence position (remaining digits)
    amino_acid = cell_str[0]
    seq_pos_str = cell_str[1:]

    try:
        seq_pos = int(seq_pos_str)
        return amino_acid, seq_pos
    except ValueError:
        return None, None


def extract_ca_coordinates_with_grn(processed_structures, grn_df, chain_id='A', use_helix_only=True,
                                     transformations=None, global_ref="7bmh"):
    """
    Extract CA coordinates with GRN mapping for structures in GRN table.

    Args:
        processed_structures: Dict of processed structure data
        grn_df: GRN table DataFrame with structure IDs as index
        chain_id: Chain to extract (default 'A')
        use_helix_only: If True, filter to helix residues only
        transformations: Dict of {struct_id: {'rotation': R, 'translation': t}} to align to global ref
        global_ref: Global reference structure ID (default '7bmh')
    """
    all_structures = {}
    grn_structures = list(grn_df.index)

    # Create case-insensitive lookup for processed structures
    proc_lower_map = {k.lower(): k for k in processed_structures.keys()}

    # Also create case-insensitive lookup for transformations if provided
    trans_lower_map = {}
    if transformations:
        trans_lower_map = {k.lower(): k for k in transformations.keys()}

    for struct_id in grn_structures:
        # Try direct match first, then case-insensitive
        if struct_id in processed_structures:
            proc_id = struct_id
        elif struct_id.lower() in proc_lower_map:
            proc_id = proc_lower_map[struct_id.lower()]
        else:
            print(f"Warning: {struct_id} has GRN assignment but not found in processed structures")
            continue

        struct_data = processed_structures[proc_id]

        # Use the same dataframe as RMSD calculation: 'df'
        if 'df' not in struct_data:
            print(f"Warning: No 'df' found for {struct_id}")
            continue

        struct_df = struct_data['df'].copy()

        # Apply the EXACT same filtering as in compute_all_vs_all_rmsd_improved()
        struct_df = struct_df[struct_df['auth_chain_id'] == chain_id]
        struct_df = struct_df[struct_df['group'] == 'ATOM']

        # Filter for helix residues if requested
        if use_helix_only and 'helix_num' in struct_df.columns:
            struct_helices = struct_df[struct_df['helix_num'] > 0]
            if len(struct_helices) > 0:
                struct_df = struct_helices
            else:
                print(f"[INFO] No TM helix residues found for {struct_id}. Using all protein residues instead.")
        elif use_helix_only:
            print(f"[INFO] helix annotations missing for {struct_id}; inferring after GRN mapping.")

        # Filter for CA atoms only
        ca_data = struct_df[struct_df['res_atom_name'] == 'CA']

        if not ca_data.empty:
            # Create GRN mapping for this structure
            grn_row = grn_df.loc[struct_id]
            seq_to_grn = {}  # Maps auth_seq_id -> GRN position
            grn_to_seq = {}  # Maps GRN position -> auth_seq_id

            for grn_pos, cell_value in grn_row.items():
                amino_acid, seq_pos = parse_grn_residue(cell_value)
                if amino_acid is not None and seq_pos is not None:
                    seq_to_grn[seq_pos] = grn_pos
                    grn_to_seq[grn_pos] = seq_pos

            # Add GRN column to CA data
            ca_data_with_grn = ca_data.copy()
            ca_data_with_grn['grn_position'] = ca_data_with_grn['auth_seq_id'].map(seq_to_grn)
            ca_data_with_grn['grn'] = ca_data_with_grn['grn_position']  # Direct GRN access

            # Enrich residue names for hover text
            if 'res_name1l' not in ca_data_with_grn.columns:
                if 'res_name3l' in ca_data_with_grn.columns:
                    ca_data_with_grn['res_name1l'] = (
                        ca_data_with_grn['res_name3l']
                        .astype(str)
                        .str.upper()
                        .map(THREE_TO_ONE_AA)
                        .fillna('X')
                    )
                else:
                    ca_data_with_grn['res_name1l'] = 'X'
            else:
                ca_data_with_grn['res_name1l'] = ca_data_with_grn['res_name1l'].fillna('X')

            if 'res_name3l' not in ca_data_with_grn.columns:
                ca_data_with_grn['res_name3l'] = ca_data_with_grn['res_name1l'].map(
                    {v: k for k, v in THREE_TO_ONE_AA.items() if len(v) == 1}
                ).fillna('UNK')

            # Determine dataset type
            dataset_type = 'predicted' if struct_id.endswith('_model_0') else 'experimental'

            # Add helix information (check if helix_num column exists in the filtered data)
            if 'helix_num' in ca_data.columns:
                ca_data_with_grn['helix_num'] = ca_data['helix_num'].fillna(0).astype(int)
            else:
                # Derive helix index from GRN label when explicit assignment is unavailable
                def _infer_helix(grn_label: str) -> int:
                    if pd.isna(grn_label):
                        return 0
                    label = str(grn_label).strip()
                    if not label:
                        return 0
                    prefix = label.split('.')[0]
                    if prefix.isdigit():
                        try:
                            helix_idx = int(prefix)
                            return helix_idx if 1 <= helix_idx <= 7 else 0
                        except ValueError:
                            return 0
                    return 0

                ca_data_with_grn['helix_num'] = ca_data_with_grn['grn_position'].apply(_infer_helix).astype(int)
                print(f"[INFO] Inferred helix numbers from GRN labels for {struct_id}")

            if use_helix_only:
                helix_mask = ca_data_with_grn['helix_num'] > 0
                if helix_mask.any():
                    ca_data_with_grn = ca_data_with_grn[helix_mask]
                else:
                    print(f"[INFO] Helix inference yielded no TM residues for {struct_id}; retaining all residues.")

            if ca_data_with_grn.empty:
                print(f"[WARN] No CA residues remain for {struct_id} after filtering; skipping structure.")
                continue

            # Capture retinal atoms if available for downstream visualization
            retinal_info = None
            df_ret = struct_data.get('df_ret')
            if isinstance(df_ret, pd.DataFrame) and not df_ret.empty:
                try:
                    retinal_df = df_ret.copy()
                    retinal_coords = retinal_df[['x', 'y', 'z']].astype(float).values
                    retinal_info = {
                        'coords': retinal_coords,
                        'atom_names': retinal_df['atom_name'].astype(str).values,
                        'res_name3l': retinal_df.get('res_name3l', pd.Series(['RET'] * len(retinal_df))).astype(str).values,
                        'chain_ids': retinal_df.get('auth_chain_id', pd.Series([''] * len(retinal_df))).astype(str).values
                    }
                except Exception as retina_exc:
                    print(f"[WARN] Failed to capture retinal coordinates for {struct_id}: {retina_exc}")
                    retinal_info = None

            # Extract raw coordinates
            raw_coords = ca_data_with_grn[['x', 'y', 'z']].astype(float).values

            # Apply alignment transformation if provided
            aligned_coords = raw_coords
            alignment_applied = False
            if transformations:
                # Try to find transformation for this structure
                trans_key = None
                if proc_id in transformations:
                    trans_key = proc_id
                elif proc_id.lower() in trans_lower_map:
                    trans_key = trans_lower_map[proc_id.lower()]
                elif struct_id in transformations:
                    trans_key = struct_id
                elif struct_id.lower() in trans_lower_map:
                    trans_key = trans_lower_map[struct_id.lower()]

                if trans_key:
                    trans = transformations[trans_key]
                    R = trans['rotation']
                    t = trans['translation']
                    # Apply transformation: coords_aligned = coords @ R.T + t
                    aligned_coords = raw_coords @ R.T + t
                    alignment_applied = True

                    # Also transform retinal if present
                    if retinal_info is not None:
                        retinal_info['coords'] = retinal_info['coords'] @ R.T + t

            # Use proc_id (lowercase) as key to match alignment paths which use lowercase IDs
            all_structures[proc_id] = {
                'coords': aligned_coords,
                'coords_raw': raw_coords,  # Keep raw coords for reference
                'residues': ca_data_with_grn['auth_seq_id'].values,
                'grn_positions': ca_data_with_grn['grn_position'].values,
                'grn': ca_data_with_grn['grn'].values,  # Direct GRN access
                'helix_numbers': ca_data_with_grn['helix_num'].values,
                'res_name1l': ca_data_with_grn['res_name1l'].values,
                'res_name3l': ca_data_with_grn['res_name3l'].values,
                'dataset': dataset_type,
                'structure_type': struct_data.get('structure_type', dataset_type),
                'seq_to_grn': seq_to_grn,
                'grn_to_seq': grn_to_seq,
                'dataframe': ca_data_with_grn,  # Store the full dataframe for efficient access
                'retinal': retinal_info,
                'grn_table_id': struct_id,  # Store original GRN table ID for reference
                'alignment_applied': alignment_applied
            }
        else:
            print(f"Warning: No CA atoms found for {struct_id} after filtering")

    print(f"Extracted CA coordinates with GRN mapping for {len(all_structures)} structures")
    return all_structures


def apply_alignment_transformations(structures, alignment_paths, reference_id='6xl3'):
    """Apply rotation matrices to align all structures to the global reference"""
    print(f"Using {reference_id} as global reference structure")

    # Also try lowercase reference_id
    if reference_id not in structures:
        if reference_id.lower() in structures:
            reference_id = reference_id.lower()
        elif reference_id.upper() in structures:
            reference_id = reference_id.upper()
        else:
            print(f"Warning: Reference {reference_id} not found in structures. Using first available.")
            reference_id = list(structures.keys())[0]
            print(f"Using {reference_id} as reference instead")

    aligned_structures = {}

    # Create case-insensitive lookup for alignment paths
    align_lower_map = {}
    for key in alignment_paths.keys():
        lower_key = (key[0].lower(), key[1].lower())
        align_lower_map[lower_key] = key

    # Reference structure (no transformation needed)
    aligned_structures[reference_id] = structures[reference_id].copy()
    aligned_structures[reference_id]['is_reference'] = True
    if structures[reference_id].get('retinal'):
        ref_retinal = structures[reference_id]['retinal']
        aligned_structures[reference_id]['retinal'] = {
            'coords': ref_retinal['coords'].copy(),
            'atom_names': ref_retinal['atom_names'].copy(),
            'res_name3l': ref_retinal['res_name3l'].copy(),
            'chain_ids': ref_retinal['chain_ids'].copy()
        }

    # Transform all other structures to align with the reference
    for struct_id in structures:
        if struct_id == reference_id:
            continue

        # Look for alignment path: we need (struct_id, reference_id) which aligns struct -> reference
        # alignment_paths[(A, B)] gives R, t such that A_coords @ R.T + t aligns A to B
        correct_key = (struct_id, reference_id)
        wrong_key = (reference_id, struct_id)
        correct_key_lower = (struct_id.lower(), reference_id.lower())
        wrong_key_lower = (reference_id.lower(), struct_id.lower())

        R = None
        t = None

        if correct_key in alignment_paths:
            # This is (struct_id, reference_id) - aligns struct -> reference (CORRECT)
            alignment_info = alignment_paths[correct_key]
            R = np.array(alignment_info['rotation'])
            t = np.array(alignment_info['translation'])
            print(f"Found alignment: {struct_id} -> {reference_id}")

        elif correct_key_lower in align_lower_map:
            orig_key = align_lower_map[correct_key_lower]
            alignment_info = alignment_paths[orig_key]
            R = np.array(alignment_info['rotation'])
            t = np.array(alignment_info['translation'])
            print(f"Found alignment (case-insensitive): {struct_id} -> {reference_id}")

        elif wrong_key in alignment_paths:
            # This is (reference_id, struct_id) - aligns reference -> struct
            # We need the inverse: stored as (struct_id, reference_id) which should exist
            # But if we're here, it means only (ref, struct) exists, so we use it inverted
            alignment_info = alignment_paths[wrong_key]
            R_orig = np.array(alignment_info['rotation'])
            t_orig = np.array(alignment_info['translation'])
            # Invert: R_inv = R.T, t_inv = -R.T @ t
            R = R_orig.T
            t = -R_orig.T @ t_orig
            print(f"Found inverse alignment: {reference_id} -> {struct_id} (inverted)")

        elif wrong_key_lower in align_lower_map:
            orig_key = align_lower_map[wrong_key_lower]
            alignment_info = alignment_paths[orig_key]
            R_orig = np.array(alignment_info['rotation'])
            t_orig = np.array(alignment_info['translation'])
            R = R_orig.T
            t = -R_orig.T @ t_orig
            print(f"Found inverse alignment (case-insensitive): {reference_id} -> {struct_id} (inverted)")

        else:
            print(f"No alignment path found for {struct_id} - EXCLUDING from visualization")
            # Skip structures without alignment paths to avoid displaying them at wrong positions
            continue

        # Apply transformation to align structure to reference frame
        # The correct formula is: (coords - struct_center) @ R.T + ref_center
        # We need to compute centroids from matched residues
        if R is not None and alignment_info is not None:
            struct_coords = structures[struct_id]['coords']
            struct_residues = structures[struct_id]['residues']
            ref_coords = structures[reference_id]['coords']
            ref_residues = structures[reference_id]['residues']

            # Get residue mapping from alignment info
            residue_mapping = alignment_info.get('residue_mapping', [])

            # Build residue -> index maps
            struct_res_to_idx = {int(res): i for i, res in enumerate(struct_residues)}
            ref_res_to_idx = {int(res): i for i, res in enumerate(ref_residues)}

            # Find matched indices
            struct_matched_idx = []
            ref_matched_idx = []
            for mapping in residue_mapping:
                s_res, r_res = mapping[0], mapping[1]
                if s_res in struct_res_to_idx and r_res in ref_res_to_idx:
                    struct_matched_idx.append(struct_res_to_idx[s_res])
                    ref_matched_idx.append(ref_res_to_idx[r_res])

            if len(struct_matched_idx) > 0:
                # Compute centroids from matched residues
                struct_center = np.mean(struct_coords[struct_matched_idx], axis=0)
                ref_center = np.mean(ref_coords[ref_matched_idx], axis=0)

                # Apply correct transformation: (coords - struct_center) @ R.T + ref_center
                coords_transformed = (struct_coords - struct_center) @ R.T + ref_center
            else:
                # Fallback: use overall centroids
                struct_center = np.mean(struct_coords, axis=0)
                ref_center = np.mean(ref_coords, axis=0)
                coords_transformed = (struct_coords - struct_center) @ R.T + ref_center

            # Copy structure data and update coordinates
            aligned_structures[struct_id] = structures[struct_id].copy()
            aligned_structures[struct_id]['coords'] = coords_transformed
            aligned_structures[struct_id]['is_reference'] = False

            # Apply same transformation to retinal coordinates if available
            if structures[struct_id].get('retinal'):
                retinal_data = structures[struct_id]['retinal']
                retinal_coords = (retinal_data['coords'] - struct_center) @ R.T + ref_center
                aligned_structures[struct_id]['retinal'] = {
                    'coords': retinal_coords,
                    'atom_names': retinal_data['atom_names'].copy(),
                    'res_name3l': retinal_data['res_name3l'].copy(),
                    'chain_ids': retinal_data['chain_ids'].copy()
                }
            else:
                aligned_structures[struct_id]['retinal'] = None

            print(f"Successfully transformed {struct_id} to reference frame")
        else:
            print(f"Failed to get transformation matrices for {struct_id}")
            aligned_structures[struct_id] = structures[struct_id].copy()
            aligned_structures[struct_id]['is_reference'] = False
            if structures[struct_id].get('retinal'):
                retinal_data = structures[struct_id]['retinal']
                aligned_structures[struct_id]['retinal'] = {
                    'coords': retinal_data['coords'].copy(),
                    'atom_names': retinal_data['atom_names'].copy(),
                    'res_name3l': retinal_data['res_name3l'].copy(),
                    'chain_ids': retinal_data['chain_ids'].copy()
                }
            else:
                aligned_structures[struct_id]['retinal'] = None

    print(f"Aligned {len(aligned_structures)-1} structures to reference {reference_id}")
    return aligned_structures


def calculate_membrane_orientation_from_helix_topology(reference_coords, reference_residues, helix_assignments):
    """
    Calculate membrane orientation using PCA component 2 as Z-axis.

    Uses PCA to find principal components, then uses the 2nd component as the transmembrane direction.

    Args:
        reference_coords: Coordinates of reference structure (N, 3)
        reference_residues: Residue sequence numbers
        helix_assignments: Dictionary with helix boundary definitions (not used in this version)

    Returns:
        rotation_matrix: 3x3 rotation matrix to orient PCA component 2 along Z-axis
    """
    print("Calculating membrane orientation using PCA component 2 as Z-axis...")

    # Center the coordinates
    centered_coords = reference_coords - np.mean(reference_coords, axis=0)

    # Perform PCA
    pca = PCA(n_components=3)
    pca.fit(centered_coords)

    # Get principal components (eigenvectors)
    principal_axes = pca.components_

    print(f"PCA explained variance ratio: {pca.explained_variance_ratio_}")
    print(f"PCA component 1: [{principal_axes[0][0]:.3f}, {principal_axes[0][1]:.3f}, {principal_axes[0][2]:.3f}]")
    print(f"PCA component 2: [{principal_axes[1][0]:.3f}, {principal_axes[1][1]:.3f}, {principal_axes[1][2]:.3f}]")
    print(f"PCA component 3: [{principal_axes[2][0]:.3f}, {principal_axes[2][1]:.3f}, {principal_axes[2][2]:.3f}]")

    # Use the 2nd PCA component as the transmembrane direction
    transmembrane_vector = principal_axes[1]  # Second component

    # Check orientation using helix topology if available
    if helix_assignments and len(helix_assignments) > 0:
        print("Checking orientation using helix topology...")

        # Define helix topology (start→end sides)
        helix_topology = {
            1: ('extracellular', 'intracellular'),
            2: ('intracellular', 'extracellular'),
            3: ('extracellular', 'intracellular'),
            4: ('intracellular', 'extracellular'),
            5: ('extracellular', 'intracellular'),
            6: ('intracellular', 'extracellular'),
            7: ('extracellular', 'intracellular')
        }

        extracellular_coords = []
        intracellular_coords = []

        # Collect coordinates based on helix topology
        for helix_num in range(1, 8):
            if helix_num in helix_assignments:
                boundaries = helix_assignments[helix_num]
                if isinstance(boundaries, list) and len(boundaries) >= 2:
                    start_pos, end_pos = boundaries[0], boundaries[1]

                    # Find residues in this helix
                    helix_mask = (reference_residues >= start_pos) & (reference_residues <= end_pos)
                    helix_coords = reference_coords[helix_mask]
                    helix_residues = reference_residues[helix_mask]

                    if len(helix_coords) > 0:
                        # Sort by residue number to get proper start→end order
                        sorted_indices = np.argsort(helix_residues)
                        sorted_coords = helix_coords[sorted_indices]

                        # Get topology for this helix
                        start_side, end_side = helix_topology[helix_num]

                        # Take first and last few residues of the helix
                        n_terminal = max(1, len(sorted_coords) // 4)  # First 25% of helix
                        c_terminal = max(1, len(sorted_coords) // 4)  # Last 25% of helix

                        start_coords = sorted_coords[:n_terminal]  # Start of helix
                        end_coords = sorted_coords[-c_terminal:]   # End of helix

                        # Add to appropriate side based on topology
                        if start_side == 'extracellular':
                            extracellular_coords.extend(start_coords)
                        else:
                            intracellular_coords.extend(start_coords)

                        if end_side == 'extracellular':
                            extracellular_coords.extend(end_coords)
                        else:
                            intracellular_coords.extend(end_coords)

        if len(extracellular_coords) > 0 and len(intracellular_coords) > 0:
            # Calculate average positions for each side
            extracellular_coords = np.array(extracellular_coords)
            intracellular_coords = np.array(intracellular_coords)

            extracellular_center = np.mean(extracellular_coords, axis=0)
            intracellular_center = np.mean(intracellular_coords, axis=0)

            # Calculate transmembrane direction from topology
            topology_transmembrane = intracellular_center - extracellular_center
            topology_transmembrane = topology_transmembrane / np.linalg.norm(topology_transmembrane)

            print(f"Topology transmembrane vector: [{topology_transmembrane[0]:.3f}, {topology_transmembrane[1]:.3f}, {topology_transmembrane[2]:.3f}]")

            # Check if PCA component 2 aligns with topology
            dot_product = np.dot(transmembrane_vector, topology_transmembrane)
            print(f"PCA component 2 vs topology alignment: {dot_product:.3f}")

            # If they point in opposite directions, flip PCA component 2
            if dot_product < 0:
                transmembrane_vector = -transmembrane_vector
                print("Flipped PCA component 2 to align with topology")

    print(f"Using transmembrane vector: [{transmembrane_vector[0]:.3f}, {transmembrane_vector[1]:.3f}, {transmembrane_vector[2]:.3f}]")

    # Target direction: transmembrane direction should point to -Z (EC at +Z, IC at -Z)
    target_direction = np.array([0, 0, -1])

    # Calculate rotation matrix using Rodrigues' formula
    v1 = transmembrane_vector / np.linalg.norm(transmembrane_vector)
    v2 = target_direction / np.linalg.norm(target_direction)

    dot_product = np.dot(v1, v2)
    if abs(dot_product) > 0.999:  # Already aligned
        if dot_product > 0:
            rotation_matrix = np.eye(3)
        else:
            # Need 180-degree rotation
            if abs(v1[0]) < 0.9:
                axis = np.cross(v1, [1, 0, 0])
            else:
                axis = np.cross(v1, [0, 1, 0])
            axis = axis / np.linalg.norm(axis)
            rotation_matrix = 2 * np.outer(axis, axis) - np.eye(3)
    else:
        # General rotation
        cross_product = np.cross(v1, v2)
        sin_angle = np.linalg.norm(cross_product)
        cos_angle = dot_product

        if sin_angle > 1e-10:
            rotation_axis = cross_product / sin_angle
            K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                         [rotation_axis[2], 0, -rotation_axis[0]],
                         [-rotation_axis[1], rotation_axis[0], 0]])
            rotation_matrix = np.eye(3) + sin_angle * K + (1 - cos_angle) * np.dot(K, K)
        else:
            rotation_matrix = np.eye(3)

    # Verify the rotation
    rotated_transmembrane = np.dot(rotation_matrix, transmembrane_vector)
    print(f"Final transmembrane vector (should point to -Z): [{rotated_transmembrane[0]:.3f}, {rotated_transmembrane[1]:.3f}, {rotated_transmembrane[2]:.3f}]")

    return rotation_matrix


def calculate_membrane_orientation_fallback(reference_coords, reference_residues):
    """Fallback PCA-based membrane orientation calculation"""
    print("Using PCA fallback method...")

    # Center the coordinates
    centered_coords = reference_coords - np.mean(reference_coords, axis=0)

    # Perform PCA
    pca = PCA(n_components=3)
    pca.fit(centered_coords)

    # The first component is the direction of maximum variance (membrane normal)
    membrane_normal = pca.components_[0]

    # Use N/C terminus to determine orientation
    min_res = np.min(reference_residues)
    max_res = np.max(reference_residues)

    n_term_mask = reference_residues <= min_res + 5
    c_term_mask = reference_residues >= max_res - 5

    if np.any(n_term_mask) and np.any(c_term_mask):
        n_term_center = np.mean(reference_coords[n_term_mask], axis=0)
        c_term_center = np.mean(reference_coords[c_term_mask], axis=0)
        nc_vector = n_term_center - c_term_center
        nc_vector = nc_vector / np.linalg.norm(nc_vector)

        if np.dot(membrane_normal, nc_vector) < 0:
            membrane_normal = -membrane_normal

    # Target direction: N-terminus to +Z
    target_direction = np.array([0, 0, 1])

    # Calculate rotation matrix
    v1 = membrane_normal / np.linalg.norm(membrane_normal)
    v2 = target_direction / np.linalg.norm(target_direction)

    dot_product = np.dot(v1, v2)
    if abs(dot_product) > 0.999:
        rotation_matrix = np.eye(3) if dot_product > 0 else -np.eye(3)
    else:
        cross_product = np.cross(v1, v2)
        sin_angle = np.linalg.norm(cross_product)
        cos_angle = dot_product

        rotation_axis = cross_product / sin_angle
        K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                     [rotation_axis[2], 0, -rotation_axis[0]],
                     [-rotation_axis[1], rotation_axis[0], 0]])
        rotation_matrix = np.eye(3) + sin_angle * K + (1 - cos_angle) * np.dot(K, K)

    return rotation_matrix


def apply_membrane_orientation(aligned_structures, reference_id='6xl3'):
    """
    Apply membrane orientation to all aligned structures based on reference structure PCA.

    Args:
        aligned_structures: Dictionary of aligned structures
        reference_id: ID of reference structure for PCA calculation

    Returns:
        oriented_structures: Dictionary of structures oriented for membrane visualization
    """
    if reference_id not in aligned_structures:
        print(f"Warning: Reference {reference_id} not found for orientation calculation")
        return aligned_structures

    print(f"\n=== Calculating Membrane Orientation from {reference_id} ===")

    # Get reference structure coordinates, residues, and helix assignments
    ref_coords = aligned_structures[reference_id]['coords']
    ref_residues = aligned_structures[reference_id]['residues']

    # Get helix assignments from the processed structures
    ref_struct_data = None
    for struct_id, processed_data in aligned_structures.items():
        if struct_id == reference_id:
            ref_struct_data = processed_data
            break

    # Try to get helix assignments from the structure data
    helix_assignments = {}
    if ref_struct_data and 'helix_definitions' in ref_struct_data:
        helix_assignments = ref_struct_data['helix_definitions']
        print(f"Found helix definitions: {helix_assignments}")
    else:
        print("Warning: No helix definitions found in reference structure")

    # Calculate rotation matrix for membrane orientation using helix topology
    membrane_rotation = calculate_membrane_orientation_from_helix_topology(
        ref_coords, ref_residues, helix_assignments
    )

    # Apply rotation to all structures
    oriented_structures = {}
    for struct_id, struct_data in aligned_structures.items():
        # Copy structure data
        oriented_structures[struct_id] = struct_data.copy()

        # Apply membrane orientation rotation
        original_coords = struct_data['coords']

        # Center coordinates for rotation
        coord_center = np.mean(original_coords, axis=0)
        centered_coords = original_coords - coord_center

        # Apply rotation (R operates on column vectors, so we need coords @ R.T)
        rotated_coords = np.dot(centered_coords, membrane_rotation.T)

        # Recenter coordinates
        oriented_coords = rotated_coords + coord_center

        # Debug: Check Z-axis spread after orientation
        if struct_id == reference_id:
            z_min, z_max = np.min(oriented_coords[:, 2]), np.max(oriented_coords[:, 2])
            z_range = z_max - z_min
            print(f"Reference structure Z-range after orientation: {z_range:.2f} (min: {z_min:.2f}, max: {z_max:.2f})")

        # Update coordinates
        oriented_structures[struct_id]['coords'] = oriented_coords

        # Rotate retinal coordinates in tandem if available
        retinal_data = struct_data.get('retinal')
        if retinal_data and retinal_data.get('coords') is not None:
            retinal_coords = retinal_data['coords']
            retinal_centered = retinal_coords - coord_center
            rotated_retinal = np.dot(retinal_centered, membrane_rotation.T) + coord_center
            oriented_structures[struct_id]['retinal'] = {
                'coords': rotated_retinal,
                'atom_names': retinal_data['atom_names'].copy(),
                'res_name3l': retinal_data['res_name3l'].copy(),
                'chain_ids': retinal_data['chain_ids'].copy()
            }
        elif retinal_data is not None:
            oriented_structures[struct_id]['retinal'] = None

        print(f"Applied membrane orientation to {struct_id}")

    print(f"Applied membrane orientation to {len(oriented_structures)} structures")
    return oriented_structures


def calculate_residue_distribution(aligned_structures, grn_df, target_grn):
    """Calculate residue distribution for a specific GRN position"""
    residue_counts = {}
    total_count = 0

    # Get the GRN row from the table
    grn_row = grn_df[target_grn] if target_grn in grn_df.columns else None
    if grn_row is None:
        return {}, 0

    # Count residues at this GRN position across all structures
    for struct_id, cell_value in grn_row.items():
        if pd.notna(cell_value) and cell_value != '-':
            cell_str = str(cell_value).strip()
            if cell_str and cell_str != '-':
                # Extract amino acid (first character)
                amino_acid = cell_str[0]
                residue_counts[amino_acid] = residue_counts.get(amino_acid, 0) + 1
                total_count += 1

    return residue_counts, total_count


def create_residue_distribution_table(residue_counts, total_count, target_grn, highlight_amino_acid=None):
    """Create a table showing residue distribution"""
    if total_count == 0:
        return go.Table(
            header=dict(values=[f"GRN {target_grn}", "Count", "%"],
                       fill_color='lightgray'),
            cells=dict(values=[["No data"], ["0"], ["0"]],
                      fill_color='white')
        )

    # Sort residues by frequency
    sorted_residues = sorted(residue_counts.items(), key=lambda x: x[1], reverse=True)

    residues = [item[0] for item in sorted_residues]
    counts = [item[1] for item in sorted_residues]
    percentages = [f"{(count/total_count)*100:.1f}%" for count in counts]

    # Color code by amino acid properties
    aa_colors = []
    for aa in residues:
        if aa in 'FWYH':  # Aromatic
            aa_colors.append('#FFB6C1')  # Light pink
        elif aa in 'AILMV':  # Hydrophobic
            aa_colors.append('#98FB98')  # Pale green
        elif aa in 'RKDE':  # Charged
            aa_colors.append('#87CEEB')  # Sky blue
        elif aa in 'STNQ':  # Polar
            aa_colors.append('#DDA0DD')  # Plum
        elif aa in 'GP':  # Special
            aa_colors.append('#F0E68C')  # Khaki
        else:
            aa_colors.append('#FFFFFF')  # White

    # Highlight selected amino acid row if requested
    highlight_colors = ['white'] * len(counts)
    if highlight_amino_acid and highlight_amino_acid in residues:
        highlight_index = residues.index(highlight_amino_acid)
        highlight_colors[highlight_index] = '#FFF2CC'  # Light yellow highlight
    elif highlight_amino_acid and highlight_amino_acid != 'ALL' and highlight_amino_acid not in residues:
        # Provide subtle cue when amino acid absent for this GRN
        highlight_colors = ['#F8F8F8'] * len(counts)

    return go.Table(
        header=dict(
            values=[f"GRN {target_grn}", "Count", "%"],
            fill_color='lightgray',
            font=dict(size=12, color='black'),
            height=25
        ),
        cells=dict(
            values=[residues, counts, percentages],
            fill_color=[aa_colors, highlight_colors, highlight_colors],
            font=dict(size=11),
            height=20
        )
    )


def create_structure_presence_table(structure_entries, target_grn, amino_acid):
    """Create a table listing structures that carry the selected residue at the GRN position."""

    header_title = f"GRN {target_grn} — {amino_acid if amino_acid and amino_acid != 'ALL' else 'Select residue'}"

    if not structure_entries:
        return go.Table(
            header=dict(
                values=[header_title],
                fill_color='lightgray',
                font=dict(size=12, color='black'),
                height=25
            ),
            cells=dict(
                values=[["No matching structures"]],
                fill_color='white',
                font=dict(size=11),
                height=20
            )
        )

    structure_entries = sorted(structure_entries, key=lambda item: (item.get('dataset', ''), item.get('structure')))
    structures = [entry.get('structure', 'Unknown') for entry in structure_entries]
    datasets = [entry.get('dataset', 'unknown') for entry in structure_entries]
    functions = [entry.get('function', 'Unknown') for entry in structure_entries]

    return go.Table(
        header=dict(
            values=[header_title, "Dataset", "Molecular Function"],
            fill_color='lightgray',
            font=dict(size=12, color='black'),
            height=25
        ),
        cells=dict(
            values=[structures, datasets, functions],
            fill_color=[['white'] * len(structures)] * 3,
            font=dict(size=10),
            height=18
        )
    )


def create_interactive_opsin_visualization(
    aligned_structures,
    grn_df,
    property_data=None,
    title="Interactive GRN-based Opsin Structure Alignment",
    width=1600,
    height=1000,
    max_structures=125,
    membrane_opacity=0.05,
    show_membrane=True
):
    """Legacy interactive visualization used for interactive_grn_alignment.html."""

    # Determine which structures will be visualized (limited by max_structures)
    visualized_struct_ids = list(aligned_structures.keys())[:max_structures]

    # Find overall center of CA coordinates for ONLY the visualized structures
    all_coords = []
    for struct_id in visualized_struct_ids:
        data = aligned_structures[struct_id]
        if 'coords' in data:
            all_coords.extend(data['coords'])

    if all_coords:
        all_coords = np.array(all_coords)
        overall_center = np.mean(all_coords, axis=0)
        print(f"Centering visualization on CA center of {len(visualized_struct_ids)} structures: {overall_center}")
    else:
        overall_center = np.array([0.0, 0.0, 0.0])

    # Apply centering offset to all structures
    for struct_id, data in aligned_structures.items():
        if 'coords' in data:
            data['coords'] = data['coords'] - overall_center
        if 'retinal' in data and data['retinal'] is not None and data['retinal'].get('coords') is not None:
            data['retinal']['coords'] = data['retinal']['coords'] - overall_center
        if 'dataframe' in data and data['dataframe'] is not None:
            df = data['dataframe']
            df[['x', 'y', 'z']] = df[['x', 'y', 'z']].values - overall_center

    existing_grns = set()
    for struct_id in visualized_struct_ids:
        data = aligned_structures[struct_id]
        df = data.get('dataframe')
        if df is not None:
            structure_grns = df['grn'].dropna().unique()
            existing_grns.update(structure_grns)
        else:
            grn_positions = data.get('grn_positions', data.get('grn', []))
            existing_grns.update([grn for grn in grn_positions if pd.notna(grn)])

    # Filter to helix GRNs (1-7) from the table that exist in at least one structure
    def is_helix_grn(grn):
        """Check if GRN is a helix position (1.XX to 7.XX)"""
        if not grn or pd.isna(grn):
            return False
        grn_str = str(grn)
        if '.' not in grn_str:
            return False
        prefix = grn_str.split('.')[0]
        if not prefix.isdigit():
            return False
        helix_num = int(prefix)
        return 1 <= helix_num <= 7

    # Get all helix GRNs from the table
    helix_grns_in_table = [grn for grn in grn_df.columns if is_helix_grn(grn)]
    # Filter to those that exist in at least one structure
    all_grn_positions = sorted([grn for grn in helix_grns_in_table if grn in existing_grns],
                                key=lambda x: (int(str(x).split('.')[0]), int(str(x).split('.')[1])))
    print(f"GRN positions in table: {len(grn_df.columns)}")
    print(f"Helix GRNs (1-7) in table: {len(helix_grns_in_table)}")
    print(f"GRN positions in structures: {len(existing_grns)}")
    print(f"Helix GRN positions for slider: {len(all_grn_positions)}")

    helix_colors = {
        1: HELIX_NUMBER_COLORS[1],
        2: HELIX_NUMBER_COLORS[2],
        3: HELIX_NUMBER_COLORS[3],
        4: HELIX_NUMBER_COLORS[4],
        5: HELIX_NUMBER_COLORS[5],
        6: HELIX_NUMBER_COLORS[6],
        7: HELIX_NUMBER_COLORS[7],
        0: '#D3D3D3'
    }

    from src.opsin_color_scheme import get_categorical_colors

    property_colors = {}
    if property_data:
        molecular_functions = set()
        for struct_id, props in property_data.items():
            if 'molecular_function' in props:
                molecular_functions.add(props['molecular_function'])

        property_colors = get_categorical_colors(list(molecular_functions), property_type='property1')
        print(f"Created property colors for {len(property_colors)} molecular functions: {list(property_colors.keys())}")

    def get_structure_color(struct_id, helix_num, coloring_mode='helix'):
        if coloring_mode == 'property' and property_data and struct_id in property_data:
            mol_func = property_data[struct_id].get('molecular_function', 'Unknown')
            return property_colors.get(mol_func, '#D3D3D3')
        return helix_colors.get(helix_num, '#D3D3D3')

    fig = make_subplots(
        rows=1, cols=2,
        column_widths=[0.75, 0.25],
        specs=[[{"type": "scatter3d"}, {"type": "table"}]],
        subplot_titles=("Opsin Structure Alignment", "Residue Distribution")
    )

    helix_traces = {}
    property_traces = {}

    structure_count = 0
    print("\n=== Creating helix-colored structure traces ===")
    for struct_id, data in aligned_structures.items():
        coords = data['coords']
        grn_positions = data.get('grn', data['grn_positions'])
        helix_numbers = data['helix_numbers']

        if structure_count < 3:
            z_min, z_max = np.min(coords[:, 2]), np.max(coords[:, 2])
            print(f"  {struct_id}: Z-range = {z_max - z_min:.2f} (min: {z_min:.2f}, max: {z_max:.2f})")

        df = data.get('dataframe')
        if df is not None:
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = coords
            grn_positions = df_aligned['grn'].values
            helix_numbers = df_aligned['helix_num'].values

        if structure_count > max_structures:
            break

        # Helper to format residue label (e.g., "F74")
        def format_res_label(res1l, seq_id):
            res = str(res1l).strip() if pd.notna(res1l) else 'X'
            seq = ''
            if pd.notna(seq_id):
                try:
                    seq = str(int(seq_id))
                except (ValueError, TypeError):
                    seq = str(seq_id).strip()
            return f"{res}{seq}"

        if df is not None:
            for helix_num, group in df_aligned.groupby('helix_num'):
                if len(group) == 0:
                    continue

                helix_coords = group[['x', 'y', 'z']].values
                helix_grn = group['grn'].values
                # Get residue info for hover text
                helix_res1l = group['res_name1l'].values if 'res_name1l' in group.columns else ['X'] * len(group)
                helix_seq_id = group['auth_seq_id'].values if 'auth_seq_id' in group.columns else [None] * len(group)
                hover_text = [
                    f"{str(grn) if pd.notna(grn) else 'No GRN'}<br>Residue: {format_res_label(r, s)}"
                    for grn, r, s in zip(helix_grn, helix_res1l, helix_seq_id)
                ]

                if helix_num == 0:
                    trace_name = "Loops/Non-helix"
                    legend_group = "Helix_0"
                else:
                    trace_name = f"Helix {int(helix_num)}"
                    legend_group = f"Helix_{int(helix_num)}"

                show_legend = legend_group not in helix_traces
                if show_legend:
                    helix_traces[legend_group] = True

                fig.add_trace(go.Scatter3d(
                    x=helix_coords[:, 0],
                    y=helix_coords[:, 1],
                    z=helix_coords[:, 2],
                    mode='markers+lines',
                    marker=dict(size=2, color=helix_colors.get(helix_num, '#D3D3D3'), opacity=0.1, line=dict(width=0)),
                    line=dict(color=helix_colors.get(helix_num, '#D3D3D3'), width=0.5),
                    name=trace_name,
                    legendgroup=legend_group,
                    showlegend=show_legend,
                    hovertemplate=f'<b>{struct_id}</b><br>' +
                                 f'Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>' +
                                 'GRN: %{text}<br>' +
                                 'X: %{x:.2f}<br>' +
                                 'Y: %{y:.2f}<br>' +
                                 'Z: %{z:.2f}<extra></extra>',
                    text=hover_text,
                    visible=True
                ), row=1, col=1)
        else:
            # Get residue info from data dict
            res_name1l = data.get('res_name1l', ['X'] * len(coords))
            residue_ids = data.get('residues', [None] * len(coords))

            for helix_num in np.unique(helix_numbers):
                helix_mask = helix_numbers == helix_num
                if not np.any(helix_mask):
                    continue

                helix_coords = coords[helix_mask]
                helix_grn = grn_positions[helix_mask]
                helix_res1l = np.array(res_name1l)[helix_mask] if res_name1l is not None else ['X'] * np.sum(helix_mask)
                helix_seq_id = np.array(residue_ids)[helix_mask] if residue_ids is not None else [None] * np.sum(helix_mask)
                hover_text = [
                    f"{str(grn) if pd.notna(grn) else 'No GRN'}<br>Residue: {format_res_label(r, s)}"
                    for grn, r, s in zip(helix_grn, helix_res1l, helix_seq_id)
                ]

                if helix_num == 0:
                    trace_name = "Loops/Non-helix"
                    legend_group = "Helix_0"
                else:
                    trace_name = f"Helix {int(helix_num)}"
                    legend_group = f"Helix_{int(helix_num)}"

                show_legend = legend_group not in helix_traces
                if show_legend:
                    helix_traces[legend_group] = True

                fig.add_trace(go.Scatter3d(
                    x=helix_coords[:, 0],
                    y=helix_coords[:, 1],
                    z=helix_coords[:, 2],
                    mode='markers+lines',
                    marker=dict(size=2, color=helix_colors.get(helix_num, '#D3D3D3'), opacity=0.1, line=dict(width=0)),
                    line=dict(color=helix_colors.get(helix_num, '#D3D3D3'), width=0.5),
                    name=trace_name,
                    legendgroup=legend_group,
                    showlegend=show_legend,
                    hovertemplate=f'<b>{struct_id}</b><br>' +
                                 f'Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>' +
                                 'GRN: %{text}<br>' +
                                 'X: %{x:.2f}<br>' +
                                 'Y: %{y:.2f}<br>' +
                                 'Z: %{z:.2f}<extra></extra>',
                    text=hover_text,
                    visible=True
                ), row=1, col=1)

        structure_count += 1

    helix_base_trace_count = len(fig.data)
    print(f"Added {helix_base_trace_count} helix-colored base traces")

    if property_data:
        print("\n=== Creating property-colored structure traces ===")
        structure_count = 0
        for struct_id, data in aligned_structures.items():
            coords = data['coords']
            grn_positions = data.get('grn', data['grn_positions'])
            helix_numbers = data['helix_numbers']

            if structure_count > max_structures:
                break

            mol_func = 'Unknown'
            if struct_id in property_data:
                mol_func = property_data[struct_id].get('molecular_function', 'Unknown')

            # Helper to format residue label (e.g., "F74")
            def format_res_label(res1l, seq_id):
                res = str(res1l).strip() if pd.notna(res1l) else 'X'
                seq = ''
                if pd.notna(seq_id):
                    try:
                        seq = str(int(seq_id))
                    except (ValueError, TypeError):
                        seq = str(seq_id).strip()
                return f"{res}{seq}"

            df = data.get('dataframe')
            if df is not None:
                df_aligned = df.copy()
                df_aligned[['x', 'y', 'z']] = coords
                grn_positions = df_aligned['grn'].values
                helix_numbers = df_aligned['helix_num'].values

                for helix_num, group in df_aligned.groupby('helix_num'):
                    if len(group) == 0:
                        continue

                    helix_coords = group[['x', 'y', 'z']].values
                    helix_grn = group['grn'].values
                    # Get residue info for hover text
                    helix_res1l = group['res_name1l'].values if 'res_name1l' in group.columns else ['X'] * len(group)
                    helix_seq_id = group['auth_seq_id'].values if 'auth_seq_id' in group.columns else [None] * len(group)
                    hover_text = [
                        f"{str(grn) if pd.notna(grn) else 'No GRN'}<br>Residue: {format_res_label(r, s)}"
                        for grn, r, s in zip(helix_grn, helix_res1l, helix_seq_id)
                    ]

                    trace_name = mol_func
                    legend_group = f"Property_{mol_func}"

                    show_legend = legend_group not in property_traces
                    if show_legend:
                        property_traces[legend_group] = True

                    property_color = get_structure_color(struct_id, helix_num, 'property')
                    fig.add_trace(go.Scatter3d(
                        x=helix_coords[:, 0],
                        y=helix_coords[:, 1],
                        z=helix_coords[:, 2],
                        mode='markers+lines',
                        marker=dict(size=2, color=property_color, opacity=0.1, line=dict(width=0)),
                        line=dict(color=property_color, width=0.5),
                        name=trace_name,
                        legendgroup=legend_group,
                        showlegend=show_legend,
                        hovertemplate=f'<b>{struct_id}</b><br>' +
                                     f'Function: {mol_func}<br>' +
                                     f'Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>' +
                                     'GRN: %{text}<br>' +
                                     'X: %{x:.2f}<br>' +
                                     'Y: %{y:.2f}<br>' +
                                     'Z: %{z:.2f}<extra></extra>',
                        text=hover_text,
                        visible=False,
                        meta={'coloring_mode': 'property'}
                    ), row=1, col=1)

            structure_count += 1

        property_base_trace_count = len(fig.data) - helix_base_trace_count
        print(f"Added {property_base_trace_count} property-colored base traces")

    base_trace_count = len(fig.data)
    print(f"Total base traces: {base_trace_count}")

    highlight_data_helix = {}
    highlight_data_property = {}

    structure_count = 0
    for struct_id, data in aligned_structures.items():
        if structure_count > max_structures:
            break

        coords = data['coords']
        df = data.get('dataframe')

        # Helper to format residue label (e.g., "F74")
        def format_res_label(res1l, seq_id):
            res = str(res1l).strip() if pd.notna(res1l) else 'X'
            seq = ''
            if pd.notna(seq_id):
                try:
                    seq = str(int(seq_id))
                except (ValueError, TypeError):
                    seq = str(seq_id).strip()
            return f"{res}{seq}"

        if df is not None:
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = coords
            grn_rows = df_aligned[df_aligned['grn'].notna()]

            for _, row in grn_rows.iterrows():
                grn_pos = row['grn']
                helix_num = row['helix_num']
                coord = [row['x'], row['y'], row['z']]
                res_label = format_res_label(row.get('res_name1l'), row.get('auth_seq_id'))

                highlight_data_helix.setdefault(grn_pos, {'coords': [], 'colors': [], 'hover_text': []})
                color_val = get_structure_color(struct_id, helix_num, 'helix')
                highlight_data_helix[grn_pos]['coords'].append(coord)
                highlight_data_helix[grn_pos]['colors'].append(color_val)
                highlight_data_helix[grn_pos]['hover_text'].append(
                    f'<b>{struct_id}</b><br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}<br>Residue: {res_label}'
                )

                if property_data:
                    highlight_data_property.setdefault(grn_pos, {'coords': [], 'colors': [], 'hover_text': []})
                    property_color = get_structure_color(struct_id, helix_num, 'property')
                    mol_func = property_data.get(struct_id, {}).get('molecular_function', 'Unknown')
                    highlight_data_property[grn_pos]['coords'].append(coord)
                    highlight_data_property[grn_pos]['colors'].append(property_color)
                    highlight_data_property[grn_pos]['hover_text'].append(
                        f'<b>{struct_id}</b><br>Function: {mol_func}<br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}<br>Residue: {res_label}'
                    )
        else:
            grn_positions = data['grn_positions']
            helix_numbers = data['helix_numbers']
            res_name1l = data.get('res_name1l', ['X'] * len(coords))
            residue_ids = data.get('residues', [None] * len(coords))

            for i, (coord, grn_pos, helix_num) in enumerate(zip(coords, grn_positions, helix_numbers)):
                if pd.notna(grn_pos):
                    res1 = res_name1l[i] if res_name1l is not None else 'X'
                    seq_id = residue_ids[i] if residue_ids is not None else None
                    res_label = format_res_label(res1, seq_id)

                    highlight_data_helix.setdefault(grn_pos, {'coords': [], 'colors': [], 'hover_text': []})
                    color_val = get_structure_color(struct_id, helix_num, 'helix')
                    highlight_data_helix[grn_pos]['coords'].append(coord)
                    highlight_data_helix[grn_pos]['colors'].append(color_val)
                    highlight_data_helix[grn_pos]['hover_text'].append(
                        f'<b>{struct_id}</b><br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}<br>Residue: {res_label}'
                    )

                    if property_data:
                        highlight_data_property.setdefault(grn_pos, {'coords': [], 'colors': [], 'hover_text': []})
                        property_color = get_structure_color(struct_id, helix_num, 'property')
                        mol_func = property_data.get(struct_id, {}).get('molecular_function', 'Unknown')
                        highlight_data_property[grn_pos]['coords'].append(coord)
                        highlight_data_property[grn_pos]['colors'].append(property_color)
                        highlight_data_property[grn_pos]['hover_text'].append(
                            f'<b>{struct_id}</b><br>Function: {mol_func}<br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}<br>Residue: {res_label}'
                        )

        structure_count += 1

    print(f"Collected helix highlight data for {len(highlight_data_helix)} GRN positions")
    if property_data:
        print(f"Collected property highlight data for {len(highlight_data_property)} GRN positions")

    helix_highlight_trace_start = len(fig.data)
    first_grn = True
    for grn_pos in all_grn_positions:
        if grn_pos in highlight_data_helix:
            coords_array = np.array(highlight_data_helix[grn_pos]['coords'])
            colors = highlight_data_helix[grn_pos]['colors']
            hover_text = highlight_data_helix[grn_pos]['hover_text']

            fig.add_trace(go.Scatter3d(
                x=coords_array[:, 0],
                y=coords_array[:, 1],
                z=coords_array[:, 2],
                mode='markers',
                marker=dict(size=4, color=colors, opacity=1.0, line=dict(width=1, color='white')),
                name=f"GRN {grn_pos}",
                showlegend=False,
                hovertemplate='%{text}<extra></extra>',
                text=hover_text,
                visible=first_grn,
                meta={'coloring_mode': 'helix'}
            ), row=1, col=1)
            first_grn = False
        else:
            fig.add_trace(go.Scatter3d(
                x=[], y=[], z=[],
                mode='markers',
                marker=dict(size=4, opacity=1.0),
                name=f"GRN {grn_pos}",
                showlegend=False,
                visible=False,
                meta={'coloring_mode': 'helix'}
            ), row=1, col=1)

    if property_data:
        for grn_pos in all_grn_positions:
            if grn_pos in highlight_data_property:
                coords_array = np.array(highlight_data_property[grn_pos]['coords'])
                colors = highlight_data_property[grn_pos]['colors']
                hover_text = highlight_data_property[grn_pos]['hover_text']

                fig.add_trace(go.Scatter3d(
                    x=coords_array[:, 0],
                    y=coords_array[:, 1],
                    z=coords_array[:, 2],
                    mode='markers',
                    marker=dict(size=4, color=colors, opacity=1.0, line=dict(width=1, color='white')),
                    name=f"GRN {grn_pos}",
                    showlegend=False,
                    hovertemplate='%{text}<extra></extra>',
                    text=hover_text,
                    visible=False,
                    meta={'coloring_mode': 'property'}
                ), row=1, col=1)
            else:
                fig.add_trace(go.Scatter3d(
                    x=[], y=[], z=[],
                    mode='markers',
                    marker=dict(size=4, opacity=1.0),
                    name=f"GRN {grn_pos}",
                    showlegend=False,
                    visible=False,
                    meta={'coloring_mode': 'property'}
                ), row=1, col=1)

    print(f"Added {len(all_grn_positions)} highlight traces")

    table_trace_start = len(fig.data)
    first_table = True
    for grn_pos in all_grn_positions:
        counts, total = calculate_residue_distribution(aligned_structures, grn_df, grn_pos)
        table = create_residue_distribution_table(counts, total, grn_pos)
        fig.add_trace(table, row=1, col=2)
        fig.data[-1].visible = first_table
        first_table = False

    print(f"Added {len(all_grn_positions)} table traces")

    print("Adding membrane reference planes...")
    all_coords = []
    for struct_id, data in list(aligned_structures.items())[:10]:
        all_coords.extend(data['coords'])
    all_coords = np.array(all_coords)

    x_range = [np.min(all_coords[:, 0]) - 5, np.max(all_coords[:, 0]) + 5]
    y_range = [np.min(all_coords[:, 1]) - 5, np.max(all_coords[:, 1]) + 5]

    membrane_trace_count = 0
    if show_membrane:
        x_vol = np.linspace(x_range[0], x_range[1], 20)
        y_vol = np.linspace(y_range[0], y_range[1], 20)
        z_vol = np.linspace(-10, 10, 15)
        X_vol, Y_vol, Z_vol = np.meshgrid(x_vol, y_vol, z_vol, indexing='ij')

        membrane_values = np.ones_like(X_vol).flatten()
        fig.add_trace(go.Volume(
            x=X_vol.flatten(),
            y=Y_vol.flatten(),
            z=Z_vol.flatten(),
            value=membrane_values,
            isomin=0.3,
            isomax=0.7,
            opacity=membrane_opacity,
            surface_count=3,
            colorscale=[[0, 'lightgray'], [0.5, 'silver'], [1, 'lightgray']],
            showscale=False,
            name='Membrane',
            legendgroup='Membrane',
            showlegend=True,
            hovertemplate='Membrane<extra></extra>'
        ), row=1, col=1)
        membrane_trace_count = 1

    # Total number of traces in fig.data for visibility array sizing
    total_traces = len(fig.data)

    # Trace indices in fig.data:
    # 0 to helix_base_trace_count-1: helix base traces
    # helix_base_trace_count to base_trace_count-1: property base traces
    # base_trace_count to base_trace_count+len(all_grn_positions)*2-1: highlight traces
    # then table traces, then membrane

    # Calculate where highlight traces start
    highlight_start_idx = base_trace_count

    def create_visibility_array(grn_index, color_mode='helix', show_grn=True):
        # Build visibility array matching exact order in fig.data
        visibility = [False] * total_traces

        # Set base traces visibility
        if color_mode == 'helix':
            for i in range(helix_base_trace_count):
                visibility[i] = True
        else:
            for i in range(helix_base_trace_count, base_trace_count):
                visibility[i] = True

        # Set membrane trace visibility (last trace)
        if membrane_trace_count > 0:
            visibility[-1] = True

        # Set highlight trace visibility - only if show_grn is True
        # Helix highlights are at indices: base_trace_count to base_trace_count + len(all_grn_positions) - 1
        # Property highlights are at indices: base_trace_count + len(all_grn_positions) to base_trace_count + 2*len(all_grn_positions) - 1
        if show_grn:
            if color_mode == 'helix':
                helix_highlight_idx = highlight_start_idx + grn_index
                if helix_highlight_idx < total_traces:
                    visibility[helix_highlight_idx] = True
            else:
                property_highlight_idx = highlight_start_idx + len(all_grn_positions) + grn_index
                if property_highlight_idx < total_traces:
                    visibility[property_highlight_idx] = True

        # Set table trace visibility - only if show_grn is True
        # Table traces are at indices: base_trace_count + 2*len(all_grn_positions) to base_trace_count + 3*len(all_grn_positions) - 1
        if show_grn:
            table_start_idx = highlight_start_idx + 2 * len(all_grn_positions)
            table_idx = table_start_idx + grn_index
            if table_idx < total_traces:
                visibility[table_idx] = True

        return visibility

    helix_steps = []
    for i, grn_pos in enumerate(all_grn_positions):
        visibility = create_visibility_array(i, 'helix')
        helix_steps.append(dict(method="restyle", args=[{"visible": visibility}], label=str(grn_pos)))

    property_steps = []
    if property_data:
        for i, grn_pos in enumerate(all_grn_positions):
            visibility = create_visibility_array(i, 'property')
            property_steps.append(dict(method="restyle", args=[{"visible": visibility}], label=str(grn_pos)))

    sliders = [dict(active=0, currentvalue={"prefix": "GRN Position: "}, pad={"t": 50}, steps=helix_steps)]

    buttons = []
    if property_data:
        buttons.append(dict(
            label="Helix Coloring",
            method="update",
            args=[
                {"visible": create_visibility_array(0, 'helix')},
                {"sliders": [dict(active=0, currentvalue={"prefix": "GRN Position: "}, pad={"t": 50}, steps=helix_steps)]}
            ]
        ))
        buttons.append(dict(
            label="Property Coloring",
            method="update",
            args=[
                {"visible": create_visibility_array(0, 'property')},
                {"sliders": [dict(active=0, currentvalue={"prefix": "GRN Position: "}, pad={"t": 50}, steps=property_steps)]}
            ]
        ))

    updatemenus = []
    if buttons:
        updatemenus.append(dict(
            type="buttons",
            direction="left",
            buttons=buttons,
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.01,
            xanchor="left",
            y=0.02,
            yanchor="bottom"
        ))

    # Add GRN toggle buttons (Show/Hide GRN highlights)
    grn_toggle_buttons = [
        dict(
            label="Show GRN",
            method="update",
            args=[
                {"visible": create_visibility_array(0, 'helix', show_grn=True)},
                {"sliders": [dict(
                    active=0,
                    currentvalue={"prefix": "GRN Position: "},
                    pad={"t": 50},
                    steps=helix_steps,
                    visible=True
                )]}
            ]
        ),
        dict(
            label="Hide GRN",
            method="update",
            args=[
                {"visible": create_visibility_array(0, 'helix', show_grn=False)},
                {"sliders": [dict(
                    active=0,
                    currentvalue={"prefix": "GRN Position: "},
                    pad={"t": 50},
                    steps=helix_steps,
                    visible=False
                )]}
            ]
        )
    ]
    updatemenus.append(dict(
        type="buttons",
        direction="left",
        buttons=grn_toggle_buttons,
        pad={"r": 10, "t": 10},
        showactive=True,
        x=0.30,  # Position to the right of color mode buttons
        xanchor="left",
        y=0.02,
        yanchor="bottom"
    ))

    layout_args = dict(
        title=f'{title} (n={structure_count})',
        paper_bgcolor='white',
        plot_bgcolor='white',
        width=width,
        height=height,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='lightgray',
            borderwidth=1
        ),
        sliders=sliders,
        margin=dict(l=0, r=0, t=80, b=50)
    )

    if updatemenus:
        layout_args['updatemenus'] = updatemenus

    fig.update_layout(**layout_args)

    fig.update_scenes(
        xaxis=dict(visible=False, showbackground=False, showgrid=False, showline=False, showticklabels=False, title=""),
        yaxis=dict(visible=False, showbackground=False, showgrid=False, showline=False, showticklabels=False, title=""),
        zaxis=dict(visible=False, showbackground=False, showgrid=False, showline=False, showticklabels=False, title=""),
        camera=dict(eye=dict(x=1.8, y=1.8, z=0.8), center=dict(x=0, y=0, z=0), up=dict(x=0, y=0, z=1)),
        aspectmode='cube',  # Force 1:1:1 aspect ratio
        bgcolor='white'
    )

    return fig


def create_interactive_opsin_visualization_extended(
    aligned_structures,
    grn_df,
    property_data=None,
    title="Interactive GRN-based Opsin Alignment (Enhanced)",
    width=1600,
    height=1000,
    max_structures=125,
    membrane_opacity=0.05,
    show_membrane=True,
    include_retinal=False,
    retinal_reference_id='6xl3',
    hover_show_residue_name=False,
    enable_amino_acid_filter=False
):
    """
    Create enhanced interactive visualization with retinal overlay and amino-acid filtering.

    Args:
        aligned_structures (dict): Dictionary of aligned structure data with coordinates
        grn_df (pd.DataFrame): GRN position table for residue distribution
        property_data (dict, optional): Property data for structures (molecular function, etc.)
        title (str): Title for the visualization
        width (int): Figure width in pixels
        height (int): Figure height in pixels
        max_structures (int): Maximum number of structures to display
        membrane_opacity (float): Opacity of membrane volume (0-1)
        show_membrane (bool): Whether to show membrane volume
        include_retinal (bool): Add retinal atoms to the visualization (reference structure only)
        retinal_reference_id (str): Structure ID whose retinal coordinates should be displayed
        hover_show_residue_name (bool): Include amino-acid identity in residue hover text
        enable_amino_acid_filter (bool): Enable per-amino-acid filtering and structure listings

    Returns:
        plotly.graph_objects.Figure: Interactive 3D visualization
    """

    def format_residue_annotation(residue_short, residue_long, seq_id):
        """Create a concise residue label for hover text (e.g., F74)."""
        short_name = str(residue_short).strip() if pd.notna(residue_short) else ''
        if not short_name:
            short_name = 'X'

        seq_label = ''
        if pd.notna(seq_id):
            try:
                seq_label = str(int(seq_id))
            except (ValueError, TypeError):
                seq_str = str(seq_id).strip()
                if seq_str:
                    seq_label = seq_str

        return f"{short_name}{seq_label}"

    def select_retinal_structure(struct_map, preferred_id):
        """Return the structure ID that should supply retinal coordinates."""

        def has_retinal(sid):
            data = struct_map.get(sid)
            if not data:
                return False
            retinal = data.get('retinal')
            if not retinal:
                return False
            coords = retinal.get('coords')
            return coords is not None and len(coords) >= 2

        if preferred_id and has_retinal(preferred_id):
            return preferred_id

        for sid in struct_map:
            if has_retinal(sid):
                return sid

        return None

    def infer_retinal_bonds(coords, min_distance=1.05, max_distance=1.95, max_bonds_per_atom=4):
        """Infer retinal bond pairs using Euclidean distances."""
        coords = np.asarray(coords, dtype=float)
        if coords.ndim != 2 or coords.shape[0] < 2:
            return []

        n_atoms = coords.shape[0]
        diff = coords[:, None, :] - coords[None, :, :]
        distances = np.linalg.norm(diff, axis=-1)
        bonds = []
        bond_counts = np.zeros(n_atoms, dtype=int)
        used_pairs = set()

        for i in range(n_atoms):
            for j in range(i + 1, n_atoms):
                dist_val = distances[i, j]
                if dist_val < min_distance or dist_val > max_distance:
                    continue
                if bond_counts[i] >= max_bonds_per_atom or bond_counts[j] >= max_bonds_per_atom:
                    continue
                pair = (i, j)
                bonds.append(pair)
                used_pairs.add(pair)
                bond_counts[i] += 1
                bond_counts[j] += 1

        # Ensure no atom is left isolated; connect to nearest neighbor if necessary
        if n_atoms > 1:
            for atom_idx in range(n_atoms):
                if bond_counts[atom_idx] == 0:
                    sorted_neighbors = np.argsort(distances[atom_idx])
                    for neighbor_idx in sorted_neighbors:
                        if neighbor_idx == atom_idx:
                            continue
                        pair = tuple(sorted((atom_idx, neighbor_idx)))
                        if pair not in used_pairs:
                            bonds.append(pair)
                            used_pairs.add(pair)
                            bond_counts[atom_idx] += 1
                            bond_counts[neighbor_idx] += 1
                        break

        return bonds

    def build_retinal_traces(struct_id, retinal_data):
        """Create bond and atom traces for a single retinal ligand."""
        coords = retinal_data.get('coords')
        if coords is None or len(coords) < 2:
            return []

        coords = np.asarray(coords, dtype=float)
        atom_names = retinal_data.get('atom_names')
        if atom_names is None:
            atom_names = [f'Atom {i + 1}' for i in range(len(coords))]

        bonds = infer_retinal_bonds(coords)
        bond_x, bond_y, bond_z = [], [], []
        for i, j in bonds:
            bond_x.extend([coords[i, 0], coords[j, 0], None])
            bond_y.extend([coords[i, 1], coords[j, 1], None])
            bond_z.extend([coords[i, 2], coords[j, 2], None])

        traces = []
        if bond_x:
            traces.append(
                go.Scatter3d(
                    x=bond_x,
                    y=bond_y,
                    z=bond_z,
                    mode='lines',
                    line=dict(color='#404040', width=4),  # Dark grey
                    name=f'{struct_id} RET bonds',
                    legendgroup='RET',
                    showlegend=True,
                    hoverinfo='skip'
                )
            )

        traces.append(
            go.Scatter3d(
                x=coords[:, 0],
                y=coords[:, 1],
                z=coords[:, 2],
                mode='markers',
                marker=dict(size=6, color='#606060', line=dict(color='#404040', width=1.5)),  # Dark grey
                name=f'{struct_id} RET atoms',
                legendgroup='RET',
                showlegend=not bool(traces),
                text=[str(name) for name in atom_names],
                hovertemplate=(
                    '<b>%{text}</b><br>'
                    f'Structure: {struct_id}<br>'
                    'X: %{x:.2f}<br>'
                    'Y: %{y:.2f}<br>'
                    'Z: %{z:.2f}<extra></extra>'
                )
            )
        )

        return traces

    # Determine which structures will be visualized (limited by max_structures)
    visualized_struct_ids = list(aligned_structures.keys())[:max_structures]

    # Find overall center of CA coordinates for ONLY the visualized structures
    all_coords = []
    for struct_id in visualized_struct_ids:
        data = aligned_structures[struct_id]
        if 'coords' in data:
            all_coords.extend(data['coords'])

    if all_coords:
        all_coords = np.array(all_coords)
        overall_center = np.mean(all_coords, axis=0)
        print(f"Centering visualization on CA center of {len(visualized_struct_ids)} structures: {overall_center}")
    else:
        overall_center = np.array([0.0, 0.0, 0.0])

    # Apply centering offset to all structures
    for struct_id, data in aligned_structures.items():
        if 'coords' in data:
            data['coords'] = data['coords'] - overall_center
        if 'retinal' in data and data['retinal'] is not None and data['retinal'].get('coords') is not None:
            data['retinal']['coords'] = data['retinal']['coords'] - overall_center
        if 'dataframe' in data and data['dataframe'] is not None:
            df = data['dataframe']
            df[['x', 'y', 'z']] = df[['x', 'y', 'z']].values - overall_center

    # Collect all GRN positions that actually exist in the visualized structures
    existing_grns = set()
    for struct_id in visualized_struct_ids:
        data = aligned_structures[struct_id]
        df = data.get('dataframe')
        if df is not None:
            # Get all non-null GRN values from this structure
            structure_grns = df['grn'].dropna().unique()
            existing_grns.update(structure_grns)
        else:
            # Fallback: use grn_positions array
            grn_positions = data.get('grn_positions', data.get('grn', []))
            structure_grns = [grn for grn in grn_positions if pd.notna(grn)]
            existing_grns.update(structure_grns)

    # Filter to helix GRNs (1-7) from the table that exist in at least one structure
    def is_helix_grn(grn):
        """Check if GRN is a helix position (1.XX to 7.XX)"""
        if not grn or pd.isna(grn):
            return False
        grn_str = str(grn)
        if '.' not in grn_str:
            return False
        prefix = grn_str.split('.')[0]
        if not prefix.isdigit():
            return False
        helix_num = int(prefix)
        return 1 <= helix_num <= 7

    # Get all helix GRNs from the table
    helix_grns_in_table = [grn for grn in grn_df.columns if is_helix_grn(grn)]
    # Filter to those that exist in at least one structure
    all_grn_positions = sorted([grn for grn in helix_grns_in_table if grn in existing_grns],
                                key=lambda x: (int(str(x).split('.')[0]), int(str(x).split('.')[1])))
    print(f"GRN positions in table: {len(grn_df.columns)}")
    print(f"Helix GRNs (1-7) in table: {len(helix_grns_in_table)}")
    print(f"GRN positions in structures: {len(existing_grns)}")
    print(f"Helix GRN positions for slider: {len(all_grn_positions)}")

    # Helix color scheme (from opsin_color_scheme.py)
    helix_colors = {
        1: HELIX_NUMBER_COLORS[1],    # '#08306B' - cold_blue_darkest
        2: HELIX_NUMBER_COLORS[2],    # '#2171B5' - cold_blue_medium
        3: HELIX_NUMBER_COLORS[3],    # '#41B6C4' - cold_cyan_medium
        4: HELIX_NUMBER_COLORS[4],    # '#FED976' - warm_yellow_medium
        5: HELIX_NUMBER_COLORS[5],    # '#FD8D3C' - warm_orange_medium
        6: HELIX_NUMBER_COLORS[6],    # '#E31A1C' - warm_red_dark
        7: HELIX_NUMBER_COLORS[7],    # '#800026' - warm_purple_dark
        0: '#D3D3D3'  # Light gray for non-helix/loop regions
    }

    # Property color scheme (from opsin_color_scheme.py)
    from src.opsin_color_scheme import get_categorical_colors

    # Initialize property colors if property data is available
    property_colors = {}
    if property_data:
        # Get unique molecular functions from property data
        molecular_functions = set()
        for struct_id, props in property_data.items():
            if 'molecular_function' in props:
                molecular_functions.add(props['molecular_function'])

        # Create color mapping for molecular functions
        property_colors = get_categorical_colors(
            list(molecular_functions),
            property_type='property1'  # Uses WARM palette
        )
        print(f"Created property colors for {len(property_colors)} molecular functions: {list(property_colors.keys())}")

    # Function to get color for a structure based on coloring mode
    def get_structure_color(struct_id, helix_num, coloring_mode='helix'):
        if coloring_mode == 'property' and property_data and struct_id in property_data:
            mol_func = property_data[struct_id].get('molecular_function', 'Unknown')
            return property_colors.get(mol_func, '#D3D3D3')  # Default to gray if unknown
        else:
            return helix_colors.get(helix_num, '#D3D3D3')  # Default helix coloring

    # Create figure with subplots - main 3D plot and residue distribution table
    fig = make_subplots(
        rows=1, cols=2,
        column_widths=[0.75, 0.25],  # 3D plot takes 75%, table takes 25%
        specs=[[{"type": "scatter3d"}, {"type": "table"}]],
        subplot_titles=("Opsin Structure Alignment", "Residue Distribution")
    )

    # Track helix legend for both coloring modes
    helix_traces = {}
    property_traces = {}

    # Add base layer traces for BOTH coloring modes (helix and property)
    # We'll create two sets of traces and toggle their visibility

    # === HELIX COLORING TRACES ===
    structure_count = 0
    print(f"\n=== Creating helix-colored structure traces ===")
    for struct_id, data in aligned_structures.items():
        # IMPORTANT: Use the aligned coordinates, not the original dataframe coordinates
        coords = data['coords']  # These are the properly aligned coordinates
        grn_positions = data.get('grn', data['grn_positions'])
        helix_numbers = data['helix_numbers']

        # Debug: Check coordinate range for first few structures
        if structure_count < 3:
            z_min, z_max = np.min(coords[:, 2]), np.max(coords[:, 2])
            z_range = z_max - z_min
            print(f"  {struct_id}: Z-range = {z_range:.2f} (min: {z_min:.2f}, max: {z_max:.2f})")

        # Get dataframe for efficient GRN access, but update it with aligned coordinates
        df = data.get('dataframe')
        if df is not None:
            # Update dataframe with aligned coordinates
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = coords
            grn_positions = df_aligned['grn'].values
            helix_numbers = df_aligned['helix_num'].values

        # Limit number of structures for performance
        if structure_count > max_structures:
            break

        # Group residues by helix for line connectivity using efficient pandas groupby
        if df is not None:
            for helix_num, group in df_aligned.groupby('helix_num'):
                if len(group) == 0:
                    continue

                helix_coords = group[['x', 'y', 'z']].values  # Now using aligned coordinates
                helix_grn = group['grn'].values

                # Create trace name for legend grouping
                if helix_num == 0:
                    trace_name = "Loops/Non-helix"
                    legend_group = "Helix_0"
                else:
                    trace_name = f"Helix {int(helix_num)}"
                    legend_group = f"Helix_{int(helix_num)}"

                # Only show legend for first structure of each helix
                show_legend = legend_group not in helix_traces
                if show_legend:
                    helix_traces[legend_group] = True

                customdata = None
                if hover_show_residue_name:
                    residue_labels = group.apply(
                        lambda row: format_residue_annotation(
                            row.get('res_name1l'),
                            row.get('res_name3l'),
                            row.get('auth_seq_id')
                        ),
                        axis=1
                    )
                    customdata = np.array(residue_labels).reshape(-1, 1)

                hovertemplate = (
                    f'<b>{struct_id}</b><br>'
                    f'Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>'
                    'GRN: %{text}<br>'
                )
                if hover_show_residue_name:
                    hovertemplate += 'Residue: %{customdata[0]}<br>'
                hovertemplate += 'X: %{x:.2f}<br>Y: %{y:.2f}<br>Z: %{z:.2f}<extra></extra>'

                trace_kwargs = dict(
                    x=helix_coords[:, 0],
                    y=helix_coords[:, 1],
                    z=helix_coords[:, 2],
                    mode='markers+lines',
                    marker=dict(
                        size=2,
                        color=helix_colors.get(helix_num, '#D3D3D3'),
                        opacity=0.1,
                        line=dict(width=0)
                    ),
                    line=dict(
                        color=helix_colors.get(helix_num, '#D3D3D3'),
                        width=0.5
                    ),
                    name=trace_name,
                    legendgroup=legend_group,
                    showlegend=show_legend,
                    hovertemplate=hovertemplate,
                    text=[str(grn) if pd.notna(grn) else 'No GRN' for grn in helix_grn],
                    visible=True
                )
                if customdata is not None:
                    trace_kwargs['customdata'] = customdata

                fig.add_trace(go.Scatter3d(**trace_kwargs), row=1, col=1)
        else:
            # Fallback to old numpy-based grouping
            residue_ids = data.get('residues')
            if residue_ids is None:
                residue_ids = [None] * len(coords)
            residue_ids = np.array(residue_ids)

            residue_names1 = data.get('res_name1l')
            if residue_names1 is None:
                residue_names1 = ['X'] * len(coords)
            residue_names1 = np.array(residue_names1)

            residue_names3 = data.get('res_name3l')
            if residue_names3 is None:
                residue_names3 = ['UNK'] * len(coords)
            residue_names3 = np.array(residue_names3)

            for helix_num in np.unique(helix_numbers):
                helix_mask = helix_numbers == helix_num
                if not np.any(helix_mask):
                    continue

                helix_coords = coords[helix_mask]
                helix_grn = grn_positions[helix_mask]
                helix_seq_ids = residue_ids[helix_mask]
                helix_res1 = residue_names1[helix_mask]
                helix_res3 = residue_names3[helix_mask]
                customdata = None
                if hover_show_residue_name:
                    residue_labels = [
                        format_residue_annotation(res1, res3, seq_id)
                        for res1, res3, seq_id in zip(helix_res1, helix_res3, helix_seq_ids)
                    ]
                    customdata = np.array(residue_labels).reshape(-1, 1)

                # Create trace name for legend grouping
                if helix_num == 0:
                    trace_name = "Loops/Non-helix"
                    legend_group = "Helix_0"
                else:
                    trace_name = f"Helix {int(helix_num)}"
                    legend_group = f"Helix_{int(helix_num)}"

                # Only show legend for first structure of each helix
                show_legend = legend_group not in helix_traces
                if show_legend:
                    helix_traces[legend_group] = True

                hovertemplate = (
                    f'<b>{struct_id}</b><br>'
                    f'Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>'
                    'GRN: %{text}<br>'
                )
                if hover_show_residue_name:
                    hovertemplate += 'Residue: %{customdata[0]}<br>'
                hovertemplate += 'X: %{x:.2f}<br>Y: %{y:.2f}<br>Z: %{z:.2f}<extra></extra>'

                trace_kwargs = dict(
                    x=helix_coords[:, 0],
                    y=helix_coords[:, 1],
                    z=helix_coords[:, 2],
                    mode='markers+lines',
                    marker=dict(
                        size=2,
                        color=helix_colors.get(helix_num, '#D3D3D3'),
                        opacity=0.1,
                        line=dict(width=0)
                    ),
                    line=dict(
                        color=helix_colors.get(helix_num, '#D3D3D3'),
                        width=0.5
                    ),
                    name=trace_name,
                    legendgroup=legend_group,
                    showlegend=show_legend,
                    hovertemplate=hovertemplate,
                    text=[str(grn) if pd.notna(grn) else 'No GRN' for grn in helix_grn],
                    visible=True
                )
                if customdata is not None:
                    trace_kwargs['customdata'] = customdata

                fig.add_trace(go.Scatter3d(**trace_kwargs), row=1, col=1)

        structure_count += 1

    helix_base_trace_count = len(fig.data)
    print(f"Added {helix_base_trace_count} helix-colored base traces")

    # === PROPERTY COLORING TRACES ===
    if property_data:
        print(f"\n=== Creating property-colored structure traces ===")
        structure_count = 0
        for struct_id, data in aligned_structures.items():
            coords = data['coords']
            grn_positions = data.get('grn', data['grn_positions'])
            helix_numbers = data['helix_numbers']

            # Limit number of structures for performance (same as helix traces)
            if structure_count > max_structures:
                break

            # Get molecular function for this structure
            mol_func = 'Unknown'
            if struct_id in property_data:
                mol_func = property_data[struct_id].get('molecular_function', 'Unknown')

            # Get dataframe for efficient processing
            df = data.get('dataframe')
            if df is not None:
                df_aligned = df.copy()
                df_aligned[['x', 'y', 'z']] = coords
                grn_positions = df_aligned['grn'].values
                helix_numbers = df_aligned['helix_num'].values

                # Group by helix but color by property
                for helix_num, group in df_aligned.groupby('helix_num'):
                    if len(group) == 0:
                        continue

                    helix_coords = group[['x', 'y', 'z']].values
                    helix_grn = group['grn'].values
                    customdata = None
                    if hover_show_residue_name:
                        residue_labels = group.apply(
                            lambda row: format_residue_annotation(
                                row.get('res_name1l'),
                                row.get('res_name3l'),
                                row.get('auth_seq_id')
                            ),
                            axis=1
                        )
                        customdata = np.array(residue_labels).reshape(-1, 1)

                    # Create trace name for property grouping
                    trace_name = mol_func
                    legend_group = f"Property_{mol_func}"

                    # Only show legend for first structure of each property
                    show_legend = legend_group not in property_traces
                    if show_legend:
                        property_traces[legend_group] = True

                    property_color = get_structure_color(struct_id, helix_num, 'property')
                    hovertemplate = (
                        f'<b>{struct_id}</b><br>'
                        f'Function: {mol_func}<br>'
                        f'Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>'
                        'GRN: %{text}<br>'
                    )
                    if hover_show_residue_name:
                        hovertemplate += 'Residue: %{customdata[0]}<br>'
                    hovertemplate += 'X: %{x:.2f}<br>Y: %{y:.2f}<br>Z: %{z:.2f}<extra></extra>'

                    trace_kwargs = dict(
                        x=helix_coords[:, 0],
                        y=helix_coords[:, 1],
                        z=helix_coords[:, 2],
                        mode='markers+lines',
                        marker=dict(
                            size=2,
                            color=property_color,
                            opacity=0.1,
                            line=dict(width=0)
                        ),
                        line=dict(
                            color=property_color,
                            width=0.5
                        ),
                        name=trace_name,
                        legendgroup=legend_group,
                        showlegend=show_legend,
                        hovertemplate=hovertemplate,
                        text=[str(grn) if pd.notna(grn) else 'No GRN' for grn in helix_grn],
                        visible=False,
                        meta={'coloring_mode': 'property'}
                    )
                    if customdata is not None:
                        trace_kwargs['customdata'] = customdata

                    fig.add_trace(go.Scatter3d(**trace_kwargs), row=1, col=1)

            structure_count += 1

        property_base_trace_count = len(fig.data) - helix_base_trace_count
        print(f"Added {property_base_trace_count} property-colored base traces")

    base_trace_count = len(fig.data)
    print(f"Total base traces: {base_trace_count}")

    retinal_trace_count = 0
    if include_retinal:
        selected_retinal_id = select_retinal_structure(aligned_structures, retinal_reference_id)
        if selected_retinal_id:
            retinal_payload = aligned_structures[selected_retinal_id].get('retinal')
            retinal_traces = build_retinal_traces(selected_retinal_id, retinal_payload)
            for trace in retinal_traces:
                fig.add_trace(trace, row=1, col=1)
            retinal_trace_count = len(retinal_traces)
            print(f"Added {retinal_trace_count} retinal traces from {selected_retinal_id}")
        else:
            print(f"[WARN] Could not find retinal coordinates for {retinal_reference_id} or any fallback structure")
    else:
        print("Retinal overlay disabled")

    # Second pass: Collect all residues with GRN positions for highlight layer (much faster)
    # Create separate highlight data for both coloring modes
    highlight_data_helix = {}  # grn_pos -> {coords, colors, hover_text} for helix coloring
    highlight_data_property = {}  # grn_pos -> {coords, colors, hover_text} for property coloring

    structure_count = 0
    for struct_id, data in aligned_structures.items():
        if structure_count > max_structures:  # Same limit
            break

        # Use aligned coordinates with efficient dataframe approach
        coords = data['coords']  # Aligned coordinates
        df = data.get('dataframe')
        if df is not None:
            # Update dataframe with aligned coordinates
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = coords

            # Super fast: filter rows with GRN values in one operation
            grn_rows = df_aligned[df_aligned['grn'].notna()]

            for _, row in grn_rows.iterrows():
                grn_pos = row['grn']
                helix_num = row['helix_num']
                coord = [row['x'], row['y'], row['z']]  # Now using aligned coordinates
                residue_label = None
                if hover_show_residue_name:
                    residue_label = format_residue_annotation(
                        row.get('res_name1l'),
                        row.get('res_name3l'),
                        row.get('auth_seq_id')
                    )

                # Initialize helix coloring data
                if grn_pos not in highlight_data_helix:
                    highlight_data_helix[grn_pos] = {
                        'coords': [],
                        'colors': [],
                        'hover_text': []
                    }

                # Add helix-colored data
                helix_color = get_structure_color(struct_id, helix_num, 'helix')
                highlight_data_helix[grn_pos]['coords'].append(coord)
                highlight_data_helix[grn_pos]['colors'].append(helix_color)
                hover_line = f'<b>{struct_id}</b><br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}'
                if hover_show_residue_name and residue_label:
                    hover_line += f'<br>Residue: {residue_label}'
                highlight_data_helix[grn_pos]['hover_text'].append(hover_line)

                # Initialize property coloring data if property data available
                if property_data:
                    if grn_pos not in highlight_data_property:
                        highlight_data_property[grn_pos] = {
                            'coords': [],
                            'colors': [],
                            'hover_text': []
                        }

                    # Add property-colored data
                    property_color = get_structure_color(struct_id, helix_num, 'property')
                    mol_func = property_data.get(struct_id, {}).get('molecular_function', 'Unknown')
                    highlight_data_property[grn_pos]['coords'].append(coord)
                    highlight_data_property[grn_pos]['colors'].append(property_color)
                    property_hover = (
                        f'<b>{struct_id}</b><br>Function: {mol_func}<br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}'
                    )
                    if hover_show_residue_name and residue_label:
                        property_hover += f'<br>Residue: {residue_label}'
                    highlight_data_property[grn_pos]['hover_text'].append(property_hover)
        else:
            # Fallback to old method
            coords = data['coords']
            grn_positions = data['grn_positions']
            helix_numbers = data['helix_numbers']
            residue_ids = data.get('residues')
            if residue_ids is None:
                residue_ids = [None] * len(coords)
            residue_names1 = data.get('res_name1l')
            if residue_names1 is None:
                residue_names1 = ['X'] * len(coords)
            residue_names3 = data.get('res_name3l')
            if residue_names3 is None:
                residue_names3 = ['UNK'] * len(coords)

            for i, (coord, grn_pos, helix_num, seq_id, res1, res3) in enumerate(
                zip(coords, grn_positions, helix_numbers, residue_ids, residue_names1, residue_names3)
            ):
                if pd.notna(grn_pos):
                    # Initialize helix coloring data
                    if grn_pos not in highlight_data_helix:
                        highlight_data_helix[grn_pos] = {
                            'coords': [],
                            'colors': [],
                            'hover_text': []
                        }

                    # Add helix-colored data
                    helix_color = get_structure_color(struct_id, helix_num, 'helix')
                    highlight_data_helix[grn_pos]['coords'].append(coord)
                    highlight_data_helix[grn_pos]['colors'].append(helix_color)
                    residue_label = None
                    if hover_show_residue_name:
                        residue_label = format_residue_annotation(res1, res3, seq_id)
                    hover_line = f'<b>{struct_id}</b><br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}'
                    if hover_show_residue_name and residue_label:
                        hover_line += f'<br>Residue: {residue_label}'
                    highlight_data_helix[grn_pos]['hover_text'].append(hover_line)

                    # Initialize property coloring data if property data available
                    if property_data:
                        if grn_pos not in highlight_data_property:
                            highlight_data_property[grn_pos] = {
                                'coords': [],
                                'colors': [],
                                'hover_text': []
                            }

                        # Add property-colored data
                        property_color = get_structure_color(struct_id, helix_num, 'property')
                        mol_func = property_data.get(struct_id, {}).get('molecular_function', 'Unknown')
                        highlight_data_property[grn_pos]['coords'].append(coord)
                        highlight_data_property[grn_pos]['colors'].append(property_color)
                        property_hover = (
                            f'<b>{struct_id}</b><br>Function: {mol_func}<br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}'
                        )
                        if hover_show_residue_name and residue_label:
                            property_hover += f'<br>Residue: {residue_label}'
                        highlight_data_property[grn_pos]['hover_text'].append(property_hover)

        structure_count += 1

    print(f"Collected helix highlight data for {len(highlight_data_helix)} GRN positions")
    if property_data:
        print(f"Collected property highlight data for {len(highlight_data_property)} GRN positions")

    # Add highlight traces for HELIX coloring
    helix_highlight_trace_start = len(fig.data)
    first_grn = True

    for grn_pos in all_grn_positions:
        if grn_pos in highlight_data_helix:
            coords_array = np.array(highlight_data_helix[grn_pos]['coords'])
            colors = highlight_data_helix[grn_pos]['colors']
            hover_text = highlight_data_helix[grn_pos]['hover_text']

            # Add helix highlight trace (high opacity)
            fig.add_trace(go.Scatter3d(
                x=coords_array[:, 0],
                y=coords_array[:, 1],
                z=coords_array[:, 2],
                mode='markers',
                marker=dict(
                    size=4,
                    color=colors,
                    opacity=1.0,  # High opacity for highlights
                    line=dict(width=1, color='white')  # White outline for visibility
                ),
                name=f"GRN {grn_pos}",
                showlegend=False,  # Don't clutter legend
                hovertemplate='%{text}<extra></extra>',
                text=hover_text,
                visible=first_grn,  # Only first GRN visible initially
                meta={'coloring_mode': 'helix'}
            ), row=1, col=1)
            first_grn = False
        else:
            # Add empty trace to maintain indexing
            fig.add_trace(go.Scatter3d(
                x=[], y=[], z=[],
                mode='markers',
                marker=dict(size=4, opacity=1.0),
                name=f"GRN {grn_pos}",
                showlegend=False,
                visible=False,
                meta={'coloring_mode': 'helix'}
            ), row=1, col=1)

    # Add highlight traces for PROPERTY coloring
    property_highlight_trace_start = len(fig.data)
    if property_data:
        first_grn = True
        for grn_pos in all_grn_positions:
            if grn_pos in highlight_data_property:
                coords_array = np.array(highlight_data_property[grn_pos]['coords'])
                colors = highlight_data_property[grn_pos]['colors']
                hover_text = highlight_data_property[grn_pos]['hover_text']

                # Add property highlight trace (high opacity, initially hidden)
                fig.add_trace(go.Scatter3d(
                    x=coords_array[:, 0],
                    y=coords_array[:, 1],
                    z=coords_array[:, 2],
                    mode='markers',
                    marker=dict(
                        size=4,
                        color=colors,
                        opacity=1.0,  # High opacity for highlights
                        line=dict(width=1, color='white')  # White outline for visibility
                    ),
                    name=f"GRN {grn_pos}",
                    showlegend=False,  # Don't clutter legend
                    hovertemplate='%{text}<extra></extra>',
                    text=hover_text,
                    visible=False,  # Property highlights hidden by default
                    meta={'coloring_mode': 'property'}
                ), row=1, col=1)
            else:
                # Add empty trace to maintain indexing
                fig.add_trace(go.Scatter3d(
                    x=[], y=[], z=[],
                    mode='markers',
                    marker=dict(size=4, opacity=1.0),
                    name=f"GRN {grn_pos}",
                    showlegend=False,
                    visible=False,
                    meta={'coloring_mode': 'property'}
                ), row=1, col=1)

    print(f"Added {len(all_grn_positions)} highlight traces")

    # Create all table traces upfront (one for each GRN position)
    table_trace_start = len(fig.data)
    first_table = True

    for grn_pos in all_grn_positions:
        counts, total = calculate_residue_distribution(aligned_structures, grn_df, grn_pos)
        table = create_residue_distribution_table(counts, total, grn_pos)
        fig.add_trace(table, row=1, col=2)

        # Only first table visible initially
        fig.data[-1].visible = first_table
        first_table = False

    print(f"Added {len(all_grn_positions)} table traces")

    # Add membrane reference planes
    print("Adding membrane reference planes...")

    # Create coordinate ranges for the planes (based on structure extent)
    all_coords = []
    for struct_id, data in list(aligned_structures.items())[:10]:  # Sample from first 10 structures
        all_coords.extend(data['coords'])
    all_coords = np.array(all_coords)

    x_range = [np.min(all_coords[:, 0]) - 5, np.max(all_coords[:, 0]) + 5]
    y_range = [np.min(all_coords[:, 1]) - 5, np.max(all_coords[:, 1]) + 5]

    # Create grid for planes
    x_plane = np.linspace(x_range[0], x_range[1], 10)
    y_plane = np.linspace(y_range[0], y_range[1], 10)
    X_plane, Y_plane = np.meshgrid(x_plane, y_plane)

    # Add membrane block as 3D volume (translucent volume between Z = -10 and +10)
    membrane_trace_count = 0
    if show_membrane:
        # Create 3D grid for volume
        x_vol = np.linspace(x_range[0], x_range[1], 20)
        y_vol = np.linspace(y_range[0], y_range[1], 20)
        z_vol = np.linspace(-10, 10, 15)  # Membrane thickness

        X_vol, Y_vol, Z_vol = np.meshgrid(x_vol, y_vol, z_vol, indexing='ij')

        # Create volume values - uniform density throughout the membrane
        membrane_values = np.ones_like(X_vol) * 0.5  # Uniform membrane density

        # Add 3D volume trace
        fig.add_trace(go.Volume(
            x=X_vol.flatten(),
            y=Y_vol.flatten(),
            z=Z_vol.flatten(),
            value=membrane_values.flatten(),
            isomin=0.3,
            isomax=0.7,
            opacity=membrane_opacity,  # Use parameter
            surface_count=3,  # Number of isosurfaces
            colorscale=[[0, 'lightgray'], [0.5, 'silver'], [1, 'lightgray']],
            showscale=False,
            name='Membrane',
            legendgroup='Membrane',
            showlegend=True,
            hovertemplate='Membrane<extra></extra>'
        ), row=1, col=1)

        membrane_trace_count = 1  # Just the volume trace
    print(f"Added {membrane_trace_count} membrane reference traces")

    # Update trace counts for complex slider system with color modes
    total_base_traces = base_trace_count + retinal_trace_count + membrane_trace_count

    # Create complex slider and button system
    # We need to handle:
    # 1. GRN position slider
    # 2. Color mode toggle (helix vs property)

    # Helper function to create visibility array for a given state
    def create_visibility_array(grn_index, color_mode='helix', show_grn=True):
        visibility = []

        # Base structure traces (helix vs property)
        if color_mode == 'helix':
            # Show helix base traces, hide property base traces
            visibility.extend([True] * helix_base_trace_count)  # Helix base traces
            if property_data:
                visibility.extend([False] * (base_trace_count - helix_base_trace_count))  # Hide property base traces
        else:  # property mode
            # Hide helix base traces, show property base traces
            visibility.extend([False] * helix_base_trace_count)  # Hide helix base traces
            if property_data:
                visibility.extend([True] * (base_trace_count - helix_base_trace_count))  # Show property base traces

        # Retinal traces (always visible when present)
        visibility.extend([True] * retinal_trace_count)

        # Membrane trace (always visible)
        visibility.extend([True] * membrane_trace_count)

        # Highlight traces (helix vs property) - only show if show_grn is True
        if color_mode == 'helix':
            # Show helix highlights, hide property highlights
            visibility.extend([False] * len(all_grn_positions))  # All helix highlights hidden
            if show_grn:
                visibility[total_base_traces + grn_index] = True  # Only current helix highlight visible
            if property_data:
                visibility.extend([False] * len(all_grn_positions))  # All property highlights hidden
        else:  # property mode
            # Hide helix highlights, show property highlights
            visibility.extend([False] * len(all_grn_positions))  # All helix highlights hidden
            if property_data:
                visibility.extend([False] * len(all_grn_positions))  # All property highlights hidden
                if show_grn:
                    visibility[total_base_traces + len(all_grn_positions) + grn_index] = True  # Only current property highlight visible

        # Table traces (same for both modes) - only show if show_grn is True
        table_start_idx = len(visibility)
        visibility.extend([False] * len(all_grn_positions))  # All table traces hidden
        if show_grn:
            visibility[table_start_idx + grn_index] = True  # Only current table visible

        return visibility

    # Create GRN slider steps for HELIX mode
    helix_steps = []
    for i, grn_pos in enumerate(all_grn_positions):
        visibility = create_visibility_array(i, 'helix')
        step = dict(
            method="restyle",
            args=[{"visible": visibility}],
            label=str(grn_pos)
        )
        helix_steps.append(step)

    # Create GRN slider steps for PROPERTY mode (if property data available)
    property_steps = []
    if property_data:
        for i, grn_pos in enumerate(all_grn_positions):
            visibility = create_visibility_array(i, 'property')
            step = dict(
                method="restyle",
                args=[{"visible": visibility}],
                label=str(grn_pos)
            )
            property_steps.append(step)

    # Add GRN slider (starts with helix mode)
    sliders = [dict(
        active=0,
        currentvalue={"prefix": "GRN Position: "},
        pad={"t": 50},
        steps=helix_steps  # Start with helix steps
    )]

    # Add color mode toggle buttons
    buttons = []
    if property_data:
        # Helix mode button
        buttons.append(dict(
            label="Helix Coloring",
            method="update",
            args=[
                {"visible": create_visibility_array(0, 'helix')},  # Set visibility for GRN position 0
                {"sliders": [dict(
                    active=0,
                    currentvalue={"prefix": "GRN Position: "},
                    pad={"t": 50},
                    steps=helix_steps
                )]}  # Update slider steps to helix mode
            ]
        ))

        # Property mode button
        buttons.append(dict(
            label="Property Coloring",
            method="update",
            args=[
                {"visible": create_visibility_array(0, 'property')},  # Set visibility for GRN position 0
                {"sliders": [dict(
                    active=0,
                    currentvalue={"prefix": "GRN Position: "},
                    pad={"t": 50},
                    steps=property_steps
                )]}  # Update slider steps to property mode
            ]
        ))

    # Create updatemenus for the color mode toggle
    updatemenus = []
    if buttons:
        updatemenus.append(dict(
            type="buttons",
            direction="left",
            buttons=buttons,
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.01,
            xanchor="left",
            y=0.02,
            yanchor="bottom"
        ))

    # Add GRN toggle buttons (Show/Hide GRN highlights)
    grn_toggle_buttons = [
        dict(
            label="Show GRN",
            method="update",
            args=[
                {"visible": create_visibility_array(0, 'helix', show_grn=True)},
                {"sliders": [dict(
                    active=0,
                    currentvalue={"prefix": "GRN Position: "},
                    pad={"t": 50},
                    steps=helix_steps,
                    visible=True
                )]}
            ]
        ),
        dict(
            label="Hide GRN",
            method="update",
            args=[
                {"visible": create_visibility_array(0, 'helix', show_grn=False)},
                {"sliders": [dict(
                    active=0,
                    currentvalue={"prefix": "GRN Position: "},
                    pad={"t": 50},
                    steps=helix_steps,
                    visible=False
                )]}
            ]
        )
    ]
    updatemenus.append(dict(
        type="buttons",
        direction="left",
        buttons=grn_toggle_buttons,
        pad={"r": 10, "t": 10},
        showactive=True,
        x=0.30,  # Position to the right of color mode buttons
        xanchor="left",
        y=0.02,
        yanchor="bottom"
    ))

    # Update layout for clean visualization with subplots
    layout_args = dict(
        title=f'{title} (n={structure_count})',
        # Clean layout
        paper_bgcolor='white',
        plot_bgcolor='white',
        width=width,
        height=height,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            bgcolor='rgba(255,255,255,0.8)',
            bordercolor='lightgray',
            borderwidth=1
        ),
        sliders=sliders,
        margin=dict(l=0, r=0, t=80, b=50)  # More bottom margin for color mode buttons
    )

    # Add updatemenus
    if updatemenus:
        layout_args['updatemenus'] = updatemenus

    fig.update_layout(**layout_args)

    # Update 3D scene (subplot 1) - CLEAN WHITE BACKGROUND
    fig.update_scenes(
        # Hide all axes and grid for clean visualization
        xaxis=dict(
            visible=False,
            showbackground=False,
            showgrid=False,
            showline=False,
            showticklabels=False,
            title=""
        ),
        yaxis=dict(
            visible=False,
            showbackground=False,
            showgrid=False,
            showline=False,
            showticklabels=False,
            title=""
        ),
        zaxis=dict(
            visible=False,
            showbackground=False,
            showgrid=False,
            showline=False,
            showticklabels=False,
            title=""
        ),
        # Set camera for proper membrane protein viewing
        camera=dict(
            eye=dict(x=1.8, y=1.8, z=0.8),  # Side view showing Z-axis as vertical
            center=dict(x=0, y=0, z=0),     # Look at center
            up=dict(x=0, y=0, z=1)          # Z-axis points up in the view
        ),
        aspectmode='cube',  # Force 1:1:1 aspect ratio
        # Clean white background
        bgcolor='white'
    )

    return fig


def create_opsin_visualization_from_workflow(
    cache_dir="opsin_output/cache",
    property_file="property/mo_exp_ST1.csv",
    grn_file="opsin_output/curated_grn_postprocessed.csv",
    output_file="opsin_output/interactive_grn_alignment_3d.html",
    reference_id='4kkh',
    **viz_kwargs
):
    """
    Convenience function to create opsin visualization from workflow cache files.

    Args:
        cache_dir (str): Path to cache directory containing workflow results
        property_file (str): Path to property CSV file
        output_file (str): Path to save HTML visualization
        reference_id (str): Reference structure ID for alignment
        **viz_kwargs: Additional arguments passed to create_interactive_opsin_visualization

    Returns:
        plotly.graph_objects.Figure: The visualization figure
    """
    print("=== Loading RMSD Cache ===")
    cache_data = load_rmsd_cache(cache_dir)
    alignment_paths = cache_data.get('alignment_paths', {})
    print(f"Found {len(alignment_paths)} alignment paths")

    print("\n=== Loading Processed Structures ===")
    processed_structures = load_processed_structures(cache_dir)

    print("\n=== Loading GRN Table ===")
    grn_df = load_grn_table(grn_file)

    print("\n=== Loading Property Data ===")
    from src.data_processing import load_opsin_property_data
    from pathlib import Path
    property_path = Path(property_file)

    property_data = None
    if property_path.exists():
        try:
            property_result = load_opsin_property_data(property_path, processed_structures)
            if property_result and 'properties' in property_result:
                property_data = property_result['properties']
                print(f"Loaded property data for {len(property_data)} structures")
            else:
                print("No property data loaded")
        except Exception as e:
            print(f"Failed to load property data: {e}")
            property_data = None
    else:
        print(f"Property file not found: {property_path}")

    print("\n=== Extracting CA Coordinates with GRN Mapping ===")
    structures = extract_ca_coordinates_with_grn(processed_structures, grn_df,
                                               chain_id='A', use_helix_only=True)

    if not structures:
        print("No structures loaded!")
        return None

    print(f"\n=== Applying Alignment Transformations ===")
    aligned_structures = apply_alignment_transformations(
        structures, alignment_paths, reference_id
    )

    print(f"\n=== Applying Membrane Orientation ===")
    oriented_structures = apply_membrane_orientation(aligned_structures, reference_id)

    # Center all structures on global center of all coordinates
    print(f"\n=== Centering on Global Center ===")
    all_coords = []
    for struct_id, data in oriented_structures.items():
        if 'coords' in data:
            all_coords.extend(data['coords'])
    if all_coords:
        global_center = np.mean(np.array(all_coords), axis=0)
        print(f"Global center: [{global_center[0]:.2f}, {global_center[1]:.2f}, {global_center[2]:.2f}]")
        for struct_id, data in oriented_structures.items():
            if 'coords' in data:
                data['coords'] = data['coords'] - global_center
            if 'retinal' in data and data['retinal'] is not None and data['retinal'].get('coords') is not None:
                data['retinal']['coords'] = data['retinal']['coords'] - global_center
        print(f"Centered {len(oriented_structures)} structures on global center")

    print(f"\n=== Creating Interactive Visualization ===")
    fig = create_interactive_opsin_visualization(oriented_structures, grn_df, property_data, **viz_kwargs)

    # Save the plot
    fig.write_html(output_file)
    print(f"Interactive visualization saved to: {output_file}")

    return fig


def create_opsin_visualization_from_workflow_b(
    cache_dir="opsin_output/cache",
    property_file="property/mo_exp_ST1.csv",
    grn_file="opsin_output/curated_grn_postprocessed.csv",
    output_file="opsin_output/interactive_grn_alignment_b.html",
    reference_id='4kkh',
    **viz_kwargs
):
    """Generate the enhanced interactive figure with retinal overlay and amino-acid filtering."""

    print("=== Loading RMSD Cache (enhanced) ===")
    cache_data = load_rmsd_cache(cache_dir)
    alignment_paths = cache_data.get('alignment_paths', {})
    print(f"Found {len(alignment_paths)} alignment paths")

    print("\n=== Loading Processed Structures ===")
    processed_structures = load_processed_structures(cache_dir)

    print("\n=== Loading GRN Table ===")
    grn_df = load_grn_table(grn_file)

    print("\n=== Loading Property Data ===")
    from src.data_processing import load_opsin_property_data
    property_path = Path(property_file)

    property_data = None
    if property_path.exists():
        try:
            property_result = load_opsin_property_data(property_path, processed_structures)
            if property_result and 'properties' in property_result:
                property_data = property_result['properties']
                print(f"Loaded property data for {len(property_data)} structures")
            else:
                print("No property data loaded")
        except Exception as exc:
            print(f"Failed to load property data: {exc}")
            property_data = None
    else:
        print(f"Property file not found: {property_path}")

    print("\n=== Extracting CA Coordinates with GRN Mapping ===")
    structures = extract_ca_coordinates_with_grn(processed_structures, grn_df, chain_id='A', use_helix_only=True)
    if not structures:
        print("No structures loaded!")
        return None

    print("\n=== Applying Alignment Transformations ===")
    aligned_structures = apply_alignment_transformations(structures, alignment_paths, reference_id)

    print("\n=== Applying Membrane Orientation ===")
    oriented_structures = apply_membrane_orientation(aligned_structures, reference_id)

    # Center all structures on global center of all coordinates
    print("\n=== Centering on Global Center ===")
    all_coords = []
    for struct_id, data in oriented_structures.items():
        if 'coords' in data:
            all_coords.extend(data['coords'])
    if all_coords:
        global_center = np.mean(np.array(all_coords), axis=0)
        print(f"Global center: [{global_center[0]:.2f}, {global_center[1]:.2f}, {global_center[2]:.2f}]")
        for struct_id, data in oriented_structures.items():
            if 'coords' in data:
                data['coords'] = data['coords'] - global_center
            if 'retinal' in data and data['retinal'] is not None and data['retinal'].get('coords') is not None:
                data['retinal']['coords'] = data['retinal']['coords'] - global_center
        print(f"Centered {len(oriented_structures)} structures on global center")

    print("\n=== Creating Enhanced Interactive Visualization ===")
    fig = create_interactive_opsin_visualization_extended(
        oriented_structures,
        grn_df,
        property_data=property_data,
        include_retinal=True,
        hover_show_residue_name=True,
        enable_amino_acid_filter=True,
        retinal_reference_id=reference_id,
        **viz_kwargs
    )

    fig.write_html(output_file)
    print(f"Enhanced interactive visualization saved to: {output_file}")

    return fig


def main():
    """Main function to create interactive GRN visualization"""
    fig = create_opsin_visualization_from_workflow()

    if fig:
        print("\n=== Visualization Complete ===")
        print("Open the HTML file in a web browser to view the interactive visualization.")
    else:
        print("\n=== Visualization Failed ===")
        print("Could not create visualization - check error messages above.")


if __name__ == "__main__":
    main()
