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

# Import protos for structure loading
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from src.data_processing import load_experimental_dataset

# Import helix color scheme
from src.opsin_color_scheme import HELIX_NUMBER_COLORS, get_categorical_colors


# Basic 3-letter to 1-letter amino acid mapping for hover text enrichment
THREE_TO_ONE_AA = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "RET": "RET", "HOH": "HOH"  # Retinal and solvent placeholders
}


def load_rmsd_cache():
    """Load alignment paths and rotation matrices from RMSD cache"""
    cache_dir = "opsin_output/cache"
    cache_files = [f for f in os.listdir(cache_dir) if f.startswith("rmsd_cache_")]
    
    if not cache_files:
        raise FileNotFoundError("No RMSD cache files found")
    
    cache_file = os.path.join(cache_dir, cache_files[0])
    print(f"Loading cache from: {cache_file}")
    
    with open(cache_file, 'rb') as f:
        cache_data = pickle.load(f)
    
    return cache_data


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
                    processed_structures = data['processed_structures']
                    print(f"Found {len(processed_structures)} processed structures in {cache_file}")
                    return processed_structures
                elif isinstance(data, dict) and len(data) > 100:
                    print(f"Found {len(data)} structures directly in {cache_file}")
                    return data
            except Exception as e:
                print(f"Error loading {cache_file}: {e}")
                continue
    
    raise FileNotFoundError("No suitable processed structures cache found")


def load_grn_table():
    """Load and parse the GRN table"""
    grn_file = "opsin_output/curated_grn.csv"
    
    print(f"Reading GRN table from: {grn_file}")
    grn_df = pd.read_csv(grn_file, index_col=0)
    
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


def extract_ca_coordinates_with_grn(processed_structures, grn_df, chain_id='A', use_helix_only=True):
    """
    Extract CA coordinates with GRN mapping for structures in GRN table
    """
    all_structures = {}
    grn_structures = list(grn_df.index)
    
    for struct_id in grn_structures:
        if struct_id in processed_structures:
            struct_data = processed_structures[struct_id]
            
            # Use the same dataframe as RMSD calculation: 'df'
            if 'df' not in struct_data:
                print(f"Warning: No 'df' found for {struct_id}")
                continue
                
            struct_df = struct_data['df'].copy()
            
            # Apply the EXACT same filtering as in compute_all_vs_all_rmsd_improved()
            struct_df = struct_df[struct_df['auth_chain_id'] == chain_id]
            struct_df = struct_df[struct_df['group'] == 'ATOM']
            
            # Filter for helix residues if requested
            if use_helix_only:
                if 'helix_num' not in struct_df.columns:
                    print(f"[WARNING] helix_num column not found for {struct_id}. Using all residues instead.")
                else:
                    struct_helices = struct_df[struct_df['helix_num'] > 0]
                    if len(struct_helices) > 0:
                        struct_df = struct_helices
                    else:
                        print(f"[WARNING] No TM helix residues found for {struct_id}. Using all protein residues instead.")
            
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
                    ca_data_with_grn['helix_num'] = 0  # Default to 0 for non-helix when column missing

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

                all_structures[struct_id] = {
                    'coords': ca_data_with_grn[['x', 'y', 'z']].astype(float).values,
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
                    'retinal': retinal_info
                }
            else:
                print(f"Warning: No CA atoms found for {struct_id} after filtering")
        else:
            print(f"Warning: {struct_id} has GRN assignment but not found in processed structures")
    
    print(f"Extracted CA coordinates with GRN mapping for {len(all_structures)} structures")
    return all_structures


def apply_alignment_transformations(structures, alignment_paths, reference_id='MerMAID1_model_0'):
    """Apply rotation matrices to align all structures to the global reference"""
    print(f"Using {reference_id} as global reference structure")
    
    if reference_id not in structures:
        print(f"Warning: Reference {reference_id} not found in structures. Using first available.")
        reference_id = list(structures.keys())[0]
        print(f"Using {reference_id} as reference instead")
    
    aligned_structures = {}
    
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
            
        # Look for alignment path between reference and this structure
        alignment_key = (reference_id, struct_id)
        reverse_key = (struct_id, reference_id)
        
        R = None
        t = None
        
        if alignment_key in alignment_paths:
            alignment_info = alignment_paths[alignment_key]
            R = np.array(alignment_info['rotation'])
            t = np.array(alignment_info['translation'])
            print(f"Found direct alignment: {reference_id} -> {struct_id}")
            
        elif reverse_key in alignment_paths:
            alignment_info = alignment_paths[reverse_key]
            R = np.array(alignment_info['rotation'])  # Already inverted
            t = np.array(alignment_info['translation'])  # Already inverted
            print(f"Found reverse alignment: {struct_id} -> {reference_id} (pre-computed inverse)")
            
        else:
            print(f"No direct alignment path found for {struct_id}")
            aligned_structures[struct_id] = structures[struct_id].copy()
            aligned_structures[struct_id]['is_reference'] = False
            continue
        
        # Apply transformation to align structure to reference frame
        if R is not None and t is not None:
            # Method 3: Apply full transformation, then recenter to match the reference
            coords_transformed = np.dot(structures[struct_id]['coords'], R) + t

            # Recenter the transformed coordinates to match the reference structure center
            ref_center = np.mean(aligned_structures[reference_id]['coords'], axis=0)
            transformed_center = np.mean(coords_transformed, axis=0)
            coords_transformed = coords_transformed - transformed_center + ref_center

            # Copy structure data and update coordinates
            aligned_structures[struct_id] = structures[struct_id].copy()
            aligned_structures[struct_id]['coords'] = coords_transformed
            aligned_structures[struct_id]['is_reference'] = False

            # Apply same transformation to retinal coordinates if available
            if structures[struct_id].get('retinal'):
                retinal_data = structures[struct_id]['retinal']
                retinal_coords = np.dot(retinal_data['coords'], R) + t
                retinal_coords = retinal_coords - transformed_center + ref_center
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


def apply_membrane_orientation(aligned_structures, reference_id='MerMAID1_model_0'):
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
    show_membrane=True,
    include_retinal=False,
    retinal_reference_id='MerMAID1_model_0',
    hover_show_residue_name=False,
    enable_amino_acid_filter=False
):
    """
    Create interactive 3D visualization of aligned opsin structures with dual coloring modes.
    
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
    
    # First, collect all GRN positions that actually exist in the structures
    existing_grns = set()
    structure_count = 0
    for struct_id, data in aligned_structures.items():
        if structure_count > max_structures:  # Use parameter instead of hardcoded value
            break
        
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
        
        structure_count += 1
    
    # Filter to only GRNs that exist in structures, and sort them
    all_grn_positions = sorted([grn for grn in grn_df.columns if grn in existing_grns])
    print(f"GRN positions in table: {len(grn_df.columns)}")
    print(f"GRN positions in structures: {len(existing_grns)}")
    print(f"Filtered GRN positions for slider: {len(all_grn_positions)}")
    
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
                
                # Base layer: static background at low opacity
                fig.add_trace(go.Scatter3d(
                    x=helix_coords[:, 0],
                    y=helix_coords[:, 1], 
                    z=helix_coords[:, 2],
                    mode='markers+lines',
                    marker=dict(
                        size=2,
                        color=helix_colors.get(helix_num, '#D3D3D3'),
                        opacity=0.1,  # Static low opacity background
                        line=dict(width=0)
                    ),
                    line=dict(
                        color=helix_colors.get(helix_num, '#D3D3D3'),
                        width=0.5
                    ),
                    name=trace_name,
                    legendgroup=legend_group,
                    showlegend=show_legend,
                    hovertemplate=f'<b>{struct_id}</b><br>' +
                                 f'Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>' +
                                 'GRN: %{text}<br>' +
                                 'X: %{x:.2f}<br>' +
                                 'Y: %{y:.2f}<br>' +
                                 'Z: %{z:.2f}<extra></extra>',
                    text=[str(grn) if pd.notna(grn) else 'No GRN' for grn in helix_grn],
                    visible=True
                ), row=1, col=1)
        else:
            # Fallback to old numpy-based grouping
            for helix_num in np.unique(helix_numbers):
                helix_mask = helix_numbers == helix_num
                if not np.any(helix_mask):
                    continue
                    
                helix_coords = coords[helix_mask]
                helix_grn = grn_positions[helix_mask]
                
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
                
                # Base layer: static background at low opacity
                fig.add_trace(go.Scatter3d(
                    x=helix_coords[:, 0],
                    y=helix_coords[:, 1], 
                    z=helix_coords[:, 2],
                    mode='markers+lines',
                    marker=dict(
                        size=2,
                        color=helix_colors.get(helix_num, '#D3D3D3'),
                        opacity=0.1,  # Static low opacity background
                        line=dict(width=0)
                    ),
                    line=dict(
                        color=helix_colors.get(helix_num, '#D3D3D3'),
                        width=0.5
                    ),
                    name=trace_name,
                    legendgroup=legend_group,
                    showlegend=show_legend,
                    hovertemplate=f'<b>{struct_id}</b><br>' +
                                 f'Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>' +
                                 'GRN: %{text}<br>' +
                                 'X: %{x:.2f}<br>' +
                                 'Y: %{y:.2f}<br>' +
                                 'Z: %{z:.2f}<extra></extra>',
                    text=[str(grn) if pd.notna(grn) else 'No GRN' for grn in helix_grn],
                    visible=True
                ), row=1, col=1)
        
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
                    
                    # Create trace name for property grouping
                    trace_name = mol_func
                    legend_group = f"Property_{mol_func}"
                    
                    # Only show legend for first structure of each property
                    show_legend = legend_group not in property_traces
                    if show_legend:
                        property_traces[legend_group] = True
                    
                    # Property-colored trace (initially hidden)
                    property_color = get_structure_color(struct_id, helix_num, 'property')
                    fig.add_trace(go.Scatter3d(
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
                        hovertemplate=f'<b>{struct_id}</b><br>' +
                                     f'Function: {mol_func}<br>' +
                                     f'Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>' +
                                     'GRN: %{text}<br>' +
                                     'X: %{x:.2f}<br>' +
                                     'Y: %{y:.2f}<br>' +
                                     'Z: %{z:.2f}<extra></extra>',
                        text=[str(grn) if pd.notna(grn) else 'No GRN' for grn in helix_grn],
                        visible=False,  # Property traces hidden by default
                        meta={'coloring_mode': 'property'}
                    ), row=1, col=1)
            
            structure_count += 1
        
        property_base_trace_count = len(fig.data) - helix_base_trace_count
        print(f"Added {property_base_trace_count} property-colored base traces")
    
    base_trace_count = len(fig.data)
    print(f"Total base traces: {base_trace_count}")
    
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
                highlight_data_helix[grn_pos]['hover_text'].append(
                    f'<b>{struct_id}</b><br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}'
                )
                
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
                    highlight_data_property[grn_pos]['hover_text'].append(
                        f'<b>{struct_id}</b><br>Function: {mol_func}<br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}'
                    )
        else:
            # Fallback to old method
            coords = data['coords']
            grn_positions = data['grn_positions'] 
            helix_numbers = data['helix_numbers']
            
            for i, (coord, grn_pos, helix_num) in enumerate(zip(coords, grn_positions, helix_numbers)):
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
                    highlight_data_helix[grn_pos]['hover_text'].append(
                        f'<b>{struct_id}</b><br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}'
                    )
                    
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
                        highlight_data_property[grn_pos]['hover_text'].append(
                            f'<b>{struct_id}</b><br>Function: {mol_func}<br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {grn_pos}'
                        )
        
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
    total_base_traces = base_trace_count + membrane_trace_count
    
    # Create complex slider and button system
    # We need to handle:
    # 1. GRN position slider
    # 2. Color mode toggle (helix vs property)
    
    # Helper function to create visibility array for a given state
    def create_visibility_array(grn_index, color_mode='helix'):
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
        
        # Membrane trace (always visible)
        visibility.extend([True] * membrane_trace_count)
        
        # Highlight traces (helix vs property)
        if color_mode == 'helix':
            # Show helix highlights, hide property highlights
            visibility.extend([False] * len(all_grn_positions))  # All helix highlights hidden
            visibility[total_base_traces + grn_index] = True  # Only current helix highlight visible
            if property_data:
                visibility.extend([False] * len(all_grn_positions))  # All property highlights hidden
        else:  # property mode
            # Hide helix highlights, show property highlights
            visibility.extend([False] * len(all_grn_positions))  # All helix highlights hidden
            if property_data:
                visibility.extend([False] * len(all_grn_positions))  # All property highlights hidden
                visibility[total_base_traces + len(all_grn_positions) + grn_index] = True  # Only current property highlight visible
        
        # Table traces (same for both modes)
        table_start_idx = len(visibility)
        visibility.extend([False] * len(all_grn_positions))  # All table traces hidden
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
    
    # Add updatemenus if we have property data
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
        aspectmode='cube',
        # Clean white background
        bgcolor='white'
    )
    
    return fig


def create_opsin_visualization_from_workflow(
    cache_dir="opsin_output/cache",
    property_file="property/mo_exp.csv", 
    output_file="opsin_output/interactive_grn_alignment_3d.html",
    reference_id='MerMAID1_model_0',
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
    cache_data = load_rmsd_cache()
    alignment_paths = cache_data.get('alignment_paths', {})
    print(f"Found {len(alignment_paths)} alignment paths")
    
    print("\n=== Loading Processed Structures ===")
    processed_structures = load_processed_structures()
    
    print("\n=== Loading GRN Table ===")
    grn_df = load_grn_table()
    
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
    
    print(f"\n=== Creating Interactive Visualization ===")
    fig = create_interactive_opsin_visualization(oriented_structures, grn_df, property_data, **viz_kwargs)
    
    # Save the plot
    fig.write_html(output_file)
    print(f"Interactive visualization saved to: {output_file}")
    
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
