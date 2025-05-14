"""
Functions for selecting reference structures and aligning other structures to them.
These functions handle the identification of reference structures and alignment operations.
"""

import numpy as np
import pandas as pd
from protos.processing.structure.struct_alignment import get_structure_alignment


def create_seq_alignment_dicts_from_paths(alignment_paths, structure_ids, global_ref, type_reference_dict=None):
    """
    Convert alignment_paths dictionary from compute_all_vs_all_rmsd_improved() to 
    seq_alignment_dicts format required by generate_grn_msa_tables()
    
    Args:
        alignment_paths: Dictionary mapping structure pairs to alignment information
            Format: {(struct1_id, struct2_id): {
                'rotation': rotation_matrix,
                'translation': translation_vector,
                'residue_mapping': [(ref_res_id, target_res_id), ...],
                'rmsd': rmsd_value
            }}
        structure_ids: List of structure IDs
        global_ref: ID of global reference structure
        type_reference_dict: Dictionary mapping types to reference structure IDs
        
    Returns:
        seq_alignment_dicts: Dictionary with global and type alignments
            Format: {
                'global': {
                    struct_id: {ref_pos: aligned_pos, ...}
                },
                'type': {
                    type_ref: {
                        struct_id: {type_pos: struct_pos, ...}
                    }
                }
            }
    """
    # Initialize seq_alignment_dicts structure
    seq_alignment_dicts = {
        'global': {},
        'type': {}
    }
    
    # Process global alignments (all structures aligned to global reference)
    for struct_id in structure_ids:
        if struct_id == global_ref:
            continue  # Skip the reference structure itself
        
        # Check if alignment exists between global reference and this structure
        if (global_ref, struct_id) in alignment_paths:
            alignment_data = alignment_paths[(global_ref, struct_id)]
            
            # Extract residue mapping and convert to seq_alignment_dicts format
            residue_mapping = alignment_data.get('residue_mapping', [])
            
            # Create mapping dictionary from reference positions to aligned positions
            mapping_dict = {}
            for ref_pos, aligned_pos in residue_mapping:
                mapping_dict[ref_pos] = aligned_pos
            
            # Add to global alignments
            seq_alignment_dicts['global'][struct_id] = mapping_dict
    
    # For type-specific alignments, use type_reference_dict if provided
    if type_reference_dict:
        for type_name, type_ref in type_reference_dict.items():
            if type_ref == global_ref:
                continue  # Skip global reference as it's already handled
            
            # Initialize type reference in the dictionary
            if type_ref not in seq_alignment_dicts['type']:
                seq_alignment_dicts['type'][type_ref] = {}
            
            # Process alignments to this type reference
            for struct_id in structure_ids:
                if struct_id == type_ref or struct_id == global_ref:
                    continue  # Skip self-alignment and global reference
                
                # Check if alignment exists for this pair
                if (type_ref, struct_id) in alignment_paths:
                    alignment_data = alignment_paths[(type_ref, struct_id)]
                    
                    # Extract residue mapping
                    residue_mapping = alignment_data.get('residue_mapping', [])
                    
                    # Create mapping dictionary
                    mapping_dict = {}
                    for ref_pos, aligned_pos in residue_mapping:
                        mapping_dict[ref_pos] = aligned_pos
                    
                    # Add to type alignments
                    seq_alignment_dicts['type'][type_ref][struct_id] = mapping_dict
    
    return seq_alignment_dicts

def find_type_references(all_structures, type_key='type'):
    """
    Find a reference structure for each group/type.
    
    Args:
        all_structures: Dictionary of all structures or Dictionary with type_key as a mapping
        type_key: Key in structure metadata for type information or dictionary of structure types
        
    Returns:
        dict: Dictionary mapping types to reference structure IDs
    """
    # Handle the case where type_key is itself a dictionary mapping structure IDs to types
    if isinstance(type_key, dict):
        group_dict = type_key
        
        # Group structures by type
        type_groups = {}
        for struct_id, struct_type in group_dict.items():
            if struct_type not in type_groups:
                type_groups[struct_type] = []
            type_groups[struct_type].append(struct_id)
            
        # Initialize type references dictionary
        type_references = {}
        
        # If all_structures is a DataFrame (RMSD matrix), use it to find best references
        if isinstance(all_structures, pd.DataFrame):
            for struct_type, struct_ids in type_groups.items():
                if not struct_ids:
                    continue
                
                # Find structure with lowest average RMSD to others in its group
                best_ref = None
                best_avg_rmsd = float('inf')
                
                for struct_id in struct_ids:
                    if struct_id in all_structures.index:
                        # Calculate average RMSD to structures in the same group
                        group_indices = [idx for idx in struct_ids if idx in all_structures.columns]
                        if group_indices:
                            avg_rmsd = all_structures.loc[struct_id, group_indices].mean()
                            if avg_rmsd < best_avg_rmsd:
                                best_avg_rmsd = avg_rmsd
                                best_ref = struct_id
                
                # If no suitable structure found, use the first one
                if best_ref is None and struct_ids:
                    for sid in struct_ids:
                        if sid in all_structures.index:
                            best_ref = sid
                            break
                    if best_ref is None and struct_ids:
                        best_ref = struct_ids[0]
                
                if best_ref is not None:
                    type_references[struct_type] = best_ref
                    
        else:
            # If all_structures is a dictionary of structures, use resolution as before
            for struct_type, struct_ids in type_groups.items():
                if not struct_ids:
                    continue
                
                # Choose the structure with best resolution if available
                best_ref = None
                best_resolution = float('inf')
                
                for struct_id in struct_ids:
                    if struct_id not in all_structures:
                        continue
                    struct_data = all_structures[struct_id]
                    if 'metadata' in struct_data and 'resolution' in struct_data['metadata']:
                        resolution = struct_data['metadata']['resolution']
                        if resolution < best_resolution:
                            best_resolution = resolution
                            best_ref = struct_id
                
                # If no resolution data available, just take the first one
                if best_ref is None and struct_ids:
                    # Find first structure that exists in all_structures
                    for sid in struct_ids:
                        if sid in all_structures:
                            best_ref = sid
                            break
                    # If none found, just take the first in the list
                    if best_ref is None and struct_ids:
                        best_ref = struct_ids[0]
                
                if best_ref is not None:
                    type_references[struct_type] = best_ref
                    
        return type_references
        
    else:
        # Original implementation for dictionary of structures with metadata
        type_groups = {}
        for struct_id, struct_data in all_structures.items():
            if isinstance(struct_data, dict) and 'metadata' in struct_data and type_key in struct_data['metadata']:
                struct_type = struct_data['metadata'][type_key]
                if struct_type not in type_groups:
                    type_groups[struct_type] = []
                type_groups[struct_type].append(struct_id)
        
        # Find a representative for each type
        type_references = {}
        for struct_type, struct_ids in type_groups.items():
            if not struct_ids:
                continue
            
            # Choose the structure with best resolution if available
            best_ref = None
            best_resolution = float('inf')
            
            for struct_id in struct_ids:
                struct_data = all_structures[struct_id]
                if 'metadata' in struct_data and 'resolution' in struct_data['metadata']:
                    resolution = struct_data['metadata']['resolution']
                    if resolution < best_resolution:
                        best_resolution = resolution
                        best_ref = struct_id
            
            # If no resolution data available, just take the first one
            if best_ref is None and struct_ids:
                best_ref = struct_ids[0]
            
            if best_ref is not None:
                type_references[struct_type] = best_ref
        
        return type_references

def find_global_reference(all_structures, type_references):
    """
    Find a global reference structure from type references.
    
    Args:
        all_structures: Dictionary of all structures or DataFrame with RMSD values
        type_references: Dictionary mapping types to reference structure IDs
        
    Returns:
        str: ID of global reference structure
    """
    # If only one type reference, use it
    if len(type_references) == 1:
        return list(type_references.values())[0]
    
    # Detect if all_structures is a DataFrame (RMSD matrix) or dict of structures
    if isinstance(all_structures, pd.DataFrame):
        # When using RMSD DataFrame, choose reference with lowest average RMSD
        best_global_ref = None
        best_avg_rmsd = float('inf')
        
        for struct_type, struct_id in type_references.items():
            if struct_id in all_structures.index:
                # Calculate average RMSD to all other structures
                avg_rmsd = all_structures.loc[struct_id].mean()
                
                # Skip structures with suspiciously low error values (<0.1)
                # which may indicate incorrectly processed structures
                if avg_rmsd < 0.1:
                    print(f"[WARNING] Skipping {struct_id} as potential reference: suspiciously low RMSD ({avg_rmsd:.4f})")
                    continue
                    
                if avg_rmsd < best_avg_rmsd:
                    best_avg_rmsd = avg_rmsd
                    best_global_ref = struct_id
        
        # If no suitable structure found, use the first one
        if best_global_ref is None and type_references:
            best_global_ref = list(type_references.values())[0]
            
        return best_global_ref
    else:
        # Original implementation for dictionary of structures
        best_global_ref = None
        best_resolution = float('inf')
        
        for struct_type, struct_id in type_references.items():
            if struct_id not in all_structures:
                continue
                
            struct_data = all_structures[struct_id]
            
            # Check for suspiciously low RMSD values in metadata
            if 'metadata' in struct_data and 'avg_rmsd' in struct_data['metadata']:
                avg_rmsd = struct_data['metadata']['avg_rmsd']
                if avg_rmsd < 0.1:
                    print(f"[WARNING] Skipping {struct_id} as potential reference: suspiciously low RMSD ({avg_rmsd:.4f})")
                    continue
                    
            if 'metadata' in struct_data and 'resolution' in struct_data['metadata']:
                resolution = struct_data['metadata']['resolution']
                if resolution < best_resolution:
                    best_resolution = resolution
                    best_global_ref = struct_id

        # BEST_GLOBAL_REF
        best_global_ref = 'CnChR2_J230_refine9'
        print("Hard set CnChR2_J230_refine9 to global ref.")

        # If no resolution data available, just take the first one
        if best_global_ref is None and type_references:
            best_global_ref = list(type_references.values())[0]
        
        return best_global_ref

def apply_alignment(coords, alignment):
    """
    Apply rotation and translation to coordinates.
    
    Args:
        coords: Coordinates to transform
        alignment: Tuple of (rotation_matrix, translation_vector, path, rmsd)
        
    Returns:
        ndarray: Transformed coordinates
    """
    R, t, _, _ = alignment
    return np.dot(coords, R.T) + t

def create_sequence_alignment_dict(ref_struct, target_struct, alignment):
    """
    Create a mapping of sequence IDs based on alignment.
    
    Args:
        ref_struct: Reference structure DataFrame
        target_struct: Target structure DataFrame
        alignment: Alignment object
        
    Returns:
        dict: Dictionary mapping target sequence IDs to reference sequence IDs
    """
    # Get unique residues in reference and target
    ref_seq_ids = sorted(ref_struct['auth_seq_id'].unique())
    target_seq_ids = sorted(target_struct['auth_seq_id'].unique())
    
    # Get aligned coordinates
    ref_coords = ref_struct[['x', 'y', 'z']].values
    target_coords = target_struct[['x', 'y', 'z']].values
    aligned_target_coords = apply_alignment(target_coords, alignment)
    
    # Calculate distance matrix between all residues
    dist_mat = np.zeros((len(target_seq_ids), len(ref_seq_ids)))
    
    for i, target_id in enumerate(target_seq_ids):
        target_residue = target_struct[target_struct['auth_seq_id'] == target_id]
        target_ca_idx = target_residue.index[0]
        target_ca_coords = aligned_target_coords[target_struct.index.get_loc(target_ca_idx)]
        
        for j, ref_id in enumerate(ref_seq_ids):
            ref_residue = ref_struct[ref_struct['auth_seq_id'] == ref_id]
            ref_ca_coords = ref_coords[ref_struct.index.get_loc(ref_residue.index[0])]
            
            # Calculate distance
            dist = np.linalg.norm(target_ca_coords - ref_ca_coords)
            dist_mat[i, j] = dist
    
    # Create mapping based on closest residue
    seq_mapping = {}
    for i, target_id in enumerate(target_seq_ids):
        closest_idx = np.argmin(dist_mat[i])
        ref_id = ref_seq_ids[closest_idx]
        seq_mapping[target_id] = ref_id
    
    return seq_mapping

def align_structure_to_reference(target_id, ref_id, structures):
    """
    Align a target structure to a reference structure.
    
    Args:
        target_id: ID of target structure
        ref_id: ID of reference structure
        structures: Dictionary of structures
        
    Returns:
        tuple: (Alignment object, sequence alignment dictionary)
    """
    if target_id not in structures or ref_id not in structures:
        return None, {}
    
    target_data = structures[target_id]
    ref_data = structures[ref_id]
    
    # Get CA atoms for alignment
    if 'df_ca_norm' in target_data:
        target_ca = target_data['df_ca_norm']
    else:
        target_df = target_data.get('df_norm', target_data.get('df'))
        target_ca = target_df[target_df['res_atom_name'] == 'CA']
    
    if 'df_ca_norm' in ref_data:
        ref_ca = ref_data['df_ca_norm']
    else:
        ref_df = ref_data.get('df_norm', ref_data.get('df'))
        ref_ca = ref_df[ref_df['res_atom_name'] == 'CA']
    
    # Ensure coordinates are numeric
    for coord in ['x', 'y', 'z']:
        target_ca[coord] = pd.to_numeric(target_ca[coord], errors='coerce')
        ref_ca[coord] = pd.to_numeric(ref_ca[coord], errors='coerce')
    
    # Get coordinates for alignment
    target_coords = target_ca[['x', 'y', 'z']].values
    ref_coords = ref_ca[['x', 'y', 'z']].values
    
    # Perform alignment
    try:
        alignment = get_structure_alignment(target_coords, ref_coords)
        
        # Create sequence mapping
        seq_mapping = create_sequence_alignment_dict(ref_ca, target_ca, alignment)
        
        return alignment, seq_mapping
    except Exception as e:
        print(f"[ERROR] Failed to align {target_id} to {ref_id}: {str(e)}")
        return None, {}

def update_structures_alignment(structures, rmsd_df=None, group_dict=None, global_ref=None, align_based_on='prot'):
    """
    Align all structures to appropriate references.
    
    Args:
        structures: Dictionary of structures
        rmsd_df: Optional DataFrame with RMSD values between structures
        group_dict: Optional dictionary mapping structure IDs to groups/types
        global_ref: Optional global reference structure ID
        align_based_on: 'prot' for protein-based alignment or 'rmsd' for RMSD-based alignment
        
    Returns:
        dict: Dictionary with alignment information
    """
    # Determine how to find type references
    if align_based_on == 'rmsd' and rmsd_df is not None and group_dict is not None:
        # Use RMSD matrix and group dictionary
        print("[INFO] Finding type references based on RMSD matrix")
        type_references = find_type_references(rmsd_df, group_dict)
    else:
        # Use structure metadata
        print("[INFO] Finding type references based on structure metadata")
        type_references = find_type_references(structures)
    
    # Find global reference if not provided
    if global_ref is None:
        if align_based_on == 'rmsd' and rmsd_df is not None:
            global_ref = find_global_reference(rmsd_df, type_references)
        else:
            global_ref = find_global_reference(structures, type_references)
    
    print(f"[INFO] Using {global_ref} as global reference")
    print(f"[INFO] Type references: {type_references}")
    
    # Initialize alignment dictionaries
    alignments = {}
    seq_alignment_dicts = {
        'global': {},
        'type': {}
    }
    
    # Align each type reference to global reference
    for struct_type, type_ref in type_references.items():
        if type_ref == global_ref:
            continue  # Skip global reference
        
        # Align type reference to global reference
        alignment, seq_mapping = align_structure_to_reference(type_ref, global_ref, structures)
        
        if alignment is not None:
            alignments[(type_ref, global_ref)] = alignment
            seq_alignment_dicts['global'][type_ref] = seq_mapping
            seq_alignment_dicts['type'][type_ref] = {}
    
    # Align each structure to its type reference
    for struct_id, struct_data in structures.items():
        if struct_id == global_ref:
            continue  # Skip global reference
        
        # Determine type reference
        struct_type = None
        if 'metadata' in struct_data and type_key in struct_data['metadata']:
            struct_type = struct_data['metadata'][type_key]
        
        if struct_type in type_references:
            type_ref = type_references[struct_type]
            
            # Skip if this is the type reference
            if struct_id == type_ref:
                continue
            
            # Align to type reference
            alignment, seq_mapping = align_structure_to_reference(struct_id, type_ref, structures)
            
            if alignment is not None:
                alignments[(struct_id, type_ref)] = alignment
                seq_alignment_dicts['type'].setdefault(type_ref, {})[struct_id] = seq_mapping
    
    # Create direct alignments to global reference for structures without type
    for struct_id, struct_data in structures.items():
        if struct_id == global_ref:
            continue  # Skip global reference
        
        # Check if already aligned through a type reference
        struct_type = None
        if 'metadata' in struct_data and type_key in struct_data['metadata']:
            struct_type = struct_data['metadata'][type_key]
        
        if struct_type in type_references and struct_id != type_references[struct_type]:
            continue  # Already aligned through type reference
        
        # Direct alignment to global reference
        alignment, seq_mapping = align_structure_to_reference(struct_id, global_ref, structures)
        
        if alignment is not None:
            alignments[(struct_id, global_ref)] = alignment
            seq_alignment_dicts['global'][struct_id] = seq_mapping
    
    return {
        'global_ref': global_ref,
        'type_references': type_references,
        'alignments': alignments,
        'seq_alignment_dicts': seq_alignment_dicts
    }



