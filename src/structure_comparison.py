"""
Functions for comparing and aligning structures.
These functions handle the calculation of RMSD between structure pairs
and comparison of binding pockets.
"""
import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from protos.analysis.structure.alignment import get_structure_alignment
from Bio.PDB.qcprot import QCPSuperimposer
from tqdm import tqdm
import pickle
import hashlib

# Import error_analysis module from the local package
from src.common_utils import compute_retinal_mean_closest_distance

from src.visualization_functions import create_and_visualize_similarity_tree, visualize_rmsd_heatmap

import matplotlib.pyplot as plt

HELIX_RET_BANDS = {
    1: (14.0, 17.8),
    2: (14.0, 14.0),
    3: (10.0, 11.0),
    4: (7.5, 14.0),
    5: (9.0, 13.0),
    6: (8.0, 17.0),
    7: (12.0, 13.0),
}

# Import retinal_utils
try:
    from protos.processing.opsin.retinal_utils import calculate_min_distances
except ImportError:
    print("[WARNING] Could not import calculate_min_distances from retinal_utils")
    # Define a simple substitute function if import fails
    def calculate_min_distances(coords1, coords2):
        """Calculate minimum distances between atoms in coords1 and coords2."""
        distances = cdist(coords1, coords2)
        return np.min(distances, axis=1)


def compute_all_vs_all_rmsd_improved(structures, align_to=None, subset='CA', chain_id='A', tm_score_threshold=0.0, 
                              speed=9, verbose=False, use_helix_only=True, cache_dir=None, force_recompute=False):
    """
    Calculate RMSD between all pairs of structures using only C-alpha atoms by default.
    Optionally align all structures to a reference structure first.
    Can use either GTalign (recommended) or basic alignment for RMSD calculation.
    
    Args:
        structures: Dictionary of structures to compare
        align_to: Optional PDB ID to use as reference for alignment
        subset: Atom subset to use ('CA' or 'backbone') - default is 'CA' for C-alpha atoms only
        chain_id: Chain ID to use
        tm_score_threshold: Minimum TM-score threshold for GTalign
        speed: Speed setting for GTalign (0-13, higher is faster but potentially less accurate)
        verbose: Whether to print verbose outputs during alignment
        use_helix_only: Whether to use only residues from helices 1-7 for alignment (requires helix_num column)
        cache_dir: Directory to store/load cached RMSD results (if None, results won't be cached)
        force_recompute: If True, ignore cached results and recompute RMSD values
        
    Returns:
        tuple: (rmsd_df, rmsd_matrix, pdb_list, alignment_paths) 
               - rmsd_df: DataFrame with RMSD values
               - rmsd_matrix: NumPy array with RMSD values
               - pdb_list: List of structure IDs
               - alignment_paths: Dictionary mapping structure pairs to alignment information
    """
    
    # Prepare the list of structures
    structure_ids = list(structures.keys())
    
    # Generate a unique cache key based on input parameters
    cache_key = f"{subset}_{chain_id}_{tm_score_threshold}_{speed}_{use_helix_only}_{sorted(structure_ids)}"
    hash_key = hashlib.md5(cache_key.encode()).hexdigest()
    cache_file = None
    
    # If cache directory provided, set up cache file path
    if cache_dir:
        if not os.path.exists(cache_dir):
            os.makedirs(cache_dir, exist_ok=True)
        cache_file = os.path.join(cache_dir, f"rmsd_cache_{hash_key}.pkl")
    
    # Try to load cached results if available and not forcing recomputation
    if cache_file and os.path.exists(cache_file) and not force_recompute:
        try:
            print(f"Loading cached RMSD results from {cache_file}")
            with open(cache_file, 'rb') as f:
                cached_data = pickle.load(f)
                
            # Verify cached data has all structures we need
            cached_ids = set(cached_data['structure_ids'])
            required_ids = set(structure_ids)
            
            if cached_ids >= required_ids:  # Check if cache contains all needed structures
                # Extract just the data we need for our current structures
                rmsd_df = cached_data['rmsd_df'].loc[structure_ids, structure_ids]
                
                # Reconstruct the matrix from the DataFrame
                rmsd_matrix = rmsd_df.values
                
                # Extract alignment paths if available
                alignment_paths = cached_data.get('alignment_paths', {})
                
                print(f"Successfully loaded cached RMSD results for {len(structure_ids)} structures")
                return rmsd_df, rmsd_matrix, structure_ids, alignment_paths
            else:
                missing_ids = required_ids - cached_ids
                print(f"Cache is missing {len(missing_ids)} structures. Recomputing RMSD matrix.")
        except Exception as e:
            print(f"Error loading cached RMSD results: {str(e)}. Recomputing RMSD matrix.")
    n_structures = len(structure_ids)
    
    # Initialize the RMSD matrix
    rmsd_matrix = np.zeros((n_structures, n_structures))
    
    # Initialize dictionary to store alignment paths
    alignment_paths = {}

    # Calculate RMSD for all pairs with progress bar
    print(f"Computing RMSD matrix using {'helix-only' if use_helix_only else 'full structure'} alignment for {len(structure_ids)} structures...")
    # Calculate total number of comparisons for the progress bar
    total_comparisons = len(structure_ids) * (len(structure_ids) - 1) // 2

    # Create progress bar
    pbar = tqdm(total=total_comparisons, desc="RMSD Calculations", unit="pairs")

    for i, struct1_id in enumerate(structure_ids):
        for j, struct2_id in enumerate(structure_ids):

            # Only compute for i < j (upper triangle) to avoid duplicate calculations
            if i < j:
                # Extract structure coordinates
                structure1 = structures[struct1_id]
                structure2 = structures[struct2_id]
                # Extract appropriate dataframes

                struct1_df = structure1['df'].copy()
                struct2_df = structure2['df'].copy()

                # Filter and process structures

                # Filter for the specified chain
                struct1_df = struct1_df[struct1_df['auth_chain_id'] == chain_id]
                struct2_df = struct2_df[struct2_df['auth_chain_id'] == chain_id]

                # Filter for protein
                struct1_df = struct1_df[struct1_df['group'] == 'ATOM']
                struct2_df = struct2_df[struct2_df['group'] == 'ATOM']

                # Filter for helix residues (1-7) if requested
                if use_helix_only:
                    if 'helix_num' not in struct1_df.columns or 'helix_num' not in struct2_df.columns:
                        print(f"[WARNING] helix_num column not found for {struct1_id} or {struct2_id}. Using all residues instead.")
                    else:
                        # Use helix-annotated residues

                        # Filter for residues with helix_num 1-7
                        struct1_helices = struct1_df[struct1_df['helix_num']>0]
                        struct2_helices = struct2_df[struct2_df['helix_num']>0]

                        # Check if we found any helix residues
                        if len(struct1_helices) > 0 and len(struct2_helices) > 0:
                            struct1_df = struct1_helices
                            struct2_df = struct2_helices
                            # Successfully filtered for helix residues
                        else:
                            print(f"[WARNING] No TM helix residues (1-7) found for {struct1_id} or {struct2_id}. Using all protein residues instead.")

                # Filter for the specified atom subset
                struct1_df = struct1_df[struct1_df['res_atom_name'] == 'CA']
                struct2_df = struct2_df[struct2_df['res_atom_name'] == 'CA']

                # Extract coordinates
                struct1_coords = struct1_df[['x', 'y', 'z']].astype(float).values
                struct2_coords = struct2_df[['x', 'y', 'z']].astype(float).values

                # Prepare coordinates for alignment
                ref_seq_ids = struct1_df['auth_seq_id'].values if 'auth_seq_id' in struct1_df.columns else list(range(len(struct1_coords)))
                target_seq_ids = struct2_df['auth_seq_id'].values if 'auth_seq_id' in struct2_df.columns else list(range(len(struct2_coords)))

                # Check if we have enough coordinates for alignment
                if len(struct1_coords) < 3 or len(struct2_coords) < 3:
                    print(f"[WARNING] Not enough atoms for alignment between {struct1_id} ({len(struct1_coords)} atoms) and {struct2_id} ({len(struct2_coords)} atoms).")
                    rmsd_matrix[i, j] = np.nan
                    rmsd_matrix[j, i] = np.nan
                    pbar.update(1)
                    continue

                # Calculate RMSD
                try:
                    # Attempt alignment of structures
                    rotation, translation, best_path, rmsd = get_structure_alignment(struct1_coords, struct2_coords)

                    # Extract indices from alignment path
                    ref_indices, target_indices = best_path

                    # Map indices to auth_seq_id values
                    ref_res_ids = [ref_seq_ids[idx] for idx in ref_indices]
                    target_res_ids = [target_seq_ids[idx] for idx in target_indices]

                    # Create mapping between reference and target residue IDs
                    residue_mapping = list(zip(ref_res_ids, target_res_ids))

                    # Store the alignment path information
                    alignment_paths[(struct1_id, struct2_id)] = {
                        'rotation': rotation.tolist(),
                        'translation': translation.tolist(),
                        'residue_mapping': residue_mapping,
                        'rmsd': rmsd
                    }

                    # Also store the reverse mapping for convenience
                    alignment_paths[(struct2_id, struct1_id)] = {
                        'rotation': np.transpose(rotation).tolist(),  # Transpose for reverse rotation
                        'translation': (-np.dot(np.transpose(rotation), translation)).tolist(),  # Reverse translation
                        'residue_mapping': [(y, x) for x, y in residue_mapping],  # Reverse the mapping
                        'rmsd': rmsd
                    }

                    rmsd_matrix[i, j] = rmsd
                    # Also set the symmetric value in the lower triangle
                    rmsd_matrix[j, i] = rmsd
                    # RMSD calculation successful
                except MemoryError:
                    print(f"[ERROR] Memory error when calculating RMSD between {struct1_id} and {struct2_id}. Skipping pair.")
                    rmsd_matrix[i, j] = np.nan
                    rmsd_matrix[j, i] = np.nan
                except Exception as e:
                    print(f"[WARNING] Failed to calculate RMSD between {struct1_id} and {struct2_id}: {str(e)}")
                    # Print more detailed error information
                    import traceback
                    traceback.print_exc()
                    rmsd_matrix[i, j] = np.nan
                    rmsd_matrix[j, i] = np.nan

            pbar.update(1)

    # Close the progress bar
    pbar.close()
    
    # Create DataFrame from RMSD matrix
    rmsd_df = pd.DataFrame(rmsd_matrix, index=structure_ids, columns=structure_ids)
    
    # Save the results to cache if a cache directory was provided
    if cache_file:
        try:
            print(f"Saving RMSD results to cache: {cache_file}")
            # Create a dictionary with all the data needed to reconstruct the results
            cache_data = {
                'rmsd_df': rmsd_df,
                'structure_ids': structure_ids,
                'alignment_paths': alignment_paths,  # Include alignment paths in cache
                'parameters': {
                    'subset': subset,
                    'chain_id': chain_id,
                    'tm_score_threshold': tm_score_threshold,
                    'speed': speed,
                    'use_helix_only': use_helix_only
                }
            }
            
            # Write to a temporary file first, then rename to avoid corruption if the process is interrupted
            tmp_file = f"{cache_file}.tmp"
            with open(tmp_file, 'wb') as f:
                pickle.dump(cache_data, f)
            
            # Use atomic rename operation
            os.replace(tmp_file, cache_file)
            print(f"Successfully saved RMSD results for {len(structure_ids)} structures")
        except Exception as e:
            print(f"Error saving RMSD results to cache: {str(e)}")
            # Continue even if saving fails - this doesn't affect the returned results
    
    # Return the DataFrame, matrix, structure IDs, and alignment paths
    return rmsd_df, rmsd_matrix, structure_ids, alignment_paths


# Assume get_structure_alignment (which uses CEalign),
# calculate_min_distances, and compute_retinal_mean_closest_distance
# are correctly defined and imported in your environment.

def apply_retinal_distance_band_filter(
    ca_df: pd.DataFrame,
    retinal_df: pd.DataFrame,
    bands: dict,
    *,
    helix_col: str = "helix_num",
    seqid_col: str = "auth_seq_id",
    atom_col: str = "res_atom_name",
    atom_name: str = "CA",
    margin: float = 1.5,
    min_per_helix: int = 4,
    min_total: int = 12
):
    """
    Keep only CA residues per helix whose distance to RET falls within a helix-specific band.
    - bands: {helix_idx: (low, high)} in Å. We'll auto-order low<=high and expand by ±margin.
    - Distances are computed as min distance from CA to any RET atom (like your pocket def).
    Returns: (filtered_df, stats_dict)
    Falls back to original ca_df if filters would drop below min_total.
    """
    stats = {"kept_per_helix": {}, "dropped_per_helix": {}, "kept_total": 0, "dropped_total": 0, "used": False}

    if ca_df.empty or retinal_df.empty or helix_col not in ca_df.columns:
        # Nothing to do / no helix annotation / no RET → return unmodified
        return ca_df, stats

    # Ensure we work on CA only
    x = ca_df[ca_df[atom_col].astype(str) == atom_name].copy()
    if x.empty:
        return ca_df, stats

    # Precompute CA→RET distances (min to any RET atom)
    ca_xyz = x[["x", "y", "z"]].astype(float).values
    ret_xyz = retinal_df[["x", "y", "z"]].astype(float).values
    dmin = calculate_min_distances(ca_xyz, ret_xyz)  # shape (N_ca,)
    if not isinstance(dmin, np.ndarray):
        dmin = np.array(dmin, dtype=float)
    x["dist_to_ret"] = dmin

    keep_mask = np.zeros(len(x), dtype=bool)
    for h in range(1, 8):
        h_mask = (x[helix_col] == h).values
        if not h_mask.any():
            continue
        lo, hi = bands.get(h, (None, None))
        if lo is None or hi is None:
            # No band for this helix → keep all of its residues
            keep_mask[h_mask] = True
            stats["kept_per_helix"][h] = int(h_mask.sum())
            stats["dropped_per_helix"][h] = 0
            continue
        # normalize band and expand with margin (handle narrow bands like 14–14)
        lo, hi = (min(lo, hi), max(lo, hi))
        lo -= margin
        hi += margin
        h_keep = h_mask & (x["dist_to_ret"].values >= lo) & (x["dist_to_ret"].values <= hi)

        # enforce a minimum per helix; if too few, relax to keep all for that helix
        if h_keep.sum() < min_per_helix:
            h_keep = h_mask  # relax: keep all residues of this helix

        keep_mask |= h_keep
        stats["kept_per_helix"][h] = int((h_keep).sum())
        stats["dropped_per_helix"][h] = int(h_mask.sum() - (h_keep).sum())

    kept = x[keep_mask]
    # If we ended up with too few residues overall, fall back
    if len(kept) < min_total:
        stats["kept_total"] = len(x)
        stats["dropped_total"] = 0
        stats["used"] = False
        return ca_df, stats  # return original (unfiltered)

    # Merge back any non-CA rows (you only align on CA, so not strictly needed here)
    stats["kept_total"] = len(kept)
    stats["dropped_total"] = int(len(x) - len(kept))
    stats["used"] = True
    return kept, stats


def calculate_binding_pocket_rmsd_for_pairs(mapping_dict, exp_processor, pred_processor, pocket_def=None, cutoff=6.0,
                                            distance_cutoff=6.0, position_tolerance=2.0, window_size=20, max_gap=4,
                                            retinal_name='RET'):
    """
    Calculate binding pocket RMSD between pairs of experimental and predicted structures.
    Uses CEalign path for pocket residue mapping. Retinal RMSD uses original logic.
    """
    effective_cutoff = distance_cutoff if distance_cutoff is not None else cutoff
    results = {}

    print(f"Processing {len(mapping_dict)} structure pairs for binding pocket RMSD...")

    for exp_id, mapping_info in tqdm(mapping_dict.items(), desc="Binding Pocket RMSD", unit="pairs"):
        if isinstance(mapping_info, dict) and 'predicted' in mapping_info:
            pred_id = mapping_info['predicted']
        else:
            pred_id = mapping_info

        print(f"[INFO] Processing {exp_id} <-> {pred_id} pair")

        exp_data_full = exp_processor.data[exp_processor.data['pdb_id'] == exp_id]
        if exp_data_full.empty:
            results[exp_id] = {'error': f'No exp data for {exp_id}'};
            continue
        pred_data_full = pred_processor.data[pred_processor.data['pdb_id'] == pred_id]
        if pred_data_full.empty:
            results[exp_id] = {'error': f'No pred data for {pred_id}'};
            continue

        exp_df = exp_data_full.copy()
        # Keep a pristine copy of pred_df for retinal RMSD, as it uses original pred retinal coords then transforms them
        pred_df_original_for_retinal_rmsd = pred_data_full.copy()
        # This copy will be transformed for pocket RMSD calculation
        pred_df_for_pocket_alignment = pred_data_full.copy()

        coord_cols = ['x', 'y', 'z']
        for df_to_clean in [exp_df, pred_df_original_for_retinal_rmsd, pred_df_for_pocket_alignment]:
            for col in coord_cols:
                df_to_clean[col] = pd.to_numeric(df_to_clean[col], errors='coerce')
            df_to_clean.dropna(subset=coord_cols, inplace=True)

        if exp_df.empty or pred_df_original_for_retinal_rmsd.empty or pred_df_for_pocket_alignment.empty:
            results[exp_id] = {'error': 'Data empty after coord cleaning'};
            continue

        # STEP 1: Extract retinal DataFrames from UNTRANSFORMED DataFrames
        # Experimental retinal (from original exp_df)
        exp_ret_df_step1 = exp_df[exp_df['res_name3l'] == retinal_name]
        if exp_ret_df_step1.empty:
            results[exp_id] = {'error': f'No exp retinal ({retinal_name})'};
            continue

        # Predicted retinal (from original, untransformed pred_df_original_for_retinal_rmsd)
        pred_ret_df_step1 = pred_df_original_for_retinal_rmsd[
            pred_df_original_for_retinal_rmsd['res_name3l'] == retinal_name]
        if pred_ret_df_step1.empty and retinal_name == 'RET':
            pred_ret_df_step1 = pred_df_original_for_retinal_rmsd[
                pred_df_original_for_retinal_rmsd['res_name3l'] == 'LIG']
            if not pred_ret_df_step1.empty:
                # If renaming, ensure all pred_df copies are consistent if 'LIG' might be used elsewhere
                pred_df_original_for_retinal_rmsd.loc[
                    pred_df_original_for_retinal_rmsd['res_name3l'] == 'LIG', 'res_name3l'] = 'RET'
                pred_df_for_pocket_alignment.loc[
                    pred_df_for_pocket_alignment['res_name3l'] == 'LIG', 'res_name3l'] = 'RET'
                pred_ret_df_step1 = pred_df_original_for_retinal_rmsd[
                    pred_df_original_for_retinal_rmsd['res_name3l'] == 'RET']

        if pred_ret_df_step1.empty:  # This is the untransformed predicted retinal DataFrame
            results[exp_id] = {'error': f'No pred retinal ({retinal_name} or LIG)'};
            continue

        # STEP 2: Define EXPERIMENTAL binding pocket residues
        exp_protein_ca_df_for_pocket = exp_df[
            (exp_df['res_atom_name'].astype(str) == 'CA') & (exp_df['res_name3l'] != retinal_name)]
        if exp_protein_ca_df_for_pocket.empty:
            results[exp_id] = {'error': 'No exp protein CAs for pocket def'};
            continue

        exp_pocket_definition_ca_coords = exp_protein_ca_df_for_pocket[['x', 'y', 'z']].values.astype(float)
        exp_ret_coords_for_pocket_def = exp_ret_df_step1[['x', 'y', 'z']].values.astype(float)

        if exp_pocket_definition_ca_coords.size == 0 or exp_ret_coords_for_pocket_def.size == 0:
            results[exp_id] = {'error': 'Empty coords for pocket def'};
            continue

        dist_to_ret = calculate_min_distances(exp_pocket_definition_ca_coords, exp_ret_coords_for_pocket_def)
        if not isinstance(dist_to_ret, np.ndarray): dist_to_ret = np.array(dist_to_ret)

        pocket_ca_indices_in_exp_protein_ca_df = np.where(dist_to_ret <= effective_cutoff)[0]
        experimental_pocket_auth_ids = exp_protein_ca_df_for_pocket.iloc[pocket_ca_indices_in_exp_protein_ca_df][
            'auth_seq_id'].unique().tolist()

        if pocket_def is not None:
            experimental_pocket_auth_ids = pocket_def

        if not experimental_pocket_auth_ids:
            print(f"[WARNING] No pocket residues defined for {exp_id}. Pocket RMSD will be NaN.")
        print(f"[INFO] Identified {len(experimental_pocket_auth_ids)} experimental pocket residues for {exp_id}")

        # STEP 3: Extract CA atoms for overall alignment
        # These are the DataFrames from which coordinates for CEalign will be extracted
        exp_ca_for_align_df = exp_df[
            (exp_df['res_atom_name'].astype(str) == 'CA') & (exp_df['res_name3l'] != retinal_name)]
        # Use the pred_df copy that will be transformed for pocket analysis
        pred_ca_for_align_df = pred_df_for_pocket_alignment[(pred_df_for_pocket_alignment['res_name3l'] != 'RET') &
                                                            (pred_df_for_pocket_alignment['res_atom_name'].astype(
                                                                str) == 'CA')]

        if exp_ca_for_align_df.empty or pred_ca_for_align_df.empty:
            results[exp_id] = {'error': 'No CAs for alignment'};
            continue

        print(
            f"[INFO] Using {len(exp_ca_for_align_df)} CAs from exp and {len(pred_ca_for_align_df)} CAs from pred for CEalign input.")

        try:
            exp_ca_for_align_df_filt, exp_band_stats = apply_retinal_distance_band_filter(
                exp_ca_for_align_df, exp_ret_df_step1, HELIX_RET_BANDS,
                margin=1.5, min_per_helix=4, min_total=12
            )
            pred_ca_for_align_df_filt, pred_band_stats = apply_retinal_distance_band_filter(
                pred_ca_for_align_df, pred_ret_df_step1, HELIX_RET_BANDS,
                margin=1.5, min_per_helix=4, min_total=12
            )

            # Adopt filtered sets if both sides passed min_total; else keep originals
            if exp_band_stats.get("used") and pred_band_stats.get("used"):
                print(f"[INFO] Retinal-band filter applied: kept EXP {exp_band_stats['kept_total']} CAs "
                      f"(dropped {exp_band_stats['dropped_total']}), "
                      f"PRED {pred_band_stats['kept_total']} CAs (dropped {pred_band_stats['dropped_total']}).")
                exp_ca_for_align_df = exp_ca_for_align_df_filt
                pred_ca_for_align_df = pred_ca_for_align_df_filt
            else:
                print(
                    "[INFO] Retinal-band filter skipped (insufficient residues after filtering on one or both sides).")
        except Exception as e_band:
            print(f"[WARN] Retinal-band filter failed; proceeding without it: {e_band}")

        """
        # STEP 4: Align structures based on CA atoms
        try:
            exp_ca_coords_for_cealign = exp_ca_for_align_df[['x', 'y', 'z']].astype(float).values
            pred_ca_coords_for_cealign = pred_ca_for_align_df[['x', 'y', 'z']].astype(float).values

            # get_structure_alignment now uses CEalign and returns path indices
            # The overall_ca_rmsd from this is based on QCP of the CEalign path.
            R, t, alignment_path_indices, overall_ca_rmsd = get_structure_alignment(
                exp_ca_coords_for_cealign, pred_ca_coords_for_cealign,
                window_size=8, max_gap=30  # CEalign specific params
            )
            print(f"[INFO] Overall CA RMSD (CEalign path + QCP): {overall_ca_rmsd:.3f}Å")

            # Apply transformation to pred_df_for_pocket_alignment ONCE. This df is used for pocket RMSD.
            pred_all_atom_coords_original_for_pocket = pred_df_for_pocket_alignment[['x', 'y', 'z']].values.astype(
                float)
            pred_df_for_pocket_alignment[['x', 'y', 'z']] = np.dot(pred_all_atom_coords_original_for_pocket, R) + t

            # Map experimental pocket residues to predicted ones using CEalign path
            pocket_residue_pairs = []
            if alignment_path_indices and len(alignment_path_indices) == 2 and \
                    len(alignment_path_indices[0]) > 0 and len(alignment_path_indices[1]) > 0:  # Check path validity
                exp_path_indices, pred_path_indices = alignment_path_indices

                # Create a set of already mapped experimental pocket auth_ids to avoid redundant processing
                # if multiple CAs of the same experimental residue map to the predicted structure.
                mapped_exp_pocket_ids_in_path = set()

                for exp_pocket_auth_id_current in experimental_pocket_auth_ids:
                    if exp_pocket_auth_id_current in mapped_exp_pocket_ids_in_path:
                        continue  # Already found a mapping for this exp pocket residue

                    for i in range(len(exp_path_indices)):
                        exp_ca_path_idx = exp_path_indices[i]
                        # Check bounds, as iloc can fail if index is out of bounds
                        if exp_ca_path_idx < len(exp_ca_for_align_df):
                            auth_id_of_exp_ca_in_path = exp_ca_for_align_df.iloc[exp_ca_path_idx]['auth_seq_id']

                            if auth_id_of_exp_ca_in_path == exp_pocket_auth_id_current:
                                pred_ca_path_idx = pred_path_indices[i]
                                if pred_ca_path_idx < len(pred_ca_for_align_df):
                                    auth_id_of_pred_ca_in_path = pred_ca_for_align_df.iloc[pred_ca_path_idx][
                                        'auth_seq_id']
                                    pocket_residue_pairs.append({
                                        'exp_auth_id': exp_pocket_auth_id_current,
                                        'pred_auth_id': auth_id_of_pred_ca_in_path
                                    })
                                    mapped_exp_pocket_ids_in_path.add(exp_pocket_auth_id_current)
                                    break  # Found mapping for this exp_pocket_auth_id_current
            else:
                print(f"[WARNING] CEalign did not return a valid or non-empty path for {exp_id} <-> {pred_id}.")

            # Ensure uniqueness of (exp_auth_id, pred_auth_id) pairs if one exp maps to multiple segments of pred
            # or vice-versa due to CEalign behavior (though less common for this simple mapping)
            if pocket_residue_pairs:
                unique_pairs_tuples = sorted(list(set(
                    (pair['exp_auth_id'], pair['pred_auth_id']) for pair in pocket_residue_pairs
                )))
                pocket_residue_pairs = [{'exp_auth_id': ep[0], 'pred_auth_id': ep[1]} for ep in unique_pairs_tuples]

            print(
                f"[INFO] Mapped {len(experimental_pocket_auth_ids)} exp pocket residues to {len(pocket_residue_pairs)} pred pocket residues via alignment path.")
        """


        # =========================
        # STEP 4: Align structures based on CA atoms  (two-pass with inlier filtering)
        # =========================
        try:
            # --- STEP 4A: initial CEalign on all protein CAs (non-RET) ---
            exp_ca_coords_for_cealign = exp_ca_for_align_df[['x', 'y', 'z']].astype(float).values
            pred_ca_coords_for_cealign = pred_ca_for_align_df[['x', 'y', 'z']].astype(float).values

            # Optional annotations (used for filtering / ordering)
            exp_ca_seqids = exp_ca_for_align_df[
                'auth_seq_id'].values if 'auth_seq_id' in exp_ca_for_align_df.columns else np.arange(
                len(exp_ca_for_align_df))
            pred_ca_seqids = pred_ca_for_align_df[
                'auth_seq_id'].values if 'auth_seq_id' in pred_ca_for_align_df.columns else np.arange(
                len(pred_ca_for_align_df))
            exp_ca_helix = exp_ca_for_align_df[
                'helix_num'].values if 'helix_num' in exp_ca_for_align_df.columns else None
            pred_ca_helix = pred_ca_for_align_df[
                'helix_num'].values if 'helix_num' in pred_ca_for_align_df.columns else None

            # Parameters for refinement
            INLIER_THRESH = 3.0  # Å, residual cutoff
            MIN_INLIERS = 12  # minimum pairs to proceed with refinement
            ENFORCE_SAME_HELIX = True  # only pair Hk↔Hk when helix_num present

            # CEalign pass 1
            R1, t1, alignment_path_indices_1, overall_ca_rmsd_1 = get_structure_alignment(
                exp_ca_coords_for_cealign, pred_ca_coords_for_cealign,
                window_size=8, max_gap=30
            )
            print(f"[INFO] Pass-1 (CEalign) CA RMSD: {overall_ca_rmsd_1:.3f} Å")

            # If path invalid/empty, fall back to single-pass logic
            if not (alignment_path_indices_1 and len(alignment_path_indices_1) == 2 and
                    len(alignment_path_indices_1[0]) > 0 and len(alignment_path_indices_1[1]) > 0):
                print(
                    f"[WARNING] CEalign (pass-1) returned empty/invalid path for {exp_id} <-> {pred_id}; using single-pass alignment.")
                R_use, t_use = R1, t1
                alignment_path_indices = alignment_path_indices_1
                overall_ca_rmsd = overall_ca_rmsd_1
            else:
                exp_idx_1, pred_idx_1 = alignment_path_indices_1

                # --- STEP 4B: analyze residuals on pass-1 path; build inlier set ---
                # IMPORTANT: keep residual computation consistent with the transform you apply later.
                # Your pipeline applies R,t to PRED to bring it into EXP frame, so compute residuals the same way:
                pred_aligned_p1 = np.dot(pred_ca_coords_for_cealign, R1) + t1  # shape (Npred, 3)

                # residual per matched pair
                residuals = []
                inlier_pairs = []
                for ii, jj in zip(exp_idx_1, pred_idx_1):
                    if ii >= len(exp_ca_coords_for_cealign) or jj >= len(pred_aligned_p1):
                        continue
                    d = np.linalg.norm(exp_ca_coords_for_cealign[ii] - pred_aligned_p1[jj])
                    # optional same-helix constraint when available
                    if ENFORCE_SAME_HELIX and (exp_ca_helix is not None) and (pred_ca_helix is not None):
                        if ii < len(exp_ca_helix) and jj < len(pred_ca_helix):
                            if not (exp_ca_helix[ii] == pred_ca_helix[jj]):
                                continue
                    residuals.append(d)
                    if d <= INLIER_THRESH:
                        inlier_pairs.append((ii, jj))

                n_total_for_coverage = max(len(exp_ca_coords_for_cealign), len(pred_ca_coords_for_cealign))
                coverage_inliers = (len(inlier_pairs) / float(n_total_for_coverage)) if n_total_for_coverage else 0.0
                print(
                    f"[INFO] Pass-1 inliers: {len(inlier_pairs)} (coverage {coverage_inliers:.2f}) with τ={INLIER_THRESH:.1f} Å; same-helix={ENFORCE_SAME_HELIX}")

                # If insufficient inliers, skip refinement and use pass-1
                if len(inlier_pairs) < MIN_INLIERS:
                    print(f"[WARN] Not enough inliers (<{MIN_INLIERS}) for refinement; using pass-1 result.")
                    R_use, t_use = R1, t1
                    alignment_path_indices = alignment_path_indices_1
                    overall_ca_rmsd = overall_ca_rmsd_1
                else:
                    # Build filtered, CONTIGUOUS arrays (sort by helix then seqid to keep CE DP happy)
                    def sort_key(pair):
                        ii, jj = pair
                        h = exp_ca_helix[ii] if exp_ca_helix is not None else 0
                        s = exp_ca_seqids[ii] if exp_ca_seqids is not None else ii
                        return (h, s)

                    inlier_pairs_sorted = sorted(inlier_pairs, key=sort_key)
                    exp_inlier_idx = np.array([ii for ii, _ in inlier_pairs_sorted], dtype=int)
                    pred_inlier_idx = np.array([jj for _, jj in inlier_pairs_sorted], dtype=int)

                    exp_coords_filt = exp_ca_coords_for_cealign[exp_inlier_idx]
                    pred_coords_filt = pred_ca_coords_for_cealign[pred_inlier_idx]

                    # --- STEP 4C: CEalign pass-2 on filtered inputs ---
                    try:
                        R2, t2, alignment_path_indices_2, overall_ca_rmsd_2 = get_structure_alignment(
                            exp_coords_filt, pred_coords_filt,
                            window_size=8, max_gap=30
                        )
                        print(
                            f"[INFO] Pass-2 (refined) CA RMSD: {overall_ca_rmsd_2:.3f} Å on {len(exp_coords_filt)} filtered CAs")

                        # Map refined path back to ORIGINAL exp/pred CA indices
                        if alignment_path_indices_2 and len(alignment_path_indices_2) == 2 and \
                                len(alignment_path_indices_2[0]) > 0 and len(alignment_path_indices_2[1]) > 0:
                            exp_idx_2_rel, pred_idx_2_rel = alignment_path_indices_2
                            # indices into filtered arrays -> map back to global
                            exp_idx_2_global = [int(exp_inlier_idx[k]) for k in exp_idx_2_rel]
                            pred_idx_2_global = [int(pred_inlier_idx[k]) for k in pred_idx_2_rel]
                            alignment_path_indices = (exp_idx_2_global, pred_idx_2_global)
                            R_use, t_use = R2, t2
                            overall_ca_rmsd = overall_ca_rmsd_2
                        else:
                            print(f"[WARN] Pass-2 returned empty/invalid path; falling back to pass-1.")
                            R_use, t_use = R1, t1
                            alignment_path_indices = alignment_path_indices_1
                            overall_ca_rmsd = overall_ca_rmsd_1
                    except Exception as e_ref:
                        print(f"[WARN] Pass-2 CEalign failed for {exp_id} <-> {pred_id}: {e_ref}; using pass-1 result.")
                        R_use, t_use = R1, t1
                        alignment_path_indices = alignment_path_indices_1
                        overall_ca_rmsd = overall_ca_rmsd_1

            # From here onward, use (R_use, t_use) and the chosen alignment_path_indices
            print(
                f"[INFO] Using CA RMSD = {overall_ca_rmsd:.3f} Å (two-pass {'refined' if (alignment_path_indices != alignment_path_indices_1) else 'initial'})")

            # Apply FINAL transform ONCE to the pocket dataframe (as in your original logic)
            pred_all_atom_coords_original_for_pocket = pred_df_for_pocket_alignment[['x', 'y', 'z']].values.astype(
                float)
            pred_df_for_pocket_alignment[['x', 'y', 'z']] = np.dot(pred_all_atom_coords_original_for_pocket,
                                                                   R_use) + t_use

            # --- Map experimental pocket residues to predicted ones using the (possibly refined) CE path ---
            pocket_residue_pairs = []
            if alignment_path_indices and len(alignment_path_indices) == 2 and \
                    len(alignment_path_indices[0]) > 0 and len(alignment_path_indices[1]) > 0:
                exp_path_indices, pred_path_indices = alignment_path_indices

                mapped_exp_pocket_ids_in_path = set()
                for exp_pocket_auth_id_current in experimental_pocket_auth_ids:
                    if exp_pocket_auth_id_current in mapped_exp_pocket_ids_in_path:
                        continue
                    # scan the CE path for this exp residue (by auth_seq_id)
                    for k in range(len(exp_path_indices)):
                        exp_ca_path_idx = exp_path_indices[k]
                        if exp_ca_path_idx < len(exp_ca_for_align_df):
                            auth_id_of_exp_ca_in_path = exp_ca_for_align_df.iloc[exp_ca_path_idx]['auth_seq_id']
                            if auth_id_of_exp_ca_in_path == exp_pocket_auth_id_current:
                                pred_ca_path_idx = pred_path_indices[k]
                                if pred_ca_path_idx < len(pred_ca_for_align_df):
                                    auth_id_of_pred_ca_in_path = pred_ca_for_align_df.iloc[pred_ca_path_idx][
                                        'auth_seq_id']
                                    pocket_residue_pairs.append({
                                        'exp_auth_id': exp_pocket_auth_id_current,
                                        'pred_auth_id': auth_id_of_pred_ca_in_path
                                    })
                                    mapped_exp_pocket_ids_in_path.add(exp_pocket_auth_id_current)
                                    break
            else:
                print(
                    f"[WARNING] CEalign did not return a valid or non-empty path for {exp_id} <-> {pred_id} (after two-pass).")

            # Uniquify residue pairs
            if pocket_residue_pairs:
                unique_pairs_tuples = sorted(list(set(
                    (pair['exp_auth_id'], pair['pred_auth_id']) for pair in pocket_residue_pairs
                )))
                pocket_residue_pairs = [{'exp_auth_id': ep[0], 'pred_auth_id': ep[1]} for ep in unique_pairs_tuples]

            print(f"[INFO] Mapped {len(experimental_pocket_auth_ids)} exp-pocket residues to "
                  f"{len(pocket_residue_pairs)} pred-pocket residues via CE path (two-pass).")

        except RuntimeError as re:
            print(f"[ERROR] Alignment runtime error for {exp_id} <-> {pred_id}: {str(re)}")
            results[exp_id] = {'error': f'Alignment failed: {str(re)}'}
            continue
        except Exception as e:
            import traceback
            print(f"[ERROR] Failed during alignment step for {exp_id} and {pred_id}: {str(e)}")
            traceback.print_exc()
            results[exp_id] = {'error': f'Processing failed (alignment): {str(e)}'}

            # STEP 5: Calculate binding pocket RMSD using mapped pairs and transformed pred_df_for_pocket_alignment
            pocket_rmsd_sum_sq = 0  # Changed from pocket_rmsd_sum to reflect squared sum
            pocket_atom_count = 0
            per_residue_rmsd = {}

            if not pocket_residue_pairs:
                print(f"[WARNING] No pocket residue pairs found after alignment path mapping for {exp_id}.")

            for pair in pocket_residue_pairs:
                exp_res_id = pair['exp_auth_id']
                pred_res_id = pair['pred_auth_id']

                exp_res_df_current = exp_df[exp_df['auth_seq_id'] == exp_res_id]
                # pred_res_df_current is from the ALREADY GLOBALLY ALIGNED pred_df_for_pocket_alignment
                pred_res_df_current = pred_df_for_pocket_alignment[
                    pred_df_for_pocket_alignment['auth_seq_id'] == pred_res_id]

                if exp_res_df_current.empty or pred_res_df_current.empty:
                    continue

                current_res_matched_atoms = []  # Use for per-residue RMSD calculation
                # Match atoms by name within this STRUCTURALLY ALIGNED pair of residues
                # This loop is from your original working binding pocket RMSD logic
                for atom_name in exp_res_df_current['res_atom_name'].unique():
                    exp_atom_df = exp_res_df_current[exp_res_df_current['res_atom_name'] == atom_name]
                    pred_atom_df = pred_res_df_current[pred_res_df_current['res_atom_name'] == atom_name]

                    if not exp_atom_df.empty and not pred_atom_df.empty:
                        exp_coords_atom = exp_atom_df[['x', 'y', 'z']].astype(float).values[0]
                        # These are from the globally aligned pred_df_for_pocket_alignment
                        pred_coords_atom_for_comparison = pred_atom_df[['x', 'y', 'z']].astype(float).values[0]

                        dist = np.sqrt(np.sum((exp_coords_atom - pred_coords_atom_for_comparison) ** 2))
                        current_res_matched_atoms.append(
                            (atom_name, dist, exp_coords_atom, pred_coords_atom_for_comparison))

                        pocket_rmsd_sum_sq += dist ** 2  # Sum of squared distances
                        pocket_atom_count += 1

                if current_res_matched_atoms:
                    atom_errors = {}
                    for atom_name_loop, dist_loop, exp_coords_loop, aligned_pred_coords_loop in current_res_matched_atoms:
                        atom_errors[atom_name_loop] = {
                            'distance': dist_loop,
                            'error_x': abs(exp_coords_loop[0] - aligned_pred_coords_loop[0]),
                            'error_y': abs(exp_coords_loop[1] - aligned_pred_coords_loop[1]),
                            'error_z': abs(exp_coords_loop[2] - aligned_pred_coords_loop[2])
                        }
                    res_rmsd_val = np.sqrt(np.mean([d ** 2 for _, d, _, _ in current_res_matched_atoms]))
                    per_residue_rmsd[f"{exp_res_id}(exp)-{pred_res_id}(pred)"] = {
                        'rmsd': res_rmsd_val,
                        'res_type_exp': exp_res_df_current['res_name3l'].iloc[0],
                        'res_type_pred': pred_res_df_current['res_name3l'].iloc[0],  # if different due to path
                        'atom_count': len(current_res_matched_atoms),
                        'atom_errors': atom_errors
                    }

            binding_pocket_rmsd = np.sqrt(pocket_rmsd_sum_sq / pocket_atom_count) if pocket_atom_count > 0 else np.nan
            print(
                f"[INFO] Binding pocket RMSD (path aligned): {binding_pocket_rmsd:.3f}Å (based on {pocket_atom_count} atoms)")

            #############################################################
            # STEP 6: Calculate retinal RMSD (Reverted to original working logic)
            #############################################################

            # Experimental retinal coordinates (from original exp_df via exp_ret_df_step1)
            exp_ret_coords = exp_ret_df_step1[['x', 'y', 'z']].astype(float).values

            # Predicted retinal coordinates (from original UNTRANSFORMED pred_df_original_for_retinal_rmsd via pred_ret_df_step1)
            pred_ret_coords_original_untransformed = pred_ret_df_step1[['x', 'y', 'z']].astype(float).values

            retinal_rmsd = np.nan  # Initialize
            if exp_ret_coords.shape[0] > 0 and pred_ret_coords_original_untransformed.shape[0] > 0:
                # Apply the global R and t (from CA alignment) to these specific original predicted retinal coordinates
                aligned_pred_ret_coords_for_rmsd = np.dot(pred_ret_coords_original_untransformed, R_use) + t_use

                try:
                    retinal_rmsd = compute_retinal_mean_closest_distance(
                        exp_ret_coords,
                        aligned_pred_ret_coords_for_rmsd  # Use the explicitly transformed retinal coords
                    )
                except Exception as e_cr_rmsd:
                    print(f"[ERROR] compute_retinal_mean_closest_distance failed for {exp_id}: {e_cr_rmsd}")
                    retinal_rmsd = np.nan
            else:
                print(
                    f"[WARNING] Could not compute retinal RMSD for {exp_id} due to empty coordinate arrays for retinal.")
            print(f"[INFO] Retinal RMSD: {retinal_rmsd:.3f}Å")

            #############################################################
            # STEP 7: Store all results
            #############################################################
            results[exp_id] = {
                'backbone_rmsd': overall_ca_rmsd,
                'overall_pocket_rmsd': binding_pocket_rmsd,
                'retinal_rmsd': retinal_rmsd,
                'per_residue_rmsd': per_residue_rmsd,
                'experimental_pocket_auth_ids': experimental_pocket_auth_ids,
                'mapped_pocket_residue_pairs': pocket_residue_pairs,
                'alignment': {
                    'rotation': R_use.tolist(),
                    'translation': t_use.tolist(),
                    'cealign_path_indices': {
                        'exp_indices': alignment_path_indices[0].tolist() if alignment_path_indices and isinstance(
                            alignment_path_indices[0], np.ndarray) else (
                            list(alignment_path_indices[0]) if alignment_path_indices and alignment_path_indices[
                                0] is not None else []),
                        'pred_indices': alignment_path_indices[1].tolist() if alignment_path_indices and isinstance(
                            alignment_path_indices[1], np.ndarray) else (
                            list(alignment_path_indices[1]) if alignment_path_indices and alignment_path_indices[
                                1] is not None else [])
                    } if alignment_path_indices and len(alignment_path_indices) == 2 else None
                }
            }

        except RuntimeError as re:
            print(f"[ERROR] Alignment runtime error for {exp_id} <-> {pred_id}: {str(re)}")
            results[exp_id] = {'error': f'Alignment failed: {str(re)}'}
        except Exception as e:
            import traceback
            print(f"[ERROR] Failed to process {exp_id} and {pred_id}: {str(e)}")
            traceback.print_exc()
            results[exp_id] = {'error': f'Processing failed: {str(e)}'}

    print(f"[INFO] Completed binding pocket RMSD calculation for {len(results)} structure pairs")
    return results


def create_exp_mapping(cp_mo_exp, cp_hide_exp):
    """
    Creates a 1:1 mapping between structures in cp_mo_exp and cp_hide_exp.
    
    Args:
        cp_mo_exp: CifBaseProcessor containing MO experimental structures
        cp_hide_exp: CifBaseProcessor containing Hideaki experimental structures
        
    Returns:
        dict: Dictionary mapping MO structure IDs to Hideaki structure IDs
    """
    mapping = {}
    
    # Get unique PDB IDs from each processor
    mo_exp_ids = cp_mo_exp.pdb_ids
    hide_exp_ids = cp_hide_exp.pdb_ids
    
    # For each MO structure, find the matching Hideaki structure
    for mo_id in mo_exp_ids:
        # Create the expected Hideaki ID by removing '_model_0' if present
        if mo_id.endswith('_model_0'):
            base_id = mo_id[:-8]  # Remove '_model_0'
        else:
            base_id = mo_id
            
        # Check if the base ID exists in Hideaki structures
        for hide_id in hide_exp_ids:
            # Remove '_model_0' from Hideaki ID if present for comparison
            if hide_id.endswith('_model_0'):
                hide_base_id = hide_id[:-8]
            else:
                hide_base_id = hide_id
                
            # If we have a match, add to mapping
            if base_id == hide_base_id:
                mapping[mo_id] = hide_id
                break
    
    return mapping


def create_exp_mapping_with_display_name(cp_mo_exp, cp_hide_exp, mo_property_file='/mnt/c/Users/hidbe/PycharmProjects/phd/projects/opsin_analysis/property/mo_exp.csv'):
    """
    Creates a 1:1 mapping between structures in cp_mo_exp and cp_hide_exp 
    using 'PDB ID' from cp_hide_exp and 'display_name + _model_0' from cp_mo_exp.
    
    Args:
        cp_mo_exp: CifBaseProcessor containing MO experimental structures
        cp_hide_exp: CifBaseProcessor containing Hideaki experimental structures
        mo_property_file: Path to the mo_exp.csv property file containing display_name data
        
    Returns:
        dict: Dictionary mapping 'PDB ID' from cp_hide_exp to corresponding 'display_name + _model_0' from cp_mo_exp
    """
    mapping = {}
    
    # Load property data from CSV file
    try:
        import pandas as pd
        import os
        
        # Check if property file exists
        if not os.path.exists(mo_property_file):
            print(f"[WARNING] Property file not found: {mo_property_file}")
            # Fall back to extracting base names from PDB IDs
            return create_exp_mapping(cp_mo_exp, cp_hide_exp)
        
        # Load the property data
        print(f"[INFO] Loading property data from {mo_property_file}")
        prop_df = pd.read_csv(mo_property_file)
        print(f"[INFO] Loaded property data with {len(prop_df)} entries")
        
        # Create dictionaries for quick lookup
        pdb_to_display = {}
        display_to_pdb = {}
        
        # Process property data to extract PDB ID and display_name mappings
        for _, row in prop_df.iterrows():
            pdb_id = row.get('PDB ID', None)
            display_name = row.get('display_name', None)
            
            # Skip if either value is missing
            if pd.isna(pdb_id) or pd.isna(display_name):
                continue
                
            # Clean up values
            pdb_id = str(pdb_id).strip()
            display_name = str(display_name).strip()
            
            # Store mappings
            if pdb_id and display_name:
                pdb_to_display[pdb_id] = display_name
                display_to_pdb[display_name] = pdb_id
        
        print(f"[INFO] Found {len(pdb_to_display)} PDB ID to display_name mappings")
        
        # Create the mapping from Hideaki experimental to MO experimental
        hide_exp_ids = cp_hide_exp.pdb_ids
        mo_exp_ids = cp_mo_exp.pdb_ids
        
        # For each Hideaki structure
        for hide_id in hide_exp_ids:
            # Remove '_model_0' if present for base comparison
            if hide_id.endswith('_model_0'):
                hide_base = hide_id[:-8]
            else:
                hide_base = hide_id
                
            # Try exact match first
            if hide_base in pdb_to_display:
                display_name = pdb_to_display[hide_base]
                mo_id = f"{display_name}_model_0"
                
                # Verify that this ID exists in mo_exp
                for mo_exp_id in mo_exp_ids:
                    if mo_exp_id == mo_id or mo_exp_id.startswith(display_name):
                        mapping[hide_id] = mo_id
                        print(f"[INFO] Mapped {hide_id} to {mo_id}")
                        break
            
            # If no match, try by similarity
            if hide_id not in mapping:
                best_match = None
                best_score = 0
                
                for display_name in display_to_pdb.keys():
                    # Simple string similarity score
                    similarity = sum(c1 == c2 for c1, c2 in zip(hide_base.lower(), display_name.lower())) / max(len(hide_base), len(display_name))
                    
                    if similarity > best_score and similarity > 0.7:  # 70% match threshold
                        best_score = similarity
                        best_match = display_name
                
                if best_match:
                    mo_id = f"{best_match}_model_0"
                    # Verify that this ID exists or is similar to IDs in mo_exp
                    for mo_exp_id in mo_exp_ids:
                        if mo_exp_id == mo_id or mo_exp_id.startswith(best_match):
                            mapping[hide_id] = mo_id
                            print(f"[INFO] Mapped {hide_id} to {mo_id} (similarity: {best_score:.2f})")
                            break
        
        print(f"[INFO] Created {len(mapping)} mappings between Hideaki and MO structures")
        return mapping
        
    except Exception as e:
        print(f"[ERROR] Error processing property data: {e}")
        import traceback
        traceback.print_exc()
        
        # Fall back to base name matching
        print("[INFO] Falling back to simple name-based mapping")
        return create_exp_mapping(cp_mo_exp, cp_hide_exp)


def create_mapping_for_rmsd_calculation(cp_mo_exp, cp_mo_pred, cp_hide_exp, cp_hide_pred):
    """
    Creates a comprehensive mapping dictionary for RMSD calculations between
    experimental and predicted structures, including both Hideaki and MO structures.
    
    Args:
        cp_mo_exp: CifBaseProcessor containing MO experimental structures
        cp_mo_pred: CifBaseProcessor containing MO predicted structures
        cp_hide_exp: CifBaseProcessor containing Hideaki experimental structures
        cp_hide_pred: CifBaseProcessor containing Hideaki predicted structures
        
    Returns:
        dict: Mapping dictionary for RMSD calculations
    """
    mapping = {}
    
    # 1. Process Hideaki experimental-predicted pairs
    print("[INFO] Creating mapping for Hideaki experimental-predicted pairs...")
    hideaki_mapping = {}
    for exp_id in cp_hide_exp.pdb_ids:
        # Standard pattern: remove '_model_0' suffix to get predicted ID
        if exp_id.endswith('_model_0'):
            base_id = exp_id[:-8]
        else:
            base_id = exp_id
            
        # Check if predicted version exists
        for pred_id in cp_hide_pred.pdb_ids:
            if pred_id.endswith('_model_0') and pred_id[:-8] == base_id:
                hideaki_mapping[exp_id] = {"predicted": pred_id}
                break
            elif pred_id == f"{base_id}_model_0":
                hideaki_mapping[exp_id] = {"predicted": pred_id}
                break
    
    print(f"[INFO] Found {len(hideaki_mapping)} Hideaki experimental-predicted pairs")
    mapping.update(hideaki_mapping)
    
    # 2. Process MO experimental-predicted pairs
    print("[INFO] Creating mapping for MO experimental-predicted pairs...")
    mo_mapping = {}
    for exp_id in cp_mo_exp.pdb_ids:
        # Check for display_name + _model_0 pattern
        base_id = exp_id
        if exp_id.endswith('_model_0'):
            base_id = exp_id[:-8]
            
        # Try different patterns for predicted structures
        for pred_id in cp_mo_pred.pdb_ids:
            # Common patterns: base_smile_model_0, base_model_0
            if (pred_id.startswith(f"{base_id}_smile") or 
                pred_id.startswith(base_id) and '_model_0' in pred_id):
                mo_mapping[exp_id] = {"predicted": pred_id}
                break
    
    print(f"[INFO] Found {len(mo_mapping)} MO experimental-predicted pairs")
    mapping.update(mo_mapping)
    
    # 3. Create cross-mappings between Hideaki and MO structures
    print("[INFO] Creating cross-mappings between Hideaki and MO structures...")
    cross_mapping = create_exp_mapping_with_display_name(cp_mo_exp, cp_hide_exp)
    
    # Convert cross-mapping to the format expected by RMSD calculation
    cross_mapping_formatted = {}
    for hide_id, mo_id in cross_mapping.items():
        cross_mapping_formatted[hide_id] = {"predicted": mo_id}
    
    print(f"[INFO] Found {len(cross_mapping_formatted)} cross-mappings between Hideaki and MO structures")
    mapping.update(cross_mapping_formatted)
    
    print(f"[INFO] Created total of {len(mapping)} structure mappings for RMSD calculation")
    return mapping


def compare_structures(data_dict, output_dir='outputs', visualize=True):
    """
    Step 5: Structure comparison

    This function calculates RMSD between pairs of unique structures and creates
    a similarity matrix without using PropertyProcessor.

    It excludes predicted structures that have experimental counterparts to avoid
    redundancy and focus on unique structural information.

    Args:
        data_dict: Dictionary with data from previous steps
        output_dir: Directory to save outputs files
        visualize: Whether to generate visualizations

    Returns:
        Dictionary with RMSD data
    """
    print("\n=== Step 5: Structure Comparison ===")
    processed_structures_complete = data_dict['processed_structures']
    structure_mapping = data_dict.get('structure_mapping', {})

    # Check if we have structures to compare
    if not processed_structures_complete:
        print("No structures available for comparison. Skipping structure comparison step.")
        return {
            'processed_structures': processed_structures_complete,
            'rmsd_df': pd.DataFrame(),
            'rmsd_matrix': np.array([]),
            'pdb_list': [],
            'group_dict': {},
            'name_dict': {}
        }

    # Make sure we have at least 2 structures to compare
    if len(processed_structures_complete) < 2:
        print(f"Only {len(processed_structures_complete)} structure(s) available. Need at least 2 for comparison.")
        pdb_id = list(processed_structures_complete.keys())[0] if processed_structures_complete else "none"

        # Return minimal data for single structure
        return {
            'processed_structures': processed_structures_complete,
            'rmsd_df': pd.DataFrame(index=[pdb_id], columns=[pdb_id], data=[[0.0]]),
            'rmsd_matrix': np.array([[0.0]]),
            'pdb_list': [pdb_id],
            'group_dict': {pdb_id: 'Unknown'},
            'name_dict': {pdb_id: pdb_id}
        }

    # Filter out predicted structures that have experimental counterparts
    # using the structure mapping
    unique_structures = processed_structures_complete.copy()
    excluded_structures = set()

    # Identify predicted structures to exclude based on the mapping
    for exp_id, pred_id in structure_mapping.items():
        # Handle both old and new mapping formats
        if isinstance(pred_id, dict) and 'predicted' in pred_id:
            pred_id = pred_id['predicted']

        # Skip if either structure is missing
        if exp_id not in processed_structures_complete or pred_id not in processed_structures_complete:
            continue

        # If both experimental and predicted structures exist, exclude the predicted one
        excluded_structures.add(pred_id)

    # Remove excluded structures from the unique set
    for struct_id in excluded_structures:
        if struct_id in unique_structures:
            del unique_structures[struct_id]

    print(f"Found {len(processed_structures_complete)} total structures")
    print(f"Excluded {len(excluded_structures)} predicted structures that have experimental counterparts")
    print(f"Computing RMSD matrix for {len(unique_structures)} unique structures using C-alpha atoms...")

    # Compute improved RMSD matrix using C-alpha atoms only for unique structures
    # Use helix residues (1-7) for alignment to focus on the transmembrane domains
    cache_dir = os.path.join(output_dir, 'cache')
    rmsd_df, rmsd_matrix, pdb_list, alignment_paths = compute_all_vs_all_rmsd_improved(
        unique_structures,
        subset='CA',  # Explicitly specify C-alpha atoms for RMSD calculation
        chain_id='A',  # Use chain A for comparison
        tm_score_threshold=0.0,  # Include all alignments regardless of TM-score
        verbose=True,  # Print main alignment results for monitoring progress
        use_helix_only=True,  # Use only helix residues (1-7) for alignment
        cache_dir=cache_dir,  # Directory to cache RMSD results
        force_recompute=False,  # Use cached results if available
    )

    # Save RMSD matrix
    rmsd_df.to_csv(os.path.join(output_dir, 'rmsd_matrix.csv'))
    print(f"Saved RMSD matrix to {os.path.join(output_dir, 'rmsd_matrix.csv')}")

    # Create group dictionary for visualization
    group_dict = {}
    name_dict = {}

    # Create groups based on metadata in processed structures
    for pdb_id in pdb_list:
        if pdb_id in processed_structures_complete:
            struct_data = processed_structures_complete[pdb_id]
            # Use property data from the loaded structure data
            if 'properties' in struct_data:
                print("available structure properties:", struct_data['properties'].keys())
                # First try to use molecular_function as the group type
                if 'molecular_function' in struct_data['properties'] and struct_data['properties'][
                    'molecular_function'] != 'Unknown':
                    group_dict[pdb_id] = struct_data['properties']['molecular_function']
                # If molecular_function is not available, fall back to domain
                elif 'domain' in struct_data['properties'] and struct_data['properties']['domain'] != 'Unknown':
                    group_dict[pdb_id] = struct_data['properties']['domain']
                else:
                    # Fallback to a default group
                    group_dict[pdb_id] = 'Unknown'
            else:
                # Fallback to a default group if no properties are available
                group_dict[pdb_id] = 'Unknown'

    # Visualize RMSD matrix only if we have structures and visualize is True
    if visualize and len(pdb_list) >= 2:
        try:
            fig = visualize_rmsd_heatmap(
                rmsd_df,
                structure_ids=pdb_list,
                group_dict=group_dict,
                name_dict=name_dict,
                group_by='molecular_function'
            )
            plt.savefig(os.path.join(output_dir, 'rmsd_matrix.png'), dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved RMSD matrix visualization to {os.path.join(output_dir, 'rmsd_matrix.png')}")
        except Exception as e:
            print(f"Warning: Could not create RMSD matrix visualization: {e}")

        # Also create a similarity tree visualization
        try:
            fig = create_and_visualize_similarity_tree(
                rmsd_data=rmsd_df,  # Pass the DataFrame directly
                group_dict=group_dict,
                name_dict=name_dict,
                group_by='molecular_function'
            )
            plt.savefig(os.path.join(output_dir, 'similarity_tree.png'), dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved similarity tree visualization to {os.path.join(output_dir, 'similarity_tree.png')}")
        except Exception as e:
            print(f"Warning: Could not create similarity tree: {e}")
            import traceback
            traceback.print_exc()

    return {
        'processed_structures': processed_structures_complete,
        'alignment_paths': alignment_paths,
        'rmsd_df': rmsd_df,
        'rmsd_matrix': rmsd_matrix,
        'pdb_list': pdb_list,
        'group_dict': group_dict,
        'name_dict': name_dict
    }


def create_unified_structure_mapping(data, property_data=None):
    """
    Create a unified structure mapping for all opsin structures

    This function creates a single unified mapping that will be used throughout the workflow.
    For MO structures (cp_mo_pred/exp), it uses mappings from property_data if available.
    For Hideaki structures (cp_hide_exp/pred), it always uses name-based mappings.

    Args:
        data: Dictionary containing dataset information and processors
        property_data: Optional dictionary with property data including structure mappings

    Returns:
        Dictionary mapping experimental structure IDs to predicted structure IDs
    """
    print("\n=== Creating Unified Structure Mapping ===")

    # Extract needed processors
    cp_mo_exp = data['cp_mo_exp']
    cp_mo_pred = data['cp_mo_pred']
    cp_hide_exp = data['cp_hide_exp']
    cp_hide_pred = data['cp_hide_pred']

    # Structure mapping will be our single source of truth
    structure_mapping = {}

    # 1. Handle MO structures using property data mapping if available
    mo_mapping = {}
    if property_data and 'structure_mapping' in property_data:
        provided_structure_mapping = property_data['structure_mapping']
        print("Using structure mapping from property data for MO structures...")
        print(f"Property data has {len(provided_structure_mapping)} total mappings")

        # Debug information to understand what's in the provided mapping
        print("First 5 mappings from property data:")
        for i, (exp_id, pred_id) in enumerate(list(provided_structure_mapping.items())[:5]):
            print(f"  - {exp_id} -> {pred_id}")

        # Debug processor IDs
        print(f"MO exp processor has {len(cp_mo_exp.pdb_ids)} structure IDs")
        print(f"MO pred processor has {len(cp_mo_pred.pdb_ids)} structure IDs")

        # Filter the mapping to only include MO structures
        # Print debug info about the PDB IDs
        property_pdb_ids = set(provided_structure_mapping.keys())

        # Print some stats about the mappings
        print(f"\nMO structure debugging:")
        print(f"Total PDB IDs in property data: {len(property_pdb_ids)}")
        print(f"PDB IDs in MO exp processor: {len(cp_mo_exp.pdb_ids)}")
        print(f"PDB IDs in MO pred processor: {len(cp_mo_pred.pdb_ids)}")

        # Find PDB IDs in dataset that are not in property data
        missing_from_property = set(cp_mo_exp.pdb_ids) - property_pdb_ids
        if missing_from_property:
            print(f"Found {len(missing_from_property)} experimental structures without property data:")
            for i, missing_id in enumerate(sorted(missing_from_property)):
                if i < 10:  # Just show a few
                    print(f"  - {missing_id}")
                elif i == 10:
                    print(f"  - ... and {len(missing_from_property) - 10} more")

        # Apply filtering to include only valid mappings
        for exp_id, pred_id in provided_structure_mapping.items():
            # Only include mappings where the experimental ID exists in our dataset
            if exp_id in cp_mo_exp.pdb_ids:
                # Check if the predicted ID exists directly
                if pred_id in cp_mo_pred.pdb_ids:
                    mo_mapping[exp_id] = pred_id
                else:
                    # If the exact predicted ID is not found, try adapting it
                    corrected_pred_id = None

                    # Fix common naming issues
                    # 1. If the ID contains dashes, replace them with underscores
                    corrected_name = pred_id.replace('-', '_')
                    if corrected_name in cp_mo_pred.pdb_ids:
                        corrected_pred_id = corrected_name

                    # 2. Try other case variations (some might be camelCase vs snake_case)
                    potential_matches = []
                    base_id = pred_id.replace('_smile_model_0', '')

                    # Be careful with short names that might be substrings of others
                    # e.g., 'Mac' vs 'MacR' - we want exact matches when possible
                    for p_id in cp_mo_pred.pdb_ids:
                        p_base = p_id.replace('_smile_model_0', '')

                        # Exact match (case insensitive)
                        if base_id.lower() == p_base.lower():
                            potential_matches = [p_id]  # This is definitely the right match
                            break

                        # Base is a substring but make sure it's not too generic
                        # (only consider substring matches for names with at least 4 characters)
                        elif len(base_id) >= 4 and base_id.lower() in p_id.lower():
                            potential_matches.append(p_id)

                    if potential_matches:
                        # Take the shortest matching ID as it's likely the most specific
                        corrected_pred_id = min(potential_matches, key=len)

                    if corrected_pred_id:
                        print(f"  - Fixed mapping: {exp_id} → {pred_id} → {corrected_pred_id}")
                        mo_mapping[exp_id] = corrected_pred_id
                    else:
                        print(f"  - Found exp ID {exp_id} but pred ID {pred_id} is missing from dataset")

        # For any experimental IDs without mappings, try a simple pattern match
        # This handles any PDB IDs that weren't included in the property file
        unmapped_exp_ids = set(cp_mo_exp.pdb_ids) - set(mo_mapping.keys())
        if unmapped_exp_ids:
            print(f"Attempting simple pattern matching for {len(unmapped_exp_ids)} unmapped experimental structures:")

            for exp_id in unmapped_exp_ids:
                # MO structures follow a simple pattern: 4-character PDB ID -> display_name_smile_model_0
                # Try to find a matching predicted structure with similar name
                for pred_id in cp_mo_pred.pdb_ids:
                    # If the experimental ID is a substring of the predicted ID
                    # (e.g., 4KLY in BPR_smile_model_0)
                    if exp_id.lower() in pred_id.lower():
                        mo_mapping[exp_id] = pred_id
                        print(f"  - Pattern match: {exp_id} -> {pred_id}")
                        break

        print(f"Found {len(mo_mapping)} MO structure pairs from property data")
    else:
        print("No property data provided for MO structures")
        print("WARNING: MO structures will not be mapped without property data")

    # 2. Handle Hideaki structures using name-based mapping
    print("Creating Hideaki structure mapping using filename patterns...")
    hideaki_mapping = {}
    for exp_id in cp_hide_exp.pdb_ids:
        # Standard pattern with model_0 as specified
        model_pred_id = f"{exp_id}_model_0"
        if model_pred_id in cp_hide_pred.pdb_ids:
            hideaki_mapping[exp_id] = model_pred_id
            continue

        # Fallback to other patterns if needed
        pred_id = f"{exp_id}_pred"
        if pred_id in cp_hide_pred.pdb_ids:
            hideaki_mapping[exp_id] = pred_id
            continue

        # Try looking at base name parts as a last resort
        exp_parts = exp_id.split('_')
        if len(exp_parts) >= 1:
            base_part = exp_parts[0]
            for pred_candidate in cp_hide_pred.pdb_ids:
                if base_part in pred_candidate and (
                        "_model_0" in pred_candidate or
                        "_pred" in pred_candidate or
                        "_smile" in pred_candidate
                ):
                    hideaki_mapping[exp_id] = pred_candidate
                    break

    # Log the mapping results for debugging
    print(f"Created {len(hideaki_mapping)} Hideaki structure pairs")
    if hideaki_mapping:
        print("Example Hideaki pairs:")
        for i, (exp_id, pred_id) in enumerate(hideaki_mapping.items()):
            if i < 3:  # Show up to 3 examples
                print(f"  - {exp_id} -> {pred_id}")

    # 3. Combine all mappings into a single unified mapping
    # MO mappings take precedence over any potential overlapping Hideaki mappings
    structure_mapping = {**hideaki_mapping, **mo_mapping}

    print(f"Created unified structure mapping with {len(structure_mapping)} experimental-predicted pairs")
    print(f"  - {len(hideaki_mapping)} Hideaki pairs based on filename patterns")
    print(f"  - {len(mo_mapping)} MO pairs from property data")

    return structure_mapping
