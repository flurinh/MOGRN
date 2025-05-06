"""
Functions for comparing and aligning structures.
These functions handle the calculation of RMSD between structure pairs
and comparison of binding pockets.
"""
import os
import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from protos.processing.structure.struct_alignment import get_structure_alignment
from Bio.PDB.qcprot import QCPSuperimposer
from tqdm import tqdm

# Import error_analysis module from the local package
from projects.opsin_analysis.common_utils import compute_retinal_mean_closest_distance

from projects.opsin_analysis.visualization_functions import create_and_visualize_similarity_tree, visualize_rmsd_heatmap

import matplotlib.pyplot as plt


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
        
# Define a fallback alignment function that doesn't use CEalign (which can cause segmentation faults)
def fallback_alignment(coords1, coords2):
    """
    A simple fallback alignment method that doesn't rely on CEalign.
    Uses QCPSuperimposer directly with no path finding, just direct superposition.
    
    Args:
        coords1: Numpy array of coordinates for the first structure (shape: n x 3)
        coords2: Numpy array of coordinates for the second structure (shape: m x 3)
        
    Returns:
        rotation matrix, translation vector, None (no path), RMSD
    """
    # Ensure we have numpy arrays
    coords1 = np.array(coords1)
    coords2 = np.array(coords2)
    
    # Get minimum number of points (we can only align as many points as are in the smaller set)
    min_points = min(len(coords1), len(coords2))
    
    # If structures have different numbers of points, use the first min_points from each
    if min_points < len(coords1):
        coords1 = coords1[:min_points]
    if min_points < len(coords2):
        coords2 = coords2[:min_points]
    
    # Create QCP superimposer instance
    qcp = QCPSuperimposer()
    
    # Set the coordinates
    qcp.set(coords1, coords2)
    
    # Perform the superposition
    qcp.run()
    
    # Return the result
    return qcp.rot, qcp.tran, None, qcp.rms

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
        verbose: Whether to print verbose output during alignment
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
    import os
    import pickle
    import hashlib
    
    # Prepare the list of structures
    structure_ids = list(structures.keys())
    n_structures = len(structure_ids)
    
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
            if (struct1_id != '2L6X') and (struct2_id != '2L6X'):
                if i == j:
                    rmsd_matrix[i, j] = 0.0
                    continue

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
                    if subset == 'CA':
                        struct1_df = struct1_df[struct1_df['res_atom_name'] == 'CA']
                        struct2_df = struct2_df[struct2_df['res_atom_name'] == 'CA']
                    elif subset == 'backbone':
                        struct1_df = struct1_df[struct1_df['res_atom_name'].isin(['CA', 'C', 'N', 'O'])]
                        struct2_df = struct2_df[struct2_df['res_atom_name'].isin(['CA', 'C', 'N', 'O'])]

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

def calculate_binding_pocket_rmsd_for_pairs(mapping_dict, exp_processor, pred_processor, pocket_def=None, cutoff=6.0,
                                            distance_cutoff=6.0, position_tolerance=2.0, window_size=20, max_gap=4,
                                            retinal_name='RET'):
    """
    Calculate binding pocket RMSD between pairs of experimental and predicted structures.
    
    Args:
        mapping_dict: Dictionary mapping experimental structure IDs to predicted structure IDs
        exp_processor: CifBaseProcessor containing experimental structures
        pred_processor: CifBaseProcessor containing predicted structures
        pocket_def: Optional predefined pocket residues
        cutoff: Distance cutoff for binding pocket definition (same as distance_cutoff)
        distance_cutoff: Distance cutoff for binding pocket definition
        position_tolerance: Tolerance for position matching
        window_size: Window size for smoothing
        max_gap: Maximum gap size
        retinal_name: Name of the retinal ligand ('RET' or 'LIG')
        
    Returns:
        dict: Dictionary of RMSD results for each structure pair
        
    Process Steps:
        1. Filter CA atoms for overall alignment
        2. Align the structures using CA atoms 
        3. Apply rotation matrix to full structure including retinal
        4. Calculate RMSD for:
           a. Overall CA atoms
           b. Binding pocket atoms (residues within cutoff distance of retinal)
           c. Retinal atoms
    """
    # Use distance_cutoff if provided, otherwise use cutoff
    effective_cutoff = distance_cutoff if distance_cutoff is not None else cutoff
    
    # Initialize results dictionary
    results = {}
    
    # Setup progress bar for processing structure pairs
    print(f"Processing {len(mapping_dict)} structure pairs for binding pocket RMSD...")
    
    # Process each experimental-predicted structure pair from the mapping with progress bar
    for exp_id, mapping_info in tqdm(mapping_dict.items(), desc="Binding Pocket RMSD", unit="pairs"):
        # Check if mapping is in new nested format or old direct format
        if isinstance(mapping_info, dict) and 'predicted' in mapping_info:
            # New format - nested dictionary
            pred_id = mapping_info['predicted']
        else:
            # Old format - direct mapping
            pred_id = mapping_info
            
        print(f"[INFO] Processing {exp_id} <-> {pred_id} pair")
        
        # Get experimental structure data
        exp_data = exp_processor.data[exp_processor.data['pdb_id'] == exp_id].copy()
        if exp_data.empty:
            print(f"[WARNING] No experimental structure data found for {exp_id}")
            continue
        
        # Get predicted structure data
        pred_data = pred_processor.data[pred_processor.data['pdb_id'] == pred_id].copy()
        if pred_data.empty:
            print(f"[WARNING] No predicted structure data found for {pred_id}")
            continue
        
        # These are already DataFrames from the processor, not dictionaries
        exp_df = exp_data
        pred_df = pred_data
        
        #############################################################
        # STEP 1: Extract retinal and validate it exists in both structures
        #############################################################
        
        # Find retinal in experimental structure
        exp_ret = exp_df[exp_df['res_name3l'] == retinal_name]
        if exp_ret.empty:
            results[exp_id] = {'error': f'No retinal ({retinal_name}) found in experimental structure'}
            continue
        
        # Find retinal in predicted structure - also check for 'LIG' if using 'RET'
        pred_ret = pred_df[pred_df['res_name3l'] == retinal_name]
        if pred_ret.empty and retinal_name == 'RET':
            # Try with 'LIG' as some prediction models use this name
            pred_ret = pred_df[pred_df['res_name3l'] == 'LIG']
            # Rename 'LIG' to 'RET' for consistency
            if not pred_ret.empty:
                pred_df.loc[pred_df['res_name3l'] == 'LIG', 'res_name3l'] = 'RET'
                pred_ret = pred_df[pred_df['res_name3l'] == 'RET']
                
        if pred_ret.empty:
            results[exp_id] = {'error': f'No retinal ({retinal_name} or LIG) found in predicted structure'}
            continue
        
        #############################################################
        # STEP 2: Define binding pocket residues
        #############################################################
        
        if pocket_def is None:
            # Extract CA atoms from experimental structure
            try:
                exp_ca = exp_df[exp_df['res_atom_name'] == 'CA']
            except TypeError:
                print("[DEBUG] Using string conversion for CA atom selection")
                exp_ca = exp_df[exp_df['res_atom_name'].astype(str) == 'CA']
            
            # Calculate distances from each CA atom to the nearest retinal atom
            dist_to_ret = calculate_min_distances(
                exp_ca[['x', 'y', 'z']].astype(float).values,
                exp_ret[['x', 'y', 'z']].astype(float).values
            )
            
            # Ensure dist_to_ret is properly handled - it may be a list or numpy array
            if isinstance(dist_to_ret, list):
                dist_to_ret = np.array(dist_to_ret)
            
            # Define pocket as residues within effective_cutoff distance
            pocket_indices = [i for i, dist in enumerate(dist_to_ret) if dist <= effective_cutoff]
            pocket_residues = exp_ca.iloc[pocket_indices]['auth_seq_id'].tolist()
            print(f"[INFO] Identified {len(pocket_residues)} pocket residues within {effective_cutoff}Å of retinal")
        else:
            pocket_residues = pocket_def
            print(f"[INFO] Using {len(pocket_residues)} predefined pocket residues")
        
        #############################################################
        # STEP 3: Extract CA atoms for overall alignment
        #############################################################
        
        # Extract CA atoms for alignment
        try:
            exp_ca = exp_df[exp_df['res_atom_name'] == 'CA']
            pred_ca = pred_df[pred_df['res_atom_name'] == 'CA']
        except TypeError:
            # Fall back to string conversion for atom selection
            print("[DEBUG] Using string conversion for CA atom selection")
            exp_ca = exp_df[exp_df['res_atom_name'].astype(str) == 'CA']
            pred_ca = pred_df[pred_df['res_atom_name'].astype(str) == 'CA']
        
        if exp_ca.empty or pred_ca.empty:
            results[exp_id] = {'error': 'Could not find CA atoms for alignment'}
            continue
            
        print(f"[INFO] Found {len(exp_ca)} CA atoms in experimental structure, {len(pred_ca)} in predicted structure")
        
        #############################################################
        # STEP 4: Align structures based on CA atoms
        #############################################################
        
        try:
            # Extract coordinates of CA atoms
            exp_ca_coords = exp_ca[['x', 'y', 'z']].astype(float).values
            pred_ca_coords = pred_ca[['x', 'y', 'z']].astype(float).values
            
            # Perform alignment and get rotation matrix and translation vector
            R, t, _, overall_ca_rmsd = get_structure_alignment(exp_ca_coords, pred_ca_coords) # (reference, query) input
            print(f"[INFO] Overall CA RMSD: {overall_ca_rmsd:.3f}Å")
            
            # Apply transformation to entire predicted structure
            pred_coords = pred_df[['x', 'y', 'z']].astype(float).values
            pred_df[['x', 'y', 'z']] = np.dot(pred_coords, R) + t
            
            #############################################################
            # STEP 5: Calculate binding pocket RMSD
            #############################################################
            
            # Calculate RMSD for binding pocket residues (all atoms)
            pocket_rmsd_sum = 0
            pocket_atom_count = 0
            per_residue_rmsd = {}
            
            for res_id in pocket_residues:
                # Extract all atoms for this residue from experimental structure
                exp_res = exp_df[exp_df['auth_seq_id'] == res_id]
                if exp_res.empty:
                    continue
                
                # Extract corresponding residue atoms from predicted structure
                pred_res = pred_df[pred_df['auth_seq_id'] == res_id]
                if pred_res.empty:
                    continue
                
                # Match atoms by name and calculate distances
                matched_atoms = []
                
                for atom_name in exp_res['res_atom_name'].unique():
                    exp_atom = exp_res[exp_res['res_atom_name'] == atom_name]
                    pred_atom = pred_res[pred_res['res_atom_name'] == atom_name]
                    
                    if not exp_atom.empty and not pred_atom.empty:
                        # Get experimental atom coordinates
                        exp_coords = exp_atom[['x', 'y', 'z']].astype(float).values[0]
                        
                        # Get predicted atom original coordinates
                        pred_atom_coords = pred_atom[['x', 'y', 'z']].astype(float).values[0]
                        
                        # Apply rotation and translation to predicted atom coordinates
                        aligned_pred_coords = np.dot(pred_atom_coords.reshape(1, 3), R.T) + t
                        aligned_pred_coords = aligned_pred_coords[0]  # Reshape to 1D array
                        
                        # Calculate distance
                        dist = np.sqrt(np.sum((exp_coords - aligned_pred_coords) ** 2))
                        matched_atoms.append((atom_name, dist, exp_coords, aligned_pred_coords))
                        
                        # Add to overall binding pocket RMSD sum
                        pocket_rmsd_sum += dist ** 2
                        pocket_atom_count += 1
                
                # Calculate RMSD for this residue
                if matched_atoms:
                    # Store error details for each atom
                    atom_errors = {}
                    for atom_name, dist, exp_coords, aligned_pred_coords in matched_atoms:
                        atom_errors[atom_name] = {
                            'distance': dist,
                            'error_x': abs(exp_coords[0] - aligned_pred_coords[0]),
                            'error_y': abs(exp_coords[1] - aligned_pred_coords[1]),
                            'error_z': abs(exp_coords[2] - aligned_pred_coords[2])
                        }
                    
                    # Calculate RMSD for this residue
                    res_rmsd = np.sqrt(np.mean([d ** 2 for _, d, _, _ in matched_atoms]))
                    res_type = exp_res['res_name3l'].iloc[0]
                    
                    # Store per-residue RMSD information
                    per_residue_rmsd[res_id] = {
                        'rmsd': res_rmsd,
                        'res_type': res_type,
                        'atom_count': len(matched_atoms),
                        'atom_errors': atom_errors
                    }
            
            # Calculate overall binding pocket RMSD
            if pocket_atom_count > 0:
                binding_pocket_rmsd = np.sqrt(pocket_rmsd_sum / pocket_atom_count)
                print(f"[INFO] Binding pocket RMSD: {binding_pocket_rmsd:.3f}Å (based on {pocket_atom_count} atoms)")
            else:
                binding_pocket_rmsd = np.nan
                print(f"[WARNING] Could not calculate binding pocket RMSD - no matching atoms")
            
            #############################################################
            # STEP 6: Calculate retinal RMSD
            #############################################################
            
            # Extract coordinates
            pred_ret_coords = pred_ret[['x', 'y', 'z']].astype(float).values
            exp_ret_coords = exp_ret[['x', 'y', 'z']].astype(float).values

            # Apply alignment to predicted retinal coordinates
            aligned_pred_ret_coords = np.dot(pred_ret_coords, R) + t

            # Calculate mean closest distance between retinal atoms
            retinal_rmsd = compute_retinal_mean_closest_distance(exp_ret_coords, aligned_pred_ret_coords)
            print(f"[INFO] Retinal RMSD: {retinal_rmsd:.3f}Å")
            
            #############################################################
            # STEP 7: Store all results
            #############################################################
            
            results[exp_id] = {
                'backbone_rmsd': overall_ca_rmsd,
                'overall_pocket_rmsd': binding_pocket_rmsd,
                'retinal_rmsd': retinal_rmsd,
                'per_residue_rmsd': per_residue_rmsd,
                'pocket_residues': pocket_residues,
                # Add rotation matrix and translation vector for reference
                'alignment': {
                    'rotation': R.tolist(),
                    'translation': t.tolist()
                }
            }
        
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


def compare_structures(data_dict, output_dir='output', visualize=True):
    """
    Step 5: Structure comparison

    This function calculates RMSD between pairs of unique structures and creates
    a similarity matrix without using PropertyProcessor.

    It excludes predicted structures that have experimental counterparts to avoid
    redundancy and focus on unique structural information.

    Args:
        data_dict: Dictionary with data from previous steps
        output_dir: Directory to save output files
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
