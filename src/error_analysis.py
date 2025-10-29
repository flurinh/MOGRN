"""
Functions for analyzing errors between experimental and predicted structures.
These functions calculate and summarize atom-level errors, distances between
structures, and create summary tables of error metrics.
"""

import numpy as np
import pandas as pd
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment
import os


# Import shared utility functions from common_utils
from src.common_utils import (
    compute_retinal_mean_closest_distance,
    find_retinal_within_cutoff
)

from src.structure_comparison import calculate_binding_pocket_rmsd_for_pairs


def calculate_atom_level_errors(exp_coords, pred_coords, atom_names):
    """
    Calculate detailed error statistics for each atom in a residue.
    
    Args:
        exp_coords: Array of experimental atom coordinates.
        pred_coords: Array of predicted atom coordinates.
        atom_names: List of atom names corresponding to coordinates.
        
    Returns:
        dict: Per-atom error statistics.
    """
    atom_errors = {}
    for i, atom_name in enumerate(atom_names):
        # Calculate Euclidean distance between corresponding atoms
        dist = np.sqrt(np.sum((exp_coords[i] - pred_coords[i]) ** 2))
        # Calculate component-wise errors
        diff_x = abs(exp_coords[i][0] - pred_coords[i][0])
        diff_y = abs(exp_coords[i][1] - pred_coords[i][1])
        diff_z = abs(exp_coords[i][2] - pred_coords[i][2])
        atom_errors[atom_name] = {
            'distance': dist,
            'error_x': diff_x,
            'error_y': diff_y,
            'error_z': diff_z
        }
    return atom_errors

def summarize_atom_errors(binding_pocket_results):
    """
    Summarize atom-level errors across all proteins and residues.
    
    Args:
        binding_pocket_results: Dictionary of results from binding pocket analysis
        
    Returns:
        DataFrame: Summarized atom-level errors
    """
    all_atom_errors = []
    for protein, results in binding_pocket_results.items():
        if 'error' in results:
            continue
        # --- Per-residue RMSD (binding pocket) ---
        per_res = results.get('per_residue_rmsd', {})
        for res_id, res_data in per_res.items():
            for atom_name, atom_error in res_data['atom_errors'].items():
                all_atom_errors.append({
                    'protein': protein,
                    'residue_id': res_id,
                    'residue_type': res_data['res_type'],
                    'atom': atom_name,
                    **atom_error
                })
        # --- If we stored RET errors separately ---
        ret_data = results.get('retinal_atom_errors')
        if ret_data:
            for atom_name, atom_error in ret_data.items():
                all_atom_errors.append({
                    'protein': protein,
                    'residue_id': 'RET',  # or some ID
                    'residue_type': 'RET',
                    'atom': atom_name,
                    **atom_error
                })
    if not all_atom_errors:
        return pd.DataFrame()  # nothing to summarize
    return pd.DataFrame(all_atom_errors)

# The compute_retinal_mean_closest_distance and find_retinal_within_cutoff functions
# have been moved to common_utils.py to avoid circular imports

def make_rmsd_table(binding_pocket_results):
    """
    Build a DataFrame with columns:
      [protein, backbone_rmsd, pocket_rmsd, retinal_rmsd]
    from the dictionary produced by calculate_binding_pocket_rmsd_for_pairs.
    
    Args:
        binding_pocket_results: Dictionary of binding pocket RMSD results
        
    Returns:
        DataFrame: Summary of RMSD values for each protein
    """
    rows = []
    for protein, data in binding_pocket_results.items():
        if 'error' in data:
            # You can either skip or fill with NaNs
            rows.append({
                'protein': protein,
                'backbone_rmsd': np.nan,
                'pocket_rmsd': np.nan,
                'retinal_rmsd': np.nan
            })
        else:
            backbone = data.get('backbone_rmsd', np.nan)
            pocket = data.get('overall_pocket_rmsd', np.nan)
            retinal = data.get('retinal_rmsd', np.nan)
            rows.append({
                'protein': protein,
                'backbone_rmsd': backbone,
                'pocket_rmsd': pocket,
                'retinal_rmsd': retinal
            })
    df = pd.DataFrame(rows, columns=['protein', 'backbone_rmsd', 'pocket_rmsd', 'retinal_rmsd'])
    return df


def calculate_structure_errors(data_dict, output_dir='outputs', visualize=True):
    """
    Step 2: Calculate errors between experimental and predicted structures

    This function calculates RMSD errors between paired experimental and
    predicted structures without using PropertyProcessor.

    Args:
        data_dict: Dictionary with data from previous step
        output_dir: Directory to save outputs files
        visualize: Whether to generate visualizations

    Returns:
        Dictionary with error data and updated processed structures
    """
    print("\n=== Step 2: Error Calculation ===")

    # Create outputs directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Unpack necessary data
    cp_mo_exp = data_dict['cp_mo_exp']
    cp_mo_pred = data_dict['cp_mo_pred']
    cp_hide_exp = data_dict['cp_hide_exp']
    cp_hide_pred = data_dict['cp_hide_pred']
    processed_structures = data_dict['processed_structures']

    # Use the GLOBALLY DEFINED structure mapping - don't recreate it
    structure_mapping = data_dict.get('structure_mapping', {})

    print(f"\n=== Using global structure mapping with {len(structure_mapping)} pairs for error calculation ===")

    # Show some example mappings to confirm we're using the right mapping
    if structure_mapping:
        print("Example structure pairs from global mapping:")
        for i, (exp_id, pred_id) in enumerate(structure_mapping.items()):
            if i < 5:  # Show up to 5 examples
                print(f"  - {exp_id} -> {pred_id}")
    else:
        print("WARNING: No structure mappings available. Error calculation may fail.")

    # Split into dataset-specific mappings for error calculation
    hideaki_mapping = {}
    mo_mapping = {}
    unmatched_mappings = {}

    for exp_id, pred_id in structure_mapping.items():
        # Determine which dataset this belongs to based on processor PDB IDs
        if exp_id in cp_hide_exp.pdb_ids and pred_id in cp_hide_pred.pdb_ids:
            hideaki_mapping[exp_id] = pred_id
        elif exp_id in cp_mo_exp.pdb_ids and pred_id in cp_mo_pred.pdb_ids:
            mo_mapping[exp_id] = pred_id
        else:
            # This pair doesn't match our expected datasets
            unmatched_mappings[exp_id] = pred_id

            # Debug the mismatch
            in_hide_exp = exp_id in cp_hide_exp.pdb_ids
            in_mo_exp = exp_id in cp_mo_exp.pdb_ids
            in_hide_pred = pred_id in cp_hide_pred.pdb_ids
            in_mo_pred = pred_id in cp_mo_pred.pdb_ids

            print(f"DEBUG: Unmatched mapping - {exp_id} -> {pred_id}")
            print(f"  Experimental in Hideaki dataset: {in_hide_exp}")
            print(f"  Experimental in MO dataset: {in_mo_exp}")
            print(f"  Predicted in Hideaki dataset: {in_hide_pred}")
            print(f"  Predicted in MO dataset: {in_mo_pred}")

    # If we have unmatched mappings, but they should work, try to fix them
    if unmatched_mappings and (len(hideaki_mapping) == 0 or len(mo_mapping) == 0):
        print(f"DEBUG: Attempting to fix {len(unmatched_mappings)} unmatched mappings")

        # Try to match structures by name similarity with fuzzy matching
        for exp_id, mapping in unmatched_mappings.items():
            pred_id = mapping  # In simplified format, mapping is just the pred_id

            # Try to match to hideaki dataset
            best_hide_exp_match = None
            best_hide_pred_match = None
            best_hide_exp_score = 0
            best_hide_pred_score = 0

            for hide_exp in cp_hide_exp.pdb_ids:
                # Simple string similarity score
                score = sum(c1 == c2 for c1, c2 in zip(exp_id, hide_exp)) / max(len(exp_id), len(hide_exp))
                if score > best_hide_exp_score:
                    best_hide_exp_score = score
                    best_hide_exp_match = hide_exp

            for hide_pred in cp_hide_pred.pdb_ids:
                score = sum(c1 == c2 for c1, c2 in zip(pred_id, hide_pred)) / max(len(pred_id), len(hide_pred))
                if score > best_hide_pred_score:
                    best_hide_pred_score = score
                    best_hide_pred_match = hide_pred

            # Try to match to MO dataset
            best_mo_exp_match = None
            best_mo_pred_match = None
            best_mo_exp_score = 0
            best_mo_pred_score = 0

            for mo_exp in cp_mo_exp.pdb_ids:
                score = sum(c1 == c2 for c1, c2 in zip(exp_id, mo_exp)) / max(len(exp_id), len(mo_exp))
                if score > best_mo_exp_score:
                    best_mo_exp_score = score
                    best_mo_exp_match = mo_exp

            for mo_pred in cp_mo_pred.pdb_ids:
                score = sum(c1 == c2 for c1, c2 in zip(pred_id, mo_pred)) / max(len(pred_id), len(mo_pred))
                if score > best_mo_pred_score:
                    best_mo_pred_score = score
                    best_mo_pred_match = mo_pred

            # Determine best dataset to assign to
            if best_hide_exp_score > 0.7 and best_hide_pred_score > 0.7:
                hideaki_mapping[best_hide_exp_match] = best_hide_pred_match
                print(f"DEBUG: Fixed mapping for Hideaki: {best_hide_exp_match} -> {best_hide_pred_match}")
            elif best_mo_exp_score > 0.7 and best_mo_pred_score > 0.7:
                mo_mapping[best_mo_exp_match] = best_mo_pred_match
                print(f"DEBUG: Fixed mapping for MO: {best_mo_exp_match} -> {best_mo_pred_match}")

    # Final report on mappings
    print(f"Final mapping counts: {len(hideaki_mapping)} Hideaki pairs, {len(mo_mapping)} MO pairs")

    # Calculate errors for Hideaki structures
    print(f"Calculating errors for {len(hideaki_mapping)} Hideaki structure pairs...")
    binding_pocket_results_hide = {}

    cp_mo_pred.format_data_types()
    cp_mo_exp.format_data_types()
    cp_hide_exp.format_data_types()
    cp_hide_pred.format_data_types()

    print(cp_hide_pred.data.dtypes)

    if hideaki_mapping:
        binding_pocket_results_hide = calculate_binding_pocket_rmsd_for_pairs(
            hideaki_mapping, cp_hide_exp, cp_hide_pred,
            retinal_name='RET',
            distance_cutoff=6.0,
            position_tolerance=2.0,
            window_size=20,
            max_gap=4
        )

        # Create RMSD table
        hide_errors_df = make_rmsd_table(binding_pocket_results_hide)
        hide_errors_df.to_csv(os.path.join(output_dir, 'hideaki_errors.csv'))
        print(f"Hideaki errors saved to {os.path.join(output_dir, 'hideaki_errors.csv')}")
    else:
        hide_errors_df = pd.DataFrame()
        print("No Hideaki structure pairs found for error calculation")

    # Calculate errors for MO structures
    print(f"Calculating errors for {len(mo_mapping)} MO structure pairs...")
    binding_results_mo_exp = {}

    if mo_mapping:
        binding_results_mo_exp = calculate_binding_pocket_rmsd_for_pairs(
            mo_mapping, cp_mo_exp, cp_mo_pred,
            retinal_name='RET',
            distance_cutoff=6.0,
            position_tolerance=2.0,
            window_size=12,
            max_gap=4
        )

        # Create RMSD table
        mo_exp_errors_df = make_rmsd_table(binding_results_mo_exp)
        mo_exp_errors_df = mo_exp_errors_df[mo_exp_errors_df['retinal_rmsd'] < 5]  # Filter out high RMSD values
        mo_exp_errors_df.to_csv(os.path.join(output_dir, 'mo_exp_errors.csv'))
        print(f"MO errors saved to {os.path.join(output_dir, 'mo_exp_errors.csv')}")
    else:
        mo_exp_errors_df = pd.DataFrame()
        print("No MO structure pairs found for error calculation")

    # Combine all error results
    all_binding_results = {**binding_pocket_results_hide, **binding_results_mo_exp}

    # Get the valid structures for further analysis
    valid_structures = set()

    # Add all experimental structures with valid error calculations
    for exp_id in list(binding_pocket_results_hide.keys()) + list(binding_results_mo_exp.keys()):
        valid_structures.add(exp_id)
        # Also add the corresponding predicted structure
        if exp_id in structure_mapping:
            valid_structures.add(structure_mapping[exp_id])

    # Add any remaining structures that we want to keep for analysis
    for pdb_id in processed_structures.keys():
        # Keep structures that have retinal
        if 'df_ret' in processed_structures[pdb_id] and not processed_structures[pdb_id]['df_ret'].empty:
            valid_structures.add(pdb_id)

    # Filter processed structures to only include valid ones
    filtered_structures = {key: processed_structures[key] for key in valid_structures
                           if key in processed_structures}

    # Add error metrics to structure metadata
    for pdb_id, data in filtered_structures.items():
        # Add error data for experimental structures
        if pdb_id in all_binding_results:
            binding_results = all_binding_results[pdb_id]
            data['binding_pocket_error'] = binding_results
            data['retinal_rmsd'] = binding_results.get('retinal_rmsd', float('nan'))
            data['pocket_rmsd'] = binding_results.get('pocket_rmsd', float('nan'))

            # Add reference to the paired structure
            if pdb_id in structure_mapping:
                paired_id = structure_mapping[pdb_id]
                data['paired_structure'] = paired_id

    # Fix some issues with the retinal data
    if hasattr(cp_mo_pred, 'data') and cp_mo_pred.data is not None and not cp_mo_pred.data.empty:
        cp_mo_pred.data.loc[cp_mo_pred.data['res_name3l'] == 'LIG', 'res_name3l'] = 'RET'
        cp_mo_pred.data.loc[cp_mo_pred.data['res_name3l'] == 'RET', 'auth_chain_id'] = 'A'

    # Print some statistics on the errors
    print("\nError Statistics:")
    if not hide_errors_df.empty:
        print(f"  Hideaki structures: {len(hide_errors_df)} pairs analyzed")
        print(f"    Mean retinal RMSD: {hide_errors_df['retinal_rmsd'].mean():.2f}Å")
        print(f"    Mean pocket RMSD: {hide_errors_df['pocket_rmsd'].mean():.2f}Å")

    if not mo_exp_errors_df.empty:
        print(f"  MO structures: {len(mo_exp_errors_df)} pairs analyzed")
        print(f"    Mean retinal RMSD: {mo_exp_errors_df['retinal_rmsd'].mean():.2f}Å")
        print(f"    Mean pocket RMSD: {mo_exp_errors_df['pocket_rmsd'].mean():.2f}Å")

    return {
        'processed_structures': filtered_structures,
        'structure_mapping': structure_mapping,
        'hide_errors_df': hide_errors_df,
        'mo_exp_errors_df': mo_exp_errors_df,
        'binding_pocket_results': binding_pocket_results_hide,
        'binding_results_mo_exp': binding_results_mo_exp
    }


def qdock():
    return

