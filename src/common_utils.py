"""
Common utility functions shared between multiple modules to avoid circular imports.
"""

import numpy as np
from scipy.spatial.distance import cdist
from scipy.optimize import linear_sum_assignment

def compute_retinal_mean_closest_distance(exp_ret_coords, pred_ret_coords):
    """
    Given two Nx3 arrays of coordinates for RET (retinal) atoms,
    finds the best 1:1 pairing (Hungarian method) that minimizes
    total distance. Returns the average matched distance.
    
    If the sets have different lengths, it matches the smaller
    set to a subset of the larger set, ignoring extra atoms.
    
    Args:
        exp_ret_coords (ndarray): shape (N,3)
        pred_ret_coords (ndarray): shape (M,3)
        
    Returns:
        float: mean distance over the matched pairs
    """
    # Ensure float arrays
    exp_ret_coords = exp_ret_coords.astype(float)
    pred_ret_coords = pred_ret_coords.astype(float)
    
    # If sets differ in size, we match whichever is smaller to a subset of the larger
    n_exp = len(exp_ret_coords)
    n_pred = len(pred_ret_coords)
    if n_exp == 0 or n_pred == 0:
        return np.nan  # or raise an error
    
    # 1) Build NxM distance matrix (for experimental vs. predicted)
    dist_mat = cdist(exp_ret_coords, pred_ret_coords)  # shape (N, M)
    
    # 2) If N \!= M, we only match min(N,M) pairs
    #    We'll do the "Hungarian method" in a square matrix:
    #    - If one set is larger, we add dummy rows or columns with 0-dist to allow linear_sum_assignment
    #      but effectively ignoring extra atoms.
    N = max(n_exp, n_pred)
    # build a square cost matrix of size NxN, fill with large values
    cost_mat = np.full((N, N), fill_value=1e6, dtype=float)
    cost_mat[:n_exp, :n_pred] = dist_mat  # fill the top-left region
    
    # 3) Solve
    row_ind, col_ind = linear_sum_assignment(cost_mat)
    
    # The matched pairs are the first min(n_exp, n_pred) among row_ind, col_ind
    # but some will be dummy rows or columns if N>n_exp or N>n_pred.
    matched_rows = []
    matched_cols = []
    for r, c in zip(row_ind, col_ind):
        # only keep if r < n_exp and c < n_pred => real, not dummy
        if r < n_exp and c < n_pred:
            matched_rows.append(r)
            matched_cols.append(c)
    
    # 4) Extract distances for those matched real pairs
    matched_dists = dist_mat[matched_rows, matched_cols]  # shape ~ min(n_exp,n_pred)
    
    # 5) Return the *average* distance
    if len(matched_dists) > 0:
        return matched_dists.mean()
    else:
        return np.nan

def find_retinal_within_cutoff(full_structure_df, chain_df, cutoff=6.0, retinal_name='RET'):
    """
    From an entire model (full_structure_df) extract all RET atoms
    that lie within 'cutoff' Å of any atom in 'chain_df'.
    
    Args:
        full_structure_df (pd.DataFrame): entire model's coordinates
            with columns ['x','y','z','res_name3l', ...].
        chain_df (pd.DataFrame): the chain subset you're aligning
            with columns ['x','y','z', ...].
        cutoff (float): distance threshold in Å
        retinal_name (str): typically 'RET'.
        
    Returns:
        pd.DataFrame: subset of 'full_structure_df' containing RET atoms
                      within 'cutoff' Å of chain_df.
    """
    # 1) Gather chain coordinates (excluding any retinal atoms in the chain)
    chain_no_ret = chain_df[chain_df['res_name3l'] != retinal_name]
    chain_no_ret = chain_no_ret[chain_no_ret['res_name3l'] != 'LIG']
    # Also exclude LYR residues as they contain retinal
    chain_no_ret = chain_no_ret[chain_no_ret['res_name3l'] != 'LYR']
    
    # If after excluding retinal the chain is empty, we can't calculate distances
    if chain_no_ret.empty:
        # Return all retinal atoms from the full structure
        ret_mask = (full_structure_df['res_name3l'] == retinal_name) | (full_structure_df['res_name3l'] == 'LIG') | (full_structure_df['res_name3l'] == 'LYR')
        return full_structure_df[ret_mask].copy()
        
    chain_coords = chain_no_ret[['x', 'y', 'z']].astype(float).values  # shape: (N,3)
    
    # 2) Subset entire model to just RET or LIG (for predicted structures)
    ret_mask = (full_structure_df['res_name3l'] == retinal_name) | (full_structure_df['res_name3l'] == 'LIG') | (full_structure_df['res_name3l'] == 'LYR')
    ret_df = full_structure_df[ret_mask].copy()
    
    # Minimal debug outputs for retinal search
    if not ret_df.empty:
        res_types = ret_df['res_name3l'].unique()
        
        # Count retinal molecules if auth_seq_id is available
        if 'auth_seq_id' in ret_df.columns:
            unique_rets = ret_df['auth_seq_id'].unique()
            if len(unique_rets) > 1:
                print(f"Found {len(unique_rets)} retinal molecules")
    
    if ret_df.empty:
        return ret_df  # no RET found at all
    
    # 3) For each RET atom, compute distance to every NON-RETINAL chain atom
    ret_coords = ret_df[['x', 'y', 'z']].astype(float).values  # shape: (M,3)
    
    # We'll compute min distance from each RET atom to the chain (excluding retinals)
    dist_mat = cdist(ret_coords, chain_coords)  # shape (M, N)
    min_dists = np.min(dist_mat, axis=1)  # shape (M,)
    
    # 4) Keep only those RET atoms where min_dist < cutoff
    keep_mask = (min_dists < cutoff)
    
    # If no retinal atoms are within cutoff, but we found retinals, use an adaptive cutoff
    if not np.any(keep_mask) and len(ret_df) > 0:
        adaptive_cutoff = min(min_dists.min() * 1.2, 20.0)  # Use slightly larger than minimum, cap at 20Å
        keep_mask = (min_dists < adaptive_cutoff)
    
    ret_in_range = ret_df.iloc[np.where(keep_mask)[0]].copy()
    
    # If we have multiple retinal molecules, we need to select just one complete molecule
    # First check if we have multiple retinal molecules in our selection
    if 'auth_seq_id' in ret_in_range.columns:
        unique_rets_in_range = ret_in_range['auth_seq_id'].unique()
        
        if len(unique_rets_in_range) > 1:
            # Count atoms for each retinal molecule
            ret_counts = []
            for ret_id in unique_rets_in_range:
                ret_atoms = ret_in_range[ret_in_range['auth_seq_id'] == ret_id]
                ret_counts.append((ret_id, len(ret_atoms)))
            
            # Sort by atom count (most atoms first)
            ret_counts.sort(key=lambda x: x[1], reverse=True)
            
            # Select the retinal with the most atoms
            best_ret_id = ret_counts[0][0]
            ret_in_range = ret_in_range[ret_in_range['auth_seq_id'] == best_ret_id].copy()
    
    # If we still couldn't find any retinal atoms within cutoff distance
    if ret_in_range.empty and not ret_df.empty:
        # Fall back to selecting the first retinal molecule if we can identify it
        if 'auth_seq_id' in ret_df.columns and len(ret_df['auth_seq_id'].unique()) > 0:
            first_ret_id = ret_df['auth_seq_id'].unique()[0]
            ret_in_range = ret_df[ret_df['auth_seq_id'] == first_ret_id].copy()
        else:
            # No way to differentiate between retinals, just return all
            ret_in_range = ret_df.copy()
    
    # 5) If we found LIG atoms and they should be renamed to RET
    if 'LIG' in ret_in_range['res_name3l'].unique():
        ret_in_range.loc[ret_in_range['res_name3l'] == 'LIG', 'res_name3l'] = retinal_name
    
    # 6) If we found LYR atoms, convert them to RET
    if 'LYR' in ret_in_range['res_name3l'].unique():
        ret_in_range.loc[ret_in_range['res_name3l'] == 'LYR', 'res_name3l'] = retinal_name
    
    return ret_in_range
