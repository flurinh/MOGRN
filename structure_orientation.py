"""
Functions for orienting and normalizing protein structures.
These functions handle structure orientation based on principal component
analysis and alignment to the z-axis for consistent membrane protein display.
"""

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

def orient_structure(df):
    """
    Orient a structure based on PCA of CA atom coordinates.
    The structure is rotated so that the first principal component
    aligns with the z-axis (membrane normal) and centered at the origin.
    
    Args:
        df: DataFrame containing structure coordinates
        
    Returns:
        DataFrame: Oriented structure
    """
    ca_df = df[df['res_atom_name'] == 'CA']
    if ca_df.empty:
        print("[WARNING] No CA atoms found for orientation; returning original.")
        return df.copy()
    
    ca_coords = ca_df[['x', 'y', 'z']].astype(float).values
    
    # Perform PCA to find principal axes
    pca = PCA(n_components=3)
    pca.fit(ca_coords)
    orientation_vec = pca.components_[0]
    
    # Set up rotation to align with z-axis
    z_axis = np.array([0, 0, 1])
    norm_vec = orientation_vec / np.linalg.norm(orientation_vec)
    
    # Calculate rotation matrix
    if np.allclose(norm_vec, z_axis):
        R = np.eye(3)
    else:
        rot_axis = np.cross(norm_vec, z_axis)
        rot_axis = rot_axis / np.linalg.norm(rot_axis)
        angle = np.arccos(np.clip(np.dot(norm_vec, z_axis), -1.0, 1.0))
        K = np.array([[0, -rot_axis[2], rot_axis[1]],
                      [rot_axis[2], 0, -rot_axis[0]],
                      [-rot_axis[1], rot_axis[0], 0]])
        R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
    
    # Apply rotation and centering
    coords = df[['x', 'y', 'z']].astype(float).values
    rotated = np.dot(coords, R.T)
    centroid = rotated.mean(axis=0)
    centered = rotated - centroid
    
    # Create new DataFrame with transformed coordinates
    df_oriented = df.copy()
    df_oriented[['x', 'y', 'z']] = centered
    
    return df_oriented

def orient_all_structures(structures):
    """
    Apply orientation to all structures in the dictionary
    
    Args:
        structures: Dictionary of structures to orient
        
    Returns:
        dict: Dictionary with oriented structures
    """
    for pdb_id, data in structures.items():
        df = data['df']
        df_norm = orient_structure(df)
        data['df_norm'] = df_norm
        data['df_ca_norm'] = df_norm[df_norm['res_atom_name'] == 'CA'].copy()
    return structures

def compute_orientation_vector_pca(ca_coords):
    """
    Compute the principal component of CA coordinates to determine orientation vector.
    Ensure the principal component aligns with the direction from N-terminus to C-terminus.
    
    Args:
        ca_coords: Numpy array of shape (N,3) with CA coordinates sorted from N-term to C-term.
        
    Returns:
        ndarray: A normalized orientation vector.
    """
    # Make sure we're working with a numpy array
    if not isinstance(ca_coords, np.ndarray):
        ca_coords = np.array(ca_coords)
    
    # Convert coordinates to numeric if they're still strings
    # We need to handle each element or column separately
    numeric_coords = np.zeros_like(ca_coords, dtype=float)
    
    # Check if we need to convert from string to numeric
    if ca_coords.dtype == object:
        print(f"[DEBUG] Converting string coordinates to numeric. Shape: {ca_coords.shape}")
        # Convert each coordinate column separately
        for i in range(ca_coords.shape[1]):
            for j in range(ca_coords.shape[0]):
                try:
                    numeric_coords[j, i] = float(ca_coords[j, i])
                except (ValueError, TypeError):
                    numeric_coords[j, i] = np.nan
    else:
        # Already numeric
        numeric_coords = ca_coords.astype(float)
    
    # Check for NaN values and report
    nan_count = np.isnan(numeric_coords).sum()
    if nan_count > 0:
        print(f"[WARNING] Found {nan_count} NaN values in coordinates after conversion")
    
    # Use clean coordinates for PCA
    clean_coords = np.nan_to_num(numeric_coords)
    
    pca = PCA(n_components=3)
    pca.fit(clean_coords)
    
    # The first principal component represents the principal axis
    principal_comp = pca.components_[0]
    
    # To ensure consistent orientation, check if it points from N to C terminus
    # by taking a reference vector from the first CA to the last CA
    ref_vec = clean_coords[-1] - clean_coords[0]
    ref_vec_norm = np.linalg.norm(ref_vec)
    if ref_vec_norm > 0:
        ref_vec = ref_vec / ref_vec_norm
    else:
        print("[WARNING] Reference vector has zero norm, using principal component directly")
        return principal_comp
    
    # If the dot product is negative, the orientation needs to be flipped
    if np.dot(principal_comp, ref_vec) < 0:
        principal_comp = -principal_comp
    
    return principal_comp

def align_vector_to_z(vec):
    """
    Compute a rotation matrix to align a vector to the z-axis.
    
    Args:
        vec: Vector to align with z-axis
        
    Returns:
        ndarray: Rotation matrix
    """
    # Normalize the vector
    vec = vec / np.linalg.norm(vec)
    
    # Define the z-axis
    z_axis = np.array([0, 0, 1])
    
    # If the vector is already aligned with z, no rotation needed
    if np.allclose(vec, z_axis) or np.allclose(vec, -z_axis):
        return np.eye(3)
    
    # Cross product gives the rotation axis
    rotation_axis = np.cross(vec, z_axis)
    rotation_axis = rotation_axis / np.linalg.norm(rotation_axis)
    
    # Dot product gives the cosine of the angle
    cos_angle = np.dot(vec, z_axis)
    angle = np.arccos(np.clip(cos_angle, -1.0, 1.0))
    
    # Rodrigues' rotation formula
    K = np.array([[0, -rotation_axis[2], rotation_axis[1]],
                 [rotation_axis[2], 0, -rotation_axis[0]],
                 [-rotation_axis[1], rotation_axis[0], 0]])
    R = np.eye(3) + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
    
    return R

def orient_structures_n_terminus_up(structures, iterations=10):
    """
    Orient structures so that the N-terminus is pointing up (positive z).
    Uses PCA on CA atoms to determine the principal axis, then aligns
    this axis with the z-axis such that the N-terminus is in the
    positive z direction.
    
    Args:
        structures: Dictionary of structures to orient
        iterations: Number of iterations for orientation refinement
        
    Returns:
        dict: Dictionary with oriented structures
    """
    for pdb_id, data in structures.items():
        # Get CA coordinates sorted by residue number
        df = data.get('df_norm', data.get('df')).copy()
        try:
            # Try filtering by res_atom_name directly
            ca_df = df[df['res_atom_name'] == 'CA'].sort_values('auth_seq_id')
        except TypeError:
            # If that fails, try string conversion
            print(f"[DEBUG] Using string conversion for CA selection in {pdb_id}")
            ca_df = df[df['res_atom_name'].astype(str) == 'CA'].sort_values('auth_seq_id')
        
        if ca_df.empty:
            print(f"[WARNING] No CA atoms found for {pdb_id}; skipping orientation.")
            continue
        
        # Make sure all coordinate values are numeric with better error handling
        for coord in ['x', 'y', 'z']:
            try:
                ca_df[coord] = pd.to_numeric(ca_df[coord], errors='coerce')
            except Exception as e:
                print(f"[WARNING] Error converting {coord} to numeric for {pdb_id}: {str(e)}")
                print(f"[DEBUG] Sample of {coord} values: {ca_df[coord].iloc[:5].tolist()}")
                # Try a more direct conversion approach
                ca_df[coord] = ca_df[coord].apply(lambda x: float(x) if isinstance(x, (str, int, float)) else np.nan)
        
        # Check that we have valid coordinates after conversion
        nan_mask = ca_df[['x', 'y', 'z']].isna().any(axis=1)
        if nan_mask.any():
            print(f"[WARNING] Found {nan_mask.sum()} CA atoms with NaN coordinates in {pdb_id}")
            ca_df = ca_df[~nan_mask].copy()
            
        if len(ca_df) < 3:
            print(f"[WARNING] Not enough valid CA atoms (only {len(ca_df)}) for {pdb_id}; skipping orientation.")
            continue
            
        print(f"[DEBUG] CA coordinates shape for {pdb_id}: {ca_df[['x', 'y', 'z']].shape}")
        ca_coords = ca_df[['x', 'y', 'z']].values
        
        # Iterative refinement of the orientation
        for _ in range(iterations):
            # Compute orientation vector
            orientation_vec = compute_orientation_vector_pca(ca_coords)
            
            # Create rotation matrix to align with z-axis
            R = align_vector_to_z(orientation_vec)
            
            # Apply rotation
            ca_coords = np.dot(ca_coords, R.T)
            
            # Check if N-terminus is pointing up (positive z)
            n_term_z = ca_coords[0, 2]
            c_term_z = ca_coords[-1, 2]
            
            if n_term_z < c_term_z:
                # Flip the structure if N-terminus is pointing down
                flip_matrix = np.diag([1, 1, -1])  # Flip around z-axis
                ca_coords = np.dot(ca_coords, flip_matrix)
        
        # Apply the final orientation to the full structure
        coords = df[['x', 'y', 'z']].astype(float).values
        rotated = np.dot(coords, R.T)
        
        # Check if we need to flip
        if n_term_z < c_term_z:
            flip_matrix = np.diag([1, 1, -1])
            rotated = np.dot(rotated, flip_matrix)
        
        # Center the structure
        centroid = rotated.mean(axis=0)
        centered = rotated - centroid
        
        # Update the structure with new coordinates
        df_oriented = df.copy()
        df_oriented[['x', 'y', 'z']] = centered
        
        # Store the oriented structure
        data['df_norm'] = df_oriented
        data['df_ca_norm'] = df_oriented[df_oriented['res_atom_name'] == 'CA'].copy()
    
    return structures


