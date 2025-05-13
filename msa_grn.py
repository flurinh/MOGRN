"""
Functions for multiple structure alignment and Generic Residue Numbering (GRN).
These functions create and analyze multiple sequence alignments with GRN position labels.
"""

import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt
from visualization_functions import visualize_msa_distances

def create_msa_table(seq_alignment_dicts, processed_structures_complete, global_ref, atom_type="all"):
    """
    Creates an MSA-like table from alignment dictionaries.
    
    Args:
        seq_alignment_dicts: Dictionary with sequence alignments
        processed_structures_complete: Dictionary with structure data
        global_ref: ID of global reference structure
        atom_type: Atom subset to use ('CA' or 'all')
        
    Returns:
        DataFrame: MSA table with residue information
    """
    # Get global alignments
    global_alignments = seq_alignment_dicts.get('global', {})
    
    # Get global reference structure
    global_ref_struct = processed_structures_complete[global_ref]
    
    # Get appropriate dataframe
    if atom_type == "CA":
        global_ref_df = global_ref_struct['df_ca_norm']
        global_ref_df = global_ref_df[global_ref_df['res_atom_name'] == 'CA']
    else:
        global_ref_df = global_ref_struct['df_norm']
        global_ref_df = global_ref_df.drop_duplicates(subset=['auth_seq_id', 'auth_chain_id'])
    
    # Get global positions (sorted auth_seq_ids)
    global_positions = sorted(global_ref_df['auth_seq_id'].unique())
    
    # Initialize sequences dictionary
    sequences = {}
    
    # Add global reference to sequences
    sequences[global_ref] = {}
    for pos in global_positions:
        residues = global_ref_df[global_ref_df['auth_seq_id'] == pos]
        if not residues.empty:
            if 'res_name1l' in residues.columns:
                aa = residues['res_name1l'].iloc[0]
            else:
                aa = '?'
            sequences[global_ref][pos] = f"{aa}{pos}"
        else:
            sequences[global_ref][pos] = '-'
    
    # Process global alignments
    for struct_id, alignment in global_alignments.items():
        struct = processed_structures_complete[struct_id]
        
        # Get appropriate dataframe
        if atom_type == "CA":
            struct_df = struct['df_ca_norm']
            struct_df = struct_df[struct_df['res_atom_name'] == 'CA']
        else:
            struct_df = struct['df_norm']
            struct_df = struct_df.drop_duplicates(subset=['auth_seq_id', 'auth_chain_id'])
        
        # Initialize sequence with gaps
        sequences[struct_id] = {pos: '-' for pos in global_positions}
        
        # Fill in aligned positions
        for ref_pos, aligned_pos in alignment.items():
            if ref_pos in global_positions:
                residues = struct_df[struct_df['auth_seq_id'] == aligned_pos]
                if not residues.empty:
                    if 'res_name1l' in residues.columns:
                        aa = residues['res_name1l'].iloc[0]
                    else:
                        aa = '?'
                    sequences[struct_id][ref_pos] = f"{aa}{aligned_pos}"
    
    # Process type-specific alignments
    type_alignments = seq_alignment_dicts.get('type', {})
    for type_ref, structs in type_alignments.items():
        if type_ref == global_ref:
            continue  # Skip if type reference is global reference
        
        # Get global alignment for this type reference
        type_to_global = global_alignments.get(type_ref, {})
        
        for struct_id, alignment in structs.items():
            struct = processed_structures_complete[struct_id]
            
            # Get appropriate dataframe
            if atom_type == "CA":
                struct_df = struct['df_ca_norm']
                struct_df = struct_df[struct_df['res_atom_name'] == 'CA']
            else:
                struct_df = struct['df_norm']
                struct_df = struct_df.drop_duplicates(subset=['auth_seq_id', 'auth_chain_id'])
            
            # Initialize sequence with gaps
            sequences[struct_id] = {pos: '-' for pos in global_positions}
            
            # Map through type reference to global reference
            for type_pos, struct_pos in alignment.items():
                # Find corresponding global position
                global_pos = None
                for g_pos, t_pos in type_to_global.items():
                    if t_pos == type_pos:
                        global_pos = g_pos
                        break
                
                if global_pos in global_positions:
                    residues = struct_df[struct_df['auth_seq_id'] == struct_pos]
                    if not residues.empty:
                        if 'res_name1l' in residues.columns:
                            aa = residues['res_name1l'].iloc[0]
                        else:
                            aa = '?'
                        sequences[struct_id][global_pos] = f"{aa}{struct_pos}"
    
    # Convert sequences dictionary to DataFrame
    msa_df = pd.DataFrame(sequences).T
    
    # Make sure columns are in the correct order (sorted global positions)
    msa_df = msa_df[global_positions]
    
    # Rename columns to sequential numbers for compatibility with old code
    msa_df.columns = range(1, len(msa_df.columns) + 1)
    
    # Store the mapping between column numbers and global positions
    # for reference when adding GRN labels
    msa_df.attrs['column_to_auth_seq'] = {i + 1: pos for i, pos in enumerate(global_positions)}
    
    return msa_df

def process_alignment_df(msa_df):
    """
    Converts alignment DataFrame to a position weight matrix.
    
    Args:
        msa_df: MSA DataFrame from create_msa_table
        
    Returns:
        DataFrame: Position weight matrix
    """
    # Extract amino acid identities from MSA
    amino_acids = {}
    
    for col in msa_df.columns:
        aa_counts = {}
        
        for val in msa_df[col]:
            if val == '-':
                continue
                
            # Extract amino acid from value (format: 'A123')
            aa = val[0] if len(val) > 0 else '?'
            
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
        
        amino_acids[col] = aa_counts
    
    # Convert to DataFrame
    aa_df = pd.DataFrame(amino_acids).fillna(0)
    
    # Calculate frequencies
    for col in aa_df.columns:
        col_sum = aa_df[col].sum()
        if col_sum > 0:
            aa_df[col] = aa_df[col] / col_sum
    
    return aa_df

def analyze_residue_composition(msa_df, positions):
    """
    Analyzes residue composition at specific positions.
    Only processes valid TM helix positions (format N.XX where N is 1-7).
    
    Args:
        msa_df: MSA DataFrame from create_msa_table
        positions: List of positions to analyze
        
    Returns:
        dict: Dictionary with position analysis
    """
    results = {}
    
    # Filter positions to only include valid TM helix positions
    filtered_positions = []
    for pos in positions:
        pos_str = str(pos)
        # Only include positions with format N.XX where N is 1-7
        if '.' in pos_str:
            try:
                helix_num = int(pos_str.split('.')[0])
                if 1 <= helix_num <= 7:
                    filtered_positions.append(pos)
            except ValueError:
                # Skip positions with non-numeric helix values
                results[pos] = {'error': 'Invalid helix number (must be 1-7)'}
                continue
        else:
            # Skip positions without dot notation
            results[pos] = {'error': 'Invalid GRN format (must be N.XX)'}
            continue
    
    for pos in filtered_positions:
        if pos not in msa_df.columns:
            results[pos] = {'error': 'Position not found'}
            continue
        
        pos_data = msa_df[pos].value_counts(dropna=False)
        
        # Extract amino acids
        aa_counts = {}
        
        for val, count in pos_data.items():
            if pd.isna(val) or val == '-':
                aa_counts['-'] = aa_counts.get('-', 0) + count
                continue
                
            # Extract amino acid from value (format: 'A123') - first character
            if isinstance(val, str) and len(val) > 0:
                aa = val[0]
                aa_counts[aa] = aa_counts.get(aa, 0) + count
            else:
                # Handle non-string or empty values
                aa_counts['?'] = aa_counts.get('?', 0) + count
        
        # Calculate frequencies
        total = sum(aa_counts.values())
        if total > 0:
            aa_freq = {aa: count / total for aa, count in aa_counts.items()}
            
            # Sort by frequency
            sorted_aa = sorted(aa_freq.items(), key=lambda x: x[1], reverse=True)
            
            results[pos] = {
                'counts': aa_counts,
                'frequencies': aa_freq,
                'sorted': sorted_aa,
                'total': total,
                'conservation': max(aa_freq.values()) if aa_freq else 0
            }
        else:
            results[pos] = {'error': 'No valid residues found'}
    
    return results

def count_residues_by_helix(df):
    """
    Counts residue positions by helix in a GRN-labeled DataFrame.
    Focuses primarily on N.YY format (where N is the helix number) to match your desired output.
    
    Args:
        df: DataFrame with GRN column labels
        
    Returns:
        dict: Dictionary with counts of positions per helix
    """
    # Initialize counters
    tm_counts = {}
    other_counts = {'n': 0, 'c': 0, 'L': 0, 'X': 0}
    
    for col in df.columns:
        col_str = str(col)
        
        # Check for TM helix positions in N.YY format (primary format)
        if '.' in col_str and len(col_str.split('.')) == 2:
            prefix, suffix = col_str.split('.')
            if prefix.isdigit() and int(prefix) >= 1 and int(prefix) <= 7:
                # This is a TM helix position in N.YY format
                tm_counts[prefix] = tm_counts.get(prefix, 0) + 1
            elif prefix == 'n':
                other_counts['n'] += 1
            elif prefix == 'c':
                other_counts['c'] += 1
            elif prefix == 'L':
                other_counts['L'] += 1
            elif prefix == 'X':
                other_counts['X'] += 1
            elif len(prefix) == 2 and prefix.isdigit():
                # This is likely a loop region in AB.CCC format
                other_counts['L'] += 1
            else:
                other_counts['X'] += 1
        
        # Check for TM helix positions in NxYY format (for backward compatibility)
        elif 'x' in col_str:
            try:
                helix, pos = col_str.split('x')
                if helix.isdigit() and int(helix) >= 1 and int(helix) <= 7:
                    tm_counts[helix] = tm_counts.get(helix, 0) + 1
                else:
                    other_counts['X'] += 1
            except:
                other_counts['X'] += 1
        
        # Any other format
        else:
            other_counts['X'] += 1
    
    result = {
        'tm_total': sum(tm_counts.values()),
        'helices': tm_counts,
        'other': other_counts
    }
    
    return result

def sort_grn_columns(df):
    """
    Sorts columns of a GRN-labeled DataFrame in the correct order based on GRN system.
    
    The order follows the GRN system:
    1. N-terminal (n.XX) 
    2. TM helices in order (1xYY through 7xYY, with positions sorted numerically)
    3. Loop regions (AB.CCC)
    4. C-terminal regions (c.XX)
    5. Any other labels
    
    Args:
        df: DataFrame with GRN column labels
        
    Returns:
        DataFrame: DataFrame with sorted columns
    """
    # Function to convert GRN label to a float value for consistent sorting
    def grn_to_sortable_value(label):
        if not isinstance(label, str):
            return 999.999  # Default high value for non-string labels
        
        # N-terminal: convert to negative value
        if label.startswith('n.'):
            try:
                distance = float(label.split('.')[1])
                return -0.01 * distance  # Negative value to sort before TM regions
            except (IndexError, ValueError):
                return -0.001  # Default for malformed n.XX
        
        # C-terminal: convert to high value
        if label.startswith('c.'):
            try:
                distance = float(label.split('.')[1])
                return 100.0 + (0.01 * distance)  # High value to sort after loops
            except (IndexError, ValueError):
                return 100.001  # Default for malformed c.XX
        
        # TM helices with dot notation (prioritize this format)
        match_dot = re.match(r'(\d+)\.(\d+)$', label)
        if match_dot and len(match_dot.groups()) == 2:
            helix, pos = match_dot.groups()
            # Make sure helices 1-7 are properly ordered
            if int(helix) >= 1 and int(helix) <= 7:
                return float(label)
        
        # TM helices with x notation: convert to helix.position format (for backward compatibility)
        match_x = re.match(r'(\d+)x(\d+)', label)
        if match_x:
            helix, pos = match_x.groups()
            return float(f"{helix}.{pos}")
        
        # Loop regions: AB.CCC format
        match_loop = re.match(r'(\d+)(\d+)\.(\d+)', label)
        if match_loop:
            closer_helix, further_helix, distance = match_loop.groups()
            # Scale the distance to be between the helix numbers
            return float(f"{closer_helix}{further_helix}.{distance}")
        
        # General loop label (L.X)
        if label.startswith('L.'):
            try:
                pos = float(label.split('.')[1])
                return 90.0 + (0.01 * pos)  # Sort between loops and C-terminal
            except (IndexError, ValueError):
                return 90.001  # Default for malformed L.X
            
        # Any other format: sort at the end
        return 999.999
    
    # Sort columns
    try:
        sorted_columns = sorted(df.columns, key=grn_to_sortable_value)
        return df[sorted_columns]
    except Exception as e:
        print(f"[WARNING] Error sorting GRN columns: {e}")
        # Fall back to original order
        return df

def create_msa_distance_table(seq_alignment_dicts, processed_structures_complete, global_ref,
                              distance_type="sidechain"):
    """
    Creates a distance table showing the closest distance from each residue to retinal.
    
    Args:
        seq_alignment_dicts: Dictionary with sequence alignments
        processed_structures_complete: Dictionary with structure data
        global_ref: ID of global reference structure
        distance_type: Either "sidechain" (closest atom) or "backbone" (CA atoms only)
        
    Returns:
        DataFrame: Distance table
    """
    from scipy.spatial.distance import cdist
    
    # First, create the alignment table to get the position mapping
    msa_df = create_msa_table(
        seq_alignment_dicts, 
        processed_structures_complete, 
        global_ref,
        atom_type="CA" if distance_type == "backbone" else "all"
    )
    
    column_to_auth_seq = msa_df.attrs.get('column_to_auth_seq', {})
    
    distances = {}
    for struct_id in msa_df.index:
        if struct_id not in processed_structures_complete or 'df_norm' not in processed_structures_complete[struct_id]:
            distances[struct_id] = {col: float('nan') for col in msa_df.columns}
            print(f"[DEBUG] Structure {struct_id} missing or lacks df_norm, setting all distances to NaN")
            continue
        
        struct = processed_structures_complete[struct_id]
        df_norm = struct['df_norm']
        
        # First try 'RET'
        ret_df = df_norm[df_norm['res_name3l'] == 'RET']
        
        # If not found, try 'LIG'
        if ret_df.empty:
            lig_df = df_norm[df_norm['res_name3l'] == 'LIG']
            if not lig_df.empty:
                # Rename 'LIG' to 'RET' for consistency
                df_norm.loc[df_norm['res_name3l'] == 'LIG', 'res_name3l'] = 'RET'
                ret_df = df_norm[df_norm['res_name3l'] == 'RET']
                print(f"[INFO] Renamed 'LIG' to 'RET' in structure {struct_id}")
        
        if ret_df.empty:
            print(f"[DEBUG] No retinal (RET or LIG) found for structure {struct_id}, setting all distances to NaN")
            distances[struct_id] = {col: float('nan') for col in msa_df.columns}
            continue
        
        # Parse auth_seq_id in the structure for consistency
        # First make an explicit copy to avoid the SettingWithCopyWarning
        if distance_type == "backbone":
            res_df = struct['df_ca_norm'].copy()
        else:
            res_df = df_norm[df_norm['res_name3l'] != 'RET'].copy()
            
        # Now modify the copy (this avoids the warning)
        res_df['numeric_seq_id'] = res_df['auth_seq_id'].apply(lambda x: int(x) if str(x).isdigit() else x)

        struct_distances = {col: float('nan') for col in msa_df.columns}
        ret_coords = ret_df[['x', 'y', 'z']].values
        
        for col in msa_df.columns:
            cell_value = msa_df.at[struct_id, col]
            if cell_value == '-':
                # Gap in the alignment
                continue
            
            try:
                # Extract auth_seq_id from cell_value
                digit_part = ''.join(c for c in cell_value if c.isdigit())
                if not digit_part:
                    continue
                
                ref_numeric_id = int(digit_part)
                
                # Find this residue in the structure using numeric_seq_id
                residue_atoms = res_df[res_df['numeric_seq_id'] == ref_numeric_id]
                if residue_atoms.empty:
                    continue
                
                res_coords = residue_atoms[['x', 'y', 'z']].values
                distances_matrix = cdist(res_coords, ret_coords)
                min_distance = float(distances_matrix.min())
                
                struct_distances[col] = min_distance
            except Exception as e:
                print(f"[WARNING] Error calculating distance for {struct_id}, column {col}: {e}")
                continue
        
        distances[struct_id] = struct_distances
    
    distance_df = pd.DataFrame(distances).T
    distance_df = distance_df[msa_df.columns]
    distance_df.attrs['column_to_auth_seq'] = column_to_auth_seq
    
    return distance_df

def find_closest_ret_residue(struct, helix_num):
    """
    Find the auth_seq_id of the residue closest to retinal in a specific helix.
    This can be used as the X.50 reference point.
    
    Args:
        struct: Structure data dictionary
        helix_num: Helix number to analyze
        
    Returns:
        auth_seq_id of closest residue or None if not found
    """
    from scipy.spatial.distance import cdist
    
    if 'df_norm' not in struct:
        return None
    
    df_norm = struct['df_norm']
    
    # Find RET
    if 'res_name3l' not in df_norm.columns:
        return None
    
    # Try 'RET'
    ret_df = df_norm[df_norm['res_name3l'] == 'RET']
    
    # If not found, try 'LIG'
    if ret_df.empty:
        lig_df = df_norm[df_norm['res_name3l'] == 'LIG']
        if not lig_df.empty:
            # Rename 'LIG' to 'RET' for consistency
            df_norm.loc[df_norm['res_name3l'] == 'LIG', 'res_name3l'] = 'RET'
            ret_df = df_norm[df_norm['res_name3l'] == 'RET']
            print(f"[INFO] Renamed 'LIG' to 'RET' in structure for closest residue calculation")
    
    if ret_df.empty:
        return None
    
    # Find residues in the specified helix
    helix_df = df_norm[(df_norm['helix_num'] == helix_num) &
                       (df_norm['tm_helix'] == True) &
                       (df_norm['res_name3l'] != 'RET')]
    
    if helix_df.empty:
        print("did not find residues with specified helix in df_norm")
        return None
    
    # Group by residue (auth_seq_id)
    min_distance = float('inf')
    closest_auth_seq = None
    
    # Calculate distances from each residue to RET
    for auth_seq, group in helix_df.groupby('auth_seq_id'):
        # Get coordinates for this residue
        res_coords = group[['x', 'y', 'z']].values
        
        # Get coordinates for RET
        ret_coords = ret_df[['x', 'y', 'z']].values
        
        # Calculate minimum distance
        distances = cdist(res_coords, ret_coords)
        min_dist = distances.min()
        
        if min_dist < min_distance:
            min_distance = min_dist
            closest_auth_seq = auth_seq
    
    return closest_auth_seq

def create_grn_column_mapping(ref_struct, original_columns):
    """
    Creates a mapping from sequential column numbers to GRN notation based on
    helix assignments in the reference structure, following the Generic Residue
    Numbering (GRN) system.
    
    Args:
        ref_struct: Reference structure data with df_ca_norm containing helix annotations
        original_columns: Original column numbers to be mapped
        
    Returns:
        dict: Mapping from original columns to GRN labels
    """
    if 'df_ca_norm' not in ref_struct:
        print("[ERROR] Reference structure lacks df_ca_norm")
        return {}
    
    ref_df = ref_struct['df_ca_norm']
    
    if 'helix_num' not in ref_df.columns or 'tm_helix' not in ref_df.columns:
        print("[ERROR] Reference structure lacks helix annotations")
        return {}
    
    # Convert original_columns to integers
    original_columns = [int(col) for col in original_columns]
    
    # Create a mapping from column number to auth_seq_id
    column_to_auth_seq = {}
    auth_seq_ids = sorted(ref_df['auth_seq_id'].unique())
    
    for i, auth_seq in enumerate(auth_seq_ids, 1):
        if i in original_columns:
            column_to_auth_seq[i] = auth_seq
    
    # Create a mapping from auth_seq_id to GRN label
    grn_mapping = {}
    
    # Process each transmembrane helix
    for helix_num in range(1, 8):  # 7TM proteins have 7 helices
        helix_df = ref_df[(ref_df['helix_num'] == helix_num) &
                          (ref_df['tm_helix'] == True)]
        
        if helix_df.empty:
            print(f"[WARNING] No TM helix {helix_num} found in reference structure")
            continue
        
        # Get auth_seq_ids for this helix in sequence order
        helix_auth_seqs = sorted(helix_df['auth_seq_id'].unique())
        
        if len(helix_auth_seqs) == 0:
            continue
        
        # Find the closest residue to RET as per GRN system (X.50 position)
        closest_to_ret = find_closest_ret_residue(ref_struct, helix_num)
        
        if closest_to_ret is not None and closest_to_ret in helix_auth_seqs:
            # Use this as the X.50 reference point as recommended in GRN system
            middle_auth_seq = closest_to_ret
            print(f"[INFO] Helix {helix_num}: Found residue {middle_auth_seq} (closest to RET) as {helix_num}.50")
        else:
            # Fallback to geometric center if retinal reference isn't available
            middle_idx = len(helix_auth_seqs) // 2
            middle_auth_seq = helix_auth_seqs[middle_idx]
            print(f"[INFO] Helix {helix_num}: Using geometric center residue {middle_auth_seq} as {helix_num}.50")
        
        # Map helix residues to their GRN labels
        try:
            middle_idx = helix_auth_seqs.index(middle_auth_seq)
        except ValueError:
            print(f"[WARNING] Could not find middle residue {middle_auth_seq} in helix {helix_num}")
            continue
        
        # Assign GRN labels using standard N.YY format (e.g., 1.50, 1.51, 1.52, etc.)
        # instead of Nx<YY> format to match your desired output
        for i, auth_seq in enumerate(helix_auth_seqs):
            offset = i - middle_idx
            # Use N.<YY> format for transmembrane helices
            grn_label = f"{helix_num}.{50 + offset}"
            
            # Find which column corresponds to this auth_seq_id
            for col, seq_id in column_to_auth_seq.items():
                if seq_id == auth_seq:
                    grn_mapping[col] = grn_label
                    break
    
    # Now handle loop regions and terminal regions
    # First, identify all TM helices and their boundaries
    helices = []
    for helix_num in range(1, 8):
        helix_df = ref_df[(ref_df['helix_num'] == helix_num) & (ref_df['tm_helix'] == True)]
        if not helix_df.empty:
            helix_seqs = sorted(helix_df['auth_seq_id'].unique())
            if helix_seqs:
                helices.append((helix_num, min(helix_seqs), max(helix_seqs)))
    
    # Sort helices by position
    helices.sort(key=lambda x: x[1])
    
    # Process remaining columns that aren't TM helices
    loop_counters = {}  # Track loop counters for each loop region
    
    for col in original_columns:
        if col in grn_mapping:
            continue  # Skip already assigned TM helix positions
            
        # Get the auth_seq_id for this column
        auth_seq = column_to_auth_seq.get(col)
        if auth_seq is None:
            continue
            
        # Find appropriate loop/terminal region
        if helices and auth_seq < helices[0][1]:
            # N-terminal region: use 'n.XX' format
            distance = helices[0][1] - auth_seq
            grn_mapping[col] = f"n.{distance}"
        elif helices and auth_seq > helices[-1][2]:
            # C-terminal region: use 'c.XX' format
            distance = auth_seq - helices[-1][2]
            grn_mapping[col] = f"c.{distance}"
        else:
            # Loop region between helices: use 'AB.CCC' format according to GRN system
            for i in range(len(helices)-1):
                prev_helix, _, prev_end = helices[i]
                next_helix, next_start, _ = helices[i+1]
                
                if prev_end < auth_seq < next_start:
                    # This is a loop between prev_helix and next_helix
                    
                    # Calculate distance from each helix
                    dist_to_prev = auth_seq - prev_end
                    dist_to_next = next_start - auth_seq
                    
                    # Determine which helix is closer
                    if dist_to_prev <= dist_to_next:
                        # Closer to prev_helix
                        closer_helix = prev_helix
                        further_helix = next_helix
                        distance = dist_to_prev
                    else:
                        # Closer to next_helix
                        closer_helix = next_helix
                        further_helix = prev_helix
                        distance = dist_to_next
                    
                    # Create loop key with closer helix first
                    loop_key = f"{closer_helix}{further_helix}"
                    
                    # Initialize counter if not present
                    if loop_key not in loop_counters:
                        loop_counters[loop_key] = 1
                    
                    # Format as AB.CCC where:
                    # - A is the closer helix
                    # - B is the further helix
                    # - CCC is the three-digit distance (with leading zeros)
                    grn_mapping[col] = f"{closer_helix}{further_helix}.{distance:03d}"
                    
                    # Increment counter for this loop region
                    loop_counters[loop_key] += 1
                    break
            else:
                # If not found in any loop, use a general loop label (fallback)
                if col not in grn_mapping:
                    grn_mapping[col] = f"L.{col}"
    
    print(f"[INFO] Created GRN mapping for {len(grn_mapping)} positions with proper GRN notation")
    return grn_mapping

def assign_helix_numbers_to_msa_tables(tables_dict, structure_data, reference_id):
    """
    Directly assigns helix numbers to MSA table columns based on the reference structure.
    
    Args:
        tables_dict: Dictionary with MSA tables
        structure_data: Dictionary with structure data
        reference_id: ID of the reference structure
        
    Returns:
        dict: Updated dictionary with tables having proper helix numbering
    """
    if reference_id not in structure_data:
        print(f"[ERROR] Reference structure {reference_id} not found")
        return tables_dict
    
    ref_struct = structure_data[reference_id]
    if 'df_ca_norm' not in ref_struct or 'helix_num' not in ref_struct['df_ca_norm'].columns:
        print(f"[ERROR] Reference structure lacks helix annotations")
        return tables_dict
    
    ref_df = ref_struct['df_ca_norm']
    
    # Create mapping from auth_seq_id to helix number
    auth_seq_to_helix = {}
    auth_seq_to_position = {}
    
    for helix_num in range(1, 8):
        tm_df = ref_df[(ref_df['helix_num'] == helix_num) & (ref_df['tm_helix'] == True)]
        auth_seqs = sorted(tm_df['auth_seq_id'].unique())
        
        if not auth_seqs:
            continue
        
        # Find the closest residue to retinal as the X.50 reference
        closest_to_ret = find_closest_ret_residue(ref_struct, helix_num)
        
        if closest_to_ret is not None and closest_to_ret in auth_seqs:
            middle_auth_seq = closest_to_ret
            print(f"[INFO] Helix {helix_num}: Found residue {middle_auth_seq} (closest to RET) as {helix_num}.50")
        else:
            middle_idx = len(auth_seqs) // 2
            middle_auth_seq = auth_seqs[middle_idx]
            print(f"[INFO] Helix {helix_num}: Using geometric center residue {middle_auth_seq} as {helix_num}.50")
        
        # Get the middle index
        try:
            middle_idx = auth_seqs.index(middle_auth_seq)
        except ValueError:
            print(f"[WARNING] Could not find middle residue {middle_auth_seq} in helix {helix_num}")
            middle_idx = len(auth_seqs) // 2
        
        # Map all residues in this helix with positions relative to middle (X.50)
        for i, auth_seq in enumerate(auth_seqs):
            offset = i - middle_idx
            position = f"{helix_num}.{50 + offset}"
            auth_seq_to_helix[auth_seq] = helix_num
            auth_seq_to_position[auth_seq] = position
    
    # Get the mapping from column positions to auth_seq_ids
    for table_name, table in tables_dict.items():
        if not isinstance(table, pd.DataFrame) or not hasattr(table, 'attrs'):
            continue
            
        column_to_auth_seq = table.attrs.get('column_to_auth_seq', {})
        if not column_to_auth_seq:
            print(f"[WARNING] Table {table_name} has no column_to_auth_seq mapping")
            continue
            
        # Create new column names
        new_columns = []
        
        for col in table.columns:
            auth_seq = column_to_auth_seq.get(int(col) if isinstance(col, (int, float)) else col)
            
            if auth_seq is not None and auth_seq in auth_seq_to_position:
                # This is a TM helix position
                new_columns.append(auth_seq_to_position[auth_seq])
            else:
                # Not in any TM helix or not mapped
                new_columns.append(f"X.{col}")
        
        # Apply new column names
        if len(new_columns) == len(table.columns):
            tables_dict[table_name].columns = new_columns
            print(f"[INFO] Applied helix numbering to {table_name} table")
    
    return tables_dict

def calculate_helix_distances(distance_table):
    """
    Calculate the mean and standard deviation of distances to retinal for each position,
    grouped by helix, using the GRN column names (N.YY format).
    
    Args:
        distance_table: DataFrame with distances (columns are GRN positions)
        
    Returns:
        dict: Dictionary with helix statistics
    """
    # Calculate column means and std deviations
    mean_distances = distance_table.mean(skipna=True)
    std_distances = distance_table.std(skipna=True)
    
    # Group by helix
    helix_stats = {}
    
    for position in mean_distances.index:
        pos_str = str(position)
        
        # Process TM helix positions in N.YY format
        if '.' in pos_str and len(pos_str.split('.')) == 2:
            prefix, suffix = pos_str.split('.')
            
            # Check if this is a TM helix position
            if prefix.isdigit() and int(prefix) >= 1 and int(prefix) <= 7:
                helix = prefix
                pos_num = float(suffix)
                
                if helix not in helix_stats:
                    helix_stats[helix] = {'positions': [], 'means': [], 'stds': []}
                
                helix_stats[helix]['positions'].append(position)
                helix_stats[helix]['means'].append(mean_distances[position])
                helix_stats[helix]['stds'].append(std_distances[position])
                
    # Sort positions within each helix and find closest position
    for helix, data in helix_stats.items():
        # Sort by position number
        sorted_idx = sorted(range(len(data['positions'])), 
                          key=lambda i: float(str(data['positions'][i]).split('.')[1]))
        
        data['positions'] = [data['positions'][i] for i in sorted_idx]
        data['means'] = [data['means'][i] for i in sorted_idx]
        data['stds'] = [data['stds'][i] for i in sorted_idx]
        
        # Find the position with minimum mean distance
        if data['means']:
            min_idx = data['means'].index(min(data['means']))
            data['closest_position'] = data['positions'][min_idx]
            data['closest_mean'] = data['means'][min_idx]
        else:
            data['closest_position'] = None
            data['closest_mean'] = None
    
    return helix_stats


import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt  # For saving plots if visualize_msa_distances is not used directly
import pickle  # For saving pickle files


# Assuming these functions are defined elsewhere (as in your provided msa_grn.py)
# from .your_module import create_msa_table, create_msa_distance_table, count_residues_by_helix, sort_grn_columns, assign_helix_numbers_to_msa_tables, calculate_helix_distances
# For standalone, you'd need these definitions or stubs
# Example stubs if not fully available:
# def create_msa_table(seq_alignment_dicts, processed_structures_complete, global_ref, atom_type="all"): return pd.DataFrame()
# def create_msa_distance_table(seq_alignment_dicts, processed_structures_complete, global_ref, distance_type="sidechain"): return pd.DataFrame()
# def count_residues_by_helix(df): return {'tm_total': 0, 'helices': {}, 'other': {'n': 0, 'c': 0, 'L': 0, 'X': 0}}
# def sort_grn_columns(df): return df
# def assign_helix_numbers_to_msa_tables(tables_dict, structure_data, reference_id): return tables_dict
# def calculate_helix_distances(distance_table): return {}
# from visualization_functions import visualize_msa_distances # Or plot_distances_with_std


def generate_grn_msa_tables(seq_alignment_dicts, processed_structures_complete, global_ref,
                            rmsd_df=None, max_rmsd_threshold=3.0, structure_mapping=None):
    """
    Creates all MSA tables with proper GRN (Generic Residue Numbering) column names.
    This function is a critical component for standardized analysis of membrane proteins:
    1. Creates tables showing sequence alignments and distances to retinal
    2. Applies Generic Residue Numbering (GRN) to enable direct comparison between structures
    3. Uses the "X.50" numbering convention where X is the helix number and 50 is the reference position
    4. Anchors each helix using the residue closest to retinal (or geometric center) as the X.50 position
    5. Filters structures with RMSD > threshold relative to global reference
    6. Prioritizes experimental structures over predicted ones when both are available

    Args:
        seq_alignment_dicts: Dictionary with sequence alignments
        processed_structures_complete: Dictionary with structure data
        global_ref: ID of global reference structure
        rmsd_df: Optional DataFrame with RMSD values between structures
        max_rmsd_threshold: Maximum RMSD to global reference for inclusion (default: 3.0)
        structure_mapping: Optional mapping from experimental to predicted structures

    Returns:
        dict: Dictionary with tables using GRN column names
    """
    # Filter structures by RMSD to global reference if rmsd_df is provided
    filtered_structures = processed_structures_complete.copy()
    filtered_seq_alignments = {
        'global': seq_alignment_dicts.get('global', {}).copy(),
        'type': seq_alignment_dicts.get('type', {}).copy()
    }

    excluded_structures = []

    # Step 1: Prioritize experimental structures over predicted ones
    if structure_mapping:
        print(f"[INFO] Prioritizing experimental structures over predicted ones using structure mapping")
        predicted_with_exp = set()
        for exp_id, pred_id_info in structure_mapping.items():
            pred_id = pred_id_info
            if isinstance(pred_id_info, dict) and 'predicted' in pred_id_info:
                pred_id = pred_id_info['predicted']

            if exp_id in filtered_structures and pred_id in filtered_structures:
                predicted_with_exp.add(pred_id)

        for pred_id in predicted_with_exp:
            if pred_id in filtered_structures:
                excluded_structures.append({'id': pred_id, 'reason': 'Experimental available'})
                del filtered_structures[pred_id]
                if pred_id in filtered_seq_alignments['global']:
                    del filtered_seq_alignments['global'][pred_id]
                for type_ref, type_aligns in filtered_seq_alignments['type'].items():
                    if pred_id in type_aligns:
                        del type_aligns[pred_id]
        if predicted_with_exp:
            print(
                f"[INFO] Excluded {len(predicted_with_exp)} predicted structures that have experimental counterparts.")

    # Step 2: Filter structures based on RMSD threshold
    if rmsd_df is not None and global_ref in rmsd_df.index:
        print(f"[INFO] Filtering structures with RMSD > {max_rmsd_threshold}Å to reference {global_ref}")
        rmsd_excluded_ids = []
        rmsd_values_for_excluded = {}
        for struct_id in list(filtered_structures.keys()):
            if struct_id == global_ref:
                continue
            if struct_id in rmsd_df.index and global_ref in rmsd_df.columns:
                rmsd = rmsd_df.loc[struct_id, global_ref]
                if pd.isna(rmsd) or rmsd > max_rmsd_threshold:
                    rmsd_excluded_ids.append(struct_id)
                    rmsd_values_for_excluded[struct_id] = rmsd if not pd.isna(rmsd) else "NaN"
                    if struct_id in filtered_structures:
                        del filtered_structures[struct_id]
                    if struct_id in filtered_seq_alignments['global']:
                        del filtered_seq_alignments['global'][struct_id]
                    for type_ref, type_aligns in filtered_seq_alignments['type'].items():
                        if struct_id in type_aligns: # <--- CORRECTED: Use struct_id
                            del type_aligns[struct_id] # <--- CORRECTED: Use struct_id

        for r_id in rmsd_excluded_ids:
            excluded_structures.append({'id': r_id,
                                        'reason': f'RMSD > {max_rmsd_threshold}Å ({rmsd_values_for_excluded[r_id]:.2f}Å)' if isinstance(
                                            rmsd_values_for_excluded[r_id],
                                            float) else f'RMSD > {max_rmsd_threshold}Å ({rmsd_values_for_excluded[r_id]})'})

        if rmsd_excluded_ids:
            print(f"[INFO] Excluded {len(rmsd_excluded_ids)} structures with RMSD > {max_rmsd_threshold}Å.")
            print(f"[DEBUG] High RMSD structures: " + ", ".join([
                                                                    f"{sid} ({rmsd_values_for_excluded[sid]:.2f}Å)" if isinstance(
                                                                        rmsd_values_for_excluded[sid],
                                                                        float) else f"{sid} ({rmsd_values_for_excluded[sid]})"
                                                                    for sid in rmsd_excluded_ids[:5]]) + (
                      "..." if len(rmsd_excluded_ids) > 5 else ""))

    if excluded_structures:
        print(f"[INFO] Total excluded structures: {len(excluded_structures)}")
        # print(f"[INFO] Excluded structures list: {excluded_structures}") # Can be very long

    if not filtered_structures or global_ref not in filtered_structures:
        print(
            f"[ERROR] No structures remaining after filtering, or global reference {global_ref} was excluded. Cannot proceed.")
        return {
            "residue_table": pd.DataFrame(), "distance_table": pd.DataFrame(),
            "ca_residue_table": pd.DataFrame(), "ca_distance_table": pd.DataFrame(),
            "excluded_structures": excluded_structures
        }
    print(f"[INFO] Using {len(filtered_structures)} structures for MSA table generation.")

    # Create initial tables with sequential numbering
    print("[INFO] Creating standard MSA table...")
    residue_table = create_msa_table(filtered_seq_alignments, filtered_structures, global_ref, atom_type="all")
    print("[INFO] Creating sidechain distance table...")
    distance_table = create_msa_distance_table(filtered_seq_alignments, filtered_structures, global_ref,
                                               distance_type="sidechain")
    print("[INFO] Creating CA-only MSA table...")
    ca_residue_table = create_msa_table(filtered_seq_alignments, filtered_structures, global_ref, atom_type="CA")
    print("[INFO] Creating CA-only distance table...")
    ca_distance_table = create_msa_distance_table(filtered_seq_alignments, filtered_structures, global_ref,
                                                  distance_type="backbone")

    if residue_table.empty:
        print("[ERROR] MSA table creation failed or resulted in an empty table.")
        return {
            "residue_table": pd.DataFrame(), "distance_table": pd.DataFrame(),
            "ca_residue_table": pd.DataFrame(), "ca_distance_table": pd.DataFrame(),
            "excluded_structures": excluded_structures
        }

    # Handle column count mismatches between residue_table and ca_residue_table
    if len(residue_table.columns) != len(ca_residue_table.columns):
        print(
            f"[WARNING] Column count mismatch: residue_table ({len(residue_table.columns)}) vs ca_residue_table ({len(ca_residue_table.columns)}). Reindexing CA tables.")
        all_cols = residue_table.columns
        ca_residue_table = ca_residue_table.reindex(columns=all_cols)
        ca_distance_table = ca_distance_table.reindex(columns=all_cols)  # Assuming distance table should also match

    residue_table_grn = residue_table.copy()
    distance_table_grn = distance_table.copy()
    ca_residue_table_grn = ca_residue_table.copy()
    ca_distance_table_grn = ca_distance_table.copy()

    print(f"[INFO] Reference structure for GRN: {global_ref}")
    ref_struct_data = filtered_structures[global_ref]

    # --- Dynamically get helix information from global_ref ---
    ref_df_for_helices = ref_struct_data.get('df_ca_norm', ref_struct_data.get('df_norm'))
    if ref_df_for_helices is None or 'helix_num' not in ref_df_for_helices.columns or 'auth_seq_id' not in ref_df_for_helices.columns:
        print(
            f"[ERROR] Helix annotation columns ('helix_num', 'auth_seq_id') not found in DataFrame of global reference {global_ref}.")
        return {
            "residue_table": residue_table, "distance_table": distance_table,
            "ca_residue_table": ca_residue_table, "ca_distance_table": ca_distance_table,
            "excluded_structures": excluded_structures
        }

    ref_helix_auth_seq_ranges = {}  # Stores {helix_num: [min_auth_id, max_auth_id]}
    ref_helix_residues = {}  # Stores {helix_num: [auth_id1, auth_id2, ...]} sorted
    for h_num in range(1, 8):
        helix_df = ref_df_for_helices[ref_df_for_helices['helix_num'] == h_num]
        if not helix_df.empty:
            # Ensure auth_seq_id are sortable (e.g. numeric or consistently formatted strings)
            # For now, assume they are directly comparable for min/max and sorting.
            # If they are like "10A", "10B", conversion might be needed for strict numeric sorting.
            unique_auth_ids = sorted(helix_df['auth_seq_id'].unique())
            if unique_auth_ids:
                ref_helix_auth_seq_ranges[h_num] = [unique_auth_ids[0], unique_auth_ids[-1]]
                ref_helix_residues[h_num] = unique_auth_ids
        else:
            print(
                f"[WARNING] No residues found for helix {h_num} in global reference {global_ref} based on 'helix_num' column.")

    if not ref_helix_auth_seq_ranges:
        print(
            f"[ERROR] No helix ranges could be determined from global reference {global_ref}. GRN assignment will be limited.")
        # Proceeding might still label N/C terms if any column_to_auth_seq exists.
    # --- End dynamic helix info ---

    column_to_auth_seq = residue_table.attrs.get('column_to_auth_seq', {})
    if not column_to_auth_seq:
        print("[ERROR] No 'column_to_auth_seq' mapping found in residue_table. Cannot proceed with GRN.")
        return {
            "residue_table": residue_table, "distance_table": distance_table,
            "ca_residue_table": ca_residue_table, "ca_distance_table": ca_distance_table,
            "excluded_structures": excluded_structures
        }
    auth_seq_to_column = {seq_id: col for col, seq_id in column_to_auth_seq.items()}

    # Map MSA columns to helix numbers using dynamic helix info
    column_to_helix = {}
    for h_num, auth_ids_in_helix in ref_helix_residues.items():
        for auth_id in auth_ids_in_helix:
            if auth_id in auth_seq_to_column:
                msa_col_num = auth_seq_to_column[auth_id]
                column_to_helix[msa_col_num] = h_num
    print(f"[INFO] Mapped {len(column_to_helix)} MSA columns to helix numbers based on global_ref annotations.")

    # Identify pivot columns (.50 position) for each helix
    helix_pivot_columns = {}
    for h_num in range(1, 8):
        msa_cols_in_this_helix = [col for col, helix_assignment in column_to_helix.items() if helix_assignment == h_num]
        if not msa_cols_in_this_helix:
            print(f"[WARNING] No MSA columns identified for Helix {h_num} for pivot calculation.")
            continue

        min_avg_dist = float('inf')
        pivot_col = None
        for msa_col_num in msa_cols_in_this_helix:
            if msa_col_num not in distance_table.columns: continue  # Should not happen if tables are consistent

            col_distances = distance_table[msa_col_num].dropna()
            if not col_distances.empty:
                avg_dist = col_distances.mean()
                if avg_dist < min_avg_dist:
                    min_avg_dist = avg_dist
                    pivot_col = msa_col_num

        if pivot_col is not None:
            helix_pivot_columns[h_num] = pivot_col
            auth_id_at_pivot = column_to_auth_seq.get(pivot_col, "N/A")
            print(
                f"[INFO] Helix {h_num}: Pivot at MSA col {pivot_col} (auth_seq_id: {auth_id_at_pivot}), avg_dist: {min_avg_dist:.2f}Å")
        else:
            # Fallback: geometric middle of the helix's MSA columns
            sorted_msa_cols_in_helix = sorted(msa_cols_in_this_helix,
                                              key=lambda c: list(residue_table.columns).index(c))
            if sorted_msa_cols_in_helix:
                middle_idx = len(sorted_msa_cols_in_helix) // 2
                pivot_col = sorted_msa_cols_in_helix[middle_idx]
                helix_pivot_columns[h_num] = pivot_col
                auth_id_at_pivot = column_to_auth_seq.get(pivot_col, "N/A")
                print(
                    f"[WARNING] Helix {h_num}: No distance-based pivot. Using geometric middle MSA col {pivot_col} (auth_seq_id: {auth_id_at_pivot}).")
            else:
                print(f"[WARNING] Helix {h_num}: Cannot determine pivot column.")

    # Rename columns to GRN format
    new_columns = []
    grn_column_mapping = {}  # original_col_num -> grn_label

    sorted_helices_by_ref_auth_seq = []
    if ref_helix_auth_seq_ranges:
        sorted_helices_by_ref_auth_seq = sorted(
            [(h_num, r[0], r[1]) for h_num, r in ref_helix_auth_seq_ranges.items() if r],  # only if range exists
            key=lambda x: x[1]  # Sort by min auth_seq_id of the helix
        )

    for original_col_num in residue_table.columns:
        grn_label = f"X.{original_col_num}"  # Default if no other rule applies
        auth_seq_id_at_col = column_to_auth_seq.get(original_col_num)
        assigned_helix_num = column_to_helix.get(original_col_num)

        if assigned_helix_num and assigned_helix_num in helix_pivot_columns:
            pivot_col_for_this_helix = helix_pivot_columns[assigned_helix_num]

            # Get all MSA columns belonging to this helix, in their MSA order
            msa_cols_for_assigned_helix = sorted(
                [col for col, h_assign in column_to_helix.items() if h_assign == assigned_helix_num],
                key=lambda c: list(residue_table.columns).index(c)
            )

            try:
                current_idx_in_helix_cols = msa_cols_for_assigned_helix.index(original_col_num)
                pivot_idx_in_helix_cols = msa_cols_for_assigned_helix.index(pivot_col_for_this_helix)
                offset = current_idx_in_helix_cols - pivot_idx_in_helix_cols
                grn_label = f"{assigned_helix_num}.{50 + offset}"
            except ValueError:
                print(
                    f"[WARNING] Error calculating offset for MSA col {original_col_num} in Helix {assigned_helix_num}. Defaulting GRN.")
                grn_label = f"X.{original_col_num}"

        elif auth_seq_id_at_col is not None and sorted_helices_by_ref_auth_seq:
            # Determine if N-term, C-term, or Loop based on global_ref helix boundaries
            first_helix_ref_start_id = sorted_helices_by_ref_auth_seq[0][1]
            last_helix_ref_end_id = sorted_helices_by_ref_auth_seq[-1][2]

            # Need robust comparison if auth_seq_id_at_col is string with insertion codes
            # For simplicity, direct comparison is used. Consider a helper for "is_before", "is_after"
            try:  # Attempt numeric conversion for comparison if possible, otherwise string compare
                num_auth_seq = int(re.sub("[^0-9]", "", str(auth_seq_id_at_col)))
                num_first_start = int(re.sub("[^0-9]", "", str(first_helix_ref_start_id)))
                num_last_end = int(re.sub("[^0-9]", "", str(last_helix_ref_end_id)))

                if num_auth_seq < num_first_start:
                    # Crude distance; better would be actual sequence distance in ref
                    distance = num_first_start - num_auth_seq
                    grn_label = f"n.{distance}"
                elif num_auth_seq > num_last_end:
                    distance = num_auth_seq - num_last_end
                    grn_label = f"c.{distance}"
                else:  # Loop
                    is_loop = False
                    for i in range(len(sorted_helices_by_ref_auth_seq) - 1):
                        prev_h_num, _, prev_h_ref_end_id = sorted_helices_by_ref_auth_seq[i]
                        next_h_num, next_h_ref_start_id, _ = sorted_helices_by_ref_auth_seq[i + 1]

                        num_prev_end = int(re.sub("[^0-9]", "", str(prev_h_ref_end_id)))
                        num_next_start = int(re.sub("[^0-9]", "", str(next_h_ref_start_id)))

                        if num_prev_end < num_auth_seq < num_next_start:
                            # distance from end of previous helix
                            distance = num_auth_seq - num_prev_end
                            grn_label = f"{prev_h_num}{next_h_num}.{distance:03d}"
                            is_loop = True
                            break
                    if not is_loop: grn_label = f"L.{original_col_num}"  # Fallback loop/inter-helix
            except ValueError:  # If auth_seq_ids are not purely numeric after stripping
                grn_label = f"L.{original_col_num}"  # Fallback for complex IDs

        grn_column_mapping[original_col_num] = grn_label
        new_columns.append(grn_label)

    print(f"[INFO] Generated {len(new_columns)} GRN column names.")
    if len(new_columns) != len(residue_table.columns):
        print(
            f"[ERROR] Mismatch between new GRN column count ({len(new_columns)}) and original table column count ({len(residue_table.columns)}). GRN assignment failed.")
        # Fallback: return tables with original sequential numbering
        return {
            "residue_table": residue_table, "distance_table": distance_table,
            "ca_residue_table": ca_residue_table, "ca_distance_table": ca_distance_table,
            "excluded_structures": excluded_structures, "grn_error": "Column count mismatch"
        }

    residue_table_grn.columns = new_columns
    distance_table_grn.columns = new_columns
    if len(ca_residue_table_grn.columns) == len(new_columns):
        ca_residue_table_grn.columns = new_columns
        ca_distance_table_grn.columns = new_columns
    else:
        print(
            f"[WARNING] CA tables column count ({len(ca_residue_table_grn.columns)}) differs from new GRN columns ({len(new_columns)}). CA tables will retain original sequential numbering or be reindexed if possible.")
        # Attempt reindex, otherwise they keep original columns
        try:
            ca_residue_table_grn = ca_residue_table_grn.rename(columns=grn_column_mapping)
            ca_distance_table_grn = ca_distance_table_grn.rename(columns=grn_column_mapping)
        except Exception as e:
            print(f"[WARNING] Failed to rename CA table columns with GRN: {e}. They may have sequential numbering.")

    residue_table_grn.attrs['grn_column_mapping'] = grn_column_mapping
    distance_table_grn.attrs['grn_column_mapping'] = grn_column_mapping
    # ... and for CA tables if successfully renamed

    output_dir_grn = "opsin_grn_tables"  # Consider passing as arg or using main output_dir
    os.makedirs(output_dir_grn, exist_ok=True)
    try:
        residue_table_grn.to_csv(os.path.join(output_dir_grn, "residue_table_grn.csv"))
        distance_table_grn.to_csv(os.path.join(output_dir_grn, "distance_table_grn.csv"))
        ca_residue_table_grn.to_csv(os.path.join(output_dir_grn, "ca_residue_table_grn.csv"))
        ca_distance_table_grn.to_csv(os.path.join(output_dir_grn, "ca_distance_table_grn.csv"))

        # Save mapping details
        mapping_records = []
        for orig_col, grn_val in grn_column_mapping.items():
            auth_s_id = column_to_auth_seq.get(orig_col, 'N/A')
            h_num = column_to_helix.get(orig_col, 'N/A')
            is_piv = "Yes" if h_num != 'N/A' and helix_pivot_columns.get(h_num) == orig_col else "No"
            mapping_records.append(
                {"Original_Column": orig_col, "GRN_Label": grn_val, "Auth_Seq_ID": auth_s_id, "Helix_in_Ref": h_num,
                 "Is_Pivot": is_piv})
        pd.DataFrame(mapping_records).to_csv(os.path.join(output_dir_grn, "grn_column_mapping_details.csv"),
                                             index=False)

        grn_pickle_data = {
            "residue_table_grn": residue_table_grn, "distance_table_grn": distance_table_grn,
            "ca_residue_table_grn": ca_residue_table_grn, "ca_distance_table_grn": ca_distance_table_grn,
            "grn_column_mapping_dict": grn_column_mapping,
            "column_to_auth_seq_map": column_to_auth_seq,
            "msa_column_to_ref_helix_map": column_to_helix,
            "helix_pivot_msa_columns": helix_pivot_columns
        }
        with open(os.path.join(output_dir_grn, "grn_tables_data.pkl"), "wb") as f:
            pickle.dump(grn_pickle_data, f)
        print(f"[INFO] GRN tables and metadata saved to {output_dir_grn}/")
    except Exception as e:
        print(f"[WARNING] Error saving GRN tables: {e}")

    # Final checks and sorting
    pos_counts = count_residues_by_helix(residue_table_grn)
    print(f"[INFO] GRN Table Stats: TM positions: {pos_counts['tm_total']}")
    for h, c in sorted(pos_counts['helices'].items(), key=lambda x: int(x[0])): print(f"  Helix {h}: {c} positions")

    if pos_counts['tm_total'] == 0:
        print(
            "[WARNING] No TM helix positions detected by GRN mapping logic. Consider fallback `assign_helix_numbers_to_msa_tables`.")
        # This is where you might call the fallback if desired, e.g.:
        # tables_dict = {"residue_table": residue_table_grn, ...}
        # updated_tables = assign_helix_numbers_to_msa_tables(tables_dict, filtered_structures, global_ref)
        # residue_table_grn = updated_tables["residue_table"] ... etc.

    residue_table_grn = sort_grn_columns(residue_table_grn)
    distance_table_grn = sort_grn_columns(distance_table_grn)
    ca_residue_table_grn = sort_grn_columns(ca_residue_table_grn)
    ca_distance_table_grn = sort_grn_columns(ca_distance_table_grn)

    helix_stats = calculate_helix_distances(distance_table_grn)
    ca_helix_stats = calculate_helix_distances(ca_distance_table_grn)
    print("[INFO] Closest residues to retinal in each helix (sidechain):")
    for h_num_str in sorted(helix_stats.keys(), key=int):  # Ensure numeric sort for helix numbers
        stats = helix_stats[h_num_str]
        if stats.get('closest_position'):
            print(f"  Helix {h_num_str}: Position {stats['closest_position']} (distance: {stats['closest_mean']:.2f}Å)")

    # Optional: Visualization call
    # try:
    #     from visualization_functions import plot_distances_with_std # Or your specific plotting function
    #     # ... plotting calls ...
    # except ImportError:
    #     print("[INFO] Visualization functions not available, skipping specific plots.")

    return {
        "residue_table": residue_table_grn,
        "distance_table": distance_table_grn,
        "ca_residue_table": ca_residue_table_grn,
        "ca_distance_table": ca_distance_table_grn,
        "helix_stats": helix_stats,
        "ca_helix_stats": ca_helix_stats,
        "excluded_structures": excluded_structures
    }


def fix_msa_table_grn_numbering(msa_tables, structure_data, reference_id, output_dir="./", prefix=""):
    """
    Fixes MSA tables with improper GRN numbering by directly assigning helix numbers,
    then visualizes the distances to retinal.

    Args:
        msa_tables: Dictionary with MSA tables from generate_grn_msa_tables function
        structure_data: Dictionary with structure data
        reference_id: ID of the reference structure
        output_dir: Directory to save plots
        prefix: Prefix for output filenames

    Returns:
        Dictionary with fixed tables and visualization paths
    """

    # Check if tables need fixing
    ca_distance_table = msa_tables.get("ca_distance_table")
    if ca_distance_table is None:
        print("[ERROR] No CA distance table found in input")
        return msa_tables

    # Count helix positions
    pos_counts = count_residues_by_helix(ca_distance_table)
    print(f"[INFO] Current TM residue positions: {pos_counts['tm_total']}")
    for helix, count in sorted(pos_counts['helices'].items(), key=lambda x: int(x[0])):
        print(f"[INFO] Helix {helix}: {count} positions")

    # If there are no TM helix positions, apply direct numbering
    if pos_counts['tm_total'] == 0:
        print("[INFO] No TM helix positions detected, applying direct helix numbering")
        fixed_tables = assign_helix_numbers_to_msa_tables(msa_tables, structure_data, reference_id)

        # Check if fixing was successful
        ca_distance_table = fixed_tables.get("ca_distance_table")
        if ca_distance_table is not None:
            pos_counts = count_residues_by_helix(ca_distance_table)
            print(f"[INFO] After fixing - TM residue positions: {pos_counts['tm_total']}")
            for helix, count in sorted(pos_counts['helices'].items(), key=lambda x: int(x[0])):
                print(f"[INFO] Helix {helix}: {count} positions")

            # Visualize the fixed tables
            vis_paths = visualize_msa_distances(fixed_tables, output_dir, prefix)

            # Add the fixed tables to the result
            result = {
                "fixed_tables": fixed_tables,
                "visualization_paths": vis_paths
            }
            return result

    # If tables already have helix numbering, just visualize them
    print("[INFO] Tables already have proper helix numbering")
    vis_paths = visualize_msa_distances(msa_tables, output_dir, prefix)

    return {
        "fixed_tables": msa_tables,
        "visualization_paths": vis_paths
    }
def complete_msa_tables(msa_tables, processed_structures, global_ref, seq_alignment_dicts, output_dir="./"):
    """
    Post-processes MSA tables to ensure all residues from each protein appear exactly once.
    This function enhances the GRN-labeled MSA tables by adding insertions that weren't
    aligned to the global reference structure.

    Args:
        msa_tables: Dictionary with MSA tables from generate_grn_msa_tables function
        processed_structures: Dictionary with structure data
        global_ref: ID of global reference structure
        seq_alignment_dicts: Dictionary with sequence alignments
        output_dir: Directory to save output files

    Returns:
        Dictionary with completed MSA tables
    """
    print("[INFO] Post-processing MSA tables to include all residues...")

    # Get the required tables
    residue_table_grn = msa_tables.get("residue_table", pd.DataFrame())
    ca_residue_table_grn = msa_tables.get("ca_residue_table", pd.DataFrame())

    if residue_table_grn.empty:
        print("[ERROR] No residue table found in input")
        return msa_tables

    # Create completed tables
    complete_residue_table = complete_msa_table_grn(processed_structures, residue_table_grn, global_ref, seq_alignment_dicts)

    # Check if CA table exists and process it
    if not ca_residue_table_grn.empty:
        complete_ca_residue_table = complete_msa_table_grn(processed_structures, ca_residue_table_grn, global_ref, seq_alignment_dicts, atom_type="CA")
    else:
        complete_ca_residue_table = pd.DataFrame()

    # Create new distance tables based on the completed residue tables
    print("[INFO] Calculating distances for the completed MSA tables...")
    complete_distance_table = pd.DataFrame()
    complete_ca_distance_table = pd.DataFrame()

    # Calculate distances for the new columns in the complete residue table
    if not complete_residue_table.empty:
        # Start with the existing distance table
        complete_distance_table = msa_tables.get("distance_table", pd.DataFrame()).copy()

        # Ensure all columns from the complete residue table exist
        for col in complete_residue_table.columns:
            if col not in complete_distance_table.columns:
                complete_distance_table[col] = np.nan

        # Fill in distances for new columns
        complete_distance_table = calculate_distances_for_missing_columns(
            complete_residue_table,
            complete_distance_table,
            processed_structures,
            atom_type="all"
        )

    # Calculate distances for the new columns in the complete CA residue table
    if not complete_ca_residue_table.empty:
        # Start with the existing CA distance table
        complete_ca_distance_table = msa_tables.get("ca_distance_table", pd.DataFrame()).copy()

        # Ensure all columns from the complete CA residue table exist
        for col in complete_ca_residue_table.columns:
            if col not in complete_ca_distance_table.columns:
                complete_ca_distance_table[col] = np.nan

        # Fill in distances for new columns
        complete_ca_distance_table = calculate_distances_for_missing_columns(
            complete_ca_residue_table,
            complete_ca_distance_table,
            processed_structures,
            atom_type="CA"
        )

    # Save the completed tables
    if not complete_residue_table.empty:
        complete_residue_table.to_csv(os.path.join(output_dir, "residue_table_grn_complete.csv"))
    if not complete_ca_residue_table.empty:
        complete_ca_residue_table.to_csv(os.path.join(output_dir, "ca_residue_table_grn_complete.csv"))
    if not complete_distance_table.empty:
        complete_distance_table.to_csv(os.path.join(output_dir, "distance_table_grn_complete.csv"))
    if not complete_ca_distance_table.empty:
        complete_ca_distance_table.to_csv(os.path.join(output_dir, "ca_distance_table_grn_complete.csv"))

    # Report on additions
    orig_cols = len(residue_table_grn.columns)
    new_cols = len(complete_residue_table.columns)
    added_cols = new_cols - orig_cols

    print(f"[INFO] Added {added_cols} new columns for insertions")
    print(f"[INFO] Original MSA table: {orig_cols} columns, Complete MSA table: {new_cols} columns")
    print(f"[INFO] Completed MSA tables saved to {output_dir}")

    # Return the completed tables in the same format as input
    completed_tables = msa_tables.copy()
    completed_tables["residue_table"] = complete_residue_table
    completed_tables["ca_residue_table"] = complete_ca_residue_table
    completed_tables["distance_table"] = complete_distance_table
    completed_tables["ca_distance_table"] = complete_ca_distance_table

    return completed_tables

def complete_msa_table_grn(processed_structures, msa_table_grn, global_ref, seq_alignment_dicts, atom_type="all"):
    """
    Create a complete MSA table that includes all residues from all proteins.

    Args:
        processed_structures: Dictionary with structure data
        msa_table_grn: Existing MSA table with GRN columns
        global_ref: ID of global reference structure
        seq_alignment_dicts: Dictionary with alignment information
        atom_type: Type of atoms to include ('all' or 'CA')

    Returns:
        DataFrame: Complete MSA table with all residues
    """
    # Return the original table if it's empty
    if msa_table_grn.empty:
        return msa_table_grn

    # 1. Start with the existing GRN-labeled MSA table
    complete_msa_table = msa_table_grn.copy()

    # Track all new columns that will be added
    new_columns = {}  # {new_grn_label: original_grn_label} mapping
    insertion_count = 0

    # 2. Extract all residues for each structure
    for struct_id, struct_data in processed_structures.items():
        if struct_id not in complete_msa_table.index:
            continue  # Skip structures not in the MSA table

        # Get appropriate dataframe based on atom_type
        if atom_type == "CA":
            if 'df_ca_norm' in struct_data:
                struct_df = struct_data['df_ca_norm']
                struct_df = struct_df[struct_df['res_atom_name'] == 'CA']
            else:
                continue  # Skip if no CA dataframe
        else:
            if 'df_norm' in struct_data:
                struct_df = struct_data['df_norm']
                struct_df = struct_df.drop_duplicates(subset=['auth_seq_id', 'auth_chain_id'])
            else:
                continue  # Skip if no normalized dataframe

        # Get all residue IDs for this structure
        all_residue_ids = sorted(struct_df['auth_seq_id'].unique())

        # Find residues already in the MSA table
        mapped_residues = {}  # {residue_id: column} mapping
        for col in complete_msa_table.columns:
            cell_value = complete_msa_table.at[struct_id, col]
            if cell_value != '-' and isinstance(cell_value, str):
                # Extract residue ID from cell value (format: 'A123')
                try:
                    digit_part = ''.join(c for c in cell_value if c.isdigit())
                    if digit_part:
                        residue_id = int(digit_part)
                        mapped_residues[residue_id] = col
                except:
                    pass

        # Identify unmapped residues (insertions relative to reference)
        unmapped_residues = [res_id for res_id in all_residue_ids if res_id not in mapped_residues]

        if not unmapped_residues:
            continue  # All residues are already in the MSA table

        print(f"[INFO] Found {len(unmapped_residues)} unmapped residues in structure {struct_id}")

        # 3. For each unmapped residue, determine where it should be inserted
        for unmapped_res in unmapped_residues:
            # Find the nearest mapped residues (before and after)
            prev_res = None
            prev_col = None
            next_res = None
            next_col = None

            for mapped_res, mapped_col in sorted(mapped_residues.items()):
                if mapped_res < unmapped_res:
                    prev_res = mapped_res
                    prev_col = mapped_col
                elif mapped_res > unmapped_res:
                    next_res = mapped_res
                    next_col = mapped_col
                    break

            # Determine which mapped residue is closer
            base_col = None
            if prev_col is not None and next_col is not None:
                # Choose the closer one
                if unmapped_res - prev_res <= next_res - unmapped_res:
                    base_col = prev_col
                    base_res = prev_res
                else:
                    base_col = next_col
                    base_res = next_res
            elif prev_col is not None:
                base_col = prev_col
                base_res = prev_res
            elif next_col is not None:
                base_col = next_col
                base_res = next_res
            else:
                # No mapped residues found, skip this unmapped residue
                continue

            # Create new GRN label with decimal extension of the base column
            try:
                if '.' in base_col:
                    # Standard GRN format (N.XX)
                    helix, pos = base_col.split('.')
                    pos_num = float(pos)

                    # Create insertion label with small decimal increment
                    # If insertion is after base position, add small increment
                    # If insertion is before, subtract small increment
                    if base_col == prev_col:
                        # Insertion after base position
                        increment = 0.001 * (unmapped_res - prev_res)
                        new_pos_num = pos_num + increment
                    else:
                        # Insertion before base position
                        increment = 0.001 * (next_res - unmapped_res)
                        new_pos_num = pos_num - increment

                    # Format with 3 decimal places
                    new_grn_label = f"{helix}.{new_pos_num:.3f}"
                else:
                    # Handle non-standard GRN format
                    new_grn_label = f"I.{unmapped_res}"
            except:
                # Fallback for any errors
                new_grn_label = f"I.{unmapped_res}"

            # Create appropriate residue name string
            res_rows = struct_df[struct_df['auth_seq_id'] == unmapped_res]
            if not res_rows.empty and 'res_name1l' in res_rows.columns:
                res_name = res_rows['res_name1l'].iloc[0]
            else:
                res_name = 'X'  # Unknown residue type

            residue_str = f"{res_name}{unmapped_res}"

            # Add column if it doesn't exist
            if new_grn_label not in complete_msa_table.columns:
                complete_msa_table[new_grn_label] = '-'
                new_columns[new_grn_label] = base_col  # Track which base column this is related to
                insertion_count += 1

            # Update this cell with the residue
            complete_msa_table.at[struct_id, new_grn_label] = residue_str

    print(f"[INFO] Added {insertion_count} new columns for insertions across all structures")

    # 4. Sort columns to maintain proper GRN order
    complete_msa_table = sort_grn_columns(complete_msa_table)

    # Return the expanded table with all residues
    return complete_msa_table

def calculate_distances_for_missing_columns(residue_table, distance_table, processed_structures, atom_type="all"):
    """
    Calculate distances to retinal for columns in the residue table that are missing from the distance table.
    
    Args:
        residue_table: MSA table with residue information
        distance_table: Distance table (potentially missing some columns)
        processed_structures: Dictionary with structure data
        atom_type: Type of atoms to include ('all' or 'CA')
        
    Returns:
        DataFrame: Updated distance table with distances for all columns
    """
    from scipy.spatial.distance import cdist
    
    # Ensure the distance table has all columns from the residue table
    for col in residue_table.columns:
        if col not in distance_table.columns:
            distance_table[col] = np.nan
    
    # Process each structure
    for struct_id in residue_table.index:
        if struct_id not in processed_structures:
            continue
            
        struct = processed_structures[struct_id]
        
        # Get appropriate dataframe
        if atom_type == "CA":
            if 'df_ca_norm' in struct:
                df_norm = struct['df_ca_norm']
                df_norm = df_norm[df_norm['res_atom_name'] == 'CA']
            else:
                continue
        else:
            if 'df_norm' in struct:
                df_norm = struct['df_norm']
            else:
                continue
        
        # Find RET (retinal) atoms
        ret_df = df_norm[df_norm['res_name3l'] == 'RET']
        if ret_df.empty:
            # Try LIG as fallback
            ret_df = df_norm[df_norm['res_name3l'] == 'LIG']
            if ret_df.empty:
                # No retinal found, skip this structure
                continue
        
        # Get retinal coordinates
        ret_coords = ret_df[['x', 'y', 'z']].values
        
        # Process each column
        for col in residue_table.columns:
            # Skip columns that already have a distance value
            if pd.notna(distance_table.at[struct_id, col]):
                continue
                
            # Get the residue ID from the residue table
            cell_value = residue_table.at[struct_id, col]
            if cell_value == '-' or not isinstance(cell_value, str):
                continue
                
            # Extract residue ID
            try:
                digit_part = ''.join(c for c in cell_value if c.isdigit())
                if not digit_part:
                    continue
                    
                residue_id = int(digit_part)
                
                # Find the residue coordinates
                res_df = df_norm[df_norm['auth_seq_id'] == residue_id]
                if res_df.empty:
                    continue
                    
                res_coords = res_df[['x', 'y', 'z']].values
                
                # Calculate distance
                distances = cdist(res_coords, ret_coords)
                min_distance = float(distances.min())
                
                # Update distance table
                distance_table.at[struct_id, col] = min_distance
            except Exception as e:
                print(f"[WARNING] Error calculating distance for {struct_id}, column {col}: {e}")
                continue
    
    return distance_table
