"""
Functions for multiple structure alignment and Generic Residue Numbering (GRN).
These functions create and analyze multiple sequence alignments with GRN position labels.
"""

import numpy as np
import pandas as pd
import re
import os
import json
from pathlib import Path
from scipy.spatial.distance import cdist
import matplotlib.pyplot as plt

from src.visualization_functions import visualize_msa_distances


def load_helix_boundaries(structure_id, helices_file='property/helices_curated.json'):
    """
    Load helix boundaries for a given structure from the helices_curated.json file.
    
    Args:
        structure_id: ID of the structure to get helix boundaries for
        helices_file: Path to the JSON file containing helix boundaries
        
    Returns:
        List of [start, end] positions for each helix, or None if not found
    """
    helices_path = Path(helices_file)
    
    # Try different possible paths
    if not helices_path.exists():
        # Try from project root
        project_root = Path(__file__).parent.parent
        helices_path = project_root / helices_file
        
    if not helices_path.exists():
        print(f"[ERROR] Helices file not found: {helices_file}")
        return None
        
    try:
        with open(helices_path, 'r') as f:
            helices_data = json.load(f)
            
        if structure_id not in helices_data:
            print(f"[ERROR] Structure {structure_id} not found in helices file")
            return None
            
        # Convert from dict format to list format expected by the code
        structure_helices = helices_data[structure_id]
        tm_ranges = []
        
        # Ensure helices are in order 1-7
        for helix_num in range(1, 8):
            helix_key = str(helix_num)
            if helix_key in structure_helices:
                tm_ranges.append(structure_helices[helix_key])
            else:
                print(f"[WARNING] Helix {helix_num} not found for structure {structure_id}")
                return None
                
        print(f"[INFO] Loaded helix boundaries for {structure_id} from {helices_file}")
        for i, (start, end) in enumerate(tm_ranges, 1):
            print(f"[INFO] Helix {i}: Auth_seq_id range {start}-{end}")
            
        return tm_ranges
        
    except Exception as e:
        print(f"[ERROR] Error loading helix boundaries: {e}")
        return None


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
    
    Uses Protos GRN sorting if available, otherwise falls back to custom sorting.
    
    The order follows the GRN system:
    1. N-terminal (n.XX) 
    2. TM helices in order (1.YY through 7.YY, with positions sorted numerically)
    3. Loop regions (AB.CCC)
    4. C-terminal regions (c.XX)
    5. Any other labels
    
    Args:
        df: DataFrame with GRN column labels
        
    Returns:
        DataFrame: DataFrame with sorted columns
    """
    # Try to use Protos GRN sorting if available
    try:
        from protos.grn.grn_utils import sort_grns_str
        
        # Convert column names to strings
        column_names = [str(col) for col in df.columns]
        
        # Sort using Protos GRN sorting
        sorted_columns = sort_grns_str(column_names, output_notation_type='dot')
        
        # Return DataFrame with sorted columns
        return df[sorted_columns]
        
    except ImportError:
        print("[INFO] Protos GRN sorting not available, using fallback sorting")
        # Fall back to custom sorting logic below
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
        match_loop = re.match(r'(\d)(\d)\.(\d+)', label)
        if match_loop:
            closer_helix, further_helix, distance = match_loop.groups()
            # Position loop after the first helix
            # e.g., 12.003 should come after helix 1 positions
            return float(closer_helix) + 0.9 + (0.0001 * float(distance))
        
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

def generate_grn_msa_tables(seq_alignment_dicts, processed_structures_complete, global_ref,
                         rmsd_df=None, max_rmsd_threshold=3.0, structure_mapping=None,
                         helices_file='property/helices_curated.json'):
    """
    Creates all MSA tables with proper GRN (Generic Residue Numbering) column names.
    This function is a critical component for standardized analysis of membrane proteins:
    1. Creates tables showing sequence alignments and distances to retinal
    2. Applies Generic Residue Numbering (GRN) to enable direct comparison between structures
    3. Uses the "X.50" numbering convention where X is the helix number and 50 is the reference position
    4. Anchors each helix using the residue closest to retinal as the X.50 position
    5. Filters structures with RMSD > threshold relative to global reference
    6. Prioritizes experimental structures over predicted ones when both are available

    Args:
        seq_alignment_dicts: Dictionary with sequence alignments
        processed_structures_complete: Dictionary with structure data
        global_ref: ID of global reference structure
        rmsd_df: Optional DataFrame with RMSD values between structures
        max_rmsd_threshold: Maximum RMSD to global reference for inclusion (default: 3.0)
        structure_mapping: Optional mapping from experimental to predicted structures
        helices_file: Path to JSON file containing helix boundaries (default: 'property/helices_curated.json')

    Returns:
        dict: Dictionary with tables using GRN column names
    """
    # Filter structures by RMSD to global reference if rmsd_df is provided
    filtered_structures = processed_structures_complete.copy()
    filtered_seq_alignments = {
        'global': seq_alignment_dicts.get('global', {}).copy(),
        'type': seq_alignment_dicts.get('type', {}).copy()
    }

    # Initialize excluded_structures for tracking filtered structures
    excluded_structures = []

    # Step 1: First prioritize experimental structures over predicted ones
    if structure_mapping:
        print(f"[INFO] Prioritizing experimental structures over predicted ones using structure mapping")

        # Identify predicted structures that have experimental counterparts
        predicted_with_exp = set()
        for exp_id, pred_id in structure_mapping.items():
            # Handle different mapping formats
            if isinstance(pred_id, dict) and 'predicted' in pred_id:
                pred_id = pred_id['predicted']

            # If both experimental and predicted structures exist, prioritize experimental
            if exp_id in filtered_structures and pred_id in filtered_structures:
                predicted_with_exp.add(pred_id)

        # Remove predicted structures that have experimental counterparts
        for pred_id in predicted_with_exp:
            if pred_id in filtered_structures:
                excluded_structures.append(pred_id)
                del filtered_structures[pred_id]
                # Also remove from sequence alignments
                if pred_id in filtered_seq_alignments['global']:
                    del filtered_seq_alignments['global'][pred_id]
                # Also check in type alignments
                for type_ref, type_aligns in filtered_seq_alignments['type'].items():
                    if pred_id in type_aligns:
                        del type_aligns[pred_id]

        print(f"[INFO] Excluded {len(predicted_with_exp)} predicted structures that have experimental counterparts")

    # Step 2: Filter structures based on RMSD threshold
    if rmsd_df is not None and global_ref in rmsd_df.index:
        print(f"[INFO] Filtering structures with RMSD > {max_rmsd_threshold}Å to reference {global_ref}")

        # Get structures with RMSD > threshold
        rmsd_excluded = []
        rmsd_values = {}
        for struct_id in list(filtered_structures.keys()):  # Use list to avoid modification during iteration
            if struct_id == global_ref:
                continue  # Always include the reference

            if struct_id in rmsd_df.index and global_ref in rmsd_df.columns:
                rmsd = rmsd_df.loc[struct_id, global_ref]
                if pd.isna(rmsd) or rmsd > max_rmsd_threshold:
                    rmsd_excluded.append(struct_id)
                    rmsd_values[struct_id] = rmsd if not pd.isna(rmsd) else "NaN"
                    # Remove from filtered structures
                    if struct_id in filtered_structures:
                        del filtered_structures[struct_id]
                    # Remove from sequence alignments
                    if struct_id in filtered_seq_alignments['global']:
                        del filtered_seq_alignments['global'][struct_id]
                    # Also check in type alignments
                    for type_ref, type_aligns in filtered_seq_alignments['type'].items():
                        if struct_id in type_aligns:
                            del type_aligns[struct_id]

        # Add RMSD-excluded structures to the overall excluded list
        excluded_structures.extend(rmsd_excluded)

        print(f"[INFO] Excluded {len(rmsd_excluded)} structures with RMSD > {max_rmsd_threshold}Å")
        # Debug: Print the excluded structures with their RMSD values
        if rmsd_excluded:
            print(f"[DEBUG] High RMSD structures: " + ", ".join([f"{sid} ({rmsd_values[sid]:.2f}Å)" if isinstance(rmsd_values[sid], float) else f"{sid} ({rmsd_values[sid]})" for sid in rmsd_excluded]))

        # Show overall exclusion statistics
        print(f"[INFO] Total excluded structures: {len(excluded_structures)}")
        if excluded_structures:
            print(f"[INFO] Excluded structures: {', '.join(excluded_structures[:10])}" +
                  (f"... and {len(excluded_structures)-10} more" if len(excluded_structures) > 10 else ""))
        print(f"[INFO] Using {len(filtered_structures)} structures for MSA table generation")

    # First create the tables with sequential numbering
    print("[INFO] Creating standard MSA table...")
    residue_table = create_msa_table(filtered_seq_alignments, filtered_structures, global_ref)

    print("[INFO] Creating sidechain distance table...")
    distance_table = create_msa_distance_table(
        filtered_seq_alignments, filtered_structures, global_ref, distance_type="sidechain"
    )

    print("[INFO] Creating CA-only MSA table...")
    ca_residue_table = create_msa_table(
        filtered_seq_alignments, filtered_structures, global_ref, atom_type="CA"
    )

    print("[INFO] Creating CA-only distance table...")
    ca_distance_table = create_msa_distance_table(
        filtered_seq_alignments, filtered_structures, global_ref, distance_type="backbone"
    )

    # Check if all tables have the same columns
    print(f"[DEBUG] Table column counts:")
    print(f"  - residue_table: {len(residue_table.columns)}")
    print(f"  - distance_table: {len(distance_table.columns)}")
    print(f"  - ca_residue_table: {len(ca_residue_table.columns)}")
    print(f"  - ca_distance_table: {len(ca_distance_table.columns)}")

    # Check for differences between regular and CA tables
    if len(residue_table.columns) != len(ca_residue_table.columns):
        print(f"[WARNING] Column count mismatch between residue_table and ca_residue_table")

        # Find which columns are different
        reg_cols = set(residue_table.columns)
        ca_cols = set(ca_residue_table.columns)

        if len(reg_cols - ca_cols) > 0:
            print(f"[DEBUG] Columns in residue_table but not in ca_residue_table: {sorted(list(reg_cols - ca_cols))}")

        if len(ca_cols - reg_cols) > 0:
            print(f"[DEBUG] Columns in ca_residue_table but not in residue_table: {sorted(list(ca_cols - reg_cols))}")

        # Ensure CA tables have the same columns as regular tables for consistent GRN mapping
        print(f"[INFO] Making CA tables consistent with regular tables")

        # Create a consistent column set - use the regular table's columns as the standard
        all_cols = residue_table.columns

        # Reindex CA tables to match the regular tables
        # This adds NaN for missing columns and removes extra columns
        ca_residue_table = ca_residue_table.reindex(columns=all_cols)
        ca_distance_table = ca_distance_table.reindex(columns=all_cols)

        print(f"[INFO] CA tables reindexed to match regular tables")
        print(f"  - ca_residue_table: {len(ca_residue_table.columns)} columns")
        print(f"  - ca_distance_table: {len(ca_distance_table.columns)} columns")

    # Create copies of the tables that will have GRN column names
    residue_table_grn = residue_table.copy()
    distance_table_grn = distance_table.copy()
    ca_residue_table_grn = ca_residue_table.copy()
    ca_distance_table_grn = ca_distance_table.copy()

    print(f"[INFO] Reference structure: {global_ref}")

    # Get the reference structure
    ref_struct = filtered_structures[global_ref]
    ref_df = ref_struct['df_ca_norm']

    # Load helix boundaries from helices_curated.json
    tm_ranges = load_helix_boundaries(global_ref, helices_file)
    
    if tm_ranges is None:
        print(f"[ERROR] Could not load helix boundaries for reference {global_ref}")
        return {
            "residue_table": pd.DataFrame(),
            "distance_table": pd.DataFrame(),
            "ca_residue_table": pd.DataFrame(),
            "ca_distance_table": pd.DataFrame(),
            "helix_stats": {},
            "ca_helix_stats": {},
            "excluded_structures": []
        }

    # Get the column to auth_seq_id mapping
    column_to_auth_seq = residue_table.attrs.get('column_to_auth_seq', {})
    if not column_to_auth_seq:
        print("[ERROR] No column_to_auth_seq mapping found in residue_table")
        return {}

    # Create a reverse mapping from auth_seq_id to column
    auth_seq_to_column = {seq_id: col for col, seq_id in column_to_auth_seq.items()}

    # Step 1: Map columns in the GRN table to the corresponding helix
    column_to_helix = {}
    for helix_num, (start, end) in enumerate(tm_ranges, 1):
        print(f"[INFO] Mapping Helix {helix_num}: Auth_seq_id range {start}-{end}")
        for auth_seq_id in range(start, end + 1):
            # Find the column that corresponds to this auth_seq_id
            if auth_seq_id in auth_seq_to_column:
                column = auth_seq_to_column[auth_seq_id]
                column_to_helix[column] = helix_num

    print(f"[INFO] Mapped {len(column_to_helix)} columns to helix numbers")

    # Step 2: Calculate the closest overall position to retinal for each helix
    helix_pivot_columns = {}

    for helix_num in range(1, 8):
        # Get all columns for this helix
        helix_columns = [col for col, h_num in column_to_helix.items() if h_num == helix_num]
        if not helix_columns:
            print(f"[WARNING] No columns found for Helix {helix_num}")
            continue

        # Calculate mean distance to retinal for each position in this helix
        min_distance = float('inf')
        pivot_column = None

        for column in helix_columns:
            # Get distances to retinal for this position
            distances = distance_table[column].dropna()
            if len(distances) == 0:
                continue

            mean_distance = distances.mean()
            if mean_distance < min_distance:
                min_distance = mean_distance
                pivot_column = column

        if pivot_column is not None:
            helix_pivot_columns[helix_num] = pivot_column
            auth_seq_id = column_to_auth_seq.get(pivot_column)
            print(f"[INFO] Helix {helix_num}: Found pivot at column {pivot_column} (auth_seq_id: {auth_seq_id}) with mean distance {min_distance:.2f}Å")
        else:
            print(f"[WARNING] Could not find pivot for Helix {helix_num}")

    # Step 3 & 4: Rename columns to format <helix>.<position> with pivot as X.50
    new_columns = []
    column_mapping = {}

    for column in residue_table.columns:
        helix_num = column_to_helix.get(column)

        if helix_num is not None and helix_num in helix_pivot_columns:
            # This is a TM helix position
            pivot_column = helix_pivot_columns[helix_num]
            offset = list(residue_table.columns).index(column) - list(residue_table.columns).index(pivot_column)
            position = 50 + offset
            grn_label = f"{helix_num}.{position}"
            column_mapping[column] = grn_label
            new_columns.append(grn_label)
        else:
            # Handle non-TM positions
            auth_seq_id = column_to_auth_seq.get(column)

            # Determine region based on auth_seq_id
            if auth_seq_id is not None:
                # Check if it's before first helix
                if auth_seq_id < tm_ranges[0][0]:
                    distance = tm_ranges[0][0] - auth_seq_id
                    grn_label = f"n.{distance}"
                # Check if it's after last helix
                elif auth_seq_id > tm_ranges[-1][1]:
                    distance = auth_seq_id - tm_ranges[-1][1]
                    grn_label = f"c.{distance}"
                else:
                    # It's in a loop region - find the closest helices
                    for i in range(len(tm_ranges) - 1):
                        prev_end = tm_ranges[i][1]
                        next_start = tm_ranges[i+1][0]

                        if prev_end < auth_seq_id < next_start:
                            # It's in a loop between helix i+1 and i+2
                            prev_helix = i + 1
                            next_helix = i + 2

                            # Determine which helix is closer
                            dist_to_prev = auth_seq_id - prev_end
                            dist_to_next = next_start - auth_seq_id

                            if dist_to_prev <= dist_to_next:
                                grn_label = f"{prev_helix}{next_helix}.{dist_to_prev:03d}"
                            else:
                                grn_label = f"{next_helix}{prev_helix}.{dist_to_next:03d}"
                            break
                    else:
                        # Fallback for positions not identified as loops
                        grn_label = f"X.{column}"
            else:
                # Columns without auth_seq_id mapping
                grn_label = f"X.{column}"

            column_mapping[column] = grn_label
            new_columns.append(grn_label)

    print(f"[INFO] Created new column names with GRN notation")

    # Verify we have the right number of column names
    if len(new_columns) != len(residue_table.columns):
        print(f"[ERROR] Column count mismatch: Generated {len(new_columns)} names for {len(residue_table.columns)} columns")
        # Fix the column count if needed
        while len(new_columns) < len(residue_table.columns):
            idx = len(new_columns)
            new_columns.append(f"X.{idx}")
            print(f"[WARNING] Added placeholder column name X.{idx}")
        if len(new_columns) > len(residue_table.columns):
            print(f"[WARNING] Trimming excess column names")
            print(f"[DEBUG] New columns length: {len(new_columns)}, residue_table columns: {len(residue_table.columns)}")
            print(f"[DEBUG] First few excess columns: {new_columns[len(residue_table.columns):len(residue_table.columns)+5]}")
            new_columns = new_columns[:len(residue_table.columns)]

    # Print column counts for debugging
    print(f"[DEBUG] Column counts before renaming:")
    print(f"  - residue_table_grn: {len(residue_table_grn.columns)}")
    print(f"  - distance_table_grn: {len(distance_table_grn.columns)}")
    print(f"  - ca_residue_table_grn: {len(ca_residue_table_grn.columns)}")
    print(f"  - ca_distance_table_grn: {len(ca_distance_table_grn.columns)}")
    print(f"  - new_columns: {len(new_columns)}")

    # Apply the new column names
    residue_table_grn.columns = new_columns
    distance_table_grn.columns = new_columns

    # Check if CA tables have the same column count
    if len(ca_residue_table_grn.columns) != len(new_columns):
        print(f"[ERROR] CA residue table has {len(ca_residue_table_grn.columns)} columns, but new_columns has {len(new_columns)}")
        print(f"[ERROR] This mismatch is causing the error. Checking columns...")

        # Check which columns might be different
        if len(ca_residue_table_grn.columns) > len(new_columns):
            print(f"[DEBUG] CA table has EXTRA columns. First extra column: {ca_residue_table_grn.columns[len(new_columns)]}")
        else:
            print(f"[DEBUG] CA table has FEWER columns than new_columns")

        # For safety, create specific column names for CA tables matching their actual column count
        ca_columns = new_columns.copy()
        if len(ca_residue_table_grn.columns) > len(ca_columns):
            # Need to add columns
            for i in range(len(ca_columns), len(ca_residue_table_grn.columns)):
                ca_columns.append(f"CA_extra_{i}")
        elif len(ca_residue_table_grn.columns) < len(ca_columns):
            # Need to trim columns
            ca_columns = ca_columns[:len(ca_residue_table_grn.columns)]

        # Apply specific column names to CA tables
        ca_residue_table_grn.columns = ca_columns
        ca_distance_table_grn.columns = ca_columns

        print(f"[WARNING] CA tables and main tables now have different column names")
    else:
        # Normal case - all tables have the same number of columns
        ca_residue_table_grn.columns = new_columns
        ca_distance_table_grn.columns = new_columns

    print(f"[DEBUG] New columns successfully applied: {len(new_columns)} columns total")

    # Debug output to check GRN column names
    print(f"[DEBUG] GRN mapping generated {len(column_mapping)} position mappings")

    # Store the column mapping in the tables for reference
    residue_table_grn.attrs['column_mapping'] = column_mapping
    distance_table_grn.attrs['column_mapping'] = column_mapping
    ca_residue_table_grn.attrs['column_mapping'] = column_mapping
    ca_distance_table_grn.attrs['column_mapping'] = column_mapping

    # Save the tables immediately after GRN column assignment (optional debugging)
    # This section saves additional debugging information
    debug_save = False  # Set to True to save debugging files
    if debug_save:
        try:
            print("[INFO] Saving MSA debugging tables with GRN column names...")
            debug_output_dir = "opsin_output/debug_grn_tables"
            os.makedirs(debug_output_dir, exist_ok=True)

            # Save the tables with GRN column names
            residue_table_grn.to_csv(f"{debug_output_dir}/residue_table_grn.csv")
            distance_table_grn.to_csv(f"{debug_output_dir}/distance_table_grn.csv")
            ca_residue_table_grn.to_csv(f"{debug_output_dir}/ca_residue_table_grn.csv")
            ca_distance_table_grn.to_csv(f"{debug_output_dir}/ca_distance_table_grn.csv")

            # Create a column mapping DataFrame
            mapping_records = []
            original_columns = list(residue_table.columns)

            for i, (orig_col, grn_label) in enumerate(zip(original_columns, new_columns)):
                auth_seq_id = column_to_auth_seq.get(orig_col, "unknown")
                helix_num = column_to_helix.get(orig_col, "")

                record = {
                    "Original_Column": orig_col,
                    "GRN_Label": grn_label,
                    "Auth_Seq_ID": auth_seq_id,
                    "Helix": helix_num if helix_num else "",
                    "Is_Pivot": "Yes" if orig_col in helix_pivot_columns.values() else "No"
                }
                mapping_records.append(record)

            # Save the column mapping as CSV
            mapping_df = pd.DataFrame(mapping_records)
            mapping_df.to_csv(f"{debug_output_dir}/column_mapping.csv", index=False)

            # Save a CSV with both original and GRN column names for each table
            for name, table, prefix in [
                ("Residue Table", residue_table, "residue"),
                ("Distance Table", distance_table, "distance"),
                ("CA Residue Table", ca_residue_table, "ca_residue"),
                ("CA Distance Table", ca_distance_table, "ca_distance")
            ]:
                # Create a copy with both original and GRN column names
                dual_table = table.copy()

                # Create new column names with both original and GRN
                dual_columns = [f"{orig}|{grn}" for orig, grn in zip(original_columns, new_columns)]
                dual_table.columns = dual_columns

                # Save with dual column names
                dual_table.to_csv(f"{debug_output_dir}/{prefix}_table_dual.csv")

            # Also save a pickle version with full metadata
            import pickle
            grn_tables = {
                "residue_table": residue_table_grn,
                "distance_table": distance_table_grn,
                "ca_residue_table": ca_residue_table_grn,
                "ca_distance_table": ca_distance_table_grn,
                "column_mapping": column_mapping,
                "column_to_auth_seq": column_to_auth_seq,
                "column_to_helix": column_to_helix,
                "helix_pivot_columns": helix_pivot_columns
            }
            with open(f"{debug_output_dir}/grn_tables.pkl", "wb") as f:
                pickle.dump(grn_tables, f)

            print(f"[INFO] Debug tables saved to {debug_output_dir}/")
            print(f"[INFO] Column mapping saved to {debug_output_dir}/column_mapping.csv")
        except Exception as e:
            print(f"[WARNING] Could not save debug GRN tables: {e}")
    
    print(f"[DEBUG] Sample column names: {', '.join(str(col) for col in list(residue_table_grn.columns)[:10])}")

    # Count columns by format
    format_counts = {'helix_format': 0, 'n_terminal': 0, 'c_terminal': 0, 'loop': 0, 'unassigned': 0}
    helix_count = {}

    for col in residue_table_grn.columns:
        col_str = str(col)

        # Primary focus on N.YY format (dot notation for helices)
        if '.' in col_str and len(col_str.split('.')) == 2:
            prefix, suffix = col_str.split('.')
            if prefix.isdigit() and int(prefix) >= 1 and int(prefix) <= 7:
                # This is a TM helix position in N.YY format
                format_counts['helix_format'] += 1
                helix_count[prefix] = helix_count.get(prefix, 0) + 1
            elif prefix == 'n':
                format_counts['n_terminal'] += 1
            elif prefix == 'c':
                format_counts['c_terminal'] += 1
            elif prefix == 'L' or (len(prefix) == 2 and prefix.isdigit()):
                format_counts['loop'] += 1
            elif prefix == 'X':
                format_counts['unassigned'] += 1
            else:
                format_counts['unassigned'] += 1
        else:
            format_counts['unassigned'] += 1

    print(f"[DEBUG] Column format counts: helix={format_counts['helix_format']}, n_terminal={format_counts['n_terminal']}, c_terminal={format_counts['c_terminal']}, loop={format_counts['loop']}, unassigned={format_counts['unassigned']}")
    print(f"[DEBUG] TM helix positions: {sum(helix_count.values())}")
    for helix, count in sorted(helix_count.items(), key=lambda x: int(x[0])):
        print(f"[DEBUG] Helix {helix}: {count} positions")

    # Count positions by helix for validation
    pos_counts = count_residues_by_helix(residue_table_grn)
    print(f"[INFO] TM residue positions: {pos_counts['tm_total']}")
    for helix, count in sorted(pos_counts['helices'].items(), key=lambda x: int(x[0])):
        print(f"[INFO] Helix {helix}: {count} positions")
    print(f"[INFO] Other positions: N-terminal={pos_counts['other']['n']}, C-terminal={pos_counts['other']['c']}, Loop={pos_counts['other']['L']}, Unassigned={pos_counts['other']['X']}")

    # Sort tables by GRN position
    residue_table_grn = sort_grn_columns(residue_table_grn)
    distance_table_grn = sort_grn_columns(distance_table_grn)
    ca_residue_table_grn = sort_grn_columns(ca_residue_table_grn)
    ca_distance_table_grn = sort_grn_columns(ca_distance_table_grn)

    # Check if the helices were properly identified - if not, use direct assignment
    pos_counts = count_residues_by_helix(residue_table_grn)
    if pos_counts['tm_total'] == 0:
        print("[WARNING] No TM helix positions detected in GRN mapping, using direct assignment")
        # Create a tables dictionary
        tables = {
            "residue_table": residue_table_grn,
            "distance_table": distance_table_grn,
            "ca_residue_table": ca_residue_table_grn,
            "ca_distance_table": ca_distance_table_grn
        }

        # Use direct assignment
        tables = assign_helix_numbers_to_msa_tables(tables, filtered_structures, global_ref)

        # Update our tables with the corrected ones
        residue_table_grn = tables["residue_table"]
        distance_table_grn = tables["distance_table"]
        ca_residue_table_grn = tables["ca_residue_table"]
        ca_distance_table_grn = tables["ca_distance_table"]

        # Verify the helix positions now
        pos_counts = count_residues_by_helix(residue_table_grn)
        print(f"[INFO] After direct assignment - TM residue positions: {pos_counts['tm_total']}")
        for helix, count in sorted(pos_counts['helices'].items(), key=lambda x: int(x[0])):
            print(f"[INFO] Helix {helix}: {count} positions")

    # Calculate helix distance statistics
    print("[INFO] Calculating helix distance statistics...")
    helix_stats = calculate_helix_distances(distance_table_grn)
    ca_helix_stats = calculate_helix_distances(ca_distance_table_grn)

    # Print the closest residue in each helix
    print("[INFO] Closest residues to retinal in each helix:")
    for helix in sorted(helix_stats.keys(), key=int):
        if helix_stats[helix]['closest_position']:
            print(f"[INFO] Helix {helix}: Position {helix_stats[helix]['closest_position']} (distance: {helix_stats[helix]['closest_mean']:.2f}Å)")

    return {
        "residue_table": residue_table_grn,
        "distance_table": distance_table_grn,
        "ca_residue_table": ca_residue_table_grn,
        "ca_distance_table": ca_distance_table_grn,
        "helix_stats": helix_stats,
        "ca_helix_stats": ca_helix_stats,
        "excluded_structures": excluded_structures if rmsd_df is not None else []
    }
