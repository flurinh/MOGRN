"""
Functions for multiple structure alignment and Generic Residue Numbering (GRN).
These functions create and analyze multiple sequence alignments (based on structural
alignment to a global reference) with GRN position labels.
"""

import numpy as np
import pandas as pd
import re
import os
import matplotlib.pyplot as plt # Keep for potential direct use or if visualize_msa_distances is used elsewhere
from tqdm import tqdm # Keep if used in generate_grn_msa_tables or other retained functions
import pickle # For generate_grn_msa_tables caching

# Assuming visualization_functions.py might contain visualize_msa_distances or other relevant functions
# If not used directly in this file after reverting, it can be removed.
try:
    from visualization_functions import visualize_msa_distances
except ImportError:
    print("[WARN] visualization_functions.py not found or visualize_msa_distances not in it.")
    def visualize_msa_distances(*args, **kwargs): # Dummy
        print("[WARN] visualize_msa_distances called but not available.")
        return {}


def create_msa_table(seq_alignment_dicts, processed_structures_complete, global_ref, atom_type="all"):
    """
    Creates an MSA-like table from structural alignment dictionaries.
    Columns are defined by the global_ref structure. Only residues from other
    structures that align to global_ref are included.

    Args:
        seq_alignment_dicts: Dictionary with sequence alignments (struct_id -> {ref_pos: aligned_pos})
        processed_structures_complete: Dictionary with structure data
        global_ref: ID of global reference structure
        atom_type: Atom subset to use for extracting residue info ('CA' or 'all')

    Returns:
        DataFrame: MSA table with residue information (AA + original auth_seq_id).
                   Columns are initially 1-based sequential integers.
                   An attribute 'column_to_auth_seq' maps these integers to
                   the auth_seq_id of the global_ref.
    """
    global_alignments = seq_alignment_dicts.get('global', {})
    if not global_ref in processed_structures_complete:
        raise ValueError(f"Global reference structure {global_ref} not found in processed_structures_complete.")
    global_ref_struct = processed_structures_complete[global_ref]

    if atom_type == "CA":
        if 'df_ca_norm' not in global_ref_struct:
            raise ValueError(f"Global reference {global_ref} missing 'df_ca_norm' for atom_type 'CA'.")
        global_ref_df = global_ref_struct['df_ca_norm']
        # Ensure we only consider CA atoms if that's the specific intent for 'df_ca_norm' usage
        global_ref_df = global_ref_df[global_ref_df['res_atom_name'] == 'CA']
    else: # "all"
        if 'df_norm' not in global_ref_struct:
            raise ValueError(f"Global reference {global_ref} missing 'df_norm' for atom_type 'all'.")
        global_ref_df = global_ref_struct['df_norm']
        # Get one representative row per residue (e.g., based on auth_seq_id and chain)
        global_ref_df = global_ref_df.drop_duplicates(subset=['auth_chain_id', 'auth_seq_id'], keep='first')

    if 'auth_seq_id' not in global_ref_df.columns:
        raise ValueError(f"Column 'auth_seq_id' not found in DataFrame for global_ref {global_ref}.")

    # Define columns based on unique, sorted auth_seq_ids of the global reference
    def sort_key_auth_id(auth_id_str_or_num):
        auth_id_str = str(auth_id_str_or_num)
        match = re.match(r"(-?\d+)([A-Za-z]*)", auth_id_str)
        if match:
            num_part, letter_part = match.groups()
            return (int(num_part), letter_part.upper())
        try: return (int(auth_id_str), "")
        except ValueError: return (float('inf'), auth_id_str.upper())

    global_ref_positions_orig_type = sorted(global_ref_df['auth_seq_id'].unique(), key=sort_key_auth_id)
    # Ensure all positions are strings for consistent dictionary keys later if needed, though not strictly for list.
    # For this function's direct use, original type is fine if comparison is consistent.

    sequences_data = {} # {struct_id: {global_ref_pos_orig_type: "AA[orig_auth_id_in_struct]"}}

    # Process global reference itself
    sequences_data[global_ref] = {}
    for ref_pos_orig in global_ref_positions_orig_type:
        res_rows = global_ref_df[global_ref_df['auth_seq_id'] == ref_pos_orig]
        if not res_rows.empty:
            aa = res_rows['res_name1l'].iloc[0] if 'res_name1l' in res_rows.columns else '?'
            sequences_data[global_ref][ref_pos_orig] = f"{aa}{ref_pos_orig}"
        else:
            sequences_data[global_ref][ref_pos_orig] = '-' # Should not happen if positions from df

    # Process other structures based on their alignment to global_ref
    all_aligned_struct_ids = set(global_alignments.keys())
    type_specific_struct_ids = set()
    if 'type' in seq_alignment_dicts:
        for type_ref_val, structs_val in seq_alignment_dicts['type'].items():
            type_specific_struct_ids.update(structs_val.keys())

    all_struct_ids_to_process = list(all_aligned_struct_ids.union(type_specific_struct_ids) - {global_ref})


    for struct_id in all_struct_ids_to_process:
        if struct_id not in processed_structures_complete:
            print(f"[WARN] create_msa_table: Structure {struct_id} from alignments not in processed_structures_complete. Skipping.")
            continue

        struct_data = processed_structures_complete[struct_id]
        current_struct_df_source = None
        if atom_type == "CA":
            if 'df_ca_norm' not in struct_data: continue
            current_struct_df_source = struct_data['df_ca_norm'][struct_data['df_ca_norm']['res_atom_name'] == 'CA']
        else: # "all"
            if 'df_norm' not in struct_data: continue
            current_struct_df_source = struct_data['df_norm'].drop_duplicates(subset=['auth_chain_id', 'auth_seq_id'], keep='first')

        if current_struct_df_source is None or current_struct_df_source.empty : continue
        if 'auth_seq_id' not in current_struct_df_source.columns : continue


        sequences_data[struct_id] = {ref_pos: '-' for ref_pos in global_ref_positions_orig_type} # Init with gaps

        alignment_to_global = global_alignments.get(struct_id)

        # If no direct global alignment, check for type-specific alignment path
        if not alignment_to_global and 'type' in seq_alignment_dicts:
            for type_ref_key, type_structs_dict in seq_alignment_dicts['type'].items():
                if struct_id in type_structs_dict:
                    type_specific_alignment = type_structs_dict[struct_id] # struct_id -> type_ref
                    type_ref_to_global_alignment = global_alignments.get(type_ref_key) # type_ref -> global_ref

                    if type_specific_alignment and type_ref_to_global_alignment:
                        # Reconstruct alignment to global_ref through type_ref
                        alignment_to_global = {}
                        for type_pos, s_pos in type_specific_alignment.items():
                            for g_pos, tr_pos in type_ref_to_global_alignment.items():
                                if tr_pos == type_pos:
                                    alignment_to_global[g_pos] = s_pos
                                    break
                        break # Found a path via a type reference

        if alignment_to_global:
            for ref_pos, aligned_s_pos in alignment_to_global.items():
                # ref_pos is an auth_seq_id from global_ref. It must be in global_ref_positions_orig_type
                if ref_pos in sequences_data[struct_id]: # Check if this ref_pos is a column
                    res_rows_struct = current_struct_df_source[current_struct_df_source['auth_seq_id'] == aligned_s_pos]
                    if not res_rows_struct.empty:
                        aa_struct = res_rows_struct['res_name1l'].iloc[0] if 'res_name1l' in res_rows_struct.columns else '?'
                        sequences_data[struct_id][ref_pos] = f"{aa_struct}{aligned_s_pos}"
        # else:
            # print(f"[DEBUG] create_msa_table: No alignment path to global_ref found for {struct_id}")

    # Ensure all PDB IDs present in processed_structures_complete (and thus potentially in msa_df.index later)
    # are also keys in sequences_data, initialized with gaps if they had no alignment data.
    # This is important if generate_grn_msa_tables passes a filtered_structures list that
    # might have entries not covered by seq_alignment_dicts.
    for sid_check in processed_structures_complete.keys():
        if sid_check not in sequences_data:
            sequences_data[sid_check] = {ref_pos: '-' for ref_pos in global_ref_positions_orig_type}
            # print(f"[DEBUG] create_msa_table: Initializing {sid_check} with all gaps as it was not in alignment dicts.")


    # Convert dictionary of dictionaries to DataFrame
    # Rows: struct_id, Columns: global_ref_positions_orig_type
    msa_df = pd.DataFrame.from_dict(sequences_data, orient='index')

    # Ensure columns are in the correct order of sorted global_ref_positions_orig_type
    # and handle cases where some structures might not have entries for all global_ref_positions if dict was sparse
    msa_df = msa_df.reindex(columns=global_ref_positions_orig_type, fill_value='-')

    # Rename columns to 1-based sequential integers for downstream compatibility (e.g., with older GRN mapping)
    # Store the mapping from these integers back to the original global_ref auth_seq_ids
    column_mapping_to_ref_auth_id = {i + 1: orig_pos for i, orig_pos in enumerate(global_ref_positions_orig_type)}
    msa_df.columns = range(1, len(msa_df.columns) + 1)
    msa_df.attrs['column_to_auth_seq'] = column_mapping_to_ref_auth_id

    return msa_df


def process_alignment_df(msa_df):
    """
    Converts alignment DataFrame (cells are "AAresnum") to a position frequency matrix.
    Args:
        msa_df: MSA DataFrame from create_msa_table (or already GRN labeled)
    Returns:
        DataFrame: Position frequency matrix (amino acids vs. positions)
    """
    amino_acid_data = {}
    for col_grn_label in msa_df.columns:
        aa_counts_at_pos = {}
        for cell_val in msa_df[col_grn_label]:
            if pd.isna(cell_val) or cell_val == '-':
                continue
            aa = str(cell_val)[0] # First char is amino acid
            aa_counts_at_pos[aa] = aa_counts_at_pos.get(aa, 0) + 1
        amino_acid_data[col_grn_label] = aa_counts_at_pos

    # Convert to DataFrame (pos as columns, aa as index) and fill NaNs with 0
    aa_df = pd.DataFrame(amino_acid_data).fillna(0).T # Transpose to get AA as columns

    # Calculate frequencies (normalize each row)
    # Ensure that we only divide by sum if sum > 0
    row_sums = aa_df.sum(axis=1)
    aa_freq_df = aa_df.apply(lambda r: r / row_sums[r.name] if row_sums[r.name] > 0 else r, axis=1)

    return aa_freq_df.T # Transpose back: AA as index, GRN pos as columns


def analyze_residue_composition(msa_df, positions_to_analyze):
    """
    Analyzes residue composition at specific GRN-labeled positions.
    Args:
        msa_df: MSA DataFrame with GRN column names.
        positions_to_analyze: List of GRN positions to analyze.
    Returns:
        dict: {position: {'counts': {}, 'frequencies': {}, 'sorted': [], 'total': int, 'conservation': float}}
    """
    results = {}
    for grn_pos in positions_to_analyze:
        if grn_pos not in msa_df.columns:
            results[grn_pos] = {'error': f'Position {grn_pos} not found in MSA table.'}
            continue

        aa_counts = {}
        total_residues_at_pos = 0
        for cell_value in msa_df[grn_pos]:
            if pd.isna(cell_value) or cell_value == '-':
                aa_counts['-'] = aa_counts.get('-', 0) + 1
                continue

            aa = str(cell_value)[0] # First character is amino acid
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
            if aa != '-': # Don't count gaps in total for conservation calculation
                total_residues_at_pos += 1

        if total_residues_at_pos > 0:
            # Exclude gap counts for frequency and conservation calculation
            aa_freq_no_gaps = {
                aa_key: count / total_residues_at_pos
                for aa_key, count in aa_counts.items() if aa_key != '-'
            }
            sorted_aa_freq = sorted(aa_freq_no_gaps.items(), key=lambda item: item[1], reverse=True)

            results[grn_pos] = {
                'counts': aa_counts, # Includes gaps
                'frequencies': aa_freq_no_gaps, # Excludes gaps from normalization
                'sorted': sorted_aa_freq,
                'total_residues_no_gaps': total_residues_at_pos,
                'conservation': max(aa_freq_no_gaps.values()) if aa_freq_no_gaps else 0.0
            }
        else: # All gaps at this position
            results[grn_pos] = {
                'counts': aa_counts,
                'frequencies': {}, 'sorted': [], 'total_residues_no_gaps': 0, 'conservation': 0.0,
                'note': 'All gaps at this position.'
            }
    return results

def sort_grn_columns(df):
    """
    Sorts columns of a GRN-labeled DataFrame (e.g., "1.50", "n.10", "1.501").
    """
    def get_sort_key(col_name_str_original):
        col_name_str = str(col_name_str_original)
        type_prefix = 5
        primary_val = float('inf')
        secondary_val = float('inf') # For TM position or insertion part

        if col_name_str.startswith('n.'):
            type_prefix = 0
            try: primary_val = float(col_name_str.split('.')[1])
            except: primary_val = str(col_name_str_original)
        elif col_name_str.startswith('c.'):
            type_prefix = 10
            try: primary_val = float(col_name_str.split('.')[1])
            except: primary_val = str(col_name_str_original)
        elif re.match(r"(\d{1,2})(\d{1,2})\.(\d+)", col_name_str): # Loop AB.CCC
            match = re.match(r"(\d{1,2})(\d{1,2})\.(\d+)", col_name_str)
            type_prefix = 8
            try:
                primary_val = int(match.group(1))
                secondary_val = float(f"{int(match.group(2))}.{int(match.group(3))}")
            except: primary_val = str(col_name_str_original)
        elif col_name_str.startswith('L.'):
            type_prefix = 9
            try: primary_val = float(col_name_str.split('.')[1])
            except: primary_val = str(col_name_str_original)
        elif re.match(r"(\d+)\.([\d\.]+)", col_name_str): # TM helix (1.50) or insertion (1.501)
            match = re.match(r"(\d+)\.(.+)", col_name_str)
            try:
                helix_num_str = match.group(1)
                pos_part_str = match.group(2)
                primary_val = int(helix_num_str)
                secondary_val = float(pos_part_str) # Handles "50" and "50.001"
                if 1 <= primary_val <= 7: type_prefix = primary_val
                else: type_prefix = 7.5 # Other X.Y sort after TM7
            except (ValueError, AttributeError):
                primary_val = str(col_name_str_original); type_prefix = 99
        elif col_name_str.startswith(('I.', 'ERR_I.', 'FMT_ERR_I.')):
            type_prefix = 98
            try: primary_val = float(re.search(r'(\d+(\.\d+)?)', col_name_str).group(1)) # Try to get number
            except: primary_val = str(col_name_str_original)
        elif col_name_str.startswith('X.'): # Fallback from assign_helix_numbers_to_msa_tables
             type_prefix = 97
             try: primary_val = int(col_name_str.split('.')[1])
             except: primary_val = str(col_name_str_original)
        else:
            type_prefix = 99
            primary_val = str(col_name_str_original)
        return (type_prefix, primary_val, secondary_val)
    try:
        string_cols = [col for col in df.columns if isinstance(col, (str, int, float))] # Allow numbers too
        # Non-string/numeric columns are problematic for this sort key, filter or handle them
        sorted_string_cols = sorted(string_cols, key=get_sort_key)
        return df[sorted_string_cols]
    except Exception as e:
        print(f"[ERROR] Sorting GRN columns failed: {e}. Returning unsorted.")
        import traceback
        traceback.print_exc()
        return df


def create_msa_distance_table(seq_alignment_dicts, processed_structures_complete, global_ref_id,
                              residue_table_with_grn_cols, # Pass the GRN labeled table
                              distance_type="sidechain"):
    """
    Creates a distance table (retinal to residue) based on an existing GRN-labeled MSA table.
    The columns of the output distance table will match the GRN columns of the input residue_table.
    """
    from scipy.spatial.distance import cdist # Import locally

    if residue_table_with_grn_cols.empty:
        print("[WARN] create_msa_distance_table: Input residue_table is empty. Returning empty DataFrame.")
        return pd.DataFrame(index=residue_table_with_grn_cols.index, columns=residue_table_with_grn_cols.columns)

    # Initialize a DataFrame for distances with NaNs, matching structure of residue_table_with_grn_cols
    distance_df_data = {
        grn_col: pd.Series(np.nan, index=residue_table_with_grn_cols.index)
        for grn_col in residue_table_with_grn_cols.columns
    }
    distance_df = pd.DataFrame(distance_df_data)

    # Iterate through each structure (row) in the residue_table
    for struct_id in tqdm(residue_table_with_grn_cols.index, desc=f"Calculating distances ({distance_type})"):
        if struct_id not in processed_structures_complete:
            # print(f"[DEBUG] {struct_id} not in processed_structures. Skipping for distance calc.")
            continue
        struct_data = processed_structures_complete[struct_id]

        # Determine source DataFrame for coordinates based on distance_type
        struct_coord_df_source = None
        if distance_type == "backbone": # CA-CA distances typically
            if 'df_ca_norm' in struct_data:
                struct_coord_df_source = struct_data['df_ca_norm'][struct_data['df_ca_norm']['res_atom_name'] == 'CA'].copy()
        else: # "sidechain" - all atoms
            if 'df_norm' in struct_data:
                struct_coord_df_source = struct_data['df_norm'].copy()

        if struct_coord_df_source is None or struct_coord_df_source.empty:
            # print(f"[DEBUG] No suitable coord df for {struct_id}, type {distance_type}. Skipping.")
            continue
        if 'res_name3l' not in struct_coord_df_source.columns or 'auth_seq_id' not in struct_coord_df_source.columns:
            # print(f"[DEBUG] Essential columns missing in coord_df for {struct_id}. Skipping.")
            continue

        # Extract retinal coordinates for this structure
        ret_df = struct_coord_df_source[struct_coord_df_source['res_name3l'] == 'RET']
        if ret_df.empty: ret_df = struct_coord_df_source[struct_coord_df_source['res_name3l'] == 'LIG'] # Fallback
        if ret_df.empty:
            # print(f"[DEBUG] No retinal found for {struct_id}. Distances will be NaN.")
            continue
        ret_coords_xyz = ret_df[['x', 'y', 'z']].values.astype(float)
        if ret_coords_xyz.size == 0: continue


        # Iterate through each GRN column of the input residue_table
        for grn_col_label in residue_table_with_grn_cols.columns:
            cell_value = residue_table_with_grn_cols.at[struct_id, grn_col_label]
            if pd.isna(cell_value) or cell_value == '-':
                continue # Skip gaps

            # Extract original auth_seq_id from the cell value (e.g., "A123" -> "123")
            # This auth_seq_id is from the *current struct_id*, not the global_ref
            match = re.match(r"[A-Z Оюн\?\*](-?[\d\w]+)", str(cell_value))
            if not match:
                # print(f"[DEBUG] Could not parse auth_seq_id from cell '{cell_value}' for {struct_id}, col {grn_col_label}")
                continue

            original_auth_seq_id_in_struct = match.group(1)
            # Attempt to convert to original type for matching if needed (int, or str if it has letters)
            try:
                # Check if it can be an int, otherwise it's a string (e.g. 10A)
                # The comparison with df['auth_seq_id'] must handle mixed types or convert consistently
                if re.fullmatch(r"-?\d+", original_auth_seq_id_in_struct):
                    key_for_lookup = int(original_auth_seq_id_in_struct)
                else:
                    key_for_lookup = original_auth_seq_id_in_struct
            except ValueError:
                 key_for_lookup = original_auth_seq_id_in_struct


            # Find atoms of this residue in the current structure's coordinate DataFrame
            # Ensure we are not selecting retinal itself if using 'all' atoms and retinal is not filtered out above
            if distance_type == "sidechain": # use 'all' atoms but exclude retinal
                residue_atoms_df = struct_coord_df_source[
                    (struct_coord_df_source['auth_seq_id'] == key_for_lookup) &
                    (struct_coord_df_source['res_name3l'] != 'RET') & # Ensure not retinal
                    (struct_coord_df_source['res_name3l'] != 'LIG')
                ]
            else: # backbone (CA only)
                 residue_atoms_df = struct_coord_df_source[struct_coord_df_source['auth_seq_id'] == key_for_lookup]


            if residue_atoms_df.empty:
                # print(f"[DEBUG] Residue {original_auth_seq_id_in_struct} not found in {struct_id}'s coord df for col {grn_col_label}")
                continue

            res_coords_xyz = residue_atoms_df[['x', 'y', 'z']].values.astype(float)
            if res_coords_xyz.size == 0: continue

            try:
                distances_matrix = cdist(res_coords_xyz, ret_coords_xyz)
                min_dist = float(np.min(distances_matrix))
                distance_df.at[struct_id, grn_col_label] = min_dist
            except Exception as e_dist:
                print(f"[WARN] cdist failed for {struct_id}, col {grn_col_label}, res {original_auth_seq_id_in_struct}: {e_dist}")

    return distance_df


def find_closest_ret_residue(struct_data, helix_num):
    """
    Finds the auth_seq_id of the residue in a specific helix that is closest to retinal.
    Uses 'df_norm' for coordinates.
    Args:
        struct_data: Dictionary for a single structure from processed_structures_complete.
        helix_num: Integer, the helix number (1-7).
    Returns:
        The auth_seq_id (original type) of the closest residue, or None.
    """
    from scipy.spatial.distance import cdist # Import locally

    if 'df_norm' not in struct_data or struct_data['df_norm'].empty: return None
    df = struct_data['df_norm']
    if not all(col in df.columns for col in ['res_name3l', 'auth_seq_id', 'helix_num', 'tm_helix', 'x', 'y', 'z']):
        return None # Missing essential columns

    ret_df = df[df['res_name3l'] == 'RET']
    if ret_df.empty: ret_df = df[df['res_name3l'] == 'LIG']
    if ret_df.empty: return None
    ret_coords = ret_df[['x', 'y', 'z']].values.astype(float)
    if ret_coords.size == 0: return None

    helix_residues_df = df[
        (df['helix_num'] == helix_num) & \
        (df['tm_helix'] == True) & \
        (~df['res_name3l'].isin(['RET', 'LIG'])) # Exclude retinal itself
    ]
    if helix_residues_df.empty: return None

    min_overall_distance = float('inf')
    closest_auth_seq_id = None

    for auth_id, res_group_df in helix_residues_df.groupby('auth_seq_id'):
        res_atoms_coords = res_group_df[['x', 'y', 'z']].values.astype(float)
        if res_atoms_coords.size == 0: continue

        dist_matrix = cdist(res_atoms_coords, ret_coords)
        min_dist_for_this_res = np.min(dist_matrix)

        if min_dist_for_this_res < min_overall_distance:
            min_overall_distance = min_dist_for_this_res
            closest_auth_seq_id = auth_id # auth_id is the original type from df

    return closest_auth_seq_id


# --- Main Orchestrating Function for GRN Table Generation ---
def generate_grn_msa_tables(seq_alignment_dicts,
                            processed_structures_complete,
                            global_ref_id,
                            rmsd_df=None,
                            max_rmsd_threshold=3.0,
                            structure_mapping=None,
                            output_dir="opsin_grn_tables"): # Added output_dir
    """
    Original functionality: Creates MSA tables (residue & distance) based on structural
    alignments to a global reference, then assigns GRN labels to columns.
    Does NOT attempt to include all residues, only structurally aligned ones.
    """
    print(f"[INFO] generate_grn_msa_tables called for global_ref: {global_ref_id}")

    # 1. Filter structures (optional, based on RMSD and experimental priority)
    filtered_structures = processed_structures_complete.copy()
    # Make a deep copy of seq_alignment_dicts to avoid modifying the original
    filtered_seq_alignments = {
        'global': {k: v.copy() for k, v in seq_alignment_dicts.get('global', {}).items()},
        'type': {
            type_ref: {struct_k: struct_v.copy() for struct_k, struct_v in type_aligns.items()}
            for type_ref, type_aligns in seq_alignment_dicts.get('type', {}).items()
        }
    }
    excluded_structures_log = []

    # Prioritize experimental if mapping provided
    if structure_mapping:
        # ... (your existing prioritization logic - ensure it modifies filtered_structures and filtered_seq_alignments) ...
        pass # Placeholder for your existing logic

    # Filter by RMSD to global_ref_id
    if rmsd_df is not None and global_ref_id in rmsd_df.index:
        # ... (your existing RMSD filtering logic - ensure it modifies filtered_structures and filtered_seq_alignments) ...
        pass # Placeholder for your existing logic

    if not filtered_structures or global_ref_id not in filtered_structures:
        print(f"[ERROR] Global reference {global_ref_id} missing or no structures left after filtering.")
        return {"residue_table": pd.DataFrame(), "distance_table": pd.DataFrame(),
                "ca_residue_table": pd.DataFrame(), "ca_distance_table": pd.DataFrame(),
                "excluded_structures": excluded_structures_log, "error": "Filtering issue"}

    print(f"[INFO] Using {len(filtered_structures)} structures after filtering for GRN table generation.")

    # 2. Create initial MSA-like tables (columns are sequential integers)
    #    These tables only contain structurally aligned residues.
    print("[INFO] Creating base MSA tables (structurally aligned residues only)...")
    base_residue_table_all = create_msa_table(filtered_seq_alignments, filtered_structures, global_ref_id, atom_type="all")
    base_residue_table_ca = create_msa_table(filtered_seq_alignments, filtered_structures, global_ref_id, atom_type="CA")

    if base_residue_table_all.empty:
        print("[ERROR] Base residue table (all atoms) is empty.")
        # Return empty structure matching expected output
        return {"residue_table": pd.DataFrame(), "distance_table": pd.DataFrame(),
                "ca_residue_table": pd.DataFrame(), "ca_distance_table": pd.DataFrame(),
                "excluded_structures": excluded_structures_log, "error": "Empty base_residue_table_all"}


    # 3. GRN Assignment Logic (copied & adapted from your `generate_grn_msa_tables` snippet)
    #    This renames the sequential integer columns of the base tables to GRN labels.
    print(f"[INFO] Assigning GRN labels to columns based on reference: {global_ref_id}")
    ref_struct_data_for_grn = filtered_structures[global_ref_id]

    # Ensure reference structure has necessary data for GRN assignment
    ref_df_for_helices = ref_struct_data_for_grn.get('df_ca_norm', ref_struct_data_for_grn.get('df_norm'))
    if ref_df_for_helices is None or not all(c in ref_df_for_helices.columns for c in ['helix_num', 'tm_helix', 'auth_seq_id']):
        print(f"[ERROR] Reference {global_ref_id} lacks helix/auth_seq_id annotations in its DataFrame. Cannot assign GRNs.")
        # Return tables with sequential columns if GRN assignment fails
        return {"residue_table": base_residue_table_all, "distance_table": pd.DataFrame(), # Distance table needs GRN
                "ca_residue_table": base_residue_table_ca, "ca_distance_table": pd.DataFrame(),
                "excluded_structures": excluded_structures_log, "grn_error": "Ref annotation missing"}

    # Get mapping from sequential column number to global_ref auth_seq_id
    # This comes from create_msa_table's attrs
    column_to_ref_auth_seq_map = base_residue_table_all.attrs.get('column_to_auth_seq', {})
    if not column_to_ref_auth_seq_map:
         print("[ERROR] 'column_to_auth_seq' mapping missing from base_residue_table_all. Cannot assign GRNs.")
         return {"residue_table": base_residue_table_all, "distance_table": pd.DataFrame(),
                 "ca_residue_table": base_residue_table_ca, "ca_distance_table": pd.DataFrame(),
                 "excluded_structures": excluded_structures_log, "grn_error": "column_to_auth_seq missing"}


    # --- GRN Naming Logic (Simplified from your generate_grn_msa_tables snippet) ---
    # This part needs to be robust and match the GRN system you expect.
    # It maps original sequential column numbers to new GRN string labels.

    # A. Determine Helix Properties from Reference Structure
    ref_helix_auth_seq_map = {} # {helix_num: [sorted auth_ids in that helix]}
    ref_helix_boundaries = {} # {helix_num: (min_auth_id, max_auth_id)}

    def sort_key_auth_id(auth_id_str_or_num): # Re-define or import if used elsewhere
        auth_id_str = str(auth_id_str_or_num)
        match = re.match(r"(-?\d+)([A-Za-z]*)", auth_id_str)
        if match: num_part, letter_part = match.groups(); return (int(num_part), letter_part.upper())
        try: return (int(auth_id_str), "")
        except ValueError: return (float('inf'), auth_id_str.upper())

    for h_num in range(1, 8):
        helix_df = ref_df_for_helices[
            (ref_df_for_helices['helix_num'] == h_num) & (ref_df_for_helices['tm_helix'] == True)
        ]
        if not helix_df.empty:
            unique_auth_ids_in_helix = sorted(helix_df['auth_seq_id'].unique(), key=sort_key_auth_id)
            if unique_auth_ids_in_helix:
                ref_helix_auth_seq_map[h_num] = unique_auth_ids_in_helix
                ref_helix_boundaries[h_num] = (unique_auth_ids_in_helix[0], unique_auth_ids_in_helix[-1])

    sorted_ref_helices_by_start = sorted(ref_helix_boundaries.items(), key=lambda item: sort_key_auth_id(item[1][0]))


    # B. Determine Pivot (X.50) for each helix in the reference
    #    For simplicity here, using find_closest_ret_residue on the ref_struct_data_for_grn
    #    Your more complex pivot logic from generate_grn_msa_tables (using average distances
    #    from a preliminary distance table) could be adapted here if needed.
    ref_helix_pivots_auth_id = {} # {h_num: pivot_auth_id}
    for h_num in ref_helix_auth_seq_map.keys():
        pivot_auth_id = find_closest_ret_residue(ref_struct_data_for_grn, h_num)
        if pivot_auth_id and pivot_auth_id in ref_helix_auth_seq_map[h_num]:
            ref_helix_pivots_auth_id[h_num] = pivot_auth_id
            # print(f"[DEBUG] GRN: Ref Helix {h_num} pivot (closest to RET): {pivot_auth_id}")
        else: # Fallback to geometric middle
            middle_idx = len(ref_helix_auth_seq_map[h_num]) // 2
            ref_helix_pivots_auth_id[h_num] = ref_helix_auth_seq_map[h_num][middle_idx]
            # print(f"[DEBUG] GRN: Ref Helix {h_num} pivot (geometric middle): {ref_helix_pivots_auth_id[h_num]}")

    # C. Generate GRN Label for each original column number
    sequential_col_to_grn_label_map = {}
    for seq_col_num, ref_auth_id_at_col in column_to_ref_auth_seq_map.items():
        grn_label = f"X.{seq_col_num}" # Default if no specific GRN rule applies

        # Is it in a TM helix?
        assigned_helix_num_for_grn = None
        for h_num_grn, auth_ids_in_h_grn in ref_helix_auth_seq_map.items():
            if ref_auth_id_at_col in auth_ids_in_h_grn:
                assigned_helix_num_for_grn = h_num_grn
                break

        if assigned_helix_num_for_grn and assigned_helix_num_for_grn in ref_helix_pivots_auth_id:
            pivot_auth_id_for_helix = ref_helix_pivots_auth_id[assigned_helix_num_for_grn]
            auth_ids_this_helix = ref_helix_auth_seq_map[assigned_helix_num_for_grn]
            try:
                current_idx_in_helix = auth_ids_this_helix.index(ref_auth_id_at_col)
                pivot_idx_in_helix = auth_ids_this_helix.index(pivot_auth_id_for_helix)
                offset = current_idx_in_helix - pivot_idx_in_helix
                grn_label = f"{assigned_helix_num_for_grn}.{50 + offset}"
            except ValueError: # Should not happen if data is consistent
                print(f"[WARN] GRN: Value error for {ref_auth_id_at_col} in helix {assigned_helix_num_for_grn}")

        # Is it N-term, C-term, or Loop? (Only if not already assigned as TM helix)
        elif sorted_ref_helices_by_start: # Need at least one defined helix
            first_helix_num, (first_helix_start_auth_id, _) = sorted_ref_helices_by_start[0]
            last_helix_num, (_, last_helix_end_auth_id) = sorted_ref_helices_by_start[-1]

            key_ref_auth = sort_key_auth_id(ref_auth_id_at_col)
            key_first_start = sort_key_auth_id(first_helix_start_auth_id)
            key_last_end = sort_key_auth_id(last_helix_end_auth_id)

            if key_ref_auth < key_first_start: # N-terminal
                # Distance needs to be calculated based on numeric part of auth_seq_id
                dist_n = key_first_start[0] - key_ref_auth[0] if isinstance(key_first_start[0], int) and isinstance(key_ref_auth[0], int) else seq_col_num
                grn_label = f"n.{dist_n}"
            elif key_ref_auth > key_last_end: # C-terminal
                dist_c = key_ref_auth[0] - key_last_end[0] if isinstance(key_last_end[0], int) and isinstance(key_ref_auth[0], int) else seq_col_num
                grn_label = f"c.{dist_c}"
            else: # Loop
                is_loop = False
                for i in range(len(sorted_ref_helices_by_start) - 1):
                    prev_h_num_loop, (_, prev_h_end_auth) = sorted_ref_helices_by_start[i]
                    next_h_num_loop, (next_h_start_auth, _) = sorted_ref_helices_by_start[i+1]

                    key_prev_end = sort_key_auth_id(prev_h_end_auth)
                    key_next_start = sort_key_auth_id(next_h_start_auth)

                    if key_prev_end < key_ref_auth < key_next_start:
                        # Distance from end of previous helix (numeric part)
                        dist_loop = key_ref_auth[0] - key_prev_end[0] if isinstance(key_ref_auth[0], int) and isinstance(key_prev_end[0], int) else seq_col_num
                        grn_label = f"{prev_h_num_loop}{next_h_num_loop}.{dist_loop:03d}" # AB.CCC format
                        is_loop = True
                        break
                if not is_loop: grn_label = f"L.{seq_col_num}" # Fallback general loop

        sequential_col_to_grn_label_map[seq_col_num] = grn_label

    # D. Apply GRN labels to table columns
    final_grn_column_labels = [sequential_col_to_grn_label_map.get(col_num, f"X.{col_num}") for col_num in base_residue_table_all.columns]

    residue_table_grn_final = base_residue_table_all.copy()
    residue_table_grn_final.columns = final_grn_column_labels

    ca_residue_table_grn_final = base_residue_table_ca.copy()
    if len(ca_residue_table_grn_final.columns) == len(final_grn_column_labels):
        ca_residue_table_grn_final.columns = final_grn_column_labels
    else:
        print(f"[WARN] CA residue table column count mismatch. Retaining sequential numbers for CA table.")
        # ca_residue_table_grn_final remains with sequential columns

    # Sort columns by GRN
    residue_table_grn_final = sort_grn_columns(residue_table_grn_final)
    if list(ca_residue_table_grn_final.columns) == list(residue_table_grn_final.columns): # Check if cols were successfully renamed
        ca_residue_table_grn_final = sort_grn_columns(ca_residue_table_grn_final)


    # 4. Create Distance Tables using the NEW GRN-labeled residue tables as template
    print("[INFO] Creating GRN-labeled distance tables...")
    distance_table_grn_final = create_msa_distance_table(
        filtered_seq_alignments, # These are used by create_msa_table *inside* create_msa_distance_table
        filtered_structures,     # if it reconstructs the base alignment.
        global_ref_id,           # Better: pass the GRN labeled table directly.
        residue_table_grn_final, # Pass the GRN labeled residue table
        distance_type="sidechain"
    )

    ca_distance_table_grn_final = create_msa_distance_table(
        filtered_seq_alignments,
        filtered_structures,
        global_ref_id,
        ca_residue_table_grn_final, # Pass GRN labeled CA residue table
        distance_type="backbone"
    )

    # Ensure distance tables are sorted consistently with residue tables
    distance_table_grn_final = sort_grn_columns(distance_table_grn_final.reindex(columns=residue_table_grn_final.columns))
    if not ca_residue_table_grn_final.empty:
         ca_distance_table_grn_final = sort_grn_columns(ca_distance_table_grn_final.reindex(columns=ca_residue_table_grn_final.columns))


    # 5. Save outputs
    os.makedirs(output_dir, exist_ok=True)
    try:
        residue_table_grn_final.to_csv(os.path.join(output_dir, "residue_table_grn.csv"))
        distance_table_grn_final.to_csv(os.path.join(output_dir, "distance_table_grn.csv"))
        if not ca_residue_table_grn_final.empty:
            ca_residue_table_grn_final.to_csv(os.path.join(output_dir, "ca_residue_table_grn.csv"))
        if not ca_distance_table_grn_final.empty:
            ca_distance_table_grn_final.to_csv(os.path.join(output_dir, "ca_distance_table_grn.csv"))

        # Save other metadata if needed (e.g., the sequential_col_to_grn_label_map for debugging)
        # with open(os.path.join(output_dir, "grn_mapping_details.pkl"), "wb") as f:
        #     pickle.dump({
        #         "sequential_col_to_ref_auth_seq": column_to_ref_auth_seq_map,
        #         "sequential_col_to_grn_label": sequential_col_to_grn_label_map
        #     }, f)
        print(f"[INFO] GRN tables saved to {output_dir}")
    except Exception as e_save:
        print(f"[WARN] Error saving GRN tables: {e_save}")

    return {
        "residue_table": residue_table_grn_final,
        "distance_table": distance_table_grn_final,
        "ca_residue_table": ca_residue_table_grn_final,
        "ca_distance_table": ca_distance_table_grn_final,
        "excluded_structures": excluded_structures_log
        # Add other relevant outputs like helix_stats if calculated here
    }


# --- calculate_helix_distances, count_residues_by_helix ---
# These helper functions can remain as they are, as they operate on GRN-labeled DataFrames.

def calculate_helix_distances(distance_table_grn): # Renamed arg for clarity
    if distance_table_grn.empty: return {}
    mean_distances = distance_table_grn.mean(skipna=True)
    std_distances = distance_table_grn.std(skipna=True)
    helix_stats = {}
    for grn_pos_label in distance_table_grn.columns:
        pos_str = str(grn_pos_label)
        match_tm = re.match(r"(\d+)\.([\d\.]+)", pos_str) # Matches "1.50" and "1.501"
        if match_tm:
            helix, pos_val_str = match_tm.groups()
            if 1 <= int(helix) <= 7:
                if helix not in helix_stats:
                    helix_stats[helix] = {'positions': [], 'means': [], 'stds': []}
                helix_stats[helix]['positions'].append(grn_pos_label) # Store original GRN label
                helix_stats[helix]['means'].append(mean_distances.get(grn_pos_label, np.nan))
                helix_stats[helix]['stds'].append(std_distances.get(grn_pos_label, np.nan))

    for helix_key, data_val in helix_stats.items():
        # Sort by the numeric part of the position
        def get_pos_float(p_str):
            try: return float(str(p_str).split('.',1)[1])
            except: return float('inf') # Error last

        sorted_indices = sorted(range(len(data_val['positions'])), key=lambda i: get_pos_float(data_val['positions'][i]))
        data_val['positions'] = [data_val['positions'][i] for i in sorted_indices]
        data_val['means'] = [data_val['means'][i] for i in sorted_indices]
        data_val['stds'] = [data_val['stds'][i] for i in sorted_indices]

        # Find closest position based on mean distance
        if data_val['means'] and not all(np.isnan(m) for m in data_val['means']):
            valid_means_with_indices = [(mean, idx) for idx, mean in enumerate(data_val['means']) if pd.notna(mean)]
            if valid_means_with_indices:
                min_mean, min_idx_in_valid = min(valid_means_with_indices, key=lambda x: x[0])
                original_min_idx = valid_means_with_indices[min_idx_in_valid][1] # map back to original list index
                data_val['closest_position'] = data_val['positions'][original_min_idx]
                data_val['closest_mean'] = min_mean
            else:
                data_val['closest_position'] = None; data_val['closest_mean'] = None
        else:
            data_val['closest_position'] = None; data_val['closest_mean'] = None
    return helix_stats

def count_residues_by_helix(df_grn_labeled): # Renamed arg
    if df_grn_labeled.empty: return {'tm_total': 0, 'helices': {}, 'other': {'n': 0, 'c': 0, 'L': 0, 'X': 0}}
    tm_counts = {str(h):0 for h in range(1,8)}
    other_counts = {'n': 0, 'c': 0, 'L': 0, 'X': 0} # L for general loops, X for unknown

    for grn_col in df_grn_labeled.columns:
        col_str = str(grn_col)
        if col_str.startswith('n.'): other_counts['n'] += 1
        elif col_str.startswith('c.'): other_counts['c'] += 1
        elif re.match(r"(\d{1,2})(\d{1,2})\.(\d+)", col_str): other_counts['L'] +=1 # Loop AB.CCC
        elif col_str.startswith('L.'): other_counts['L'] +=1 # General Loop L.X
        elif re.match(r"(\d+)\.([\d\.]+)", col_str): # TM Helix like 1.50 or 1.501
            helix_part = col_str.split('.')[0]
            if helix_part.isdigit() and 1 <= int(helix_part) <= 7:
                tm_counts[helix_part] = tm_counts.get(helix_part, 0) + 1
            else:
                other_counts['X'] += 1 # Non-standard helix number
        else: # Fallback for other formats like I.X, ERR_I.X, X.sequential
            other_counts['X'] += 1

    return {
        'tm_total': sum(tm_counts.values()),
        'helices': tm_counts,
        'other': other_counts
    }

