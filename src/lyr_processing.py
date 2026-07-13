# lyr_processing.py
"""
LYR (Lysine-Retinal Schiff base) processing utilities.

This module provides functions to:
1. Split LYR residues into LYS (lysine) + RET (retinal) components
2. Rename LIG (generic ligand) to RET for predicted structures
3. Standardize all retinal naming at the beginning of the workflow

IMPORTANT: These functions should be called EARLY in the workflow,
before any caching or downstream processing, to ensure consistent
residue naming throughout the pipeline.
"""

import pandas as pd
import numpy as np # For NaN handling if needed

# Define canonical atom names for the Lysine part when it's within an LYR residue
# These are standard Lysine atom names.
LYS_ATOM_NAMES_IN_LYR = {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'}
# Any other atom found in an LYR entry will be assumed to be part of the Retinal moiety.

def convert_single_lyr_entry(lyr_df: pd.DataFrame, retinal_res_name: str = 'RET') -> pd.DataFrame:
    """
    Converts a DataFrame containing atoms of a single LYR residue instance
    into two sets of atoms: one for LYS (ATOM group) and one for RET (HETATM group).

    Args:
        lyr_df (pd.DataFrame): DataFrame containing all atoms for one LYR residue instance.
                               It's assumed all rows belong to the same original LYR residue
                               (i.e., same pdb_id, auth_chain_id, auth_seq_id,
                               potentially different res_atom_name or alt_id).
        retinal_res_name (str): The res_name3l to assign to the retinal moiety.

    Returns:
        pd.DataFrame: A DataFrame containing the separated LYS and RET atoms.
                      Returns an empty DataFrame if input is empty or processing fails.
    """
    if lyr_df.empty:
        return pd.DataFrame()

    # Ensure 'res_atom_name' column exists
    if 'res_atom_name' not in lyr_df.columns:
        print("Warning: 'res_atom_name' column missing in LYR DataFrame. Cannot process.")
        return lyr_df # Return as is, or an empty DF

    lys_part_rows = []
    ret_part_rows = []

    # Preserve original columns, make copies for modification
    # This ensures all original data (atom_id, coords, b_factor, alt_id, etc.) is kept
    for _, atom_row_series in lyr_df.iterrows():
        atom_row = atom_row_series.copy() # Make a mutable copy of the row

        if atom_row['res_atom_name'] in LYS_ATOM_NAMES_IN_LYR:
            atom_row['res_name3l'] = 'LYS'
            atom_row['res_name1l'] = 'K'
            atom_row['group'] = 'ATOM'
            lys_part_rows.append(atom_row)
        else:
            # Assume all other atoms are part of the retinal moiety
            atom_row['res_name3l'] = retinal_res_name
            atom_row['group'] = 'HETATM'  # Retinal is a HETATM
            ret_part_rows.append(atom_row)

    # Combine the processed parts
    processed_parts = []
    if lys_part_rows:
        processed_parts.append(pd.DataFrame(lys_part_rows))
    if ret_part_rows:
        processed_parts.append(pd.DataFrame(ret_part_rows))

    if not processed_parts:
        print("Warning: LYR entry resulted in no atoms after splitting. Original LYR atoms:")
        return pd.DataFrame() # Or return original lyr_df if preferred

    return pd.concat(processed_parts, ignore_index=True)


def process_lyr_in_dataframe(structure_df: pd.DataFrame, retinal_res_name: str = 'RET') -> pd.DataFrame:
    """
    Processes all LYR residues within a single structure's DataFrame.
    Identifies LYR residues, converts them to LYS and RET components,
    and reconstructs the DataFrame.

    Args:
        structure_df (pd.DataFrame): The DataFrame for a single PDB structure.
        retinal_res_name (str): The name to assign to the retinal moiety (e.g., 'RET').

    Returns:
        pd.DataFrame: The modified DataFrame with LYR processed.
    """
    if 'res_name3l' not in structure_df.columns or not (structure_df['res_name3l'] == 'LYR').any():
        return structure_df  # No LYR residues to process or column missing

    # Separate LYR atoms from the rest of the structure
    non_lyr_df = structure_df[structure_df['res_name3l'] != 'LYR'].copy()
    all_lyr_atoms_df = structure_df[structure_df['res_name3l'] == 'LYR'].copy()

    if all_lyr_atoms_df.empty:
        return structure_df # Should be caught by .any() above, but defensive

    processed_individual_lyr_dfs = []
    residue_instance_id_cols = ['auth_chain_id', 'auth_seq_id']

    # Check for insertion code column (e.g., 'pdbx_PDB_ins_code' or 'ins_code')
    # Use the one present in the DataFrame.
    ins_code_col = None
    if 'pdbx_PDB_ins_code' in all_lyr_atoms_df.columns:
        ins_code_col = 'pdbx_PDB_ins_code'
    elif 'ins_code' in all_lyr_atoms_df.columns:
        ins_code_col = 'ins_code'

    if ins_code_col:
        if all_lyr_atoms_df[ins_code_col].notna().any() and \
           (all_lyr_atoms_df[ins_code_col].astype(str).str.strip().replace('', 'nan') != 'nan').any():
            residue_instance_id_cols.append(ins_code_col)

    # Ensure all grouping columns actually exist
    valid_group_cols = [col for col in residue_instance_id_cols if col in all_lyr_atoms_df.columns]
    if 'auth_seq_id' not in valid_group_cols: # auth_seq_id is critical
        print("Warning: 'auth_seq_id' missing from LYR data. Cannot group LYR instances. Skipping LYR processing.")
        return structure_df

    # Group by each unique LYR residue instance
    for _, lyr_instance_df in all_lyr_atoms_df.groupby(valid_group_cols, dropna=False):
        converted_df = convert_single_lyr_entry(lyr_instance_df, retinal_res_name=retinal_res_name)
        if not converted_df.empty:
            processed_individual_lyr_dfs.append(converted_df)

    if not processed_individual_lyr_dfs: # No LYR residues were successfully converted
        # This means all_lyr_atoms_df will be added back, effectively no change for LYRs
        print("Warning: No LYR residues were successfully converted. Original LYR data will be kept.")
        return pd.concat([non_lyr_df, all_lyr_atoms_df], ignore_index=True).sort_values(by='atom_id').reset_index(drop=True)

    # Combine non-LYR parts with newly processed LYS/RET parts
    final_df_parts = [non_lyr_df] + processed_individual_lyr_dfs
    final_df = pd.concat(final_df_parts, ignore_index=True)

    # It's good practice to re-sort the DataFrame, e.g., by atom_id,
    # to maintain an order similar to original PDB files if desired.
    if 'atom_id' in final_df.columns:
        # Ensure atom_id is numeric for sorting, coercing errors
        final_df['atom_id'] = pd.to_numeric(final_df['atom_id'], errors='coerce')
        final_df = final_df.sort_values(by='atom_id').reset_index(drop=True)

    return final_df


def process_lyr_in_structures_dict(structures_dict: dict, retinal_res_name: str = 'RET') -> dict:
    """
    Processes LYR residues in all structures within the 'processed_structures' dictionary.
    Modifies the 'df' (and 'df_norm' if present and different) in each structure's data.
    """
    print(f"\n--- Starting LYR processing for {len(structures_dict)} structures ---")
    processed_global_dict = {} # Create a new dict to avoid modifying input dict during iteration issues

    for pdb_id, structure_data_entry in structures_dict.items():
        # Create a shallow copy of the structure_data_entry dictionary.
        # DataFrames within it will be replaced, not modified in place initially.
        current_struc_data_copy = structure_data_entry.copy()

        if 'df' in current_struc_data_copy and isinstance(current_struc_data_copy['df'], pd.DataFrame):
            original_df = current_struc_data_copy['df']
            num_lyr_before = (original_df['res_name3l'] == 'LYR').sum()

            if num_lyr_before > 0:
                print(f"  Processing LYR in {pdb_id} (found {num_lyr_before} LYR atoms)...")
                processed_df = process_lyr_in_dataframe(original_df, retinal_res_name=retinal_res_name)
                current_struc_data_copy['df'] = processed_df

                num_lyr_after = (processed_df['res_name3l'] == 'LYR').sum()
                num_lys_added = (processed_df['res_name3l'] == 'LYS').sum() - (original_df['res_name3l'] == 'LYS').sum()
                num_ret_added = (processed_df['res_name3l'] == retinal_res_name).sum() - (original_df['res_name3l'] == retinal_res_name).sum()

                if num_lyr_after == 0:
                    print(f"    Successfully processed LYR for {pdb_id}. LYS atoms added: {num_lys_added}, {retinal_res_name} atoms added: {num_ret_added}.")
                else:
                    print(f"    Warning: LYR atoms still present in {pdb_id} after processing: {num_lyr_after} atoms.")

                if len(processed_df) != len(original_df):
                    print(f"    Warning: Atom count changed for {pdb_id} during LYR processing. "
                          f"Before: {len(original_df)}, After: {len(processed_df)}")

            # Optionally, process 'df_norm' if it exists and is a distinct DataFrame
            if 'df_norm' in current_struc_data_copy and \
               isinstance(current_struc_data_copy['df_norm'], pd.DataFrame) and \
               current_struc_data_copy['df_norm'] is not original_df: # Check it's not the same object as 'df'

                original_df_norm = current_struc_data_copy['df_norm']
                if (original_df_norm['res_name3l'] == 'LYR').any():
                    print(f"  Processing LYR in df_norm for {pdb_id}...")
                    current_struc_data_copy['df_norm'] = process_lyr_in_dataframe(original_df_norm, retinal_res_name=retinal_res_name)

        processed_global_dict[pdb_id] = current_struc_data_copy

    print("--- LYR processing finished ---")
    return processed_global_dict


def process_lyr_in_processor_data(processor, retinal_res_name: str = 'RET') -> None:
    """
    Processes LYR residues in a CifBaseProcessor's main `data` DataFrame (if it exists).
    Modifies `processor.data` in place.
    """
    if hasattr(processor, 'data') and isinstance(processor.data, pd.DataFrame) and not processor.data.empty:
        print(f"\n--- Starting LYR processing for processor: {getattr(processor, 'name', 'Unnamed Processor')} ---")
        original_data_df = processor.data
        num_lyr_before = (original_data_df['res_name3l'] == 'LYR').sum()

        if num_lyr_before > 0:
            print(f"  Found {num_lyr_before} LYR atoms in processor data. Processing...")

            pdb_ids_in_processor = original_data_df['pdb_id'].unique()
            processed_dfs_for_processor = []

            for pdb_id in pdb_ids_in_processor:
                df_single_pdb = original_data_df[original_data_df['pdb_id'] == pdb_id].copy()
                if (df_single_pdb['res_name3l'] == 'LYR').any():
                     print(f"    Processing LYR for {pdb_id} within processor...")
                     processed_dfs_for_processor.append(
                         process_lyr_in_dataframe(df_single_pdb, retinal_res_name=retinal_res_name))
                else:
                     processed_dfs_for_processor.append(df_single_pdb)

            if processed_dfs_for_processor:
                processor.data = pd.concat(processed_dfs_for_processor, ignore_index=True)
                # Re-sort the entire processor.data if needed
                if 'atom_id' in processor.data.columns and 'pdb_id' in processor.data.columns:
                     processor.data['atom_id'] = pd.to_numeric(processor.data['atom_id'], errors='coerce')
                     processor.data = processor.data.sort_values(by=['pdb_id', 'atom_id']).reset_index(drop=True)

                num_lyr_after = (processor.data['res_name3l'] == 'LYR').sum()
                if num_lyr_after == 0:
                    print(f"  Successfully processed LYR for processor data.")
                else:
                    print(f"  Warning: LYR atoms still present in processor data after processing: {num_lyr_after} atoms.")
            else: # Should not happen if there were LYRs
                print(f"  No data resulted from LYR processing in processor.")

        print(f"--- LYR processing for processor finished ---")
    else:
        print(f"Processor {getattr(processor, 'name', 'Unnamed Processor')} has no data or data is not a DataFrame. Skipping LYR processing.")


def standardize_retinal_naming(df: pd.DataFrame, retinal_res_name: str = 'RET', verbose: bool = False) -> pd.DataFrame:
    """
    Standardize all retinal-related residue naming in a structure DataFrame.

    This function should be called EARLY in the workflow (before caching)
    to ensure consistent residue naming throughout the pipeline.

    Processing order:
    1. LYR → LYS + RET (split covalently-bound retinal-lysine)
    2. LIG → RET (rename generic ligand to retinal)

    Args:
        df: Structure DataFrame with res_name3l column
        retinal_res_name: Name to use for retinal (default: 'RET')
        verbose: Print processing details

    Returns:
        DataFrame with standardized retinal naming
    """
    if df is None or df.empty:
        return df

    if 'res_name3l' not in df.columns:
        return df

    df = df.copy()

    # Step 1: Process LYR → LYS + RET
    has_lyr = (df['res_name3l'] == 'LYR').any()
    if has_lyr:
        num_lyr_before = (df['res_name3l'] == 'LYR').sum()
        if verbose:
            print(f"    Processing LYR: {num_lyr_before} atoms")
        df = process_lyr_in_dataframe(df, retinal_res_name=retinal_res_name)
        num_lyr_after = (df['res_name3l'] == 'LYR').sum()
        if verbose and num_lyr_after == 0:
            num_lys = (df['res_name3l'] == 'LYS').sum()
            num_ret = (df['res_name3l'] == retinal_res_name).sum()
            print(f"    LYR split: LYS atoms, RET atoms added")

    # Step 2: Rename LIG → RET
    has_lig = (df['res_name3l'] == 'LIG').any()
    if has_lig:
        num_lig = (df['res_name3l'] == 'LIG').sum()
        if verbose:
            print(f"    Renaming LIG → {retinal_res_name}: {num_lig} atoms")
        df.loc[df['res_name3l'] == 'LIG', 'res_name3l'] = retinal_res_name

    return df


def standardize_retinal_in_structures_dict(structures_dict: dict, retinal_res_name: str = 'RET', verbose: bool = True) -> dict:
    """
    Standardize retinal naming in all structures in a dictionary.

    This processes both 'df' and 'df_norm' if present.

    Args:
        structures_dict: Dictionary of {struct_id: {'df': DataFrame, ...}}
        retinal_res_name: Name to use for retinal
        verbose: Print processing summary

    Returns:
        Dictionary with standardized structures
    """
    if verbose:
        print(f"\n--- Standardizing retinal naming for {len(structures_dict)} structures ---")

    lyr_processed = 0
    lig_processed = 0

    for struct_id, struct_data in structures_dict.items():
        # Process main DataFrame
        if 'df' in struct_data and isinstance(struct_data['df'], pd.DataFrame):
            df = struct_data['df']
            had_lyr = (df['res_name3l'] == 'LYR').any() if 'res_name3l' in df.columns else False
            had_lig = (df['res_name3l'] == 'LIG').any() if 'res_name3l' in df.columns else False

            struct_data['df'] = standardize_retinal_naming(df, retinal_res_name=retinal_res_name, verbose=False)

            if had_lyr:
                lyr_processed += 1
            if had_lig:
                lig_processed += 1

        # Process normalized DataFrame if different from main
        if 'df_norm' in struct_data and isinstance(struct_data['df_norm'], pd.DataFrame):
            if struct_data['df_norm'] is not struct_data.get('df'):
                struct_data['df_norm'] = standardize_retinal_naming(
                    struct_data['df_norm'], retinal_res_name=retinal_res_name, verbose=False
                )

    if verbose:
        print(f"  LYR→LYS+RET processed: {lyr_processed} structures")
        print(f"  LIG→RET processed: {lig_processed} structures")
        print("--- Retinal standardization complete ---")

    return structures_dict