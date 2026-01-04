"""
Functions for loading, filtering, and preprocessing opsin structure data.

NOTE: Most dataset loading functions in this module are DEPRECATED.
Use the new protos StructureProcessor API instead. See opsin_analysis_workflow.py.

This module retains utility functions like:
- format_cif_columns
- ensure_structure_dtypes
"""

import os
import pandas as pd
import numpy as np
from typing import Dict, List, Any, Optional, Union
from pathlib import Path
import pickle
import json
from tqdm import tqdm

# New protos API
from protos.processing.structure import StructureProcessor
from protos.io.paths.path_config import ProtosPaths


try:
    from src.lyr_processing import process_lyr_in_processor_data
    print("[INFO] Successfully imported LYR processing utilities")
except ImportError:
    print("[WARNING] Could not import LYR processing utilities. LYR residues will not be processed correctly.")


def load_experimental_dataset(dataset_name='mo_exp_A'):
    """
    Load a dataset using the new StructureProcessor API.

    Available datasets (created by prepare_data.py):
    - mo_exp_A: Experimental structures Set A (pre-Sept 2021, within Boltz training)
    - mo_exp_B: Experimental structures Set B (post-Sept 2021) + Hideaki exp
    - mo_pred_exp: Boltz predictions of experimental structures (63 + 8 Hideaki = 71)
    - mo_pred_novel: Boltz predictions of novel opsins (~58)

    Args:
        dataset_name: Name of the dataset to load

    Returns:
        StructureProcessor: Loaded processor with dataset
    """
    processor = StructureProcessor("dataset_loader")
    processor.load_dataset(dataset_name)
    print(f"[INFO] {dataset_name} loaded with {len(processor.structure_ids)} structures.")
    return processor

def filter_structures_by_chain(cp, chain='A'):
    """
    Extract structures for a specific chain from a CifProcessor dataset
    
    Args:
        cp: CifBaseProcessor with loaded dataset
        chain: Chain identifier to extract
        
    Returns:
        dict: Dictionary of structures with the specified chain
    """
    structures = {}
    for pdb_id in cp.pdb_ids:
        df = cp.data[cp.data['pdb_id'] == pdb_id].copy()
        df_chain = df[df['auth_chain_id'] == chain].copy()
        if df_chain.empty:
            print(f"[WARNING] {pdb_id}: Chain {chain} not found; skipping.")
            continue
        structures[pdb_id] = {'df': df_chain}
    print(f"[INFO] {len(structures)} structures retained for chain {chain}.")
    return structures

def filter_structures_by_chain_and_retinal(processor, chain='A', retinal_name='RET', cutoff=6.0):
    """
    Filter structures by chain and only include the retinal closest to that chain.
    Also handles the case where retinal is named 'LIG' instead of 'RET'.
    
    Args:
        processor: CifBaseProcessor instance
        chain: Chain ID to filter by (default: 'A')
        retinal_name: Name of retinal residue (default: 'RET')
        cutoff: Distance cutoff in Angstroms for retinal selection (default: 6.0)
    
    Returns:
        Dictionary of filtered dataframes by PDB ID, with only chain A and one retinal
        Also updates the processor's data with the filtered structures.
    """
    from src.error_analysis import find_retinal_within_cutoff
    
    filtered_structures = {}
    filtered_dfs = []
    
    for pdb_id in processor.pdb_ids:
        # Get all atoms for this PDB ID
        df_pdb = processor.data[processor.data['pdb_id'] == pdb_id]
        
        # Check for and remove duplicate atoms
        if 'atom_id' in df_pdb.columns:
            initial_count = len(df_pdb)
            df_pdb = df_pdb.drop_duplicates(subset=['atom_id'], keep='first')
            if len(df_pdb) < initial_count:
                print(f"{pdb_id}: Removed {initial_count - len(df_pdb)} duplicate atoms")
        
        # Extract only the specified chain
        df_chain = df_pdb[df_pdb['auth_chain_id'] == chain].copy()
        
        if df_chain.empty:
            continue
        
        # Find retinal within cutoff distance of the chain
        retinal_df = find_retinal_within_cutoff(df_pdb, df_chain, cutoff=cutoff, retinal_name=retinal_name)
        
        # If we found retinal, combine it with the chain data
        if not retinal_df.empty:
            # Check for duplicates when combining
            combined_df = pd.concat([df_chain, retinal_df], ignore_index=True)
            if 'atom_id' in combined_df.columns:
                combined_df = combined_df.drop_duplicates(subset=['atom_id'], keep='first')
            
            filtered_structures[pdb_id] = combined_df
            filtered_dfs.append(combined_df)
            print(f"{pdb_id}: {len(df_chain)} chain atoms, {len(retinal_df)} retinal atoms")
        else:
            # If no retinal found within cutoff, just keep the chain
            filtered_structures[pdb_id] = df_chain
            filtered_dfs.append(df_chain)
            print(f"{pdb_id}: {len(df_chain)} chain atoms, 0 retinal atoms (none found within {cutoff}Å)")
    
    # Update processor's data with the filtered data
    if filtered_dfs:
        # Stack all filtered dataframes
        stacked_filtered_df = pd.concat(filtered_dfs, ignore_index=True)
        print(f"Updating processor data: {len(processor.data)} atoms → {len(stacked_filtered_df)} atoms")
        
        # Replace the processor's data with the filtered data
        processor.data = stacked_filtered_df
        
        # Update PDB IDs to reflect only those that have filtered data
        processor.pdb_ids = list(filtered_structures.keys())
        print(f"Updated processor with filtered data: {len(processor.pdb_ids)} structures")
    
    return filtered_structures


def load_opsin_property_data(property_file, processed_structures):
    """
    Load property data from CSV file and create structure mappings
    
    Args:
        property_file: Path to CSV file with property data
        processed_structures: Dictionary of processed structures
        
    Returns:
        Dictionary with properties and structure mappings
    """
    print(f"\n=== Loading Property Data from {property_file}===")
    
    try:
        # Read the property CSV file
        df_properties = pd.read_csv(property_file)
        print(f"Loaded property data with {len(df_properties)} entries")
        
        # Clean up property data
        # Replace question marks in function annotations
        if 'molecular_function' in df_properties.columns:
            df_properties['molecular_function'] = df_properties['molecular_function'].apply(
                lambda x: str(x).replace('?', '').strip() if pd.notna(x) else ''
            )
        
        # Create structure property dictionary
        properties = {}
        structure_mapping = {}
        
        # Helper function to clean structure names for predicted structures
        def clean_structure_name(name):
            if not name:
                return name
            # Remove special characters (. or +) and replace - with _
            cleaned = str(name).strip()
            cleaned = cleaned.replace('.', '')
            cleaned = cleaned.replace('+', '')
            cleaned = cleaned.replace('-', '_')
            return cleaned
        
        # Process each property row
        for _, row in df_properties.iterrows():
            row_dict = dict(row)
            
            # Get domain/taxonomic information
            # Primary keys to check for domain (in order of preference)
            domain_keys = ['Rhodopsin Type (Microbial)']
            domain = 'Unknown'
            for key in domain_keys:
                if key in row_dict and pd.notna(row_dict[key]):
                    domain = str(row_dict[key]).strip()
                    if domain:  # Only break if we found a non-empty value
                        break
            
            # Get functional classification
            # Primary keys to check for function (in order of preference)
            function_keys = ['molecular_function']
            molecular_function = 'Unknown'
            for key in function_keys:
                if key in row_dict and pd.notna(row_dict[key]):
                    molecular_function = str(row_dict[key]).strip()
                    if molecular_function:  # Only break if we found a non-empty value
                        break
            
            # Clean up strings
            if isinstance(domain, str):
                domain = domain.replace('?', '').strip()
            if isinstance(molecular_function, str):
                molecular_function = molecular_function.replace('?', '').strip()
            
            # Extract property data
            property_data = {
                'domain': domain,
                'molecular_function': molecular_function
            }
            
            # Process additional properties (e.g., experimental status)
            if 'experimental' in row_dict and pd.notna(row_dict['experimental']):
                property_data['experimentally_determined'] = bool(row_dict['experimental'])
            
            # Check for pdb_id and short_name to build the mapping
            pdb_id = None
            short_name = None
            
            if 'pdb_id' in row_dict and pd.notna(row_dict['pdb_id']):
                pdb_id = str(row_dict['pdb_id']).strip()
            
            if 'short_name' in row_dict and pd.notna(row_dict['short_name']):
                short_name = str(row_dict['short_name']).strip()
                # Clean short_name for predicted structure and append _smile_model_0
                predicted_id = clean_structure_name(short_name)
                # ALWAYS append _smile_model_0 to predicted IDs for consistent naming
                predicted_id_with_suffix = f"{predicted_id}_smile_model_0"
            
            # Create mapping if both pdb_id and short_name exist
            if pdb_id and short_name:
                # Map pdb_id to predicted_id_with_suffix (short_name with _smile_model_0 suffix)
                structure_mapping[pdb_id] = predicted_id_with_suffix
                
                # Add properties for both the pdb_id and the predicted_id_with_suffix
                properties[pdb_id] = property_data.copy()
                properties[predicted_id_with_suffix] = property_data.copy()
                
                # Add properties to processed structures if they exist
                if pdb_id in processed_structures:
                    if 'properties' not in processed_structures[pdb_id]:
                        processed_structures[pdb_id]['properties'] = {}
                    processed_structures[pdb_id]['properties'].update(property_data)
                
                # Check if the predicted structure with suffix exists in processed_structures
                if predicted_id_with_suffix in processed_structures:
                    if 'properties' not in processed_structures[predicted_id_with_suffix]:
                        processed_structures[predicted_id_with_suffix]['properties'] = {}
                    processed_structures[predicted_id_with_suffix]['properties'].update(property_data)
            
            # For entries that have only pdb_id or only short_name, still add properties
            elif pdb_id and pdb_id in processed_structures:
                properties[pdb_id] = property_data.copy()
                if 'properties' not in processed_structures[pdb_id]:
                    processed_structures[pdb_id]['properties'] = {}
                processed_structures[pdb_id]['properties'].update(property_data)
            
            elif short_name:
                # Always use the name with _smile_model_0 suffix for predicted structures
                predicted_id = clean_structure_name(short_name)
                predicted_id_with_suffix = f"{predicted_id}_smile_model_0"
                
                properties[predicted_id_with_suffix] = property_data.copy()
                
                if predicted_id_with_suffix in processed_structures:
                    if 'properties' not in processed_structures[predicted_id_with_suffix]:
                        processed_structures[predicted_id_with_suffix]['properties'] = {}
                    processed_structures[predicted_id_with_suffix]['properties'].update(property_data)
        
        # Filter mapping to include only pairs where both keys and values are in processed_structures
        filtered_mapping = {}
        for exp_id, pred_id in structure_mapping.items():
            if exp_id in processed_structures and (
                pred_id in processed_structures or 
                # Also check for the predicted ID with other common suffixes
                pred_id.replace("_smile_model_0", "_model_0") in processed_structures
            ):
                # If the exact suffix isn't found but a variant is, use the variant
                if pred_id not in processed_structures:
                    variant = pred_id.replace("_smile_model_0", "_model_0")
                    if variant in processed_structures:
                        filtered_mapping[exp_id] = variant
                else:
                    filtered_mapping[exp_id] = pred_id
        
        print(f"Added properties for {len(properties)} structures")
        print(f"Created {len(filtered_mapping)} experimental-predicted structure mappings")
        
        return {
            'properties': properties,
            'structure_mapping': filtered_mapping
        }
    
    except Exception as e:
        print(f"Error loading property data: {e}")
        import traceback
        traceback.print_exc()
        return {'properties': {}, 'structure_mapping': {}}

def rebuild_hideaki_predicted_dataset(model_files, dataset):
    """
    Load predicted structures from a list of prediction file paths,
    create a new dataset and save it as 'mo_hide_pred'.
    
    This function expects each file path to point to a CIF file.
    If the file path already ends with '.cif', it is removed before calling load_structure,
    because load_structure appends the extension.
    
    Args:
        model_files: List of file paths to predicted structure CIF files
        dataset: Name for the new dataset
        
    Returns:
        tuple: (CifBaseProcessor with loaded structures, list of error files)
    """
    cp_pred = CifBaseProcessor()
    all_structs = []
    pdb_ids = []
    errors = []
    for file in model_files:
        try:
            # Remove the extension if it exists.
            if file.lower().endswith('.cif'):
                file_no_ext = file[:-4]
            else:
                file_no_ext = file
            struct = load_structure(file_no_ext, folder="")  # adjust folder if needed
            # Rename any ligand labeled 'LIG' to 'RET'
            struct.loc[struct['res_name3l'].str.upper() == 'LIG', 'res_name3l'] = 'RET'
            pdb_id = extract_pdb_id(file) + '_pred'
            struct['pdb_id'] = pdb_id
            for coord in ['x', 'y', 'z']:
                struct[coord] = struct[coord].astype(float)
            all_structs.append(struct)
            pdb_ids.append(pdb_id)
        except:
            errors.append(file)
    cp_pred.dfl = all_structs
    cp_pred.pdb_ids = pdb_ids
    cp_pred.concat_data()  # Combine list of DataFrames into one
    cp_pred.to_pkl(overwrite=True, folder='data/structures/')
    cp_pred.save_dataset(dataset)
    print(
        f"[INFO] Predicted dataset saved as '{dataset}' with {len(pdb_ids)} structures. Ran into {len(errors)} errors!")
    return cp_pred, errors

def extract_pdb_id(model_file):
    """
    Extract PDB ID from a model file path
    
    Args:
        model_file: Path to the model file
        
    Returns:
        str: Extracted PDB ID
    """
    # Get the last part of the path
    pdb_id = os.path.basename(model_file)
    # If it ends with '.cif', remove the extension
    if pdb_id.lower().endswith('.cif'):
        pdb_id = pdb_id[:-4]
    return pdb_id

def collect_cif_files(root_folder):
    """
    Recursively search for all CIF files in root_folder and its subdirectories.
    
    Args:
        root_folder (str): Path to the root folder.
        
    Returns:
        List[str]: A list of full file paths to CIF files.
    """
    cif_files = []
    for dirpath, dirnames, filenames in os.walk(root_folder):
        for filename in filenames:
            if filename.lower().endswith('.cif'):
                full_path = os.path.join(dirpath, filename)
                cif_files.append(full_path)
    return cif_files


def format_cif_columns(df: pd.DataFrame) -> pd.DataFrame:
    """
    Applies proper data type formatting to CIF DataFrame columns.
    
    This function ensures correct data types for all common columns in CIF files:
    - Converts coordinates (x, y, z, model_x, model_y, model_z) to float
    - Converts sequence and atom IDs to integers (auth_seq_id, label_seq_id, atom_id)
    - Ensures proper string formatting for identifiers (pdb_id, chain_id, etc.)
    - Handles B-factors, occupancy and other numeric fields
    
    Args:
        df: DataFrame containing CIF data with standard column names

    Returns:
        DataFrame with proper data types for all columns
    """
    # Make a copy to avoid modifying the original
    formatted_df = df.copy()
    
    # Define column groups with their target data types
    column_types = {
        # Coordinate columns (float)
        'float_columns': [
            'x', 'y', 'z',                          # Standard coordinates 
            'model_x', 'model_y', 'model_z',        # Model coordinates
            'cartn_x', 'cartn_y', 'cartn_z',        # Alternative names
            'B_iso_or_equiv', 'B_factor',           # B-factors/temperature factors
            'occupancy', 'probability',             # Occupancy-related
            'auth_seq_id', 'label_seq_id',          # Sequence positions as float for calculations
            'atom_id', 'label_entity_id',           # Atom and entity IDs
            'pdbx_formal_charge',                   # Charge
        ],
        # Integer columns
        'int_columns': [
            'auth_seq_id', 'label_seq_id',          # Sequence IDs (cast to int after float)
            'atom_id', 'label_entity_id',           # Atom and entity IDs
            'model_id',                             # Model number
            'pdbx_PDB_model_num',                   # PDB model number
        ],
        # String columns
        'str_columns': [
            'pdb_id', 'struct_id',                  # Structure identifiers
            'auth_asym_id', 'label_asym_id',        # Chain identifiers
            'auth_chain_id', 'label_chain_id',      # Alternative chain identifiers
            'auth_comp_id', 'label_comp_id',        # Component identifiers
            'auth_atom_id', 'label_atom_id',        # Atom identifiers
            'type_symbol', 'element',               # Element symbols
            'group',                                # Atom group (ATOM/HETATM)
            'res_name', 'res_name3l',               # Residue names
            'ins_code', 'pdbx_PDB_ins_code',        # Insertion codes
            'alt_id', 'label_alt_id',               # Alternate location identifiers
            'res_atom_name',                        # Residue atom name
        ]
    }
    
    # Process float columns first (important for those that may later be integers)
    for column in column_types['float_columns']:
        if column in formatted_df.columns:
            try:
                formatted_df[column] = pd.to_numeric(formatted_df[column], errors='coerce').astype('float64')
            except Exception as e:
                print(f"Warning: Could not convert {column} to float: {e}")
    
    # Process integer columns (after float conversion)
    for column in column_types['int_columns']:
        if column in formatted_df.columns:
            try:
                # First ensure it's a float (to handle NaN values properly)
                if not pd.api.types.is_numeric_dtype(formatted_df[column]):
                    formatted_df[column] = pd.to_numeric(formatted_df[column], errors='coerce')
                    
                # Convert non-NaN values to int, keeping NaN as is
                # First create a mask for non-NaN values
                mask = ~pd.isna(formatted_df[column])
                
                # Create a new Series with integers where possible and NaN elsewhere
                int_values = pd.Series(index=formatted_df.index, dtype='Int64')  # nullable integer type
                int_values.loc[mask] = formatted_df.loc[mask, column].astype('int64')
                
                # Assign back to DataFrame
                formatted_df[column] = int_values
            except Exception as e:
                print(f"Warning: Could not convert {column} to integer: {e}")
    
    # Process string columns
    for column in column_types['str_columns']:
        if column in formatted_df.columns:
            try:
                # Convert to string, handling NaN values properly
                formatted_df[column] = formatted_df[column].astype(str)
                # Replace 'nan' with empty string
                formatted_df[column] = formatted_df[column].replace('nan', '')
            except Exception as e:
                print(f"Warning: Could not convert {column} to string: {e}")
    
    # Handle boolean columns if present
    bool_columns = ['is_hetatm', 'is_ligand', 'is_water']
    for column in bool_columns:
        if column in formatted_df.columns:
            try:
                formatted_df[column] = formatted_df[column].astype(bool)
            except Exception as e:
                print(f"Warning: Could not convert {column} to boolean: {e}")
    
    return formatted_df


def ensure_structure_dtypes(structures_dict: Dict[str, Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
    """
    Ensures all structure DataFrames have proper data types.
    
    Processes all structures in a dictionary of structure data, applying
    data type formatting to each DataFrame.
    
    Args:
        structures_dict: Dictionary of structure data with DataFrames
        
    Returns:
        Dictionary with properly formatted structure data
    """
    result = {}
    
    for pdb_id, data in structures_dict.items():
        # Create a copy to avoid modifying original
        result[pdb_id] = data.copy()
        
        # Format the main structure DataFrame
        if 'df' in data and isinstance(data['df'], pd.DataFrame) and not data['df'].empty:
            result[pdb_id]['df'] = format_cif_columns(data['df'])
        
        # Format normalized DataFrames if present
        if 'df_norm' in data and isinstance(data['df_norm'], pd.DataFrame) and not data['df_norm'].empty:
            result[pdb_id]['df_norm'] = format_cif_columns(data['df_norm'])
        
        if 'df_ca_norm' in data and isinstance(data['df_ca_norm'], pd.DataFrame) and not data['df_ca_norm'].empty:
            result[pdb_id]['df_ca_norm'] = format_cif_columns(data['df_ca_norm'])
        
        # Format retinal DataFrame if present
        if 'df_ret' in data and isinstance(data['df_ret'], pd.DataFrame) and not data['df_ret'].empty:
            result[pdb_id]['df_ret'] = format_cif_columns(data['df_ret'])
    
    return result


def format_cif_processor_data(processor) -> None:
    """
    Applies proper data type formatting to CifProcessor's data.
    
    Args:
        processor: CifProcessor or CifBaseProcessor instance
        
    Returns:
        None (modifies processor.data in place)
    """
    if not hasattr(processor, 'data') or processor.data is None or processor.data.empty:
        print("Warning: Processor has no data to format")
        return
    
    # Format the data DataFrame
    processor.data = format_cif_columns(processor.data)
    
    print(f"Formatted data types for {len(processor.pdb_ids)} structures in processor")


def load_opsin_structures(data_dir, output_dir='outputs', chain_id='A', visualize=True, use_cache=True, cache_raw=True,
                          retinal_name='RET', retinal_cutoff=6.0):
    """
    Step 1: Load and process opsin structures from prepared datasets

    This function loads the opsin structures from the prepared datasets:
    - mo_exp: Microbial opsins experimental structures
    - mo_pred: Microbial opsins predicted structures
    - hideaki_exp: Hideaki experimental structures
    - hideaki_pred: Hideaki predicted structures

    It doesn't use PropertyProcessor, relying purely on the CIF files and dataset information.
    Uses a two-stage caching system:
    1. Raw data cache: Stores unfiltered structure data (independent of chain_id)
    2. Processed data cache: Stores chain-specific processed structures

    Args:
        data_dir: Directory where PROTOS data is stored (mmcif files, datasets, etc.)
        output_dir: Path to save analysis results, figures, and cache files
        chain_id: Chain ID to use for processing (default: 'A')
        visualize: Whether to generate visualizations
        use_cache: Whether to use cached structure data (default: True)
        cache_raw: Whether to cache raw unfiltered data (default: True)
        retinal_name: Name of retinal residue (default: 'RET')
        retinal_cutoff: Distance cutoff in Angstroms for retinal selection (default: 6.0)

    Returns:
        Dictionary containing processed structures and processor objects
    """
    print("=== Step 1: Loading Experimental and Predicted Structures ===")

    # Ensure paths are Path objects
    data_dir = Path(data_dir)
    output_dir = Path(output_dir)

    # Set environment variables to force ProtosPaths to use the data_dir for PROTOS data
    os.environ["PROTOS_DATA_ROOT"] = str(data_dir.absolute())
    os.environ["PROTOS_REF_DATA_ROOT"] = str(data_dir.absolute())  # Force reference data to match user data

    print(f"Using PROTOS data directory: {data_dir}")
    print(f"Using analysis outputs directory: {output_dir}")

    # Define cache file paths - store cache in output_dir for better organization
    cache_dir = output_dir / "cache"
    os.makedirs(cache_dir, exist_ok=True)
    print(f"Using cache directory: {cache_dir}")

    # Create necessary directories in the PROTOS data path
    structure_dir = data_dir / "structure"
    mmcif_dir = structure_dir / "mmcif"
    dataset_dir = structure_dir / "structure_dataset"
    standard_dir = dataset_dir / "standard"

    # Create all directories
    for directory in [structure_dir, mmcif_dir, dataset_dir, standard_dir]:
        os.makedirs(directory, exist_ok=True)
        print(f"Ensured PROTOS data directory exists: {directory}")

    # Two-stage caching: raw data and processed data
    raw_cache_file = cache_dir / "raw_structures.pkl"
    processed_cache_file = cache_dir / f"processed_structures_{chain_id}.pkl"
    print(f"Raw cache file: {raw_cache_file}")
    print(f"Processed cache file: {processed_cache_file}")

    # STAGE 1: Load raw datasets (from cache or from source)
    # Check if raw cache exists and should be used
    raw_data = None
    if use_cache and cache_raw and raw_cache_file.exists():
        print(f"Found cached raw structure data at {raw_cache_file}")
        try:
            # Load cached raw data
            with open(raw_cache_file, 'rb') as f:
                raw_data = pickle.load(f)
                print(
                    f"Loaded raw structure data from cache with {sum(1 for p in raw_data.values() if isinstance(p, CifBaseProcessor))} processors")
        except Exception as e:
            print(f"Error loading raw cache file: {e}")
            print("Will load structures from original sources")
            raw_data = None

    # Branch based on whether we're using cached data or loading from scratch
    if raw_data is None:
        # Loading from scratch - initialize paths handler with explicit paths
        print("Creating new processor instances and loading from files")
        paths = ProtosPaths(
            user_data_root=str(data_dir.absolute()),
            ref_data_root=str(data_dir.absolute()),  # Force reference data to be the same as user data
            create_dirs=True
        )

        print(f"Path resolver user_data_root: {paths.user_data_root}")
        print(f"Path resolver ref_data_root: {paths.ref_data_root}")

        # Initialize processors with explicit paths
        cp_mo_exp = CifBaseProcessor(
            name="mo_exp_processor",
            data_root=str(data_dir.absolute()),
            processor_data_dir="structure"
        )

        cp_hide_exp = CifBaseProcessor(
            name="hide_exp_processor",
            data_root=str(data_dir.absolute()),
            processor_data_dir="structure"
        )

        cp_mo_pred = CifBaseProcessor(
            name="mo_pred_processor",
            data_root=str(data_dir.absolute()),
            processor_data_dir="structure"
        )

        cp_hide_pred = CifBaseProcessor(
            name="hide_pred_processor",
            data_root=str(data_dir.absolute()),
            processor_data_dir="structure"
        )

        # Override key paths to ensure consistency
        for processor in [cp_mo_exp, cp_hide_exp, cp_mo_pred, cp_hide_pred]:
            processor.path_structure_dir = os.path.join(str(data_dir.absolute()), "structure", "mmcif")
            processor.path_dataset_dir = os.path.join(str(data_dir.absolute()), "structure", "structure_dataset")
            processor.path_alignment_dir = os.path.join(str(data_dir.absolute()), "structure", "alignments")

        # Set empty datasets initially
        datasets = {}
    else:
        # Using cached data - retrieve processors and datasets
        print("Using processors and datasets from cache")
        cp_mo_exp = raw_data.get('cp_mo_exp')
        cp_mo_pred = raw_data.get('cp_mo_pred')
        cp_hide_exp = raw_data.get('cp_hide_exp')
        cp_hide_pred = raw_data.get('cp_hide_pred')
        datasets = raw_data.get('datasets', {})

        # Verify processors loaded correctly
        for name, processor in [
            ('cp_mo_exp', cp_mo_exp),
            ('cp_mo_pred', cp_mo_pred),
            ('cp_hide_exp', cp_hide_exp),
            ('cp_hide_pred', cp_hide_pred)
        ]:
            if processor is None:
                print(f"WARNING: {name} not found in cache, creating new instance")
                processor = CifBaseProcessor(
                    name=name.replace('cp_', '') + "_processor",
                    data_root=data_dir
                )

    # Create processor dictionary for convenience
    processors = {
        "mo_exp": cp_mo_exp,
        "hideaki_exp": cp_hide_exp,
        "mo_pred": cp_mo_pred,
        "hideaki_pred": cp_hide_pred
    }

    # Check if structure_dataset/standard directory exists and list files
    dataset_standard_dir = structure_dir / "structure_dataset" / "standard"
    if dataset_standard_dir.exists():
        print(f"Standard dataset directory exists: {dataset_standard_dir}")
        json_files = list(dataset_standard_dir.glob("*.json"))
        print(f"Found {len(json_files)} JSON files in standard dataset directory:")
        for json_file in json_files:
            print(f"  - {json_file.name}")
    else:
        print(f"Standard dataset directory not found: {dataset_standard_dir}")
        os.makedirs(dataset_standard_dir, exist_ok=True)
        print(f"Created standard dataset directory")

    # Check if we need to load datasets and structures
    if raw_data is None:
        # Check registry files
        registry_path = structure_dir / "registry.json"
        if registry_path.exists():
            try:
                with open(registry_path, 'r') as f:
                    registry_data = json.load(f)
                    print(f"Registry file exists with {len(registry_data.get('datasets', []))} datasets")
            except Exception as e:
                print(f"Error reading registry file: {e}")
        else:
            print(f"Registry file not found: {registry_path}")
            with open(registry_path, 'w') as f:
                json.dump({"datasets": []}, f)
            print(f"Created empty registry file")

        # Try to load each dataset when loading from scratch
        for dataset_name, processor in processors.items():
            print(f"Loading {dataset_name} dataset...")
            try:
                # First check if the dataset is already in the registry
                available_datasets = []
                if hasattr(processor, 'list_datasets'):
                    available_datasets = processor.list_datasets()

                print(f"  Available datasets in registry: {available_datasets}")

                if dataset_name in available_datasets:
                    # Dataset already in registry, load it with proper data type formatting
                    print(f"  Loading {dataset_name} from registry with data type formatting...")
                    processor.load_dataset(dataset_name, apply_dtypes=True, debug=True)
                    print(f"  Loaded {dataset_name} from registry with {len(processor.pdb_ids)} structures")

                    # Store dataset info
                    datasets[dataset_name] = {
                        "id": dataset_name,
                        "name": dataset_name,
                        "description": "",
                        "pdb_ids": processor.pdb_ids
                    }
                    continue

                # If not in registry, check in multiple locations:
                # 1. First in the current output_dir
                # 2. Then in the default opsin_output location
                # 3. Finally in the protos reference data

                # Option 1: Check in current outputs directory
                dataset_path = structure_dir / "structure_dataset" / "standard" / f"{dataset_name}.json"

                # Option 2: Check in default opsin_output location
                if not dataset_path.exists():
                    default_path = Path(__file__).resolve().parent / "opsin_output" / "structure" / "structure_dataset" / "standard" / f"{dataset_name}.json"
                    if default_path.exists():
                        dataset_path = default_path

                # Option 3: Check in projects/opsin_analysis/data
                if not dataset_path.exists():
                    # Check in data directory within opsin_analysis
                    data_path = Path(__file__).resolve().parent / "data" / "structure" / "structure_dataset" / "standard" / f"{dataset_name}.json"
                    if data_path.exists():
                        dataset_path = data_path

                print(f"  Looking for dataset file at: {dataset_path}")
                print(f"  File exists: {dataset_path.exists()}")

                if dataset_path.exists():
                    print(f"  Loading dataset from {dataset_path}")
                    # Load content from JSON file directly
                    with open(dataset_path, 'r') as f:
                        dataset_data = json.load(f)

                    pdb_ids = dataset_data.get("pdb_ids", [])
                    print(f"  Found {len(pdb_ids)} PDB IDs in dataset file")

                    # Create the dataset properly in the processor to register it
                    if hasattr(processor, 'create_dataset'):
                        processor.create_dataset(
                            dataset_id=dataset_name,
                            name=dataset_data.get("name", dataset_name),
                            description=dataset_data.get("description", ""),
                            content=pdb_ids
                        )
                        print(f"  Created dataset {dataset_name} with {len(pdb_ids)} structures in registry")

                        # Also update the registry manually if needed
                        if hasattr(processor, 'registry'):
                            if 'datasets' not in processor.registry:
                                processor.registry['datasets'] = []

                            # Add dataset to registry if not already there
                            dataset_ids = [d.get('id') for d in processor.registry['datasets'] if 'id' in d]
                            if dataset_name not in dataset_ids:
                                processor.registry['datasets'].append({
                                    'id': dataset_name,
                                    'path': str(dataset_path),
                                    'type': 'json'
                                })
                                print(f"  Added {dataset_name} to registry manually")

                                # Save the updated registry
                                registry_path = structure_dir / "registry.json"
                                try:
                                    with open(registry_path, 'w') as f:
                                        json.dump(processor.registry, f, indent=2)
                                    print(f"  Updated registry file at: {registry_path}")
                                except Exception as e:
                                    print(f"  Error saving registry: {e}")

                    # Load the structures - check mmcif directory first
                    mmcif_dir = structure_dir / "mmcif"
                    if mmcif_dir.exists():
                        cif_files = list(mmcif_dir.glob("*.cif"))
                        print(f"  Found {len(cif_files)} CIF files in mmcif directory")

                    # If the dataset comes from another location, make sure paths are properly set
                    if hasattr(processor, 'paths') and hasattr(processor.paths, 'mmcif_dir'):
                        # Set the correct mmcif directory in processor
                        processor.paths.mmcif_dir = str(mmcif_dir)
                        print(f"  Updated mmcif directory path in processor: {processor.paths.mmcif_dir}")

                    # Check if we need to copy any missing CIF files from alternate locations
                    # 1. Check if we have all CIF files or if some are missing
                    missing_cifs = []
                    for pdb_id in pdb_ids:
                        cif_path = mmcif_dir / f"{pdb_id.lower()}.cif"
                        if not cif_path.exists():
                            missing_cifs.append(pdb_id)

                    if missing_cifs:
                        print(f"  Missing {len(missing_cifs)} CIF files: {missing_cifs[:5]}...")

                        # Check reference data in projects/opsin_analysis/data
                        ref_mmcif_dir = Path(__file__).resolve().parent / "data" / "structure" / "mmcif"
                        if ref_mmcif_dir.exists():
                            print(f"  Checking reference data at: {ref_mmcif_dir}")
                            copied_count = 0
                            for pdb_id in missing_cifs:
                                ref_cif_path = ref_mmcif_dir / f"{pdb_id.lower()}.cif"
                                if ref_cif_path.exists():
                                    target_path = mmcif_dir / f"{pdb_id.lower()}.cif"
                                    try:
                                        with open(ref_cif_path, 'r') as src, open(target_path, 'w') as dst:
                                            dst.write(src.read())
                                        copied_count += 1
                                    except Exception as e:
                                        print(f"  Warning: Failed to copy {pdb_id}: {e}")

                            if copied_count > 0:
                                print(f"  Copied {copied_count} CIF files from reference data")

                    # Load the structures with proper data type formatting
                    print(f"  Loading structures with data type formatting...")
                    processor.load_structures(pdb_ids, apply_dtypes=True, debug=True)
                    print(f"  Loaded {len(processor.pdb_ids)}/{len(pdb_ids)} structures for {dataset_name}")

                    # Store dataset info
                    datasets[dataset_name] = {
                        "id": dataset_name,
                        "name": dataset_data.get("name", dataset_name),
                        "description": dataset_data.get("description", ""),
                        "pdb_ids": pdb_ids
                    }

                    # If the dataset file doesn't exist in the outputs directory, copy it there
                    output_dataset_path = structure_dir / "structure_dataset" / "standard" / f"{dataset_name}.json"
                    if dataset_path != output_dataset_path and not output_dataset_path.exists():
                        try:
                            print(f"  Copying dataset file to outputs directory: {output_dataset_path}")
                            os.makedirs(output_dataset_path.parent, exist_ok=True)
                            with open(dataset_path, 'r') as src, open(output_dataset_path, 'w') as dst:
                                dst.write(src.read())
                        except Exception as e:
                            print(f"  Warning: Could not copy dataset file: {e}")
                else:
                    print(f"  No standard dataset file found for {dataset_name}")
                    datasets[dataset_name] = {
                        "id": dataset_name,
                        "name": dataset_name,
                        "description": "",
                        "pdb_ids": []
                    }

            except Exception as e:
                print(f"  Error loading {dataset_name}: {str(e)}")
                # Create empty dataset information
                datasets[dataset_name] = {
                    "id": dataset_name,
                    "name": dataset_name,
                    "description": "",
                    "pdb_ids": []
                }
    else:
        # Using cached data - processors and datasets are already loaded
        for name, processor in processors.items():
            if processor and hasattr(processor, 'pdb_ids'):
                print(f"Using cached {name} processor with {len(processor.pdb_ids)} structures")
            else:
                print(f"WARNING: Cached {name} processor has no structures")

    # Fix LYR
    process_lyr_in_processor_data(cp_mo_exp, retinal_res_name='RET')

    # Filter structures by chain and retinal
    filtered_structures = {}

    # Rename 'LIG' to 'RET'
    cp_mo_pred.data.loc[cp_mo_pred.data['res_name3l'] == 'LIG', 'res_name3l'] = 'RET'
    cp_hide_pred.data.loc[cp_hide_pred.data['res_name3l'] == 'LIG', 'res_name3l'] = 'RET'

    # set retinal chain to A
    cp_mo_pred.data.loc[cp_mo_pred.data['res_name3l'] == retinal_name, 'auth_chain_id'] = chain_id
    cp_hide_pred.data.loc[cp_hide_pred.data['res_name3l'] == retinal_name, 'auth_chain_id'] = chain_id

    for dataset_name, processor in processors.items():
        filtered_structures[dataset_name] = filter_structures_by_chain_and_retinal(processor, chain=chain_id,
                                                                                   retinal_name='RET',
                                                                                   cutoff=retinal_cutoff)



    # If we loaded from scratch (not from cache), save the processors to cache
    if raw_data is None and cache_raw:
        # Create fresh raw data dictionary for saving to cache
        raw_data_to_save = {
            'cp_mo_exp': cp_mo_exp,
            'cp_mo_pred': cp_mo_pred,
            'cp_hide_exp': cp_hide_exp,
            'cp_hide_pred': cp_hide_pred,
            'datasets': datasets
        }

        # Save raw data to cache - this is independent of chain_id and can be reused
        try:
            print(f"Saving raw structure data to cache: {raw_cache_file}")
            with open(raw_cache_file, 'wb') as f:
                pickle.dump(raw_data_to_save, f)
            print(f"Successfully saved raw structure data to cache")
        except Exception as e:
            print(f"Warning: Failed to save raw cache file: {e}")

    # Print summary of loaded structures
    print("\nStructure Loading Summary:")
    for name, processor in processors.items():
        if processor and hasattr(processor, 'pdb_ids'):
            print(f"  {name}: {len(processor.pdb_ids)} structures")
        else:
            print(f"  {name}: No structures")

    # STAGE 2: Process structures for this specific chain ID
    # Check if processed cache exists and should be used
    if use_cache and processed_cache_file.exists():
        print(f"Found cached processed structure data at {processed_cache_file}")
        try:
            # Load cached processed data
            with open(processed_cache_file, 'rb') as f:
                processed_data = pickle.load(f)
                print(f"Loaded {len(processed_data.get('processed_structures', {}))} processed structures from cache")

                # Ensure all structures have coordinates as float
                if 'processed_structures' in processed_data:
                    # Use the imported format_cif_columns function from data_processing.py
                    processed_data['processed_structures'] = ensure_structure_dtypes(
                        processed_data['processed_structures'])
                    print("Converted coordinate data types to float for all structures")

                # Combine with the raw processors and datasets
                result = {
                    'processed_structures': processed_data.get('processed_structures', {}),
                    'structure_mapping': processed_data.get('structure_mapping', {}),
                    'cp_mo_exp': cp_mo_exp,
                    'cp_mo_pred': cp_mo_pred,
                    'cp_hide_exp': cp_hide_exp,
                    'cp_hide_pred': cp_hide_pred,
                    'datasets': datasets
                }

                return result
        except Exception as e:
            print(f"Error loading processed cache file: {e}")
            print("Will process structures from raw data")

    # Create processor dictionary for filtering
    processors = {
        "mo_exp": cp_mo_exp,
        "hideaki_exp": cp_hide_exp,
        "mo_pred": cp_mo_pred,
        "hideaki_pred": cp_hide_pred
    }

    # Rename 'LIG' to 'RET' in predicted structures
    if hasattr(cp_mo_pred, 'data') and cp_mo_pred.data is not None and not cp_mo_pred.data.empty:
        cp_mo_pred.data.loc[cp_mo_pred.data['res_name3l'] == 'LIG', 'res_name3l'] = retinal_name
        print(f"Renamed LIG to {retinal_name} in mo_pred dataset")

    if hasattr(cp_hide_pred, 'data') and cp_hide_pred.data is not None and not cp_hide_pred.data.empty:
        cp_hide_pred.data.loc[cp_hide_pred.data['res_name3l'] == 'LIG', 'res_name3l'] = retinal_name
        print(f"Renamed LIG to {retinal_name} in hideaki_pred dataset")

    # Set retinal chain to match the specified chain_id
    for processor in [cp_mo_pred, cp_hide_pred, cp_mo_exp, cp_hide_exp]:
        if hasattr(processor, 'data') and processor.data is not None and not processor.data.empty:
            processor.data.loc[processor.data['res_name3l'] == retinal_name, 'auth_chain_id'] = chain_id

    # Print retinal counts for debugging
    for name, processor in processors.items():
        if hasattr(processor, 'data') and processor.data is not None and not processor.data.empty:
            ret_count = len(processor.data[processor.data['res_name3l'] == retinal_name])
            print(f"Found {ret_count} {retinal_name} atoms in {name} dataset")

    # Filter structures by chain and retinal
    filtered_structures = {}

    for dataset_name, processor in processors.items():
        filtered_structures[dataset_name] = filter_structures_by_chain_and_retinal(processor,
                                                                                   chain=chain_id,
                                                                                   retinal_name=retinal_name,
                                                                                   cutoff=retinal_cutoff)

    # Combine all structures for processing
    all_data = []
    all_pdb_ids = []

    # Add data from all processors if they have data
    for processor in processors.values():
        if hasattr(processor, 'data') and processor.data is not None and not processor.data.empty:
            all_data.append(processor.data)
            all_pdb_ids.extend(processor.pdb_ids)

    # Create combined processor
    cp_complete = CifBaseProcessor(
        name="combined_processor",
        data_root=data_dir
    )

    if all_data:
        cp_complete.data = pd.concat(all_data)
        cp_complete.pdb_ids = all_pdb_ids

    # Print statistics
    print(f"\nLoaded structures:")
    print(f"  - Experimental MO: {len(cp_mo_exp.pdb_ids)}")
    print(f"  - Predicted MO: {len(cp_mo_pred.pdb_ids)}")
    print(f"  - Experimental Hideaki: {len(cp_hide_exp.pdb_ids)}")
    print(f"  - Predicted Hideaki: {len(cp_hide_pred.pdb_ids)}")
    print(f"  - Total: {len(cp_complete.pdb_ids)}")

    # Process all structures
    processed_structures_complete = {}

    # Extract structures
    for pdb_id in tqdm(cp_complete.pdb_ids, desc="Processing structures"):
        df_pdb = cp_complete.data[cp_complete.data['pdb_id'] == pdb_id]
        df_chain = df_pdb[df_pdb['auth_chain_id'] == chain_id].copy()

        if df_chain.empty:
            print(f"  [WARNING] {pdb_id}: Chain {chain_id} not found, skipping.")
            continue

        # Select retinal (RET) by finding RET residues
        try:
            # Create a safer implementation of retinal selection
            df_ret = df_chain[df_chain['res_name3l'] == 'RET'].copy()

            if not df_ret.empty:
                # Count unique RET residues (by auth_seq_id)
                ret_residues = df_ret['auth_seq_id'].unique() if 'auth_seq_id' in df_ret.columns else []

                if len(ret_residues) == 1:
                    print(f"[DEBUG] Exactly one RET with {len(df_ret)} atoms found in the entire structure.")
                    chosen_ret = df_ret
                elif len(ret_residues) > 1:
                    # Multiple RET molecules, find closest to protein chain
                    print(f"[DEBUG] Multiple RET molecules found, selecting closest to chain {chain_id}")

                    # For each RET residue, calculate the average distance to protein chain
                    distances = []
                    protein_atoms = df_chain[df_chain['res_name3l'] != 'RET']

                    if not protein_atoms.empty:
                        for ret_id in ret_residues:
                            ret_atoms = df_ret[df_ret['auth_seq_id'] == ret_id]

                            # Simple distance calculation: pick center of RET
                            ret_center = ret_atoms[['x', 'y', 'z']].mean()

                            # Get distances to protein atoms
                            dists = np.sqrt(((protein_atoms[['x', 'y', 'z']] - ret_center) ** 2).sum(axis=1))
                            min_dist = dists.min()
                            distances.append((ret_id, min_dist))

                        # Sort by distance
                        distances.sort(key=lambda x: x[1])
                        closest_ret_id = distances[0][0]
                        chosen_ret = df_ret[df_ret['auth_seq_id'] == closest_ret_id].copy()
                    else:
                        # No protein atoms? Just pick the first RET
                        chosen_ret = df_ret[df_ret['auth_seq_id'] == ret_residues[0]].copy()
                else:
                    chosen_ret = pd.DataFrame()
            else:
                chosen_ret = pd.DataFrame()
        except Exception as e:
            print(f"  [ERROR] Error selecting retinal for {pdb_id}: {str(e)}")
            import traceback
            traceback.print_exc()
            chosen_ret = pd.DataFrame()

        processed_structures_complete[pdb_id] = {
            'chain_id': chain_id,
            'df': df_chain,
            'df_ret': chosen_ret if chosen_ret is not None else pd.DataFrame()
        }

    # Filter ligands other than RET
    for pdb_id, data in processed_structures_complete.items():
        df = data['df']
        is_atom_record = (df['group'] != 'HETATM')
        is_ret_hetatm = (df['group'] == 'HETATM') & (df['res_name3l'] == 'RET')
        is_lyr_hetatm = (df['group'] == 'HETATM') & (df['res_name3l'] == 'LYR')
        keep_condition = is_atom_record | is_ret_hetatm | is_lyr_hetatm
        filtered_df = df[keep_condition]
        processed_structures_complete[pdb_id]['df'] = filtered_df

    # Set up df_norm and df_ca_norm with proper type conversion
    for pdb_id, data in processed_structures_complete.items():
        if 'df' not in data or data['df'].empty:
            continue
        df = data['df']
        df_norm = df.copy()

        # Ensure coordinates are numeric
        for coord in ['x', 'y', 'z']:
            if coord in df_norm.columns:
                try:
                    df_norm[coord] = pd.to_numeric(df_norm[coord], errors='coerce')
                except Exception as e:
                    print(f"Warning: Could not convert {coord} to numeric for {pdb_id}: {str(e)}")

        data['df_norm'] = df_norm

        # Extract CA atoms with proper type conversion
        df_ca_norm = df_norm[df_norm['res_atom_name'] == 'CA'].copy()

        # Verify coordinate types and report any issues
        for coord in ['x', 'y', 'z']:
            if coord in df_ca_norm.columns:
                if not pd.api.types.is_numeric_dtype(df_ca_norm[coord]):
                    print(f"Warning: {coord} is not numeric in CA atoms for {pdb_id}")
                    print(f"Sample values: {df_ca_norm[coord].head()}")
                    # Force conversion
                    df_ca_norm[coord] = pd.to_numeric(df_ca_norm[coord], errors='coerce')

        data['df_ca_norm'] = df_ca_norm

        # Add metadata from filename to enable proper classification
        if '_pred' in pdb_id:
            data['structure_type'] = 'predicted'
            base_name = pdb_id.replace('_pred', '')
            data['base_name'] = base_name
        else:
            data['structure_type'] = 'experimental'
            data['base_name'] = pdb_id

    # Create a mapping for structure pairs (experimental-predicted)
    structure_mapping = {}

    # For structure selection in step 5 (compare_structures), we'll use:
    # 1. Experimental structures for hideaki's set
    # 2. Experimental MO structures where available
    # 3. Predicted MO structures for those lacking experimental counterparts

    # Example: If we have "A1ACR1_J318_refine3" and "A1ACR1_J318_refine3_pred" or "A1ACR1_J318_refine3_model_0",
    # we'll map them as a pair

    # First, organize structures into predicted vs experimental
    experimental_structures = {}
    predicted_structures = {}

    for pdb_id, data in processed_structures_complete.items():
        is_predicted = False
        base_name = pdb_id

        # Check if this is a predicted structure by filename patterns
        if any(pattern in pdb_id for pattern in ['_pred', '_model_0', '_smile']):
            is_predicted = True
            # Extract base name - try different patterns
            if '_pred' in pdb_id:
                base_name = pdb_id.replace('_pred', '')
            elif '_model_0' in pdb_id:
                base_name = pdb_id.replace('_model_0', '')
                if '_smile' in base_name:
                    base_name = base_name.replace('_smile', '')

        # Store in appropriate dictionary
        if is_predicted:
            predicted_structures[pdb_id] = {
                'base_name': base_name,
                'data': data
            }
        else:
            experimental_structures[pdb_id] = {
                'base_name': base_name,
                'data': data
            }

    print(
        f"DEBUG: Found {len(experimental_structures)} experimental and {len(predicted_structures)} predicted structures")

    # Create mappings between experimental and predicted structures
    pairs_found = 0
    for exp_id, exp_info in experimental_structures.items():
        exp_base = exp_info['base_name']

        # Direct matches first
        matched = False
        for pred_id, pred_info in predicted_structures.items():
            pred_base = pred_info['base_name']

            # Check if bases match exactly or are very similar
            if exp_base == pred_base or (
                    len(exp_base) > 3 and len(pred_base) > 3 and
                    (exp_base in pred_base or pred_base in exp_base)
            ):
                structure_mapping[exp_id] = pred_id
                pairs_found += 1
                matched = True
                break

        # If no direct match, try matching base parts
        if not matched:
            exp_parts = exp_base.split('_')
            if len(exp_parts) >= 1:
                exp_first_part = exp_parts[0]
                for pred_id, pred_info in predicted_structures.items():
                    if exp_first_part in pred_id:
                        structure_mapping[exp_id] = pred_id
                        pairs_found += 1
                        break

    print(f"DEBUG: Created {pairs_found} experimental-predicted structure pairs by base name matching")

    # Print some example pairs for debugging
    if structure_mapping:
        print("DEBUG: Example structure pairs:")
        for i, (exp_id, pred_id) in enumerate(structure_mapping.items()):
            if i < 5:  # Show up to 5 examples
                print(f"  - {exp_id} -> {pred_id}")

    # Create processed result dictionary (chain-specific)
    processed_data = {
        'processed_structures': processed_structures_complete,
        'structure_mapping': structure_mapping
    }

    # Save processed result to cache file
    try:
        print(f"Saving processed structures to cache: {processed_cache_file}")
        with open(processed_cache_file, 'wb') as f:
            pickle.dump(processed_data, f)
        print(f"Successfully saved {len(processed_structures_complete)} processed structures to cache")
    except Exception as e:
        print(f"Warning: Failed to save processed cache file: {e}")

    # Create the complete result by combining raw and processed data
    result = {
        'processed_structures': processed_structures_complete,
        'structure_mapping': structure_mapping,
        'cp_mo_exp': cp_mo_exp,
        'cp_mo_pred': cp_mo_pred,
        'cp_hide_exp': cp_hide_exp,
        'cp_hide_pred': cp_hide_pred,
        'datasets': datasets
    }

    return result
