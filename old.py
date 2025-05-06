import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from scipy.spatial.distance import cdist, squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import tempfile
import pickle
import json
from pathlib import Path
import csv

# Import from protos package for structure processing
from protos.io.paths.path_config import ProtosPaths, DataSource
from protos.processing.structure.struct_base_processor import CifBaseProcessor

# Import from protos package for opsin analysis
try:
    from protos.processing.opsin.retinal_utils import calculate_min_distances, scale_sizes
except ImportError:
    print("Warning: Could not import retinal_utils, some functions may not work")

# Import visualization functions
try:
    from projects.opsin_analysis.visualization_functions import (
        plot_distances_with_std,
        create_residue_conservation_plot,
        plot_average_distances_by_helix,
        plot_distance_heatmap,
        # Structure Visualization
        visualize_single_7tm_bundle,

        # Enhanced RMSD Visualizations
        visualize_rmsd_heatmap,
        create_and_visualize_similarity_tree,
        visualize_rmsd_matrix_improved,

        # Sequence Analysis Visualizations
        print_residue_composition,
        plot_average_distances_by_helix,
        plot_distance_heatmap,
        create_residue_conservation_plot,
        create_sequence_logo,

        visualize_msa_distances
    )
except ImportError:
    print("Warning: Could not import visualization_functions, some visualizations may not work")

try:
    from protos.visualization.structure_vis import (
        structure_vis,
        visualize_protein_with_ret,
        plot_aligned_structures
    )
except ImportError:
    print("Warning: Could not import structure_vis, visualizations may not work")

# Import FoldMason helper functions
try:
    from projects.opsin_analysis.foldmason_helpers import (
        run_easy_msa,
        run_structuremsa,
        run_refinemsa,
        run_msa2lddtreport,
        run_structuremsacluster
    )
except ImportError:
    print("Warning: Could not import foldmason_helpers, FoldMason integration may not work")

# Import opsin analysis functions from modular structure
# Error Analysis
from projects.opsin_analysis.error_analysis import (
    make_rmsd_table
)

# Structure Comparison
from projects.opsin_analysis.structure_comparison import (
    compute_all_vs_all_rmsd_improved,
    calculate_binding_pocket_rmsd_for_pairs
)

# Reference Alignment
from projects.opsin_analysis.reference_alignment import (
    find_type_references,
    find_global_reference,
    create_seq_alignment_dicts_from_paths
)

# MSA and GRN
from projects.opsin_analysis.msa_grn import (
    analyze_residue_composition,
    generate_grn_msa_tables,
    calculate_helix_distances
)


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
    from projects.opsin_analysis.error_analysis import find_retinal_within_cutoff

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


def load_opsin_structures(output_dir='output', chain_id='A', visualize=True, use_cache=True, cache_raw=True,
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
        output_dir: Path to save output files
        chain_id: Chain ID to use for processing (default: 'A')
        visualize: Whether to generate visualizations
        use_cache: Whether to use cached structure data (default: True)
        cache_raw: Whether to cache raw unfiltered data (default: True)

    Returns:
        Dictionary containing processed structures and processor objects
    """
    # Create output directory if it doesn't exist

    print("=== Step 1: Loading Experimental and Predicted Structures ===")

    # Set up data paths to match prepare_data.py
    # Use the same directory where prepare_data.py stored the datasets
    data_dir = Path(__file__).resolve().parent / "opsin_output"
    output_dir = data_dir
    os.makedirs(data_dir, exist_ok=True)
    print(f"Using data root: {data_dir}")

    # Define cache file paths - store cache in output_dir for better organization
    cache_dir = Path(output_dir) / "cache"
    os.makedirs(cache_dir, exist_ok=True)
    print(f"Using cache directory: {cache_dir}")

    # Create necessary directories
    structure_dir = data_dir / "structure"
    mmcif_dir = structure_dir / "mmcif"
    dataset_dir = structure_dir / "structure_dataset"
    standard_dir = dataset_dir / "standard"

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

    # Set up common directory paths
    structure_dir = data_dir / "structure"
    mmcif_dir = structure_dir / "mmcif"
    dataset_dir = structure_dir / "structure_dataset"
    standard_dir = dataset_dir / "standard"

    # Create necessary directories
    for directory in [structure_dir, mmcif_dir, dataset_dir, standard_dir]:
        os.makedirs(directory, exist_ok=True)

    # Branch based on whether we're using cached data or loading from scratch
    if raw_data is None:
        # Loading from scratch - initialize paths handler
        print("Creating new processor instances and loading from files")
        paths = ProtosPaths(
            user_data_root=data_dir,
            create_dirs=True
        )

        # Initialize processors for different datasets
        cp_mo_exp = CifBaseProcessor(
            name="mo_exp_processor",
            data_root=data_dir
        )

        cp_hide_exp = CifBaseProcessor(
            name="hide_exp_processor",
            data_root=data_dir
        )

        cp_mo_pred = CifBaseProcessor(
            name="mo_pred_processor",
            data_root=data_dir
        )

        cp_hide_pred = CifBaseProcessor(
            name="hide_pred_processor",
            data_root=data_dir
        )

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

                # If not in registry, try to load directly from JSON file
                dataset_path = structure_dir / "structure_dataset" / "standard" / f"{dataset_name}.json"

                print(f"  Looking for dataset file: {dataset_path}")
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

                    # Load the structures - check mmcif directory first
                    mmcif_dir = structure_dir / "mmcif"
                    if mmcif_dir.exists():
                        cif_files = list(mmcif_dir.glob("*.cif"))
                        print(f"  Found {len(cif_files)} CIF files in mmcif directory")

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
        filtered_df = df[~((df['group'] == 'HETATM') & (df['res_name3l'] != 'RET'))]
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


def calculate_structure_errors(data_dict, output_dir='output', visualize=True):
    """
    Step 2: Calculate errors between experimental and predicted structures

    This function calculates RMSD errors between paired experimental and
    predicted structures without using PropertyProcessor.

    Args:
        data_dict: Dictionary with data from previous step
        output_dir: Directory to save output files
        visualize: Whether to generate visualizations

    Returns:
        Dictionary with error data and updated processed structures
    """
    print("\n=== Step 2: Error Calculation ===")

    # Create output directory if it doesn't exist
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


def ensure_structure_dtypes(structures):
    """
    Ensures all structure dataframes have correct data types.

    Args:
        structures: Dictionary of structures

    Returns:
        dict: Updated dictionary with correct data types
    """
    for pdb_id, data in structures.items():
        print(f"[DEBUG] Converting data types for {pdb_id}...")
        for struct_key in ['df', 'df_norm', 'df_ca_norm']:
            if struct_key in data:
                df = data[struct_key]
                # Convert coordinate columns to float
                for col in ['x', 'y', 'z']:
                    if col in df.columns:
                        try:
                            df[col] = pd.to_numeric(df[col], errors='coerce')
                        except Exception as e:
                            print(f"[WARNING] Error converting {col} in {pdb_id}: {e}")
                # Convert auth_seq_id to appropriate type (could be int or string)
                if 'auth_seq_id' in df.columns:
                    # Some auth_seq_id values might have insertion codes (e.g., '123A')
                    # so we'll keep those as strings, but convert numeric ones to int when possible
                    try:
                        # First check if all values can be converted to numeric
                        df['auth_seq_id_numeric'] = pd.to_numeric(df['auth_seq_id'], errors='coerce')
                        # If no NaNs, then all values are numeric and can be integers
                        if not df['auth_seq_id_numeric'].isna().any():
                            df['auth_seq_id'] = df['auth_seq_id_numeric'].astype(int)
                        # Drop the temporary column
                        df = df.drop(columns=['auth_seq_id_numeric'])
                    except:
                        # If conversion fails, keep as is
                        pass
                # Apply changes back to the structure data
                data[struct_key] = df
    return structures


def define_reference_helices(reference_structure, helix_ref_file=None):
    """
    Define the helix boundaries for the reference structure.
    This can either load from a file or use hardcoded values.

    Args:
        reference_structure: The reference structure data
        helix_ref_file: Optional path to a JSON file containing helix definitions

    Returns:
        Dictionary of helix definitions {helix_num: {'start': pos, 'end': pos}}
    """
    import json
    from pathlib import Path

    # First, try to load from file if provided
    if helix_ref_file and Path(helix_ref_file).exists():
        try:
            with open(helix_ref_file, 'r') as f:
                helix_data = json.load(f)

            # Check if the reference structure ID is in the helix data
            ref_id = next(iter(helix_data.keys()))

            # Convert to standard format
            helix_defs = {}
            for helix_num, bounds in helix_data[ref_id].items():
                if isinstance(bounds, list) and len(bounds) == 2:
                    helix_defs[helix_num] = {'start': bounds[0], 'end': bounds[1]}

            if helix_defs:
                print(f"[INFO] Loaded helix definitions from {helix_ref_file}: {helix_defs}")
                return helix_defs
        except Exception as e:
            print(f"[WARNING] Error loading helix definitions from file: {e}")

    # Fallback: Use hardcoded definitions for CnChR2_J230_refine9
    # These values are taken from the JSON file we fixed earlier
    ref_id = reference_structure.get('pdb_id', 'unknown')
    print(f"[INFO] Using hardcoded helix definitions for {ref_id}")

    if ref_id == 'CnChR2_J230_refine9':
        return {
            '1': {'start': 88, 'end': 111},
            '2': {'start': 120, 'end': 136},
            '3': {'start': 157, 'end': 174},
            '4': {'start': 188, 'end': 206},
            '5': {'start': 214, 'end': 232},
            '6': {'start': 250, 'end': 269},
            '7': {'start': 285, 'end': 304}
        }

    # If no valid definitions, return empty dict
    print(f"[WARNING] No helix definitions available for {ref_id}")
    return {}


def annotate_helices_from_alignments(processed_structures, reference_id, helix_definitions,
                                     alignment_paths, chain_id='A'):
    """
    Annotate helices in all structures based on alignments to reference structure.

    Args:
        processed_structures: Dictionary of processed structures
        reference_id: ID of reference structure with defined helices
        helix_definitions: Dictionary of helix boundaries in reference
        alignment_paths: Paths from structure_comparison.get_structure_alignment
        chain_id: Chain ID to use

    Returns:
        Updated processed_structures with helix annotations
    """
    print(f"[INFO] Annotating helices based on reference structure: {reference_id}")

    if reference_id not in processed_structures:
        print(f"[ERROR] Reference structure {reference_id} not found in processed structures")
        return processed_structures

    # Get reference structure data
    ref_structure = processed_structures[reference_id]

    # Extract mapping from auth_seq_id to position for reference structure
    ref_df = None
    if 'df_norm' in ref_structure:
        ref_df = ref_structure['df_norm']
    elif 'df' in ref_structure:
        ref_df = ref_structure['df']

    if ref_df is None or ref_df.empty:
        print(f"[ERROR] No dataframe found for reference structure {reference_id}")
        return processed_structures

    # Filter for the specified chain
    ref_df = ref_df[ref_df['auth_chain_id'] == chain_id]

    # Filter for CA atoms
    ref_ca_df = ref_df[ref_df['res_atom_name'] == 'CA']
    if ref_ca_df.empty:
        try:
            # Try string comparison
            ref_ca_df = ref_df[ref_df['res_atom_name'].astype(str) == 'CA']
        except:
            pass

    if ref_ca_df.empty:
        print(f"[ERROR] No CA atoms found in reference structure {reference_id}")
        return processed_structures

    # Create mapping from auth_seq_id to indices for reference
    ref_seq_ids = ref_ca_df['auth_seq_id'].values
    ref_seq_map = {seq_id: i for i, seq_id in enumerate(ref_seq_ids)}

    # Create helix assignments for the reference structure
    helices_by_residue = {}
    for helix_num, helix_info in helix_definitions.items():
        start_pos = helix_info['start']
        end_pos = helix_info['end']

        # Assign all residues in range to this helix
        for res_id in range(start_pos, end_pos + 1):
            helices_by_residue[res_id] = int(helix_num)

    # Store helix assignments in the reference structure
    ref_structure['helix_assignments'] = helices_by_residue
    ref_structure['helix_definitions'] = helix_definitions

    # Now annotate all other structures using alignment paths
    for struct_id, structure in processed_structures.items():
        if struct_id == reference_id:
            continue  # Skip reference, already annotated

        # Get alignment path between reference and this structure
        if (reference_id, struct_id) in alignment_paths:
            path_info = alignment_paths[(reference_id, struct_id)]
        elif (struct_id, reference_id) in alignment_paths:
            # If found in reverse order, we need to flip the mapping
            path_info = alignment_paths[(struct_id, reference_id)]
            # Flip the residue mapping (B->A becomes A->B)
            if 'residue_mapping' in path_info:
                path_info['residue_mapping'] = [(b, a) for a, b in path_info['residue_mapping']]
        else:
            print(f"[WARNING] No alignment path found between {reference_id} and {struct_id}")
            continue

        # Extract residue mapping
        if 'residue_mapping' not in path_info:
            print(f"[WARNING] No residue mapping in alignment path for {struct_id}")
            continue

        residue_mapping = path_info['residue_mapping']

        # Create reverse mapping (reference -> structure)
        target_helix_assignments = {}

        for ref_res_id, target_res_id in residue_mapping:
            # Check if the reference residue has a helix assignment
            if ref_res_id in helices_by_residue:
                helix_num = helices_by_residue[ref_res_id]
                target_helix_assignments[target_res_id] = helix_num

        # Store the helix assignments in the structure
        structure['helix_assignments'] = target_helix_assignments

        # Create helix definitions for the target structure
        target_helix_defs = {}
        for helix_num in range(1, 8):  # Assume 7TM helices
            helix_residues = [res_id for res_id, h_num in target_helix_assignments.items()
                              if h_num == helix_num]
            if helix_residues:
                target_helix_defs[str(helix_num)] = {
                    'start': min(helix_residues),
                    'end': max(helix_residues)
                }

        structure['helix_definitions'] = target_helix_defs
        print(f"[INFO] Annotated {len(target_helix_assignments)} residues with helix assignments in {struct_id}")

        # Add a tm_helices list for compatibility with visualization
        structure['tm_helices'] = list(range(1, 8))

    return processed_structures


def orient_and_annotate_structures(data_dict, output_dir='output', visualize=True):
    """
    Step 3 & 4: Orient structures and annotate helices using alignments

    Instead of using PCA-based orientation and algorithmic helix finding,
    this function uses structure alignments to transfer helix annotations
    from a reference structure to all other structures.

    Args:
        data_dict: Dictionary with data from previous steps
        output_dir: Directory to save output files
        visualize: Whether to generate visualizations

    Returns:
        Dictionary with oriented and annotated structures
    """
    print("\n=== Step 3 & 4: Annotating Helices Using Alignments ===")
    processed_structures = data_dict['processed_structures']
    alignment_paths = data_dict.get('alignment_paths', {})

    # Check if we have structures and alignment paths
    if not processed_structures:
        print("[ERROR] No structures available for annotation.")
        return {
            'processed_structures': processed_structures
        }

    if not alignment_paths:
        print("[ERROR] No alignment paths available. Make sure to run structure comparison first.")
        print("[ERROR] Falling back to old approach.")

        # Fallback to old approach (not recommended)
        # This is just a placeholder, the old code is removed as requested
        return {
            'processed_structures': processed_structures
        }

    # Ensure all structures have correct data types
    print("[INFO] Ensuring correct data types for all structures...")
    processed_structures = ensure_structure_dtypes(processed_structures)

    # Extract CA atoms for each structure (needed for helix annotation)
    for pdb_id, data in processed_structures.items():
        if 'df_norm' not in data or data['df_norm'].empty:
            if 'df' in data and not data['df'].empty:
                df_norm = data['df'].copy()
                data['df_norm'] = df_norm
            else:
                print(f"[WARNING] {pdb_id}: No structure data available.")
                continue

        df_norm = data['df_norm']
        try:
            # Try direct comparison first
            df_ca_norm = df_norm[df_norm['res_atom_name'] == 'CA'].copy()
            # If we didn't get any CA atoms, try string conversion
            if df_ca_norm.empty:
                print(f"[DEBUG] No CA atoms found with direct comparison for {pdb_id}, trying string conversion")
                df_ca_norm = df_norm[df_norm['res_atom_name'].astype(str) == 'CA'].copy()
        except TypeError:
            # If TypeError occurs, use string conversion
            print(f"[DEBUG] TypeError when extracting CA atoms for {pdb_id}, using string conversion")
            df_ca_norm = df_norm[df_norm['res_atom_name'].astype(str) == 'CA'].copy()

        if df_ca_norm.empty:
            print(f"[WARNING] No CA atoms found for {pdb_id} even after string conversion!")

        data['df_ca_norm'] = df_ca_norm

    # Choose reference structure for helix definitions
    # 1. Try to use the global reference if available
    reference_id = data_dict.get('global_ref')

    # 2. If not available, use CnChR2_J230_refine9 if present
    if not reference_id or reference_id not in processed_structures:
        if 'CnChR2_J230_refine9' in processed_structures:
            reference_id = 'CnChR2_J230_refine9'
        else:
            # 3. Otherwise use the first structure
            reference_id = next(iter(processed_structures.keys()))

    print(f"[INFO] Using {reference_id} as reference structure for helix annotation")

    # Get helix definitions for reference structure
    helix_ref_file = os.path.join(os.path.dirname(__file__), 'property', 'helix_ref_CnChR2_J230_refine9.json')
    helix_definitions = define_reference_helices(
        processed_structures[reference_id],
        helix_ref_file=helix_ref_file
    )

    if not helix_definitions:
        print("[ERROR] No helix definitions available. Cannot proceed with annotation.")
        return {
            'processed_structures': processed_structures,
            'reference_structure': reference_id
        }

    # Annotate helices in all structures based on reference
    processed_structures = annotate_helices_from_alignments(
        processed_structures,
        reference_id,
        helix_definitions,
        alignment_paths
    )

    # Optionally visualize the annotated structures
    if visualize:
        for sid, data in processed_structures.items():
            if 'tm_helices' in data:
                try:
                    fig = visualize_single_7tm_bundle(processed_structures, sid)
                    fig.write_html(os.path.join(output_dir, f'{sid}_7tm_bundle.html'))
                    print(
                        f"[INFO] Saved 7TM bundle visualization to {os.path.join(output_dir, f'{sid}_7tm_bundle.html')}")
                except Exception as e:
                    print(f"[WARNING] Failed to visualize 7TM bundle for {sid}: {e}")

    return {
        'processed_structures': processed_structures,
        'reference_structure': reference_id,
        'helix_definitions': helix_definitions
    }


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


def align_and_assign_grn(data_dict, output_dir='output', visualize=True):
    """
    Step 6: Structure alignment and GRN assignment

    This function assigns Generic Residue Numbers (GRN) using cached alignment paths
    without recalculating alignments.

    Args:
        data_dict: Dictionary with data from previous steps
        output_dir: Directory to save output files
        visualize: Whether to generate visualizations

    Returns:
        Dictionary with alignment data
    """
    print("\n=== Step 6: Structure Alignment & GRN Assignment ===")

    processed_structures_complete = data_dict['processed_structures']
    alignment_paths = data_dict['alignment_paths']

    rmsd_df = data_dict.get('rmsd_df', pd.DataFrame())
    pdb_list = data_dict.get('pdb_list', [])
    group_dict = data_dict.get('group_dict', {})

    # Check if we have structures to align
    if not processed_structures_complete or len(processed_structures_complete) < 2:
        print(f"Need at least 2 structures for alignment. Found {len(processed_structures_complete)}.")
        # Return empty data structures
        return {
            'processed_structures': processed_structures_complete,
            'seq_alignment_dicts': {},
            'msa_df': pd.DataFrame(),
            'distance_table': pd.DataFrame(),
            'ca_msa_df': pd.DataFrame(),
            'ca_distance_table': pd.DataFrame(),
            'global_ref': None,
            'type_reference_dict': {}
        }

    # Check if RMSD matrix is valid
    if rmsd_df.empty or len(pdb_list) < 2:
        print("No valid RMSD matrix available for alignment.")
        return {
            'processed_structures': processed_structures_complete,
            'seq_alignment_dicts': {},
            'msa_df': pd.DataFrame(),
            'distance_table': pd.DataFrame(),
            'ca_msa_df': pd.DataFrame(),
            'ca_distance_table': pd.DataFrame(),
            'global_ref': None,
            'type_reference_dict': {}
        }

    # Find type references and global reference
    type_reference_dict = find_type_references(rmsd_df.loc[pdb_list, pdb_list], group_dict)

    print(type_reference_dict)

    global_ref = find_global_reference(rmsd_df.loc[pdb_list, pdb_list], type_reference_dict)
    print(f"Global reference structure: {global_ref}")

    # Use cached alignment paths instead of recalculating alignments
    print("Using cached alignment paths from step 5 instead of recalculating alignments...")
    seq_alignment_dicts = create_seq_alignment_dicts_from_paths(
        alignment_paths=alignment_paths,
        structure_ids=pdb_list,
        global_ref=global_ref,
        type_reference_dict=type_reference_dict
    )

    # Generate all MSA tables with GRN labeling
    print("Generating MSA and distance tables with GRN labeling...")

    try:
        # Pass the RMSD matrix to filter structures based on RMSD to reference
        tables = generate_grn_msa_tables(
            seq_alignment_dicts,
            processed_structures_complete,
            global_ref,
            rmsd_df=rmsd_df,  # Pass RMSD matrix for filtering
            max_rmsd_threshold=3.0  # Filter structures with RMSD > 3.0 to reference
        )

        # Extract tables
        msa_df = tables["residue_table"]
        distance_table = tables["distance_table"]
        ca_msa_df = tables["ca_residue_table"]
        ca_distance_table = tables["ca_distance_table"]

        # Report on any excluded structures
        if "excluded_structures" in tables and tables["excluded_structures"]:
            excluded_count = len(tables["excluded_structures"])
            print(f"\n[INFO] {excluded_count} structures were excluded from MSA due to high RMSD (>3.0Å)")
            print(f"[INFO] Final MSA includes {len(msa_df)} structures and {len(msa_df.columns)} positions")

        # Save the tables using direct file operations
        msa_df.to_csv(os.path.join(output_dir, "msa_table_grn.csv"))
        distance_table.to_csv(os.path.join(output_dir, "distance_table_grn.csv"))
        ca_msa_df.to_csv(os.path.join(output_dir, "ca_msa_table_grn.csv"))
        ca_distance_table.to_csv(os.path.join(output_dir, "ca_distance_table_grn.csv"))

        _50_positions = []

        print("MSA and distance tables generated and saved.")

        # Display statistics about the tables
        print("\nMSA and Distance Table Statistics:")
        print(f"Number of structures: {len(msa_df)}")
        print(f"Number of aligned positions: {len(msa_df.columns)}")

        # Count TM helices in GRN positions
        tm_positions = [col for col in msa_df.columns if '.' in col and not col.startswith('L.')]
        print(f"TM residue positions: {len(tm_positions)}")

        # Count positions by helix
        for helix in range(1, 8):
            helix_positions = [col for col in msa_df.columns if col.startswith(f"{helix}.")]
            print(f"  Helix {helix}: {len(helix_positions)} positions")

        # Distance statistics - handle NaN values properly
        print(f"\nDistance statistics:")
        if not distance_table.empty and not ca_distance_table.empty:
            avg_sidechain = distance_table.mean(skipna=True).mean(skipna=True)
            avg_backbone = ca_distance_table.mean(skipna=True).mean(skipna=True)
            print(f"  Average distance to RET (sidechain): {avg_sidechain:.2f}Å")
            print(f"  Average distance to RET (backbone): {avg_backbone:.2f}Å")

            # Find closest residues to RET (handling NaN values)
            closest_residues = []
            for col in distance_table.columns:
                if '.' in col and not col.startswith('L.'):
                    avg_dist = distance_table[col].mean(skipna=True)
                    if not pd.isna(avg_dist):  # Skip columns with all NaN
                        closest_residues.append((col, avg_dist))

            closest_residues.sort(key=lambda x: x[1])
            print("\nTop 10 closest residues to RET (across all structures):")
            for pos, dist in closest_residues[:10]:
                print(f"  {pos}: {dist:.2f}Å")

        # Visualizations
        if visualize:
            try:
                # Plot distance heatmap and line plots
                if not distance_table.empty:
                    fig1 = plot_average_distances_by_helix(distance_table)
                    fig1.savefig(os.path.join(output_dir, 'distances_by_helix.png'), dpi=300, bbox_inches='tight')
                    plt.close(fig1)
                    print(f"Saved distance by helix plot to {os.path.join(output_dir, 'distances_by_helix.png')}")

                    fig2 = plot_distance_heatmap(distance_table)
                    fig2.savefig(os.path.join(output_dir, 'distance_heatmap.png'), dpi=300, bbox_inches='tight')
                    plt.close(fig2)
                    print(f"Saved distance heatmap to {os.path.join(output_dir, 'distance_heatmap.png')}")

                # Analyze residue composition at key positions
                if not msa_df.empty:
                    positions_to_analyze = [f"{helix}.50" for helix in range(1, 8)]
                    positions_in_df = [pos for pos in positions_to_analyze if pos in msa_df.columns]

                    if positions_in_df:
                        residue_composition = analyze_residue_composition(msa_df, positions_in_df)
                        print_residue_composition(residue_composition)

                    # Also analyze binding pocket positions
                    binding_pocket_positions = ["3.37", "3.40", "6.48", "6.51", "7.43", "7.46", "7.51"]
                    positions_in_df = [pos for pos in binding_pocket_positions if pos in msa_df.columns]

                    if positions_in_df:
                        binding_pocket_composition = analyze_residue_composition(msa_df, positions_in_df)
                        print_residue_composition(binding_pocket_composition)

                    # Create sequence logo using the fixed function
                    try:
                        fig = create_sequence_logo(msa_df)
                        fig.savefig(os.path.join(output_dir, "sequence_logo.png"), dpi=300, bbox_inches='tight')
                        plt.close(fig)
                        print(f"Saved sequence logo to {os.path.join(output_dir, 'sequence_logo.png')}")
                    except Exception as e:
                        print(f"[WARNING] Error creating sequence logo: {str(e)}")

                    # Create residue conservation plot
                    fig = create_residue_conservation_plot(msa_df)
                    fig.savefig(os.path.join(output_dir, "residue_conservation.png"), dpi=300, bbox_inches='tight')
                    plt.close(fig)
                    print(f"Saved residue conservation plot to {os.path.join(output_dir, 'residue_conservation.png')}")

            except Exception as e:
                print(f"Warning: Error during visualization: {e}")

        return {
            'processed_structures': processed_structures_complete,
            'seq_alignment_dicts': seq_alignment_dicts,
            'msa_df': msa_df,
            'distance_table': distance_table,
            'ca_msa_df': ca_msa_df,
            'ca_distance_table': ca_distance_table,
            'global_ref': global_ref,
            'type_reference_dict': type_reference_dict
        }

    except Exception as e:
        print(f"Error during GRN assignment: {e}")
        import traceback
        traceback.print_exc()

        return {
            'processed_structures': processed_structures_complete,
            'seq_alignment_dicts': {},
            'msa_df': pd.DataFrame(),
            'distance_table': pd.DataFrame(),
            'ca_msa_df': pd.DataFrame(),
            'ca_distance_table': pd.DataFrame(),
            'global_ref': global_ref if 'global_ref' in locals() else None,
            'type_reference_dict': type_reference_dict if 'type_reference_dict' in locals() else {}
        }


def prepare_structures_for_foldmason(processed_structures_complete, tmp_dir='tmp_foldmason'):
    """
    Export processed structures to temporary CIF files for FoldMason using direct file operations
    """
    print("[INFO] Preparing structures for FoldMason alignment...")
    os.makedirs(tmp_dir, exist_ok=True)
    structure_files = []

    # Create a temporary CifProcessor for writing files
    from protos.processing.structure.struct_base_processor import CifBaseProcessor
    data_dir = os.environ.get("PROTOS_DATA_ROOT", os.path.abspath("opsin_output"))
    cp_temp = CifBaseProcessor(name="tmp_foldmason", data_root=data_dir)

    for pdb_id, data in tqdm(processed_structures_complete.items(), desc="Creating temporary CIF files"):
        if 'df_norm' not in data or data['df_norm'].empty:
            print(f"[WARNING] {pdb_id}: No df_norm available, skipping.")
            continue

        # Create a temporary dataframe with the oriented, processed structure
        temp_df = data['df_norm'].copy()

        # Ensure all necessary columns are present for CIF output
        if 'pdb_id' not in temp_df.columns:
            temp_df['pdb_id'] = pdb_id

        # Ensure atom names are properly set if missing
        required_columns = ['group', 'type_symbol', 'auth_asym_id', 'auth_atom_id']
        for col in required_columns:
            if col not in temp_df.columns:
                if col == 'group':
                    temp_df[col] = 'ATOM'
                elif col == 'type_symbol':
                    # Derive from res_atom_name if available
                    if 'res_atom_name' in temp_df.columns:
                        temp_df[col] = temp_df['res_atom_name'].str[0]
                    else:
                        temp_df[col] = 'C'  # Default to carbon
                elif col == 'auth_asym_id':
                    if 'auth_chain_id' in temp_df.columns:
                        temp_df[col] = temp_df['auth_chain_id']
                    else:
                        temp_df[col] = 'A'
                elif col == 'auth_atom_id':
                    if 'res_atom_name' in temp_df.columns:
                        temp_df[col] = temp_df['res_atom_name']
                    else:
                        temp_df[col] = 'CA'  # Default to CA atoms

        # Write directly to a temporary CIF file
        temp_file = os.path.join(tmp_dir, f"{pdb_id}.cif")

        try:
            # Use cif_handler directly for more reliable output
            from protos.io.cif_handler import write_cif_file
            write_cif_file(temp_df, temp_file)
            structure_files.append(temp_file)
        except Exception as e:
            print(f"[ERROR] Failed to create temporary CIF for {pdb_id}: {str(e)}")

    print(f"[INFO] Created {len(structure_files)} temporary CIF files for FoldMason.")
    return structure_files


def align_with_foldmason(data_dict, output_dir='output', tmp_dir='tmp_foldmason', visualize=True):
    """
    Use FoldMason for multiple structure alignment instead of pairwise alignment
    """
    print("\n=== Step 6: FoldMason Multiple Structure Alignment & GRN Assignment ===")

    processed_structures_complete = data_dict['processed_structures']

    # Make sure output directories exist
    os.makedirs(output_dir, exist_ok=True)
    tmp_path = os.path.join(output_dir, tmp_dir)
    os.makedirs(tmp_path, exist_ok=True)

    # 1. Export structures to temporary files
    structure_files = prepare_structures_for_foldmason(processed_structures_complete, tmp_path)

    if not structure_files:
        print("[ERROR] No structure files were created. Cannot proceed with FoldMason alignment.")
        return data_dict

    try:
        # 2. Run FoldMason's easy-msa for initial alignment
        output_prefix = os.path.join(output_dir, "opsin_msa")
        print(f"[INFO] Running FoldMason easy-msa on {len(structure_files)} structures...")
        run_easy_msa(structure_files, output_prefix, tmp_path, report_mode=2)

        # 3. Generate refined structural MSA
        structure_db = f"{output_prefix}_strucDB"
        msa_prefix = os.path.join(output_dir, "opsin_structural_msa")

        print("[INFO] Running structural MSA with FoldMason...")
        run_structuremsa(structure_db, msa_prefix)

        # 4. Refine the MSA to improve quality
        input_msa = f"{msa_prefix}.fasta"
        refined_msa = os.path.join(output_dir, "opsin_refined.fasta")

        print("[INFO] Refining MSA with FoldMason...")
        run_refinemsa(structure_db, input_msa, refined_msa, refine_iters=1000)

        # 5. Generate LDDT data (avoid HTML format)
        # Instead of using msa2lddtreport, use msa2lddt to get LDDT scores directly
        from projects.opsin_analysis.foldmason_helpers import run_msa2lddt

        lddt_output = os.path.join(output_dir, "alignment_lddt_scores.txt")
        print("[INFO] Calculating LDDT scores...")
        try:
            # Run FoldMason's LDDT calculation and capture the output
            lddt_result = run_msa2lddt(structure_db, refined_msa)

            # Save the LDDT scores to a text file
            with open(lddt_output, 'w') as f:
                f.write(lddt_result)

            print("[INFO] LDDT scores saved to", lddt_output)
        except Exception as e:
            print(f"[WARNING] Error calculating LDDT scores: {str(e)}")
            lddt_output = None
            lddt_result = None

        # Process results as in original function...
        # Rest of the function implementation remains the same

        # Parse and analyze FoldMason results here...

        return {
            'processed_structures': processed_structures_complete,
            'foldmason_alignment': refined_msa,
            'foldmason_lddt_scores': lddt_output,
            'foldmason_lddt_result': lddt_result
        }

    except Exception as e:
        print(f"[ERROR] FoldMason alignment failed: {str(e)}")
        import traceback
        traceback.print_exc()
        return data_dict


def load_opsin_property_data(property_file, processed_structures):
    """
    Load property data from CSV file and create structure mappings

    Args:
        property_file: Path to CSVFailed to visualize 7TM bundle file with property data
        processed_structures: Dictionary of processed structures

    Returns:
        Dictionary with properties and structure mappings
    """
    print("\n=== Loading Property Data ===")

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

        # Process each property row
        for _, row in df_properties.iterrows():
            row_dict = dict(row)

            # Get domain/taxonomic information
            # Primary keys to check for domain (in order of preference)
            domain_keys = ['Domain', 'domain', 'Rhodopsin Type (Microbial)', 'rhodopsin_type']
            domain = 'Unknown'
            for key in domain_keys:
                if key in row_dict and pd.notna(row_dict[key]):
                    domain = str(row_dict[key]).strip()
                    if domain:  # Only break if we found a non-empty value
                        break

            # Get functional classification
            # Primary keys to check for function (in order of preference)
            function_keys = ['Function', 'function', 'molecular_function_normalized', 'molecular_function']
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

            # Create mappings for MO structures
            if pd.notna(row.get('pdb_id')) and pd.notna(row.get('short_name')):
                # Experimental MO structure ID (pdb_id)
                exp_id = str(row['pdb_id']).strip()

                # Get the short name for predicted structures
                short_name = str(row['short_name']).strip()

                # Ensure short_name is a string (to handle potential numerical values like 1E12)
                if not isinstance(short_name, str):
                    short_name = str(short_name)

                # Clean the short_name - convert dashes to underscores, remove plus signs
                clean_name = short_name.replace('-', '_').replace('+', '')

                # Handle special cases that need exact matching
                special_cases = {
                    '7BMH': 'MacR',  # Be specific about MacR (not Mac)
                    '3UG9': 'CrChR1',
                    '4PXK': 'HmBRI_D94N',
                    '1E12': '1E12',  # Ensure this is treated as a string
                }

                if exp_id in special_cases:
                    clean_name = special_cases[exp_id]
                    print(f"  Using special case mapping for {exp_id}: {clean_name}")
                else:
                    print(f"  Processing MO mapping for exp ID {exp_id}: short_name = '{short_name}' → '{clean_name}'")

                # MO predicted structures consistently use '_smile_model_0' suffix
                pred_id = f"{clean_name}_smile_model_0"

                # Add to properties dictionary - store by all relevant keys
                properties[exp_id] = property_data
                properties[pred_id] = property_data
                properties[clean_name] = property_data  # Add clean name as a key too

                # Also store clean_name in the property data itself for easier lookup
                property_data['clean_name'] = clean_name

                # Add the mapping - this is the standard pattern for all MO structures
                structure_mapping[exp_id] = pred_id

                # Also add a mapping from clean_name to both IDs for easier lookup
                structure_mapping[clean_name] = {
                    'experimental': exp_id,
                    'predicted': pred_id
                }

            # Create mappings for Hideaki structures
            # Only create Hideaki mapping if we haven't created a MO mapping for this row
            if pd.notna(row.get('short_name')) and not pd.notna(row.get('pdb_id')):
                # Experimental Hideaki structure ID (short_name)
                short_name = row['short_name'].strip()
                exp_id = short_name

                # For Hideaki structures, just convert dashes to underscores
                clean_name = short_name.replace('-', '_')
                print(f"  Processing Hideaki mapping: {short_name} → {clean_name}")

                # Predicted Hideaki structure ID (short_name + '_model_0')
                pred_id = f"{clean_name}_model_0"

                # Add to properties dictionary - store by all relevant keys
                properties[exp_id] = property_data
                properties[pred_id] = property_data
                properties[clean_name] = property_data  # Add clean name as a key too

                # Also store clean_name in the property data itself for easier lookup
                property_data['clean_name'] = clean_name

                # Add to structure mapping
                structure_mapping[exp_id] = pred_id

                # Also add a mapping from clean_name to both IDs for easier lookup
                structure_mapping[clean_name] = {
                    'experimental': exp_id,
                    'predicted': pred_id
                }

        # Associate properties with processed structures
        matched_count = 0
        for pdb_id, data in processed_structures.items():
            # Try direct lookup by PDB ID
            if pdb_id in properties:
                data['properties'] = properties[pdb_id]
                matched_count += 1
                continue

            # If not found, try extracting and looking up the clean name
            # Remove common suffixes to get a potential clean name
            potential_clean_name = pdb_id
            for suffix in ['_smile_model_0', '_model_0', '_pred']:
                if pdb_id.endswith(suffix):
                    potential_clean_name = pdb_id[:-len(suffix)]
                    break

            # Look up by potential clean name
            if potential_clean_name in properties:
                data['properties'] = properties[potential_clean_name]
                matched_count += 1

        print(f"Created property mappings for {len(properties)} structure IDs")
        print(f"Created structure mappings for {len(structure_mapping)} experimental-predicted pairs")
        print(f"Associated property data with {matched_count}/{len(processed_structures)} structures")

        # Create classification groups based on property data
        domain_groups = {}
        function_groups = {}
        for pdb_id, props in properties.items():
            # Classification by domain (taxonomic group)
            domain_type = props.get('domain', 'Unknown')
            if domain_type not in domain_groups:
                domain_groups[domain_type] = []
            domain_groups[domain_type].append(pdb_id)

            # Classification by function
            mol_function = props.get('molecular_function', 'Unknown')
            if mol_function not in function_groups:
                function_groups[mol_function] = []
            function_groups[mol_function].append(pdb_id)

        print("\nClassification by Domain:")
        for group, ids in domain_groups.items():
            print(f"  - {group}: {len(ids)} structures")

        print("\nClassification by Function:")
        for group, ids in function_groups.items():
            print(f"  - {group}: {len(ids)} structures")

        # For backward compatibility, also update processed_structures
        return {
            'processed_structures': processed_structures,
            'property_data': df_properties,
            'function_groups': function_groups,
            'domain_groups': domain_groups,
            'classification_groups': {**domain_groups, **function_groups},  # Combine both group types for compatibility
            'structure_mapping': structure_mapping,
            'properties': properties
        }

    except Exception as e:
        print(f"Error loading property data: {e}")
        import traceback
        traceback.print_exc()
        return {
            'processed_structures': processed_structures,
            'property_data': None,
            'classification_groups': {},
            'structure_mapping': {},
            'properties': {}
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


def align_to_reference_and_annotate_helices(data_dict, output_dir='output', visualize=True):
    """
    Step 4: Custom step that aligns all structures to a reference structure
    (from helix_ref_CnChR2_J230_refine9.json) and annotates all structures
    with helix numbers based on the reference helices.

    If helices.json already exists in property directory, it will load helix
    definitions from there instead of recalculating alignments.

    Args:
        data_dict: Dictionary with data from previous steps
        output_dir: Directory to save output files
        visualize: Whether to generate visualizations

    Returns:
        Dictionary with aligned structures and helix annotations
    """
    print("\n=== Step 4: Aligning to Reference and Annotating Helices ===")

    import json
    import os
    import numpy as np
    import pandas as pd
    from pathlib import Path
    from protos.processing.structure.struct_alignment import get_structure_alignment

    processed_structures = data_dict['processed_structures']

    # Check if we have structures to align
    if not processed_structures:
        print("[ERROR] No structures available for alignment.")
        return {
            'processed_structures': processed_structures
        }

    # Define paths for helix files
    property_dir = os.path.join(os.path.dirname(__file__), 'property')
    helix_ref_file = os.path.join(property_dir, 'helix_ref_CnChR2_J230_refine9.json')
    helix_cache_file = os.path.join(property_dir, 'helices.json')

    # Initialize variables
    global_helix_annotations = {}
    alignment_paths = {}
    formatted_helix_defs = {}
    ref_id = None

    # First, check if we already have helices.json with annotations for existing structures
    if os.path.exists(helix_cache_file):
        try:
            print(f"[INFO] Found existing helix definitions at {helix_cache_file}")

            # Check file size to detect potential truncation
            file_size = os.path.getsize(helix_cache_file)
            if file_size < 100:  # Suspiciously small for a JSON with structure definitions
                print(f"[WARNING] Helix cache file is suspiciously small ({file_size} bytes)")
                print(f"[WARNING] This may indicate a truncated or corrupted file")
                raise ValueError("Suspicious file size indicates potential corruption")

            # Try to read the first few lines to check format
            with open(helix_cache_file, 'r') as f:
                # Read first 1000 characters to check structure
                file_start = f.read(1000)
                if not file_start.strip().startswith('{') or '}' not in file_start:
                    print(f"[WARNING] Helix cache file does not appear to be valid JSON")
                    raise ValueError("Invalid JSON format detected")
                # Reset file pointer
                f.seek(0)

                # Try to load the full JSON
                try:
                    global_helix_annotations = json.load(f)
                except json.JSONDecodeError as json_err:
                    print(f"[WARNING] JSON parsing error: {json_err}")

                    # Try to read the whole file to get more info
                    f.seek(0)
                    file_content = f.read()
                    last_complete_brace = file_content.rfind('}')

                    if last_complete_brace > 0:
                        print("[INFO] Attempting to recover partial JSON data...")
                        try:
                            # Try to repair by closing the JSON properly
                            repaired_json = file_content[:last_complete_brace + 1]
                            global_helix_annotations = json.loads(repaired_json)
                            print(
                                f"[INFO] Successfully recovered {len(global_helix_annotations)} structure definitions")
                        except:
                            print("[WARNING] Failed to recover JSON data")
                            raise ValueError("Could not repair corrupted JSON file")
                    else:
                        raise ValueError("JSON file appears severely corrupted")

            # Extract reference structure ID - assume it's the first one
            existing_structures = list(global_helix_annotations.keys())
            if existing_structures:
                # Try to find the reference structure (CnChR2_J230_refine9 or similar)
                ref_candidates = [s for s in existing_structures if 'CnChR2_J230_refine9' in s]
                if ref_candidates:
                    ref_id = ref_candidates[0]
                else:
                    # Just use the first one
                    ref_id = existing_structures[0]

                print(f"[INFO] Using {ref_id} as reference structure from cached helix definitions")

                # Format helix definitions
                helix_definitions = global_helix_annotations[ref_id]
                for helix_num, bounds in helix_definitions.items():
                    if isinstance(bounds, list) and len(bounds) == 2:
                        formatted_helix_defs[helix_num] = {'start': bounds[0], 'end': bounds[1]}

                # Validate that we got all 7 helices for the reference structure
                if len(formatted_helix_defs) < 7:
                    print(
                        f"[WARNING] Reference structure only has {len(formatted_helix_defs)} helices defined, expected 7")
                    print(f"[WARNING] This may indicate incomplete or corrupted data")

                    # If we have fewer than 4 helices, consider the data unreliable
                    if len(formatted_helix_defs) < 4:
                        print(f"[WARNING] Too few helices defined. Will recalculate all helix definitions")
                        raise ValueError("Insufficient helix definitions in reference structure")
            else:
                print(f"[WARNING] No structures found in helix cache file")
                raise ValueError("Empty structure list in helix cache")

            # Check if we have definitions for all current structures
            current_struct_ids = set(processed_structures.keys())
            cached_struct_ids = set(global_helix_annotations.keys())
            missing_structs = current_struct_ids - cached_struct_ids

            # Check if any cached structures have incomplete definitions
            incomplete_structs = []
            for sid, struct_def in global_helix_annotations.items():
                if sid not in processed_structures:
                    continue  # Skip structures we don't need

                # Check if this structure has a reasonable number of helices
                if len(struct_def) < 4:  # Expect at least 4 of 7 helices
                    incomplete_structs.append(sid)

            if incomplete_structs:
                print(f"[WARNING] Found {len(incomplete_structs)} structures with incomplete helix definitions")
                print(f"[WARNING] These will be recalculated: {', '.join(incomplete_structs[:5])}")
                if len(incomplete_structs) > 5:
                    print(f"[WARNING] ... and {len(incomplete_structs) - 5} more")
                missing_structs.update(incomplete_structs)

            if missing_structs:
                print(f"[INFO] Found {len(missing_structs)} structures that need helix definitions")
                print(f"[INFO] Will calculate helix definitions for: {', '.join(list(missing_structs)[:5])}")
                if len(missing_structs) > 5:
                    print(f"[INFO] ... and {len(missing_structs) - 5} more")

                # We'll need to load the reference helix definitions and calculate alignments
                # for the missing structures
                need_alignment = True
            else:
                print(f"[INFO] Found helix definitions for all {len(current_struct_ids)} structures")
                need_alignment = False

        except Exception as e:
            print(f"[WARNING] Error loading cached helix definitions: {e}")
            print(f"[INFO] Will recalculate helix definitions for all structures")

            # Make a backup of the problematic file
            if os.path.exists(helix_cache_file):
                backup_file = f"{helix_cache_file}.bak"
                try:
                    import shutil
                    shutil.copy2(helix_cache_file, backup_file)
                    print(f"[INFO] Created backup of problematic helix cache at {backup_file}")
                except Exception as backup_err:
                    print(f"[WARNING] Failed to create backup: {backup_err}")

            # Reset to empty dictionary and recalculate
            global_helix_annotations = {}
            need_alignment = True
    else:
        print(f"[INFO] No cached helix definitions found at {helix_cache_file}")
        print(f"[INFO] Will calculate helix definitions for all structures")
        need_alignment = True

    # If we need to calculate alignments for all or some structures
    if need_alignment:
        # Load reference helix definitions if not already loaded
        if not ref_id or not formatted_helix_defs:
            if not os.path.exists(helix_ref_file):
                print(f"[ERROR] Reference helix file not found: {helix_ref_file}")
                return {
                    'processed_structures': processed_structures
                }

            try:
                with open(helix_ref_file, 'r') as f:
                    helix_data = json.load(f)

                # Get the reference structure ID
                ref_id = next(iter(helix_data.keys()))
                helix_definitions = helix_data[ref_id]

                print(f"[INFO] Loaded helix definitions for reference structure {ref_id}")

                # Format helix definitions as dictionary of {helix_num: {'start': pos, 'end': pos}}
                formatted_helix_defs = {}
                for helix_num, bounds in helix_definitions.items():
                    if isinstance(bounds, list) and len(bounds) == 2:
                        formatted_helix_defs[helix_num] = {'start': bounds[0], 'end': bounds[1]}

                if not formatted_helix_defs:
                    print("[ERROR] Invalid helix definition format.")
                    return {
                        'processed_structures': processed_structures
                    }

                # Add reference structure helices to global annotations
                global_helix_annotations[ref_id] = helix_definitions

            except Exception as e:
                print(f"[ERROR] Failed to load helix reference file: {e}")
                return {
                    'processed_structures': processed_structures
                }

        # Check if reference structure exists in our dataset
        if ref_id not in processed_structures:
            print(f"[WARNING] Reference structure {ref_id} not found in processed structures.")
            print(f"[WARNING] Will try to find a structure with similar name.")

            # Try to find a structure with a similar name
            potential_ref_ids = [sid for sid in processed_structures.keys() if ref_id in sid]
            if potential_ref_ids:
                # Use the first match
                ref_id = potential_ref_ids[0]
                print(f"[INFO] Using {ref_id} as reference structure.")
            else:
                print(f"[ERROR] Cannot find suitable reference structure. Aborting helix annotation.")
                return {
                    'processed_structures': processed_structures
                }

        # Ensure all structures have correct data types
        print("[INFO] Ensuring correct data types for all structures...")
        for pdb_id, data in processed_structures.items():
            if 'df' in data:
                df = data['df']
                for col in ['x', 'y', 'z']:
                    if col in df.columns:
                        try:
                            df[col] = pd.to_numeric(df[col], errors='coerce')
                        except Exception as e:
                            print(f"[WARNING] Error converting {col} in {pdb_id}: {e}")

                # Create df_norm if it doesn't exist
                if 'df_norm' not in data or data['df_norm'].empty:
                    data['df_norm'] = df.copy()

                # Extract CA atoms for alignment
                try:
                    df_ca = df[df['res_atom_name'] == 'CA'].copy()
                    if df_ca.empty:
                        # Try string conversion
                        df_ca = df[df['res_atom_name'].astype(str) == 'CA'].copy()
                    data['df_ca'] = df_ca
                except Exception as e:
                    print(f"[WARNING] Error extracting CA atoms in {pdb_id}: {e}")

        # Get reference structure data
        ref_structure = processed_structures[ref_id]

        # Extract reference CA atoms
        if 'df_ca' not in ref_structure or ref_structure['df_ca'].empty:
            print(f"[ERROR] No CA atoms found in reference structure {ref_id}")
            return {
                'processed_structures': processed_structures
            }

        ref_ca_df = ref_structure['df_ca']
        ref_ca_coords = ref_ca_df[['x', 'y', 'z']].astype(float).values
        ref_seq_ids = ref_ca_df['auth_seq_id'].values if 'auth_seq_id' in ref_ca_df.columns else list(
            range(len(ref_ca_coords)))

        # Determine which structures need alignment
        structures_to_align = []
        for struct_id in processed_structures.keys():
            if struct_id == ref_id:
                continue  # Skip reference structure

            # Check if structure already has helix definitions in the cache
            if struct_id in global_helix_annotations:
                # Verify the definitions are complete (all 7 helices)
                helix_defs = global_helix_annotations[struct_id]
                if len(helix_defs) < 7:
                    print(
                        f"[INFO] Structure {struct_id} has incomplete helix definitions ({len(helix_defs)}/7 helices)")
                    structures_to_align.append(struct_id)
                else:
                    print(f"[INFO] Using cached helix definitions for {struct_id}")
            else:
                # Structure not in cache, needs alignment
                structures_to_align.append(struct_id)

        # Align structures that don't have helix definitions yet
        if structures_to_align:
            print(f"[INFO] Aligning {len(structures_to_align)} structures to reference {ref_id} and mapping helices...")

            for struct_id in structures_to_align:
                structure = processed_structures[struct_id]

                # Skip structures without CA atoms
                if 'df_ca' not in structure or structure['df_ca'].empty:
                    print(f"[WARNING] Skipping {struct_id} - No CA atoms found.")
                    continue

                # Get target structure CA atoms
                target_ca_df = structure['df_ca']
                target_ca_coords = target_ca_df[['x', 'y', 'z']].astype(float).values
                target_seq_ids = target_ca_df['auth_seq_id'].values if 'auth_seq_id' in target_ca_df.columns else list(
                    range(len(target_ca_coords)))

                try:
                    # Perform structure alignment
                    rotation, translation, best_path, rmsd = get_structure_alignment(ref_ca_coords, target_ca_coords)

                    # Extract indices from alignment path
                    ref_indices, target_indices = best_path

                    # Map indices to auth_seq_id values
                    ref_res_ids = [ref_seq_ids[idx] for idx in ref_indices]
                    target_res_ids = [target_seq_ids[idx] for idx in target_indices]

                    # Create mapping between reference and target residue IDs
                    residue_mapping = list(zip(ref_res_ids, target_res_ids))

                    # Store alignment path
                    alignment_paths[(ref_id, struct_id)] = {
                        'rotation': rotation.tolist(),
                        'translation': translation.tolist(),
                        'coord_indices': best_path,
                        'residue_mapping': residue_mapping,
                        'rmsd': rmsd
                    }

                    # Now map helix definitions from reference to target
                    target_helix_defs = {}
                    target_residue_helices = {}

                    # For each helix in reference, find corresponding residues in target
                    for helix_num, bounds in formatted_helix_defs.items():
                        ref_start = bounds['start']
                        ref_end = bounds['end']

                        # Find all aligned residues in this helix
                        helix_mappings = []
                        for ref_res, target_res in residue_mapping:
                            if ref_start <= ref_res <= ref_end:
                                helix_mappings.append((ref_res, target_res))
                                target_residue_helices[target_res] = int(helix_num)

                        # If we found mappings for this helix, create target helix definition
                        if helix_mappings:
                            target_res_in_helix = [t for _, t in helix_mappings]
                            target_helix_defs[helix_num] = {
                                'start': min(target_res_in_helix),
                                'end': max(target_res_in_helix)
                            }

                    # Add to global annotations
                    global_helix_annotations[struct_id] = {
                        str(h_num): [bounds['start'], bounds['end']]
                        for h_num, bounds in target_helix_defs.items()
                    }

                    print(f"[INFO] Successfully aligned and annotated {struct_id} (RMSD: {rmsd:.3f}Å)")

                except Exception as e:
                    print(f"[ERROR] Failed to align and annotate {struct_id}: {e}")

        # Save global helix annotations to property directory
        try:
            # Make sure all helix definitions are serializable
            clean_helix_annotations = {}
            for struct_id, helix_defs in global_helix_annotations.items():
                # Create a clean copy with only strings and simple types
                clean_struct_defs = {}
                for helix_id, bounds in helix_defs.items():
                    # Convert bounds to simple list with integers if needed
                    if isinstance(bounds, list) and len(bounds) == 2:
                        # Ensure bounds are integers
                        clean_bounds = [int(bounds[0]), int(bounds[1])]
                        clean_struct_defs[str(helix_id)] = clean_bounds
                # Only add structures with at least one helix definition
                if clean_struct_defs:
                    clean_helix_annotations[struct_id] = clean_struct_defs

            # First, write to a temporary file to avoid corruption
            temp_file = f"{helix_cache_file}.tmp"
            with open(temp_file, 'w') as f:
                json.dump(clean_helix_annotations, f, indent=2)

            # Now rename the temporary file to the actual file (atomic operation)
            import os
            if os.path.exists(helix_cache_file):
                # Create a backup of the current file first
                backup_file = f"{helix_cache_file}.bak"
                try:
                    import shutil
                    shutil.copy2(helix_cache_file, backup_file)
                except Exception as backup_err:
                    print(f"[WARNING] Failed to create backup: {backup_err}")

            # Rename temp file to actual file name
            import shutil
            shutil.move(temp_file, helix_cache_file)

            print(f"[INFO] Saved helix annotations for {len(clean_helix_annotations)} structures to {helix_cache_file}")
        except Exception as e:
            print(f"[WARNING] Failed to save helix annotations: {e}")
            import traceback
            traceback.print_exc()

    # Now apply helix annotations to all structures
    print("[INFO] Applying helix annotations to all structures...")

    for struct_id, structure in processed_structures.items():
        # Skip structures without helix definitions
        if struct_id not in global_helix_annotations:
            print(f"[WARNING] No helix definition found for {struct_id}")

            # Instead of skipping, initialize with empty helix data
            # This ensures all structures have the expected fields
            structure['helix_definitions'] = {}
            structure['residue_to_helix'] = {}
            structure['tm_helices'] = []

            # Initialize dataframes with helix_num = 0
            for df_key in ['df', 'df_norm', 'df_ca', 'df_ca_norm']:
                if df_key in structure and not structure[df_key].empty:
                    if 'helix_num' not in structure[df_key].columns:
                        structure[df_key]['helix_num'] = 0

            continue

        # Get helix definitions for this structure
        struct_helix_defs = global_helix_annotations[struct_id]

        # Convert to internal format
        formatted_struct_helix_defs = {}
        target_residue_helices = {}

        for helix_num, bounds in struct_helix_defs.items():
            if isinstance(bounds, list) and len(bounds) == 2:
                formatted_struct_helix_defs[helix_num] = {'start': bounds[0], 'end': bounds[1]}

                # Create mapping from residue ID to helix number
                start_pos = bounds[0]
                end_pos = bounds[1]
                for res_id in range(start_pos, end_pos + 1):
                    target_residue_helices[res_id] = int(helix_num)

        # Store helix definitions and mapping in structure data
        structure['helix_definitions'] = formatted_struct_helix_defs
        structure['residue_to_helix'] = target_residue_helices

        # Also store as list for compatibility with visualization functions
        structure['tm_helices'] = list(range(1, 8))

        # Process every dataframe in the structure to ensure consistent helix annotations
        for df_key in ['df', 'df_norm', 'df_ca', 'df_ca_norm']:
            if df_key in structure and not structure[df_key].empty:
                df = structure[df_key]

                # Initialize helix_num column to 0 (no helix)
                if 'helix_num' in df.columns:
                    # Reset existing column
                    df['helix_num'] = 0
                else:
                    # Create new column
                    df['helix_num'] = 0

                # Update helix_num based on mapped residues
                if 'auth_seq_id' in df.columns:
                    # Use vectorized operations when possible for efficiency
                    if len(target_residue_helices) > 0:
                        # Create a mapping series for faster lookup
                        import pandas as pd
                        import numpy as np

                        # Extract all residue IDs and convert to numeric if needed
                        residue_ids = df['auth_seq_id'].unique()

                        # Update each residue's helix number
                        for res_id, helix_num in target_residue_helices.items():
                            # Apply the helix number to all atoms in this residue
                            mask = df['auth_seq_id'] == res_id
                            df.loc[mask, 'helix_num'] = int(helix_num)

                # Store updated dataframe back in the structure
                structure[df_key] = df

        # Also update df_ret if it exists to maintain consistency
        if 'df_ret' in structure and not structure['df_ret'].empty:
            # Retinal is not part of helices, but should have helix_num column for consistency
            if 'helix_num' not in structure['df_ret'].columns:
                structure['df_ret']['helix_num'] = 0

        # Debugging print statement showing helix details for this structure
        unique_helix_numbers = sorted(list(set([int(h) for h in formatted_struct_helix_defs.keys()])))
        helix_info = f"{struct_id}: {len(formatted_struct_helix_defs)} helices => {unique_helix_numbers}"

        # Add details for each helix
        for h_num in range(1, 8):
            h_str = str(h_num)
            if h_str in formatted_struct_helix_defs:
                start = formatted_struct_helix_defs[h_str]['start']
                end = formatted_struct_helix_defs[h_str]['end']
                helix_info += f" | Helix {h_num}: {start}-{end}"
            else:
                helix_info += f" | Helix {h_num}: missing"

        print(f"[DEBUG] {helix_info}")

    return {
        'processed_structures': processed_structures,
        'reference_structure': ref_id,
        'helix_definitions': formatted_helix_defs,
        'alignment_paths': alignment_paths,
        'helix_annotations_file': helix_cache_file
    }


def run_opsin_analysis_workflow(output_dir='output', visualize=True, use_foldmason=True, use_cache=True, cache_raw=True,
                                property_file=None, chain_id='A', retinal_name='RET', retinal_cutoff=6.0):
    """
    Run the full opsin analysis workflow with standardized dataset handling
    without using PropertyProcessor

    Args:
        output_dir: Directory to save output files
        visualize: Whether to generate visualizations
        use_foldmason: Whether to use FoldMason for structure alignment
        use_cache: Whether to use cached structure data (default: True)
        cache_raw: Whether to cache raw unfiltered data (default: True)
        property_file: Path to CSV file with property data (optional)
        chain_id: Chain ID to use for analysis (default: 'A')
        retinal_name: Name of retinal residue (default: 'RET')
        retinal_cutoff: Distance cutoff in Angstroms for retinal selection (default: 6.0)

    Returns:
        Dictionary with analysis results
    """
    # Step 1: Load and process structures
    print("\n" + "=" * 80)
    print("RUNNING OPSIN ANALYSIS WORKFLOW")
    print("=" * 80 + "\n")

    # Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # Step 1: Load opsin structures from datasets (with two-stage caching)
    data = load_opsin_structures(output_dir, chain_id=chain_id, visualize=visualize,
                                 use_cache=use_cache, cache_raw=cache_raw,
                                 retinal_name=retinal_name, retinal_cutoff=retinal_cutoff)

    # Step 2: Load property data if provided
    property_data = None
    if property_file:
        default_property_path = Path(__file__).resolve().parent / "property" / "mo_exp.csv"
        property_path = property_file if os.path.exists(property_file) else default_property_path

        if os.path.exists(property_path):
            print(f"Using property data from: {property_path}")
            property_data = load_opsin_property_data(property_path, data['processed_structures'])
            data.update(property_data)
        else:
            print(f"Property file not found: {property_path}")

    # Step 3: Create a single unified structure mapping that will be used globally
    # This is the ONLY place where structure mappings are created in the workflow
    structure_mapping = create_unified_structure_mapping(data, property_data)

    # Replace any existing mapping with our unified mapping
    data['structure_mapping'] = structure_mapping

    print(f"\nUnified structure mapping created with {len(structure_mapping)} pairs")
    print(f"This mapping will be used throughout the workflow for all structure comparisons")

    # Step 4: Calculate errors between experimental and predicted structures
    # Pass the unified structure mapping to all subsequent functions
    data.update(calculate_structure_errors(data, output_dir=output_dir, visualize=visualize))

    # Step 4a: Align all structures to reference and annotate helices
    # This is our new custom step that replaces the old orient_and_annotate_structures
    data.update(align_to_reference_and_annotate_helices(data, output_dir, visualize=visualize))

    # Step 5: Structure comparison
    data.update(compare_structures(data, output_dir, visualize=visualize))

    # Step 6: Structure alignment and GRN assignment
    data.update(align_and_assign_grn(data, output_dir, visualize=visualize))

    """
    # Use FoldMason if requested (only if it exists in the environment)
    if use_foldmason:
        try:
            from projects.opsin_analysis.foldmason_helpers import run_easy_msa
            print("\n[INFO] Running additional FoldMason multiple structure alignment...")
            foldmason_data = align_with_foldmason(data, output_dir, visualize=visualize)

            # Only update if FoldMason provided results
            if foldmason_data and len(foldmason_data) > 1:
                data.update(foldmason_data)
                print("[INFO] FoldMason alignment successfully added to results")
            else:
                print("[WARNING] FoldMason alignment didn't produce usable results")
        except ImportError:
            print("[WARNING] FoldMason integration not available, skipping FoldMason alignment")
            print("[INFO] To use FoldMason, make sure the foldmason_helpers.py module is available")
    """

    # Save a summary of the results
    try:
        summary = {
            "datasets": list(data.get('datasets', {}).keys()),
            "structures_count": len(data.get('processed_structures', {})),
            "exp_pred_pairs": len(data.get('structure_mapping', {})),
            "timestamp": pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')
        }

        with open(os.path.join(output_dir, 'analysis_summary.json'), 'w') as f:
            json.dump(summary, f, indent=2)
    except Exception as e:
        print(f"Warning: Could not save summary: {e}")

    print("\n" + "=" * 80)
    print(f"ANALYSIS COMPLETE. Results saved to {output_dir}")
    print("=" * 80)

    return data


if __name__ == '__main__':
    import argparse

    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Run opsin analysis workflow')
    parser.add_argument('--output-dir', type=str, default='opsin_analysis_results',
                        help='Directory to save output files')
    parser.add_argument('--no-visualize', action='store_false', dest='visualize',
                        help='Disable visualization generation')
    parser.add_argument('--no-foldmason', action='store_false', dest='use_foldmason',
                        help='Disable FoldMason alignment')
    parser.add_argument('--no-cache', action='store_false', dest='use_cache',
                        help='Disable structure caching')
    parser.add_argument('--no-raw-cache', action='store_false', dest='cache_raw',
                        help='Disable raw data caching')
    parser.add_argument('--property-file', type=str, default=None,
                        help='Path to CSV file with property data')
    parser.set_defaults(visualize=True, use_foldmason=True, use_cache=True, cache_raw=True)

    # Parse arguments
    args = parser.parse_args()

    # Configuration parameters
    output_directory = args.output_dir
    enable_visualization = args.visualize
    use_foldmason = args.use_foldmason
    use_cache = args.use_cache
    cache_raw = args.cache_raw
    property_file = args.property_file

    # Run the workflow
    results = run_opsin_analysis_workflow(
        output_dir=output_directory,
        visualize=enable_visualization,
        use_foldmason=use_foldmason,
        use_cache=use_cache,
        cache_raw=cache_raw,
        property_file=property_file
    )


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
    from projects.opsin_analysis.msa_grn import assign_helix_numbers_to_msa_tables, count_residues_by_helix

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
