#!/usr/bin/env python
"""
Comprehensive script to load precomputed opsin analysis data and generate a
curated set of visualizations.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import argparse
from pathlib import Path
from scipy.spatial.distance import squareform  # Used in create_and_visualize_similarity_tree
from scipy.cluster.hierarchy import linkage, dendrogram  # Used in create_and_visualize_similarity_tree
import json  # For summary saving example

# --- Import visualization functions ---
try:
    from visualization_functions import (
        create_opsin_overview_plot, create_rmsd_color_scale_figure,
        create_and_visualize_similarity_tree, visualize_rmsd_matrix_improved,
        plot_similarity_tree, plot_rmsd_heatmap,
        plot_distances_with_std, plot_average_distances_by_helix, plot_distance_heatmap,
        create_residue_conservation_plot, plot_helix_logo_plots, plot_conservation_around_x50,
        create_combined_distance_logo_plot, visualize_binding_pocket
    )
except ImportError as e1:
    print(f"Error importing from visualization_functions: {e1}")
    print("Attempting to import from current directory as a fallback for development.")
    try:
        # This is a common structure if plot.py and visualization_functions.py are siblings
        from visualization_functions import (  # type: ignore
            create_opsin_overview_plot, create_rmsd_color_scale_figure,
            create_and_visualize_similarity_tree, visualize_rmsd_matrix_improved,
            plot_similarity_tree, plot_rmsd_heatmap,
            plot_distances_with_std, plot_average_distances_by_helix, plot_distance_heatmap,
            create_residue_conservation_plot, plot_helix_logo_plots, plot_conservation_around_x50,
            create_combined_distance_logo_plot, visualize_binding_pocket
        )
    except ImportError as e2:
        print(f"Fallback import failed: {e2}. Ensure visualization_functions.py is accessible.")


        # Define dummy functions if import fails completely, to allow script to try running
        def _dummy_plot_func(*args, **kwargs):
            print(
                f"Warning: Plot function called but not loaded due to import error: {kwargs.get('title', 'Unknown Plot')}")
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, "Plotting function unavailable", ha='center', va='center')
            return fig


        create_opsin_overview_plot = create_rmsd_color_scale_figure = _dummy_plot_func
        create_and_visualize_similarity_tree = visualize_rmsd_matrix_improved = _dummy_plot_func
        plot_similarity_tree = plot_rmsd_heatmap = _dummy_plot_func
        plot_distances_with_std = plot_average_distances_by_helix = plot_distance_heatmap = _dummy_plot_func
        create_residue_conservation_plot = plot_helix_logo_plots = plot_conservation_around_x50 = _dummy_plot_func
        create_combined_distance_logo_plot = visualize_binding_pocket = _dummy_plot_func

# Import the color scheme tools
try:
    from opsin_color_scheme import get_group_colors
except ImportError:
    print("Error: opsin_color_scheme.py not found.")


    def get_group_colors(items_list_or_dict, palette_name=None):  # Adjusted dummy
        if isinstance(items_list_or_dict, dict):
            items = list(items_list_or_dict.keys())
        else:
            items = list(items_list_or_dict)
        return {item: "#CCCCCC" for item in items}


def load_data(input_dir, output_dir_ref, chain_id='A'):
    """
    Load precomputed data from workflow cache or CSV files.
    Args:
        input_dir: Directory containing the data files (often the root of the project output)
        output_dir_ref: Directory where CSVs might be found if not in cache (typically input_dir or a subfolder like 'opsin_output')
        chain_id: Chain ID to use (default: 'A')
    Returns:
        Dictionary with loaded data
    """
    data = {}
    cache_dir = Path(input_dir) / 'cache'

    print(f"DEBUG: Attempting to load data. Cache directory: {cache_dir}")
    print(f"DEBUG: output_dir_ref for CSVs: {Path(output_dir_ref).resolve()}")

    processed_structures_loaded_from_cache = False
    if cache_dir.exists():
        print(f"DEBUG: Found cache directory: {cache_dir}")

        # Prioritize processed_structures for debugging properties
        proc_struct_fname = f'processed_structures_{chain_id}.pkl'
        cache_path_proc_struct = cache_dir / proc_struct_fname
        if cache_path_proc_struct.exists():
            try:
                print(f"DEBUG: Loading from cache: {proc_struct_fname}...")
                with open(cache_path_proc_struct, 'rb') as f:
                    file_data = pickle.load(f)
                    # Ensure 'processed_structures' key is correctly handled
                    if isinstance(file_data, dict) and 'processed_structures' in file_data:
                        data['processed_structures'] = file_data['processed_structures']
                    elif isinstance(file_data, dict):  # If the pkl itself is the dict of structures
                        data['processed_structures'] = file_data
                    else:
                        print(f"  DEBUG: Unexpected data type in {proc_struct_fname}")

                    processed_structures_loaded_from_cache = True

                print(
                    f"DEBUG: Successfully loaded {proc_struct_fname}. It contains {len(data.get('processed_structures', {}))} structures.")
                if 'processed_structures' in data:
                    count = 0
                    for sid_sample, s_data_sample in data['processed_structures'].items():
                        if count < 3:
                            print(f"  DEBUG_SAMPLE_CACHE {sid_sample}: type={type(s_data_sample)}")
                            if isinstance(s_data_sample, dict):
                                print(f"    DEBUG_SAMPLE_CACHE {sid_sample} keys: {list(s_data_sample.keys())}")
                                if 'properties' in s_data_sample:
                                    props_sample = s_data_sample['properties']
                                    print(f"    DEBUG_SAMPLE_CACHE {sid_sample} properties type: {type(props_sample)}")
                                    if isinstance(props_sample, dict):
                                        print(
                                            f"      DEBUG_SAMPLE_CACHE {sid_sample} properties: molecular_function='{props_sample.get('molecular_function', 'MISSING_KEY_MF')}', domain='{props_sample.get('domain', 'MISSING_KEY_DOM')}'")
                                    else:
                                        print(
                                            f"      DEBUG_SAMPLE_CACHE {sid_sample} 'properties' is not a dict. Value: {props_sample}")
                                else:
                                    print(f"    DEBUG_SAMPLE_CACHE {sid_sample} NO 'properties' key found.")
                            else:
                                print(
                                    f"    DEBUG_SAMPLE_CACHE {sid_sample} s_data_sample is not a dict. Value: {s_data_sample}")
                        count += 1
                    if len(data['processed_structures']) > 3:
                        print(
                            f"  DEBUG_SAMPLE_CACHE ... and {len(data['processed_structures']) - 3} more structures in cache.")
            except Exception as e:
                print(f"  DEBUG: Error loading {proc_struct_fname}: {e}")
        else:
            print(f"DEBUG: Cache file {proc_struct_fname} not found at {cache_path_proc_struct}")

        # Load other cache files
        all_other_cache_files = [
            'raw_structures.pkl',
            f'helix_annotations_{chain_id}.pkl',
            f'structure_comparison_{chain_id}.pkl',
            f'structure_errors_{chain_id}.pkl',
            f'grn_assignment_{chain_id}.pkl'
        ]
        for file_name in all_other_cache_files:
            cache_path = cache_dir / file_name
            if cache_path.exists():
                try:
                    print(f"DEBUG: Loading from cache: {file_name}...")
                    with open(cache_path, 'rb') as f:
                        file_data_other = pickle.load(f)  # Use a different variable name
                        if file_name == f'structure_comparison_{chain_id}.pkl':
                            print(f"  DEBUG: Keys in structure_comparison: {list(file_data_other.keys())}")
                            if 'rmsd_df' in file_data_other and isinstance(file_data_other['rmsd_df'], pd.DataFrame):
                                data['rmsd_df'] = file_data_other['rmsd_df']
                                print(
                                    f"    DEBUG: Loaded 'rmsd_df' from structure_comparison. Shape: {data['rmsd_df'].shape}")
                            elif 'rmsd_matrix' in file_data_other and isinstance(file_data_other['rmsd_matrix'],
                                                                                 pd.DataFrame):
                                data['rmsd_df'] = file_data_other[
                                    'rmsd_matrix']  # Fallback if 'rmsd_df' key was 'rmsd_matrix'
                                print(
                                    f"    DEBUG: Loaded 'rmsd_df' (from 'rmsd_matrix' key) from structure_comparison. Shape: {data['rmsd_df'].shape}")

                            # Load other potential keys from this file too, avoiding overwrite if 'rmsd_df' was main dict key
                            for key, value in file_data_other.items():
                                if key not in ['rmsd_df',
                                               'rmsd_matrix'] or key not in data:  # don't overwrite if already set
                                    data[key] = value
                            continue  # Skip generic update for this file if rmsd_df handled

                        # For grn_assignment, expect specific keys
                        if file_name == f'grn_assignment_{chain_id}.pkl' and isinstance(file_data_other, dict):
                            expected_grn_keys = ["residue_table", "distance_table", "ca_residue_table",
                                                 "ca_distance_table"]
                            for grn_key in expected_grn_keys:
                                if grn_key in file_data_other:
                                    data[grn_key] = file_data_other[grn_key]
                                    print(
                                        f"    DEBUG: Loaded '{grn_key}' from grn_assignment. Shape: {data[grn_key].shape if isinstance(data[grn_key], pd.DataFrame) else type(data[grn_key])}")
                            # Load other keys from grn_assignment if they don't clash
                            for k, v in file_data_other.items():
                                if k not in expected_grn_keys and k not in data:
                                    data[k] = v
                            continue

                        # Generic update for other files
                        if isinstance(file_data_other, dict):
                            data.update(file_data_other)
                        else:
                            # Store non-dict data under a key derived from filename if not clashing
                            data_key_name = file_name.replace(f'_{chain_id}.pkl', '').replace('.pkl', '')
                            if data_key_name not in data:
                                data[data_key_name] = file_data_other

                    print(f"  DEBUG: Successfully loaded {file_name}")
                except Exception as e:
                    print(f"  DEBUG: Error loading {file_name}: {e}")
            else:
                print(f"  DEBUG: Warning: Cache file not found: {cache_path}")
    else:
        print(f"DEBUG: Cache directory not found: {cache_dir}. Will rely on CSVs.")

    # 2. Load CSV files
    csv_files_to_load = [
        ('rmsd_matrix.csv', 'rmsd_df'),
        ('molecular_functions.csv', 'molecular_functions_df'),
        ('ca_distance_table_grn.csv', 'ca_distance_table'),
        ('distance_table_grn.csv', 'distance_table'),
        ('ca_msa_table_grn.csv', 'msa_table'),  # Often preferred for CA-based analysis
        ('residue_table_grn.csv', 'residue_table'),  # Alternative or full atom
        ('mo_exp_errors.csv', 'mo_exp_errors_df'),
        ('hideaki_errors.csv', 'hideaki_errors_df')
    ]

    molecular_functions_df_loaded_from_csv = False
    for file_name, key_name in csv_files_to_load:
        if key_name in data and data[key_name] is not None and not (
                isinstance(data[key_name], pd.DataFrame) and data[key_name].empty):
            print(f"  DEBUG: Data for '{key_name}' already loaded (likely from cache). Skipping CSV {file_name}.")
            if key_name == 'molecular_functions_df': molecular_functions_df_loaded_from_csv = True  # Mark as available even if from cache
            continue

        csv_path = Path(output_dir_ref) / file_name
        if csv_path.exists():
            print(f"DEBUG: Loading from CSV: {file_name}...")
            try:
                if key_name in ['rmsd_df', 'ca_distance_table', 'distance_table', 'msa_table', 'residue_table']:
                    df = pd.read_csv(csv_path, index_col=0)
                else:
                    df = pd.read_csv(csv_path)
                data[key_name] = df
                print(f"  DEBUG: Successfully loaded {file_name} into '{key_name}'. Shape: {df.shape}")
                if key_name == 'molecular_functions_df':
                    molecular_functions_df_loaded_from_csv = True
                    print(f"    DEBUG_CSV: molecular_functions_df columns: {df.columns.tolist()}")
                    if not df.empty:
                        print(f"    DEBUG_CSV: molecular_functions_df head:\n{df.head()}")
            except Exception as e:
                print(f"  DEBUG: Error loading {file_name}: {e}")
        else:
            print(f"  DEBUG: Warning: CSV file not found: {csv_path}")

    # 3. Post-processing and data structuring
    print("DEBUG: --- Populating group_dict (Molecular Function) ---")
    data['group_dict'] = {}
    populated_mf_from_processed_structures = False
    if 'processed_structures' in data and isinstance(data['processed_structures'], dict):
        temp_group_dict = {}
        unknown_mf_count_cache = 0
        known_mf_count_cache = 0
        for sid, s_data in data['processed_structures'].items():
            mf = "Unknown"
            if isinstance(s_data, dict) and 'properties' in s_data and isinstance(s_data['properties'], dict):
                mf_val = s_data['properties'].get('molecular_function')  # Use .get for safety
                if pd.notna(mf_val) and mf_val != "Unknown" and str(mf_val).strip() != "":
                    mf = str(mf_val)
                    known_mf_count_cache += 1
                else:
                    unknown_mf_count_cache += 1
            else:
                unknown_mf_count_cache += 1
            temp_group_dict[sid] = mf

        if temp_group_dict:
            data['group_dict'] = temp_group_dict
            populated_mf_from_processed_structures = True
            print(f"DEBUG: Populated 'group_dict' from 'processed_structures' cache.")
            print(
                f"DEBUG:   Known MF from cache: {known_mf_count_cache}, Unknown MF from cache: {unknown_mf_count_cache}")
        else:
            print(f"DEBUG: 'processed_structures' cache present but yielded no 'group_dict' entries.")

    if not populated_mf_from_processed_structures and 'molecular_functions_df' in data:
        df_mol_func = data['molecular_functions_df']
        if isinstance(df_mol_func,
                      pd.DataFrame) and 'structure_id' in df_mol_func.columns and 'molecular_function' in df_mol_func.columns:
            temp_group_dict_csv = {}
            unknown_mf_count_csv, known_mf_count_csv = 0, 0
            for _, row in df_mol_func.iterrows():
                sid_csv = str(row['structure_id'])
                mf_csv = row['molecular_function']
                if pd.notna(mf_csv) and mf_csv != "Unknown" and str(mf_csv).strip() != "":
                    temp_group_dict_csv[sid_csv] = str(mf_csv)
                    known_mf_count_csv += 1
                else:
                    temp_group_dict_csv[sid_csv] = "Unknown"
                    unknown_mf_count_csv += 1

            if temp_group_dict_csv:
                data['group_dict'] = temp_group_dict_csv
                print(f"DEBUG: Populated 'group_dict' from 'molecular_functions_df' (CSV).")
                print(f"DEBUG:   Known MF from CSV: {known_mf_count_csv}, Unknown MF from CSV: {unknown_mf_count_csv}")
            else:
                print(f"DEBUG: 'molecular_functions_df' CSV present but yielded no 'group_dict' entries.")

    elif not populated_mf_from_processed_structures and (
            'molecular_functions_df' not in data or not molecular_functions_df_loaded_from_csv):
        print(
            f"DEBUG: Neither 'processed_structures' cache nor 'molecular_functions.csv' yielded 'group_dict'. Attempting PDB ID fallback.")
        if 'rmsd_df' in data and isinstance(data['rmsd_df'], pd.DataFrame):
            print("  DEBUG: Deriving 'group_dict' from PDB IDs in 'rmsd_df' as a fallback.")
            temp_group_dict_fallback = {}
            for sid_fallback in data['rmsd_df'].index:
                if 'ChR' in sid_fallback or 'channel' in sid_fallback.lower():
                    temp_group_dict_fallback[sid_fallback] = "Cation channel"
                elif 'HR' in sid_fallback or 'pump' in sid_fallback.lower() or 'PR' in sid_fallback:
                    temp_group_dict_fallback[sid_fallback] = "Proton pump"
                elif 'ACR' in sid_fallback or 'chloride' in sid_fallback.lower():
                    temp_group_dict_fallback[sid_fallback] = "Chloride pump"
                else:
                    temp_group_dict_fallback[sid_fallback] = "Unknown"
            if temp_group_dict_fallback: data['group_dict'] = temp_group_dict_fallback

    if not data['group_dict']:
        print("DEBUG: 'group_dict' is empty after all attempts. Defaulting all SIDs in rmsd_df to Unknown.")
        if 'rmsd_df' in data and isinstance(data['rmsd_df'], pd.DataFrame):
            data['group_dict'] = {sid: "Unknown" for sid in data['rmsd_df'].index}

    unknown_final_mf_count = sum(1 for v in data['group_dict'].values() if v == "Unknown" or pd.isna(v))
    total_final_mf_count = len(data['group_dict'])
    print(
        f"DEBUG: Final 'group_dict' for MF has {total_final_mf_count} entries. {unknown_final_mf_count} are 'Unknown'.")
    if 0 < total_final_mf_count < 10: print(f"DEBUG: Final group_dict content: {data['group_dict']}")

    print("DEBUG: --- Populating domain_dict_for_plots (Domain & Error) ---")
    data['domain_dict_for_plots'] = {}
    all_sids_for_domain_dict = []
    if 'rmsd_df' in data and isinstance(data['rmsd_df'], pd.DataFrame):
        all_sids_for_domain_dict = data['rmsd_df'].index.tolist()

    unknown_domain_count, known_domain_count = 0, 0
    for sid in all_sids_for_domain_dict:
        domain_info = {'domain': "Unknown", 'average_error': None}
        domain_source = "default_unknown"

        if 'processed_structures' in data and sid in data['processed_structures']:
            s_proc_data = data['processed_structures'][sid]
            if isinstance(s_proc_data, dict) and 'properties' in s_proc_data and isinstance(s_proc_data['properties'],
                                                                                            dict):
                props = s_proc_data['properties']
                domain_val = props.get('domain')
                if pd.notna(domain_val) and domain_val != "Unknown" and str(domain_val).strip() != "":
                    domain_info['domain'] = str(domain_val)
                    domain_source = "cache_properties"

        # Get average_error (prioritizing structure_errors cache, then specific error CSVs)
        error_source = "N/A"
        if 'structure_errors' in data and isinstance(data['structure_errors'], dict) and sid in data[
            'structure_errors']:
            s_error_data = data['structure_errors'][sid]
            if isinstance(s_error_data, dict) and 'average_error' in s_error_data:
                avg_err_val = s_error_data['average_error']
                if isinstance(avg_err_val, (int, float)) and pd.notna(avg_err_val):
                    domain_info['average_error'] = float(avg_err_val)
                    error_source = "cache_structure_errors"

        if domain_info['average_error'] is None:  # Fallback to error CSVs
            for err_df_key in ['mo_exp_errors_df', 'hideaki_errors_df']:
                if err_df_key in data and isinstance(data[err_df_key], pd.DataFrame):
                    err_df = data[err_df_key]
                    if 'structure_id' in err_df.columns and sid in err_df['structure_id'].values:
                        row = err_df[err_df['structure_id'] == sid].iloc[0]
                        error_cols = ['backbone_rmsd', 'pocket_rmsd', 'retinal_rmsd',
                                      'average_error']  # Check 'average_error' too
                        valid_errors = [row[col] for col in error_cols if
                                        col in row and pd.notna(row[col]) and isinstance(row[col], (int, float))]
                        if valid_errors:
                            domain_info['average_error'] = float(np.mean(valid_errors))
                            error_source = f"csv_{err_df_key}"
                            # Update domain from error CSV if cache didn't provide it and CSV has it
                            if domain_info['domain'] == "Unknown" and 'domain' in row and pd.notna(row['domain']):
                                domain_info['domain'] = str(row['domain'])
                                domain_source = f"csv_err_domain_{err_df_key}"
                            break

        if domain_info['domain'] == "Unknown":
            unknown_domain_count += 1
        else:
            known_domain_count += 1

        if sid in ['7D77_A', '7BZ2_R', 'AF-P0C6S6-F1-model_v4', 'INSERT_A_KNOWN_PDB_ID_HERE']:  # Add specific PDBs
            print(
                f"  DEBUG_DOMAIN {sid}: final domain='{domain_info['domain']}' (source: {domain_source}), error='{domain_info['average_error']}' (source: {error_source})")

        data['domain_dict_for_plots'][sid] = domain_info

    print(f"DEBUG: Populated 'domain_dict_for_plots' for {len(data['domain_dict_for_plots'])} structures.")
    print(f"DEBUG:   Known domains: {known_domain_count}, Unknown domains: {unknown_domain_count}")

    if not data['domain_dict_for_plots'] and 'rmsd_df' in data and isinstance(data['rmsd_df'], pd.DataFrame):
        print("  DEBUG: Warning: 'domain_dict_for_plots' is empty. Filling with Unknown domains and no errors.")
        data['domain_dict_for_plots'] = {
            sid: {'domain': "Unknown", 'average_error': None} for sid in data['rmsd_df'].index
        }

    # Ensure 'msa_table' and 'residue_table' are available
    if 'msa_table' not in data and 'ca_msa_table_grn' in data:
        data['msa_table'] = data['ca_msa_table_grn']
        print("  DEBUG: Used 'ca_msa_table_grn' as 'msa_table'.")
    if 'residue_table' not in data and 'ca_residue_table_grn' in data:  # Prefer ca_residue_table if specific residue_table is missing
        data['residue_table'] = data['ca_residue_table_grn']
        print("  DEBUG: Used 'ca_residue_table_grn' as 'residue_table'.")
    elif 'residue_table' not in data and 'msa_table' in data:  # Fallback to msa_table
        data['residue_table'] = data['msa_table']
        print("  DEBUG: Used 'msa_table' as 'residue_table'.")

    # For opsin overview plot
    if 'processed_structures' in data and 'rmsd_df' in data and isinstance(data['rmsd_df'], pd.DataFrame):
        overview_data_list = []
        for sid in data['rmsd_df'].index:
            s_data = data.get('processed_structures', {}).get(sid, {})  # .get for safety
            props = s_data.get('properties', {}) if isinstance(s_data, dict) else {}

            mf_overview = data.get('group_dict', {}).get(sid, "Unknown")  # Use already populated group_dict
            domain_overview = data.get('domain_dict_for_plots', {}).get(sid, {}).get('domain',
                                                                                     "Unknown")  # Use already populated

            overview_data_list.append({
                'short_name': props.get('display_name', sid) if isinstance(props, dict) else sid,
                'molecular_function_normalized': mf_overview,
                'domain': domain_overview,
                'experimentally_determined': props.get('is_experimental', False) if isinstance(props, dict) else False
            })
        data['overview_df'] = pd.DataFrame(overview_data_list)
        print(f"  DEBUG: Prepared 'overview_df' for opsin overview plot with {len(data['overview_df'])} entries.")
        if not data['overview_df'].empty:
            print(f"    DEBUG: overview_df head:\n{data['overview_df'].head()}")

    print("DEBUG: Data loading and initial structuring complete.")
    return data


def generate_summary_csv(data, output_dir):
    """
    Generate a CSV file with protein info: name, avg RMSD, domain, function, error.
    """
    summary_data = []
    if 'rmsd_df' not in data or not isinstance(data['rmsd_df'], pd.DataFrame):
        print("Skipping summary CSV: 'rmsd_df' not available or not a DataFrame.")
        return None

    rmsd_df = data['rmsd_df']
    processed_structures = data.get('processed_structures', {})
    domain_dict_for_plots = data.get('domain_dict_for_plots', {})
    group_dict_for_summary = data.get('group_dict', {})  # Use the final group_dict

    print("DEBUG: --- Generating Summary CSV ---")
    sids_in_rmsd = rmsd_df.index.tolist()
    print(f"DEBUG: Number of SIDs in rmsd_df for summary: {len(sids_in_rmsd)}")

    for protein_id in sids_in_rmsd:
        rmsd_values = rmsd_df.loc[protein_id].values
        valid_rmsd_values = rmsd_values[~np.isnan(rmsd_values) & (rmsd_values > 1e-6)]  # Exclude self (0) and NaNs
        avg_rmsd_val = np.mean(valid_rmsd_values) if len(valid_rmsd_values) > 0 else np.nan

        s_data = processed_structures.get(protein_id, {})
        props = s_data.get('properties', {}) if isinstance(s_data, dict) else {}
        display_name = props.get('display_name', protein_id) if isinstance(props, dict) else protein_id

        domain_plot_info = domain_dict_for_plots.get(protein_id, {'domain': "Unknown", 'average_error': np.nan})
        domain = domain_plot_info['domain']
        error_val = domain_plot_info['average_error']

        function = group_dict_for_summary.get(protein_id, "Unknown")

        if protein_id in ['7D77_A', '7BZ2_R', 'AF-P0C6S6-F1-model_v4', 'INSERT_A_KNOWN_PDB_ID_HERE']:  # Add PDBs
            print(
                f"  DEBUG_SUMMARY {protein_id}: Name='{display_name}', AvgRMSD={avg_rmsd_val:.2f if pd.notna(avg_rmsd_val) else 'N/A'}, Domain='{domain}', Function='{function}', Error='{error_val:.2f if pd.notna(error_val) else 'N/A'}'")

        summary_data.append([
            display_name,
            round(avg_rmsd_val, 2) if pd.notna(avg_rmsd_val) else 'N/A',
            domain,
            function,
            round(error_val, 2) if pd.notna(error_val) else 'N/A'
        ])

    summary_df = pd.DataFrame(summary_data,
                              columns=['Protein', 'Average RMSD', 'Domain', 'Molecular Function', 'Error'])
    summary_df = summary_df.sort_values(by=['Molecular Function', 'Domain', 'Average RMSD'],
                                        na_position='last')  # Improved sort

    csv_path = Path(output_dir) / 'protein_summary.csv'
    summary_df.to_csv(csv_path, index=False)
    print(f"Saved protein summary to {csv_path}")

    # Save color assignments
    all_functions = sorted(list(set(summary_df['Molecular Function'].unique()) - {'Unknown'}))
    all_domains = sorted(list(set(summary_df['Domain'].unique()) - {'Unknown'}))

    function_colors = get_group_colors(all_functions)  # Uses opsin_color_scheme.py
    domain_colors = get_group_colors(all_domains)

    color_data = []
    for func, color in function_colors.items(): color_data.append(['Function', func, color])
    for dom, color in domain_colors.items(): color_data.append(['Domain', dom, color])

    if color_data:  # Only save if we have colors
        color_df = pd.DataFrame(color_data, columns=['Type', 'Value', 'Color'])
        color_path = Path(output_dir) / 'color_assignments.csv'
        color_df.to_csv(color_path, index=False)
        print(f"Saved color assignments to {color_path}")
    else:
        print("DEBUG: No color assignments to save (all functions/domains might be 'Unknown').")

    return csv_path


def generate_plots(data, output_dir_figures, error_threshold=3.0):
    """
    Generate all curated plots and save them to the output directory.
    """
    print(f"\n--- Generating plots in {output_dir_figures} ---")
    Path(output_dir_figures).mkdir(parents=True, exist_ok=True)

    # --- A. Overview & Reference ---
    print("\nDEBUG: 1. Generating Opsins Overview Plot...")
    if 'overview_df' in data and isinstance(data['overview_df'], pd.DataFrame) and not data['overview_df'].empty:
        try:
            fig = create_opsin_overview_plot(data['overview_df'])
            fig.savefig(Path(output_dir_figures) / 'A1_opsin_overview.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("  DEBUG: Saved A1_opsin_overview.png")
        except Exception as e:
            print(f"  DEBUG: Error generating opsin overview plot: {e}")
    else:
        print("  DEBUG: Skipped Opsins Overview: 'overview_df' not available, not a DataFrame, or empty.")
        if 'overview_df' in data: print(
            f"    DEBUG: overview_df type: {type(data['overview_df'])}, empty: {data['overview_df'].empty if isinstance(data['overview_df'], pd.DataFrame) else 'N/A'}")

    print("\nDEBUG: 2. Generating RMSD Color Scale Reference...")
    try:
        fig = create_rmsd_color_scale_figure()
        fig.savefig(Path(output_dir_figures) / 'A2_rmsd_color_scale.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("  DEBUG: Saved A2_rmsd_color_scale.png")
    except Exception as e:
        print(f"  DEBUG: Error generating RMSD color scale: {e}")

    if 'rmsd_df' not in data or not isinstance(data['rmsd_df'], pd.DataFrame) or data['rmsd_df'].empty:
        print("\nDEBUG: CRITICAL: 'rmsd_df' is missing or empty. Skipping all RMSD-based plots (B, C).")
    else:
        rmsd_df_orig = data['rmsd_df']
        print(f"DEBUG: Using rmsd_df_orig with shape {rmsd_df_orig.shape} for plots B and C.")

        # --- C. RMSD & Structural Similarity (Simple/Unfiltered) ---
        print("\nDEBUG: 3. Generating Simple Unfiltered Structural Similarity Tree...")
        try:
            fig_simple_tree = plot_similarity_tree(
                rmsd_df_orig,
                title="Unfiltered Structural Similarity Tree"
            )
            fig_simple_tree.savefig(Path(output_dir_figures) / 'C3_similarity_tree_unfiltered.png', dpi=300,
                                    bbox_inches='tight')
            plt.close(fig_simple_tree)
            print("  DEBUG: Saved C3_similarity_tree_unfiltered.png")
        except Exception as e:
            print(f"  DEBUG: Error generating simple unfiltered similarity tree: {e}")

        # --- B. RMSD & Structural Similarity (Linked & Filtered) ---
        print(f"\nDEBUG: --- Preparing for Linked & Filtered RMSD Plots (Threshold: {error_threshold} Å) ---")
        pdb_list_orig = rmsd_df_orig.index.tolist()
        group_dict_orig = data.get('group_dict', {sid: "Unknown" for sid in pdb_list_orig})
        domain_dict_for_plots = data.get('domain_dict_for_plots',
                                         {sid: {'domain': "Unknown", 'average_error': None} for sid in pdb_list_orig})

        avg_rmsds_from_matrix = rmsd_df_orig.mean(axis=1)  # Mean RMSD of a structure to all others

        kept_ids = []
        filtered_out_details = []
        print("  DEBUG: Filtering structures for linked plots...")
        for pdb_id in pdb_list_orig:
            error_val = None
            source = "N/A"
            # Try explicit average_error from domain_dict_for_plots first
            if pdb_id in domain_dict_for_plots and domain_dict_for_plots[pdb_id].get(
                    'average_error') is not None and pd.notna(domain_dict_for_plots[pdb_id]['average_error']):
                error_val = domain_dict_for_plots[pdb_id]['average_error']
                source = "domain_dict_avg_error"
            # Fallback to average RMSD from matrix if explicit error not available
            elif pdb_id in avg_rmsds_from_matrix and pd.notna(avg_rmsds_from_matrix[pdb_id]):
                error_val = avg_rmsds_from_matrix[pdb_id]
                source = "avg_rmsd_from_matrix"

            if error_val is not None and error_val <= error_threshold:
                kept_ids.append(pdb_id)
            elif error_val is not None:
                filtered_out_details.append(
                    f"    - {pdb_id}: value={error_val:.2f} Å (source: {source}) > {error_threshold}")
            else:
                kept_ids.append(pdb_id)
                # print(f"    DEBUG: {pdb_id}: No error data or avg RMSD, kept by default.")

        if filtered_out_details:
            print(f"  DEBUG: Filtered out {len(filtered_out_details)} structures:")
            for detail in filtered_out_details[:5]: print(detail)
            if len(filtered_out_details) > 5: print(f"    ... and {len(filtered_out_details) - 5} more.")
        else:
            print("  DEBUG: No structures filtered based on error/RMSD threshold for linked plots.")

        if len(kept_ids) < 2:
            print("  DEBUG: Skipped Linked Plots: Less than 2 structures remaining after filtering.")
        else:
            print(f"  DEBUG: Proceeding with {len(kept_ids)} structures for linked visualizations.")
            filtered_rmsd_df = rmsd_df_orig.loc[kept_ids, kept_ids].copy()

            filtered_group_dict = {k: v for k, v in group_dict_orig.items() if k in kept_ids}
            filtered_domain_dict_for_plots = {k: v for k, v in domain_dict_for_plots.items() if k in kept_ids}

            Z_linkage = None
            try:
                print("  DEBUG: Calculating linkage matrix Z for filtered data...")
                matrix_for_linkage = filtered_rmsd_df.values.copy()
                # Handle NaNs and Infs before squareform
                if np.any(~np.isfinite(matrix_for_linkage)):
                    mean_finite = np.nanmean(matrix_for_linkage[np.isfinite(matrix_for_linkage)])
                    fill_val = mean_finite if pd.notna(mean_finite) else max(10.0, np.nanmax(
                        matrix_for_linkage[np.isfinite(matrix_for_linkage)] * 1.1 if np.any(
                            np.isfinite(matrix_for_linkage)) else 10.0))  # Robust fallback
                    matrix_for_linkage = np.nan_to_num(matrix_for_linkage, nan=fill_val, posinf=fill_val, neginf=0.0)
                np.fill_diagonal(matrix_for_linkage, 0.0)
                if not np.allclose(matrix_for_linkage, matrix_for_linkage.T, atol=1e-5):
                    print("    DEBUG: Warning: RMSD matrix for linkage is not perfectly symmetric. Forcing symmetry.")
                    matrix_for_linkage = (matrix_for_linkage + matrix_for_linkage.T) / 2.0
                    np.fill_diagonal(matrix_for_linkage, 0.0)

                condensed_matrix = squareform(matrix_for_linkage, checks=True)
                Z_linkage = linkage(condensed_matrix, method='average')
                print("    DEBUG: Successfully calculated linkage matrix Z.")
            except Exception as e_link:
                print(f"    DEBUG: ERROR calculating linkage matrix: {e_link}.")
                # import traceback; traceback.print_exc() # Uncomment for full trace

            if Z_linkage is not None:
                print("\nDEBUG: 4. Generating Filtered Structural Similarity Tree (Linked)...")
                try:
                    tree_fig, ordered_tree_ids = create_and_visualize_similarity_tree(
                        rmsd_data=filtered_rmsd_df,  # Pass DataFrame
                        linkage_matrix=Z_linkage,
                        group_dict=filtered_group_dict,
                        domain_dict=filtered_domain_dict_for_plots
                    )
                    tree_fig.savefig(Path(output_dir_figures) / 'B4_similarity_tree_filtered_linked.png', dpi=300,
                                     bbox_inches='tight')
                    plt.close(tree_fig)
                    print("  DEBUG: Saved B4_similarity_tree_filtered_linked.png")
                except Exception as e_tree:
                    print(f"  DEBUG: Error generating filtered similarity tree: {e_tree}")
                    # import traceback; traceback.print_exc()

                print("\nDEBUG: 5. Generating Filtered RMSD Heatmap (Linked & Clustered)...")
                try:
                    # visualize_rmsd_matrix_improved is expected to save its own figure
                    heatmap_clustermap = visualize_rmsd_matrix_improved(
                        rmsd_df=filtered_rmsd_df,
                        linkage_matrix=Z_linkage,
                        group_dict=filtered_group_dict,
                        domain_dict=filtered_domain_dict_for_plots,  # This is the complex dict
                        output_file=Path(output_dir_figures) / 'B5_rmsd_heatmap_filtered_linked.png'
                    )
                    if heatmap_clustermap:
                        plt.close(heatmap_clustermap.fig)  # Close the figure associated with clustermap
                        print("  DEBUG: Saved B5_rmsd_heatmap_filtered_linked.png (by visualize_rmsd_matrix_improved)")
                    else:
                        print(
                            "  DEBUG: Clustermap generation failed or returned None from visualize_rmsd_matrix_improved.")
                except Exception as e_heatmap:
                    print(f"  DEBUG: Error generating filtered RMSD heatmap: {e_heatmap}")
                    # import traceback; traceback.print_exc()
            else:
                print("  DEBUG: Skipped Linked Tree and Heatmap due to Z_linkage calculation failure.")

    # --- D. Distance-to-Retinal Plots ---
    # ... (Your D plots logic - add DEBUG prints similarly if needed) ...
    plot_configs_dist = [
        ('ca_distance_table', 'D6_ca_distance_retinal_std.png', "CA Distance to Retinal by Position", True),
        ('distance_table', 'D7_sidechain_distance_retinal_std.png', "Sidechain Distance to Retinal by Position", False),
    ]
    for key, fname, title, use_ca_flag in plot_configs_dist:
        print(f"\nDEBUG: {fname.split('_')[0]} {title}...")
        if key in data and isinstance(data[key], pd.DataFrame) and not data[key].empty:
            try:
                fig = plot_distances_with_std(data[key], title=title, use_ca=use_ca_flag)
                fig.savefig(Path(output_dir_figures) / fname, dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"  DEBUG: Saved {fname}")
            except Exception as e:
                print(f"  DEBUG: Error generating {title}: {e}")
        else:
            print(f"  DEBUG: Skipped {title}: Data key '{key}' not available, not DataFrame, or empty.")

    # --- E. Sequence Conservation & Composition ---
    # ... (Your E plots logic - add DEBUG prints similarly if needed) ...
    print("\nDEBUG: 12. Residue Conservation Bar Plot...")
    if 'msa_table' in data and isinstance(data['msa_table'], pd.DataFrame) and not data['msa_table'].empty:
        try:
            fig = create_residue_conservation_plot(data['msa_table'], helix_highlighting=True)
            fig.savefig(Path(output_dir_figures) / 'E12_residue_conservation.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("  DEBUG: Saved E12_residue_conservation.png")
        except Exception as e:
            print(f"  DEBUG: Error generating residue conservation plot: {e}")
    else:
        print("  DEBUG: Skipped Residue Conservation: 'msa_table' not available, not DataFrame, or empty.")

    # --- F. Combined Plots (Optional) ---
    # ... (Your F plots logic - add DEBUG prints similarly if needed) ...
    print("\nDEBUG: 15. Combined Distance Line Plot & Sequence Logo...")
    if 'distance_table' in data and isinstance(data['distance_table'], pd.DataFrame) and not data[
        'distance_table'].empty and \
            'msa_table' in data and isinstance(data['msa_table'], pd.DataFrame) and not data['msa_table'].empty:
        try:
            fig = create_combined_distance_logo_plot(data['distance_table'], data['msa_table'])
            fig.savefig(Path(output_dir_figures) / 'F15_combined_distance_logo.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("  DEBUG: Saved F15_combined_distance_logo.png")
        except Exception as e:
            print(
                f"  DEBUG: Error generating combined distance logo plot: {e}")  # import traceback; traceback.print_exc()
    else:
        print(
            "  DEBUG: Skipped Combined Distance Logo: 'distance_table' or 'msa_table' not available, not DataFrame, or empty.")

    print(f"\n--- All plot generation attempts complete. Check {output_dir_figures}. ---")


def verify_filtering_in_summary(summary_csv_path, error_threshold=3.0):
    if not summary_csv_path or not Path(summary_csv_path).exists():
        print("\nDEBUG: Verification: Protein summary CSV not found for verification.")
        return

    print(f"\n--- Verifying Filtering in Protein Summary ({summary_csv_path}) ---")
    df = pd.read_csv(summary_csv_path)
    violations = 0
    if 'Average RMSD' in df.columns:
        df['Average RMSD Numeric'] = pd.to_numeric(df['Average RMSD'], errors='coerce')
        high_rmsd_entries = df[df['Average RMSD Numeric'] > error_threshold]
        if not high_rmsd_entries.empty:
            print(f"  DEBUG: WARNING: Found {len(high_rmsd_entries)} entries with Average RMSD > {error_threshold}:")
            violations += len(high_rmsd_entries)
    if 'Error' in df.columns:
        df['Error Numeric'] = pd.to_numeric(df['Error'], errors='coerce')
        high_error_entries = df[df['Error Numeric'] > error_threshold]
        if not high_error_entries.empty:
            print(f"  DEBUG: WARNING: Found {len(high_error_entries)} entries with Error > {error_threshold}:")
            violations += len(high_error_entries)
    if violations == 0:
        print(
            f"  DEBUG: Verification PASSED (summary check): No entries found exceeding threshold {error_threshold} Å.")
    else:
        print(f"  DEBUG: Verification FAILED (summary check): {violations} potential violations.")


def main():
    parser = argparse.ArgumentParser(description='Generate plots for opsin analysis.')
    parser.add_argument('--input-dir', '-i', type=str, default='opsin_output/',
                        help='Directory containing input data files (cache subdir, CSVs).')
    parser.add_argument('--output-dir', '-o', type=str, default='opsin_plots_output/',  # Changed default
                        help='Directory to save output plots and summaries.')
    parser.add_argument('--chain-id', '-c', type=str, default='A',
                        help='Chain ID used in the analysis (default: A)')
    parser.add_argument('--quality', '-q', type=str, choices=['low', 'medium', 'high'], default='high',
                        help='Figure quality (affects DPI). Default: high')
    parser.add_argument('--error-threshold', '-et', type=float, default=3.0,
                        help='RMSD/Error threshold in Angstrom for filtering linked plots. Default: 3.0')
    args = parser.parse_args()

    output_dir_figures = Path(args.output_dir) / 'figures'
    output_dir_summaries = Path(args.output_dir) / 'summaries'
    output_dir_figures.mkdir(parents=True, exist_ok=True)
    output_dir_summaries.mkdir(parents=True, exist_ok=True)

    print(f"Opsin Plot Generation Script - DEBUG ENABLED")
    print(f"Input Data Directory: {Path(args.input_dir).resolve()}")
    print(f"Output Directory (Plots): {output_dir_figures.resolve()}")
    print(f"Output Directory (Summaries): {output_dir_summaries.resolve()}")

    dpi_map = {'low': 100, 'medium': 200, 'high': 300}
    plt.rcParams['figure.dpi'] = dpi_map[args.quality]
    plt.rcParams['savefig.dpi'] = dpi_map[args.quality]

    loaded_data = load_data(args.input_dir, args.input_dir, args.chain_id)

    if not loaded_data or not isinstance(loaded_data, dict):  # Added check for dict
        print("No data loaded or data is not a dictionary. Exiting.")
        return

    generate_plots(loaded_data, output_dir_figures, error_threshold=args.error_threshold)
    summary_csv_path = generate_summary_csv(loaded_data, output_dir_summaries)
    if summary_csv_path:
        verify_filtering_in_summary(summary_csv_path, error_threshold=args.error_threshold)

    print(f"\nAll operations complete. Outputs are in {Path(args.output_dir).resolve()}")


if __name__ == "__main__":
    main()