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
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import json  # For summary saving
from src.data_processing import load_opsin_property_data  # Import property loading function

# --- Import visualization functions ---
try:
    from src.visualization_functions import (
        create_opsin_overview_plot, create_rmsd_color_scale_figure,
        create_and_visualize_similarity_tree, visualize_rmsd_matrix_improved,
        plot_similarity_tree, plot_rmsd_heatmap,
        plot_distances_with_std, plot_average_distances_by_helix, plot_distance_heatmap,
        create_residue_conservation_plot, plot_helix_logo_plots, plot_conservation_around_x50,
        create_combined_distance_logo_plot, visualize_binding_pocket
    )

    VIS_FUNCS_LOADED = True
except ImportError as e1:
    print(f"Error importing from visualization_functions: {e1}")
    print("Attempting to import from current directory as a fallback for development.")
    try:
        from visualization_functions import (  # type: ignore # Fallback for local dev
            create_opsin_overview_plot, create_rmsd_color_scale_figure,
            create_and_visualize_similarity_tree, visualize_rmsd_matrix_improved,
            plot_similarity_tree, plot_rmsd_heatmap,
            plot_distances_with_std, plot_average_distances_by_helix, plot_distance_heatmap,
            create_residue_conservation_plot, plot_helix_logo_plots, plot_conservation_around_x50,
            create_combined_distance_logo_plot, visualize_binding_pocket
        )

        VIS_FUNCS_LOADED = True
    except ImportError as e2:
        print(f"Fallback import failed: {e2}. Plotting functions will be dummied.")
        VIS_FUNCS_LOADED = False
        import inspect  # Keep inspect import local to this block


        # Define dummy functions
        def _dummy_plot_func(*args, **kwargs):
            fn_name = inspect.stack()[1].function if len(inspect.stack()) > 1 else "UnknownFunction"
            print(
                f"Warning: Plot function '{fn_name}' called but not loaded due to import error. Args: {args}, Kwargs: {kwargs}")
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, f"Plotting function\n{fn_name}\nunavailable", ha='center', va='center', color='red')
            # Adjust for functions returning a tuple (fig, other_data)
            if fn_name in ["create_and_visualize_similarity_tree"]:
                return fig, None
            return fig


        create_opsin_overview_plot = _dummy_plot_func
        create_rmsd_color_scale_figure = _dummy_plot_func
        create_and_visualize_similarity_tree = _dummy_plot_func
        visualize_rmsd_matrix_improved = _dummy_plot_func
        plot_similarity_tree = _dummy_plot_func
        plot_rmsd_heatmap = _dummy_plot_func
        plot_distances_with_std = _dummy_plot_func
        plot_average_distances_by_helix = _dummy_plot_func
        plot_distance_heatmap = _dummy_plot_func
        create_residue_conservation_plot = _dummy_plot_func
        plot_helix_logo_plots = _dummy_plot_func
        plot_conservation_around_x50 = _dummy_plot_func
        create_combined_distance_logo_plot = _dummy_plot_func
        visualize_binding_pocket = _dummy_plot_func

# Import the color scheme tools
try:
    from src.opsin_color_scheme import get_group_colors, RMSD_COMPACT_CMAP

    COLOR_SCHEME_LOADED = True
except ImportError:
    print("Error: opsin_color_scheme.py not found. Using fallback colors.")
    COLOR_SCHEME_LOADED = False


    def get_group_colors(items_list_or_dict, palette_name=None):  # type: ignore
        if isinstance(items_list_or_dict, dict):
            items = list(items_list_or_dict.keys())
        else:
            items = list(items_list_or_dict)
        num_items = len(items)
        # Generate a simple list of distinct colors for fallback
        fallback_palette = plt.cm.get_cmap('viridis', num_items if num_items > 0 else 1)
        return {item: fallback_palette(i) for i, item in enumerate(items)}


    RMSD_COMPACT_CMAP = "viridis"  # Fallback colormap


def load_data(input_dir_path_str, output_dir_ref_str, chain_id='A'):  # Renamed args for clarity
    """
    Load precomputed data from workflow cache or CSV files, and consolidate properties.
    """
    data = {}
    input_dir = Path(input_dir_path_str)  # Convert to Path
    output_dir_ref = Path(output_dir_ref_str)  # Convert to Path
    cache_dir = input_dir / 'cache'
    print(f"DEBUG: PlotFigures - Loading data. Cache: {cache_dir}, CSV Ref: {output_dir_ref.resolve()}")

    # --- Stage 1: Load Core Data from Cache (Primary) & CSVs (Fallback) ---
    proc_struct_fname = f'processed_structures_{chain_id}.pkl'
    cache_path_proc_struct = cache_dir / proc_struct_fname
    if cache_path_proc_struct.exists():
        try:
            with open(cache_path_proc_struct, 'rb') as f:
                loaded_proc_data = pickle.load(f)
                if isinstance(loaded_proc_data, dict) and 'processed_structures' in loaded_proc_data and isinstance(
                        loaded_proc_data['processed_structures'], dict):
                    data['processed_structures'] = loaded_proc_data['processed_structures']
                elif isinstance(loaded_proc_data, dict):  # if pkl directly saves the dict
                    data['processed_structures'] = loaded_proc_data
                else:
                    print(f"  DEBUG: Unexpected data type in {proc_struct_fname}: {type(loaded_proc_data)}")
            print(
                f"DEBUG: PlotFigures - Loaded {proc_struct_fname} ({len(data.get('processed_structures', {}))} entries)")
        except Exception as e:
            print(f"DEBUG: PlotFigures - Error loading {proc_struct_fname}: {e}")
    else:
        print(f"DEBUG: PlotFigures - Cache file {proc_struct_fname} not found at {cache_path_proc_struct}.")

    rmsd_cache_fname = f'structure_comparison_{chain_id}.pkl'
    cache_path_rmsd = cache_dir / rmsd_cache_fname
    if cache_path_rmsd.exists():
        try:
            with open(cache_path_rmsd, 'rb') as f:
                comp_data = pickle.load(f)
                # Prioritize 'rmsd_df' then 'rmsd_matrix'
                if 'rmsd_df' in comp_data and isinstance(comp_data['rmsd_df'], pd.DataFrame):
                    data['rmsd_df'] = comp_data['rmsd_df']
                elif 'rmsd_matrix' in comp_data and isinstance(comp_data['rmsd_matrix'], pd.DataFrame):
                    data['rmsd_df'] = comp_data['rmsd_matrix']  # Store as 'rmsd_df' for consistency

                if 'group_dict' in comp_data: data['group_dict_from_comparison_cache'] = comp_data['group_dict']
                if 'domain_dict' in comp_data: data['domain_dict_from_comparison_cache'] = comp_data['domain_dict']
                if 'pdb_list' in comp_data: data['pdb_list_from_comparison_cache'] = comp_data['pdb_list']

            print(
                f"DEBUG: PlotFigures - Loaded {rmsd_cache_fname}. RMSD_DF shape: {data.get('rmsd_df', pd.DataFrame()).shape}")
        except Exception as e:
            print(f"DEBUG: PlotFigures - Error loading {rmsd_cache_fname}: {e}")

    if 'rmsd_df' not in data or data['rmsd_df'].empty:
        csv_path_rmsd = output_dir_ref / 'rmsd_matrix.csv'
        if csv_path_rmsd.exists():
            try:
                data['rmsd_df'] = pd.read_csv(csv_path_rmsd, index_col=0)
                print(f"DEBUG: PlotFigures - Loaded rmsd_matrix.csv. Shape: {data['rmsd_df'].shape}")
            except Exception as e:
                print(f"DEBUG: PlotFigures - Error loading rmsd_matrix.csv: {e}")
        else:
            print("DEBUG: PlotFigures - CRITICAL: RMSD matrix not found in cache or CSV. Many plots might fail.")

    grn_cache_fname = f'grn_assignment_{chain_id}.pkl'
    cache_path_grn = cache_dir / grn_cache_fname
    if cache_path_grn.exists():
        try:
            with open(cache_path_grn, 'rb') as f:
                grn_cache_data = pickle.load(f)
                if isinstance(grn_cache_data, dict):
                    # Define keys expected from grn_assignment.pkl
                    grn_keys_from_cache = ["residue_table", "distance_table",
                                           "ca_residue_table", "ca_distance_table", "msa_table"]
                    for key in grn_keys_from_cache:
                        if key in grn_cache_data and isinstance(grn_cache_data[key], pd.DataFrame):
                            data[key] = grn_cache_data[key]
                            print(
                                f"  DEBUG: PlotFigures - Loaded '{key}' from {grn_cache_fname}. Shape: {data[key].shape}")
        except Exception as e:
            print(f"DEBUG: PlotFigures - Error loading {grn_cache_fname}: {e}")

    # CSV Fallbacks: Define mapping from CSV filename to the desired key in `data` dict
    # Ensure keys here match what downstream plotting functions expect.
    csv_data_map = {
        'ca_distance_table_grn.csv': 'ca_distance_table',
        'distance_table_grn.csv': 'distance_table',  # Use 'distance_table' to match cache key
        'msa_table_grn.csv': 'msa_table',  # Correct source for MSA data
        'ca_msa_table_grn.csv': 'ca_msa_table',  # Correct source for CA-specific MSA data
        'residue_table_grn_complete.csv': 'residue_table',  # Use the file that exists for residue info
    }

    print(f"\nDEBUG: PlotFigures - Checking CSV files in {output_dir_ref}")
    try:
        if output_dir_ref.exists() and output_dir_ref.is_dir():
            files = list(output_dir_ref.glob('*.csv'))
            print(f"DEBUG: PlotFigures - Found {len(files)} CSV files: {[f.name for f in files]}")
        else:
            print(
                f"DEBUG: PlotFigures - CSV reference directory {output_dir_ref} does not exist or is not a directory!")
    except Exception as e_glob:
        print(f"DEBUG: PlotFigures - Error listing CSV files: {e_glob}")

    for csv_file, data_key in csv_data_map.items():
        # Load from CSV if key not in data, or if it's an empty DataFrame from cache
        if data_key not in data or (isinstance(data.get(data_key), pd.DataFrame) and data[data_key].empty):
            csv_path = output_dir_ref / csv_file
            if csv_path.exists():
                try:
                    # Attempt to load, ensuring index_col=0 for DataFrames that need it
                    loaded_df_csv = pd.read_csv(csv_path, index_col=0)
                    data[data_key] = loaded_df_csv
                    print(f"DEBUG: PlotFigures - Loaded {csv_file} into '{data_key}'. Shape: {data[data_key].shape}")
                except Exception as e:
                    print(f"DEBUG: PlotFigures - Error loading {csv_file}: {e}")
            else:
                print(f"DEBUG: PlotFigures - CSV file for '{data_key}' not found: {csv_path}")
        else:
            print(
                f"  DEBUG: PlotFigures - Data for '{data_key}' already loaded (likely from cache). Skipping CSV load for {csv_file}.")

    errors_cache_fname = f'structure_errors_{chain_id}.pkl'
    cache_path_errors = cache_dir / errors_cache_fname
    if cache_path_errors.exists():
        try:
            with open(cache_path_errors, 'rb') as f:
                errors_data_cache = pickle.load(f)
                if isinstance(errors_data_cache, dict) and 'structure_errors' in errors_data_cache:
                    data['structure_errors'] = errors_data_cache['structure_errors']
                elif isinstance(errors_data_cache, dict):
                    data['structure_errors'] = errors_data_cache
            print(
                f"DEBUG: PlotFigures - Loaded {errors_cache_fname} ({len(data.get('structure_errors', {}))} entries).")
        except Exception as e:
            print(f"DEBUG: PlotFigures - Error loading {errors_cache_fname}: {e}")

    # --- Load Opsin Property Data ---
    property_file_path = Path("property/mo_exp.csv")  # Use Path for robustness
    if property_file_path.exists():
        print(f"\nDEBUG: PlotFigures - Loading property data from {property_file_path}")
        try:
            # Pass resolved input_dir if load_opsin_property_data needs it for file paths
            property_data_loaded = load_opsin_property_data(str(property_file_path),
                                                            data.get('processed_structures', {}))
            if property_data_loaded and 'properties' in property_data_loaded:
                data['property_data'] = property_data_loaded['properties']
                print(
                    f"DEBUG: PlotFigures - Loaded properties for {len(data['property_data'])} structures from property file.")
            else:
                print("DEBUG: PlotFigures - No 'properties' key found in data returned by load_opsin_property_data.")
        except Exception as e:
            print(f"DEBUG: PlotFigures - Error loading property data from {property_file_path}: {e}")
            import traceback
            traceback.print_exc()
    else:
        print(f"DEBUG: PlotFigures - Property file not found: {property_file_path}")

    # --- Stage 2: Consolidate Properties (Molecular Function, Domain, Error) ---
    sids_in_analysis = []
    if 'rmsd_df' in data and isinstance(data['rmsd_df'], pd.DataFrame) and not data['rmsd_df'].empty:
        sids_in_analysis = data['rmsd_df'].index.tolist()
        print(f"DEBUG: PlotFigures - Using {len(sids_in_analysis)} SIDs from rmsd_df for property consolidation.")
    elif 'processed_structures' in data and data['processed_structures']:
        sids_in_analysis = list(data['processed_structures'].keys())
        print(
            f"DEBUG: PlotFigures - Using {len(sids_in_analysis)} SIDs from processed_structures for property consolidation (rmsd_df was missing/empty).")
    else:
        print(
            "DEBUG: PlotFigures - CRITICAL: No SIDs for analysis (rmsd_df and processed_structures missing/empty). Property mapping will be empty.")
        data['group_dict'] = {}
        data['domain_dict_for_plots'] = {}
        data['overview_df'] = pd.DataFrame()
        return data

    final_group_dict = {}
    final_domain_dict_for_plots = {}
    overview_data_list_for_df = []
    processed_structures_map = data.get('processed_structures', {})
    property_data_map = data.get('property_data', {})  # From mo_exp.csv

    for sid in sids_in_analysis:
        mf, domain, avg_error, is_experimental, display_name = "Unknown", "Unknown", None, False, sid

        # Prioritize properties from property_data_map (mo_exp.csv) as it's often more curated
        # Try direct match for sid
        prop_entry = property_data_map.get(sid)

        # If direct match fails, try base ID for predicted structures
        if not prop_entry and ('_model_' in sid or '_smile_model_' in sid):
            base_id = sid.split('_model_')[0] if '_model_' in sid else sid.split('_smile_model_')[0]
            # Search for any key in property_data_map starting with base_id
            for p_sid, p_data_val in property_data_map.items():
                if p_sid.startswith(base_id):
                    prop_entry = p_data_val
                    # print(f"DEBUG: Matched property for predicted {sid} using base ID {base_id} (key: {p_sid})")
                    break

        if prop_entry and isinstance(prop_entry, dict):
            mf_val = prop_entry.get('molecular_function')
            if pd.notna(mf_val) and str(mf_val).strip() and str(mf_val).lower() != "unknown": mf = str(mf_val)

            domain_val = prop_entry.get('domain')
            if pd.notna(domain_val) and str(domain_val).strip() and str(domain_val).lower() != "unknown": domain = str(
                domain_val)

            is_experimental_prop = prop_entry.get('experimentally_determined', False)
            # Convert potential string "True"/"False" to boolean
            if isinstance(is_experimental_prop, str):
                is_experimental = is_experimental_prop.lower() == 'true'
            elif isinstance(is_experimental_prop, bool):
                is_experimental = is_experimental_prop

            display_name_prop = prop_entry.get('short_name', prop_entry.get('display_name'))  # Check both
            if pd.notna(display_name_prop) and str(display_name_prop).strip(): display_name = str(display_name_prop)

        # Fallback or supplement with processed_structures data if properties still Unknown
        struct_info = processed_structures_map.get(sid, {})
        if isinstance(struct_info, dict):
            # Update MF/Domain if still Unknown and available in struct_info.properties
            struct_props = struct_info.get('properties', {})
            if isinstance(struct_props, dict):
                if mf == "Unknown":
                    mf_val_struct = struct_props.get('molecular_function')
                    if pd.notna(mf_val_struct) and str(mf_val_struct).strip() and str(
                        mf_val_struct).lower() != "unknown": mf = str(mf_val_struct)
                if domain == "Unknown":
                    domain_val_struct = struct_props.get('domain')
                    if pd.notna(domain_val_struct) and str(domain_val_struct).strip() and str(
                        domain_val_struct).lower() != "unknown": domain = str(domain_val_struct)

                # is_experimental from struct_info.properties if not set by property_data_map
                if not prop_entry:  # Only if prop_entry was not found/used
                    is_experimental_struct = struct_props.get('is_experimental', False)
                    if isinstance(is_experimental_struct, str):
                        is_experimental = is_experimental_struct.lower() == 'true'
                    elif isinstance(is_experimental_struct, bool):
                        is_experimental = is_experimental_struct

                if display_name == sid:  # If display_name is still default PDB ID
                    display_name_struct = struct_props.get('display_name')
                    if pd.notna(display_name_struct) and str(display_name_struct).strip(): display_name = str(
                        display_name_struct)

            # Get average_error (prioritize structure_errors, then struct_info)
            if 'structure_errors' in data and isinstance(data['structure_errors'], dict) and sid in data[
                'structure_errors']:
                s_error_data = data['structure_errors'][sid]
                if isinstance(s_error_data, dict) and 'average_error' in s_error_data and pd.notna(
                        s_error_data['average_error']):
                    avg_error = float(s_error_data['average_error'])
            elif 'average_error' in struct_info and pd.notna(struct_info['average_error']):
                avg_error = float(struct_info['average_error'])

        final_group_dict[sid] = mf
        final_domain_dict_for_plots[sid] = {'domain': domain, 'average_error': avg_error}
        overview_data_list_for_df.append({
            'id': sid,  # Add original ID for easier joining/debugging
            'short_name': display_name, 'molecular_function_normalized': mf,
            'domain': domain, 'experimentally_determined': is_experimental
        })

    data['group_dict'] = final_group_dict
    data['domain_dict_for_plots'] = final_domain_dict_for_plots
    data['overview_df'] = pd.DataFrame(overview_data_list_for_df)

    unknown_mf_final = sum(1 for v in data['group_dict'].values() if v == "Unknown")
    unknown_dom_final = sum(
        1 for v_dict in data['domain_dict_for_plots'].values() if v_dict.get('domain', "Unknown") == "Unknown")
    print(f"DEBUG: PlotFigures - Final group_dict (MF): {len(data['group_dict'])} entries, {unknown_mf_final} Unknown.")
    print(
        f"DEBUG: PlotFigures - Final domain_dict_for_plots (Domain): {len(data['domain_dict_for_plots'])} entries, {unknown_dom_final} Unknown.")
    if not data['overview_df'].empty:
        print(f"DEBUG: PlotFigures - overview_df created with {len(data['overview_df'])} entries.")
        # print(data['overview_df'].head()) # For debugging overview_df content
    else:
        print("DEBUG: PlotFigures - overview_df is empty.")

    # Fallbacks for msa_table and residue_table (ensure msa_table is actual MSA data)
    if 'msa_table' not in data or (isinstance(data['msa_table'], pd.DataFrame) and data['msa_table'].empty):
        if data.get('ca_msa_table') is not None and not data['ca_msa_table'].empty:
            data['msa_table'] = data['ca_msa_table']  # Use ca_msa_table if msa_table from msa_table_grn.csv failed
            print("  DEBUG: PlotFigures - Used 'ca_msa_table' as 'msa_table'.")
        else:
            print("  DEBUG: PlotFigures - CRITICAL: 'msa_table' could not be populated with MSA data.")
            # Potentially create an empty DataFrame to prevent crashes, but plots will be empty/wrong
            data['msa_table'] = pd.DataFrame()

    if 'residue_table' not in data or (isinstance(data['residue_table'], pd.DataFrame) and data['residue_table'].empty):
        # If 'residue_table' from residue_table_grn_complete.csv wasn't loaded,
        # use the (hopefully correctly loaded) msa_table as fallback,
        # as plotting functions (logo, conservation) can extract AAs from MSA format.
        if data.get('msa_table') is not None and not data['msa_table'].empty:
            data['residue_table'] = data['msa_table']
            print("  DEBUG: PlotFigures - Used correctly loaded 'msa_table' as 'residue_table' (fallback).")
        else:
            print("  DEBUG: PlotFigures - 'residue_table' could not be populated as 'msa_table' is also missing/empty.")
            data['residue_table'] = pd.DataFrame()

    print("DEBUG: PlotFigures - Data loading and property consolidation complete.")
    return data


def generate_summary_csv(data, output_dir):
    summary_data = []
    if 'rmsd_df' not in data or not isinstance(data['rmsd_df'], pd.DataFrame) or data['rmsd_df'].empty:
        print("DEBUG: Skipping summary CSV: 'rmsd_df' not available or empty.")
        return None

    rmsd_df = data['rmsd_df']
    domain_dict_for_plots = data.get('domain_dict_for_plots', {})
    group_dict_for_summary = data.get('group_dict', {})

    # Use overview_df for display names, matching by 'id' column
    overview_df_for_names = data.get('overview_df', pd.DataFrame())
    if not overview_df_for_names.empty and 'id' in overview_df_for_names.columns:
        name_map = pd.Series(overview_df_for_names.short_name.values, index=overview_df_for_names.id).to_dict()
    else:
        name_map = {}

    print("DEBUG: --- Generating Summary CSV ---")
    sids_in_rmsd = rmsd_df.index.tolist()

    for protein_id in sids_in_rmsd:
        rmsd_values = rmsd_df.loc[protein_id].values
        # Exclude diagonal (self-RMSD = 0) and NaNs for average calculation
        valid_rmsd_values = rmsd_values[~np.isnan(rmsd_values) & (np.abs(rmsd_values) > 1e-9)]  # Small epsilon for zero
        avg_rmsd_val = np.mean(valid_rmsd_values) if len(valid_rmsd_values) > 0 else np.nan

        display_name = name_map.get(protein_id, protein_id)  # Get from map, fallback to ID

        domain_plot_info = domain_dict_for_plots.get(protein_id, {'domain': "Unknown", 'average_error': np.nan})
        domain = domain_plot_info.get('domain', "Unknown")  # Use .get for safety
        error_val = domain_plot_info.get('average_error', np.nan)
        function = group_dict_for_summary.get(protein_id, "Unknown")

        summary_data.append([
            display_name,
            round(avg_rmsd_val, 2) if pd.notna(avg_rmsd_val) else 'N/A',
            domain, function,
            round(error_val, 2) if pd.notna(error_val) else 'N/A'
        ])

    summary_df = pd.DataFrame(summary_data,
                              columns=['Protein', 'Average RMSD', 'Domain', 'Molecular Function', 'Error'])
    # Sort, ensuring Unknown comes last for string columns if desired (pandas default is usually fine)
    summary_df = summary_df.sort_values(by=['Molecular Function', 'Domain', 'Average RMSD'],
                                        na_position='last',
                                        key=lambda x: x.replace('Unknown', 'ZZZUnknown') if x.name in [
                                            'Molecular Function', 'Domain'] else x)

    csv_path = Path(output_dir) / 'protein_summary.csv'
    summary_df.to_csv(csv_path, index=False)
    print(f"Saved protein summary to {csv_path}")

    all_functions = sorted(
        list(set(s for s in summary_df['Molecular Function'].unique() if
                 s not in ["Unknown", "ZZZUnknown"] and pd.notna(s))))
    all_domains = sorted(
        list(set(s for s in summary_df['Domain'].unique() if s not in ["Unknown", "ZZZUnknown"] and pd.notna(s))))

    if COLOR_SCHEME_LOADED:
        function_colors = get_group_colors(all_functions)
        domain_colors = get_group_colors(all_domains)
        color_data = []
        for func, color in function_colors.items(): color_data.append(['Function', func, color])
        for dom, color in domain_colors.items(): color_data.append(['Domain', dom, color])
        if color_data:
            color_df = pd.DataFrame(color_data, columns=['Type', 'Value', 'Color'])
            color_path = Path(output_dir) / 'color_assignments.csv'
            color_df.to_csv(color_path, index=False)
            print(f"Saved color assignments to {color_path}")
    else:
        print("DEBUG: Color scheme not loaded, skipping color_assignments.csv.")
    return csv_path


def generate_plots(data, output_dir_figures, error_threshold=3.0):  # Removed input_dir from args
    """
    Generate all curated plots and save them to the output directory.
    """
    if not VIS_FUNCS_LOADED:
        print("CRITICAL DEBUG: Visualization functions not loaded. Most plots will be skipped or dummied.")
        # return # Optionally exit early

    print(f"\n--- Generating plots in {output_dir_figures} ---")
    Path(output_dir_figures).mkdir(parents=True, exist_ok=True)

    # --- A. Overview & Reference ---
    print("\nDEBUG: 1. Generating Opsins Overview Plot...")
    if 'overview_df' in data and isinstance(data['overview_df'], pd.DataFrame) and not data['overview_df'].empty:
        try:
            fig_overview = create_opsin_overview_plot(data['overview_df'])
            if fig_overview:
                fig_overview.savefig(Path(output_dir_figures) / 'A1_opsin_overview.png', dpi=300, bbox_inches='tight')
                plt.close(fig_overview)
                print("  DEBUG: Saved A1_opsin_overview.png")
        except Exception as e:
            print(f"  DEBUG: Error generating opsin overview plot: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("  DEBUG: Skipped Opsins Overview: 'overview_df' not valid or empty.")

    print("\nDEBUG: 2. Generating RMSD Color Scale Reference...")
    try:
        fig_rmsd_scale = create_rmsd_color_scale_figure()
        if fig_rmsd_scale:
            fig_rmsd_scale.savefig(Path(output_dir_figures) / 'A2_rmsd_color_scale.png', dpi=300, bbox_inches='tight')
            plt.close(fig_rmsd_scale)
            print("  DEBUG: Saved A2_rmsd_color_scale.png")
    except Exception as e:
        print(f"  DEBUG: Error generating RMSD color scale: {e}")
        import traceback
        traceback.print_exc()

    if 'rmsd_df' not in data or not isinstance(data['rmsd_df'], pd.DataFrame) or data['rmsd_df'].empty:
        print("\nDEBUG: CRITICAL: 'rmsd_df' is missing or empty. Skipping RMSD-based plots (B, C).")
    else:
        rmsd_df_orig = data['rmsd_df']
        # --- C. RMSD & Structural Similarity (Simple/Unfiltered) ---
        print("\nDEBUG: 3. Generating Simple Unfiltered Structural Similarity Tree...")
        try:
            fig_simple_tree = plot_similarity_tree(rmsd_df_orig,
                                                   title="Unfiltered Structural Similarity Tree")
            if fig_simple_tree:
                fig_simple_tree.savefig(Path(output_dir_figures) / 'C3_similarity_tree_unfiltered.png', dpi=300,
                                        bbox_inches='tight')
                plt.close(fig_simple_tree)
                print("  DEBUG: Saved C3_similarity_tree_unfiltered.png")
        except Exception as e:
            print(f"  DEBUG: Error generating simple unfiltered similarity tree: {e}")
            import traceback
            traceback.print_exc()

        # --- B. RMSD & Structural Similarity (Linked & Filtered) ---
        print(f"\nDEBUG: --- Preparing for Linked & Filtered RMSD Plots (Threshold: {error_threshold} Å) ---")
        pdb_list_orig = rmsd_df_orig.index.tolist()
        group_dict_for_plots = data.get('group_dict', {})
        domain_dict_for_plots_with_error = data.get('domain_dict_for_plots', {})

        kept_ids = []
        # avg_rmsds_from_matrix = rmsd_df_orig.mean(axis=1) # Fallback error metric not ideal

        for pdb_id in pdb_list_orig:
            error_val = np.nan  # Default to NaN if no error info
            domain_info = domain_dict_for_plots_with_error.get(pdb_id, {})
            if isinstance(domain_info, dict) and 'average_error' in domain_info and pd.notna(
                    domain_info['average_error']):
                error_val = domain_info['average_error']

            # Keep if error is below threshold OR if error_val is NaN (meaning no error data, assume keep)
            if pd.isna(error_val) or error_val <= error_threshold:
                kept_ids.append(pdb_id)
            else:
                print(f"  DEBUG: Filtering out {pdb_id} due to error {error_val:.2f} > {error_threshold:.2f}")

        print(f"  DEBUG: Filtered to {len(kept_ids)} structures for linked plots.")

        if len(kept_ids) >= 2:
            filtered_rmsd_df = rmsd_df_orig.loc[kept_ids, kept_ids].copy()
            # Ensure dicts are also filtered for the kept_ids
            filtered_group_dict = {k: v for k, v in group_dict_for_plots.items() if k in kept_ids}
            filtered_domain_dict = {k: v for k, v in domain_dict_for_plots_with_error.items() if k in kept_ids}

            Z_linkage = None
            try:
                # Robust Z_linkage calculation
                matrix_for_linkage = filtered_rmsd_df.copy()
                # Fill NaN with a large value (or mean/median of non-NaNs) before squareform
                if matrix_for_linkage.isnull().values.any():
                    fill_val = np.nanmax(matrix_for_linkage.values) * 2 if pd.notna(
                        np.nanmax(matrix_for_linkage.values)) else 10.0
                    matrix_for_linkage.fillna(fill_val, inplace=True)

                np.fill_diagonal(matrix_for_linkage.values, 0.0)  # Ensure diagonal is zero
                condensed_matrix = squareform(matrix_for_linkage.values, checks=False)
                Z_linkage = linkage(condensed_matrix, method='average')
            except Exception as e_link:
                print(f"    DEBUG: ERROR calculating Z_linkage: {e_link}")
                import traceback;
                traceback.print_exc()

            if Z_linkage is not None:
                print("\nDEBUG: 4. Generating Filtered Structural Similarity Tree (Linked)...")
                try:
                    tree_fig, ordered_tree_ids = create_and_visualize_similarity_tree(
                        rmsd_data=filtered_rmsd_df, linkage_matrix=Z_linkage,
                        group_dict=filtered_group_dict, domain_dict=filtered_domain_dict
                    )
                    if tree_fig:
                        tree_fig.savefig(Path(output_dir_figures) / 'B4_similarity_tree_filtered_linked.png', dpi=300,
                                         bbox_inches='tight')
                        plt.close(tree_fig)
                        print("  DEBUG: Saved B4_similarity_tree_filtered_linked.png")
                except Exception as e_tree:
                    print(f"  DEBUG: Error generating filtered similarity tree: {e_tree}")
                    import traceback;
                    traceback.print_exc()

                print("\nDEBUG: 5. Generating Filtered RMSD Heatmap (Linked & Clustered)...")
                try:
                    heatmap_clustermap = visualize_rmsd_matrix_improved(
                        rmsd_df=filtered_rmsd_df, linkage_matrix=Z_linkage,
                        group_dict=filtered_group_dict, domain_dict=filtered_domain_dict,
                        output_file=Path(output_dir_figures) / 'B5_rmsd_heatmap_filtered_linked.png'
                    )
                    if heatmap_clustermap and hasattr(heatmap_clustermap, 'fig'): plt.close(heatmap_clustermap.fig)
                    # print("  DEBUG: Saved B5_rmsd_heatmap_filtered_linked.png (by visualize_rmsd_matrix_improved)") # Saved internally
                except Exception as e_heatmap:
                    print(f"  DEBUG: Error generating filtered RMSD heatmap: {e_heatmap}")
                    import traceback;
                    traceback.print_exc()

        else:
            print(
                "  DEBUG: Skipped Linked Tree and Heatmap (not enough structures after filtering or Z_linkage failed).")

    # --- D. Distance-to-Retinal Plots ---
    print("\nDEBUG: Generating Distance to Retinal with STD (Plot D6 - All-Atom)...")
    distance_table_for_plot_d6 = None
    if 'distance_table' in data and isinstance(data['distance_table'], pd.DataFrame) and not data[
        'distance_table'].empty:
        distance_table_for_plot_d6 = data['distance_table']
        print("  DEBUG: Using data['distance_table'] for D6 plot.")

    if distance_table_for_plot_d6 is not None:
        try:
            fig_d6 = plot_distances_with_std(distance_table_for_plot_d6,
                                             title="All-Atom Distance to Retinal by Position",
                                             use_ca=False)
            if fig_d6:
                fig_d6.savefig(Path(output_dir_figures) / 'D6_all_atom_distance_retinal_std.png', dpi=300,
                               bbox_inches='tight')
                plt.close(fig_d6)
                print("  DEBUG: Saved D6_all_atom_distance_retinal_std.png")
        except Exception as e:
            print(f"  DEBUG: Error generating D6 (All-Atom Distance STD) plot: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("  DEBUG: Skipped D6 (All-Atom Distance STD) plot: 'distance_table' not found or empty in loaded_data.")

    print("\nDEBUG: Generating CA Atom Distance to Retinal with STD (Plot D_CA)...")
    ca_distance_table_for_plot = None
    if 'ca_distance_table' in data and isinstance(data['ca_distance_table'], pd.DataFrame) and not data[
        'ca_distance_table'].empty:
        ca_distance_table_for_plot = data['ca_distance_table']
        print("  DEBUG: Using data['ca_distance_table'] for CA Distance STD plot.")

    if ca_distance_table_for_plot is not None:
        try:
            fig_ca = plot_distances_with_std(ca_distance_table_for_plot,
                                             title="CA Atom Distance to Retinal by Position",
                                             use_ca=True)
            if fig_ca:
                fig_ca.savefig(Path(output_dir_figures) / 'D_ca_atom_distance_retinal_std.png', dpi=300,
                               bbox_inches='tight')
                plt.close(fig_ca)
                print("  DEBUG: Saved D_ca_atom_distance_retinal_std.png")
        except Exception as e:
            print(f"  DEBUG: Error generating CA Atom Distance STD plot: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("  DEBUG: Skipped CA Atom Distance STD plot: 'ca_distance_table' not found or empty.")

    # --- E. Sequence Conservation & Composition ---
    print("\nDEBUG: Generating Residue Conservation Bar Plot (Plot E12)...")
    if 'msa_table' in data and isinstance(data['msa_table'], pd.DataFrame) and not data['msa_table'].empty:
        try:
            fig_e12 = create_residue_conservation_plot(data['msa_table'], helix_highlighting=True)
            if fig_e12:
                fig_e12.savefig(Path(output_dir_figures) / 'E12_residue_conservation.png', dpi=300, bbox_inches='tight')
                plt.close(fig_e12)
                print("  DEBUG: Saved E12_residue_conservation.png")
        except Exception as e:
            print(f"  DEBUG: Error generating residue conservation plot (E12): {e}")
            import traceback
            traceback.print_exc()
    else:
        print("  DEBUG: Skipped Residue Conservation plot (E12): 'msa_table' invalid or empty.")

    print("\nDEBUG: Generating Helix Logo Plots around X.50 (Plot E_LogoX50)...")
    # plot_helix_logo_plots expects 'residue_table' to be MSA-like
    residue_table_for_logo = data.get('residue_table', data.get('msa_table'))  # Fallback to msa_table
    if residue_table_for_logo is not None and isinstance(residue_table_for_logo,
                                                         pd.DataFrame) and not residue_table_for_logo.empty:
        try:
            fig_logo_x50 = plot_helix_logo_plots(residue_table_for_logo)
            if fig_logo_x50:
                fig_logo_x50.savefig(Path(output_dir_figures) / 'E_helix_logos_x50.png', dpi=300, bbox_inches='tight')
                plt.close(fig_logo_x50)
                print("  DEBUG: Saved E_helix_logos_x50.png")
        except Exception as e:
            print(f"  DEBUG: Error generating Helix Logo plots (X.50): {e}")
            import traceback;
            traceback.print_exc()
    else:
        print("  DEBUG: Skipped Helix Logo plots (X.50): 'residue_table' (or msa_table fallback) invalid or empty.")

    # --- F. Combined Plots (Optional) ---
    print("\nDEBUG: Generating Combined Distance/Logo Plot (Plot F_Combined)...")
    dist_table_for_combined = data.get('distance_table')  # All-atom distance
    msa_table_for_combined = data.get('msa_table')
    if dist_table_for_combined is not None and not dist_table_for_combined.empty and \
            msa_table_for_combined is not None and not msa_table_for_combined.empty:
        try:
            fig_combined = create_combined_distance_logo_plot(dist_table_for_combined, msa_table_for_combined)
            if fig_combined:
                fig_combined.savefig(Path(output_dir_figures) / 'F_combined_distance_logo.png', dpi=300,
                                     bbox_inches='tight')
                plt.close(fig_combined)
                print("  DEBUG: Saved F_combined_distance_logo.png")
        except Exception as e:
            print(f"  DEBUG: Error generating Combined Distance/Logo plot: {e}")
            import traceback;
            traceback.print_exc()
    else:
        print("  DEBUG: Skipped Combined Distance/Logo plot: Missing 'distance_table' or 'msa_table'.")

    print(f"\n--- All plot generation attempts complete. Check {output_dir_figures}. ---")


def verify_filtering_in_summary(summary_csv_path, error_threshold=3.0):
    # This function was not part of the error, so it remains as is.
    pass


def main():
    parser = argparse.ArgumentParser(description='Generate plots for opsin analysis.')
    parser.add_argument('--input-dir', '-i', type=str, default='opsin_output/', help='Directory for input data.')
    parser.add_argument('--output-dir', '-o', type=str, default='opsin_output/',
                        help='Directory for output plots.')
    parser.add_argument('--chain-id', '-c', type=str, default='A', help='Chain ID (default: A)')
    parser.add_argument('--quality', '-q', type=str, choices=['low', 'medium', 'high'], default='medium',
                        help='Figure quality (default: medium)')
    parser.add_argument('--error-threshold', '-et', type=float, default=3.0,
                        help='Error threshold for filtering (default: 3.0)')
    args = parser.parse_args()

    output_dir_figures = Path(args.output_dir) / 'figures'
    output_dir_summaries = Path(args.output_dir) / 'summaries'
    output_dir_figures.mkdir(parents=True, exist_ok=True)
    output_dir_summaries.mkdir(parents=True, exist_ok=True)

    print(f"Opsin Plot Generation Script - DEBUG ENABLED")
    print(f"Input Dir: {Path(args.input_dir).resolve()}")
    print(f"Output Plots Dir: {output_dir_figures.resolve()}")
    print(f"Output Summaries Dir: {output_dir_summaries.resolve()}")

    dpi_map = {'low': 100, 'medium': 200, 'high': 300}
    plt.rcParams['figure.dpi'] = dpi_map[args.quality]
    plt.rcParams['savefig.dpi'] = dpi_map[args.quality]

    # Pass args.input_dir as both input_dir_path_str and output_dir_ref_str
    # if CSV fallbacks are expected in the main input directory.
    loaded_data = load_data(args.input_dir, args.input_dir, args.chain_id)

    if not loaded_data or not isinstance(loaded_data, dict):
        print("CRITICAL DEBUG: No data loaded or data is not a dictionary. Exiting plot generation.")
        return

    debug_property_data(loaded_data)  # Keep this for insights

    generate_plots(loaded_data, output_dir_figures, error_threshold=args.error_threshold)
    summary_csv_path = generate_summary_csv(loaded_data, output_dir_summaries)
    # if summary_csv_path:
    #     verify_filtering_in_summary(summary_csv_path, error_threshold=args.error_threshold)

    print(f"\nAll operations complete. Outputs are in {Path(args.output_dir).resolve()}")


def debug_property_data(loaded_data):
    """
    Debug the property data to identify missing data patterns
    """
    print("\n=== DEBUG: PROPERTY DATA ANALYSIS ===")

    group_dict = loaded_data.get('group_dict', {})
    domain_dict_for_plots = loaded_data.get('domain_dict_for_plots', {})  # Corrected key
    property_data_map = loaded_data.get('property_data', {})  # Renamed for clarity

    unknown_mf_sids = [sid for sid, val in group_dict.items() if val == "Unknown"]
    unknown_domain_sids = [sid for sid, d_info in domain_dict_for_plots.items() if
                           d_info.get('domain') == "Unknown"]  # Use .get

    print(f"Total structures with properties in property_data_map: {len(property_data_map)} structures")
    print(
        f"Structures with 'Unknown' molecular function after consolidation: {len(unknown_mf_sids)}/{len(group_dict)} structures")
    print(
        f"Structures with 'Unknown' domain after consolidation: {len(unknown_domain_sids)}/{len(domain_dict_for_plots)} structures")

    if property_data_map and unknown_mf_sids:
        exp_missing_mf_in_prop_map = 0
        pred_missing_mf_in_prop_map = 0  # Structures that *are* in property_data_map but MF is still Unknown

        actually_missing_from_prop_map = 0  # SIDs not found in property_data_map at all

        for sid in unknown_mf_sids:
            if sid in property_data_map:
                is_exp = False  # Default
                exp_status = property_data_map[sid].get('experimentally_determined')
                if isinstance(exp_status, str):
                    is_exp = exp_status.lower() == 'true'
                elif isinstance(exp_status, bool):
                    is_exp = exp_status

                if property_data_map[sid].get('molecular_function', "Unknown").lower() == "unknown":
                    if is_exp:
                        exp_missing_mf_in_prop_map += 1
                    else:
                        pred_missing_mf_in_prop_map += 1
            else:  # SID itself not found in property_data_map
                actually_missing_from_prop_map += 1

        print(f"\nAnalysis of {len(unknown_mf_sids)} structures with final 'Unknown' MF:")
        print(
            f"  - SIDs where MF is 'Unknown' *within* their property_data_map entry (Experimental): {exp_missing_mf_in_prop_map}")
        print(
            f"  - SIDs where MF is 'Unknown' *within* their property_data_map entry (Predicted/Other): {pred_missing_mf_in_prop_map}")
        print(
            f"  - SIDs completely missing from property_data_map (thus MF couldn't be sourced there): {actually_missing_from_prop_map}")

        if unknown_mf_sids:
            print("\nExample structures with final 'Unknown' MF (max 5 shown):")
            for sid in unknown_mf_sids[:5]:
                prop_entry_debug = property_data_map.get(sid)
                struct_info_debug = loaded_data.get('processed_structures', {}).get(sid, {}).get('properties', {})
                print(f"  - SID: {sid}")
                print(f"    Property_data_map entry: {prop_entry_debug}")
                print(f"    Processed_structures.properties entry: {struct_info_debug}")
                print(f"    Final Group Dict MF: {group_dict.get(sid)}")
                print(f"    Final Domain Dict Info: {domain_dict_for_plots.get(sid)}")


if __name__ == "__main__":
    main()