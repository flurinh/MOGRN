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
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import json  # For summary saving

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
        from src.visualization_functions import (  # type: ignore
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


        # Define dummy functions
        def _dummy_plot_func(*args, **kwargs):
            fn_name = inspect.stack()[1].function if len(inspect.stack()) > 1 else "UnknownFunction"
            print(
                f"Warning: Plot function '{fn_name}' called but not loaded due to import error. Args: {args}, Kwargs: {kwargs}")
            fig, ax = plt.subplots()
            ax.text(0.5, 0.5, f"Plotting function\n{fn_name}\nunavailable", ha='center', va='center', color='red')
            return fig, None if fn_name == "create_and_visualize_similarity_tree" else fig  # Adjust for return tuple


        import inspect

        create_opsin_overview_plot = create_rmsd_color_scale_figure = _dummy_plot_func
        create_and_visualize_similarity_tree = _dummy_plot_func  # Note: this one returns a tuple
        visualize_rmsd_matrix_improved = _dummy_plot_func
        plot_similarity_tree = plot_rmsd_heatmap = _dummy_plot_func
        plot_distances_with_std = plot_average_distances_by_helix = plot_distance_heatmap = _dummy_plot_func
        create_residue_conservation_plot = plot_helix_logo_plots = plot_conservation_around_x50 = _dummy_plot_func
        create_combined_distance_logo_plot = visualize_binding_pocket = _dummy_plot_func

# Import the color scheme tools
try:
    from src.opsin_color_scheme import get_group_colors, RMSD_COMPACT_CMAP  # Import specific items needed

    COLOR_SCHEME_LOADED = True
except ImportError:
    print("Error: opsin_color_scheme.py not found. Using fallback colors.")
    COLOR_SCHEME_LOADED = False


    def get_group_colors(items_list_or_dict, palette_name=None):
        if isinstance(items_list_or_dict, dict):
            items = list(items_list_or_dict.keys())
        else:
            items = list(items_list_or_dict)
        num_items = len(items)
        # Generate a simple list of distinct colors for fallback
        fallback_palette = plt.cm.get_cmap('viridis', num_items if num_items > 0 else 1)
        return {item: fallback_palette(i) for i, item in enumerate(items)}


    RMSD_COMPACT_CMAP = "viridis"  # Fallback colormap


def load_data(input_dir, output_dir_ref, chain_id='A'):
    """
    Load precomputed data from workflow cache or CSV files, and consolidate properties.
    """
    data = {}
    cache_dir = Path(input_dir) / 'cache'
    print(f"DEBUG: PlotFigures - Loading data. Cache: {cache_dir}, CSV Ref: {Path(output_dir_ref).resolve()}")

    # --- Stage 1: Load Core Data from Cache (Primary) & CSVs (Fallback) ---
    # Load processed_structures.pkl
    proc_struct_fname = f'processed_structures_{chain_id}.pkl'
    cache_path_proc_struct = cache_dir / proc_struct_fname
    if cache_path_proc_struct.exists():
        try:
            with open(cache_path_proc_struct, 'rb') as f:
                loaded_proc_data = pickle.load(f)
                if isinstance(loaded_proc_data, dict) and 'processed_structures' in loaded_proc_data and isinstance(
                        loaded_proc_data['processed_structures'], dict):
                    data['processed_structures'] = loaded_proc_data['processed_structures']
                elif isinstance(loaded_proc_data, dict):
                    data['processed_structures'] = loaded_proc_data
                else:
                    print(f"  DEBUG: Unexpected data type in {proc_struct_fname}: {type(loaded_proc_data)}")
            print(
                f"DEBUG: PlotFigures - Loaded {proc_struct_fname} ({len(data.get('processed_structures', {}))} entries)")
        except Exception as e:
            print(f"DEBUG: PlotFigures - Error loading {proc_struct_fname}: {e}")
    else:
        print(f"DEBUG: PlotFigures - Cache file {proc_struct_fname} not found at {cache_path_proc_struct}.")

    # Load RMSD data (structure_comparison.pkl or rmsd_matrix.csv)
    rmsd_cache_fname = f'structure_comparison_{chain_id}.pkl'
    cache_path_rmsd = cache_dir / rmsd_cache_fname
    if cache_path_rmsd.exists():
        try:
            with open(cache_path_rmsd, 'rb') as f:
                comp_data = pickle.load(f)
                if 'rmsd_df' in comp_data and isinstance(comp_data['rmsd_df'], pd.DataFrame):
                    data['rmsd_df'] = comp_data['rmsd_df']
                elif 'rmsd_matrix' in comp_data and isinstance(comp_data['rmsd_matrix'], pd.DataFrame):
                    data['rmsd_df'] = comp_data['rmsd_matrix']
                # Load other useful items from comparison data if needed by plots directly
                if 'group_dict' in comp_data: data['group_dict_from_comparison_cache'] = comp_data['group_dict']
                if 'domain_dict' in comp_data: data['domain_dict_from_comparison_cache'] = comp_data['domain_dict']
                if 'pdb_list' in comp_data: data['pdb_list_from_comparison_cache'] = comp_data['pdb_list']

            print(
                f"DEBUG: PlotFigures - Loaded {rmsd_cache_fname}. RMSD_DF shape: {data.get('rmsd_df', pd.DataFrame()).shape}")
        except Exception as e:
            print(f"DEBUG: PlotFigures - Error loading {rmsd_cache_fname}: {e}")

    if 'rmsd_df' not in data:  # Fallback to CSV for RMSD
        csv_path_rmsd = Path(output_dir_ref) / 'rmsd_matrix.csv'
        if csv_path_rmsd.exists():
            try:
                data['rmsd_df'] = pd.read_csv(csv_path_rmsd, index_col=0)
                print(f"DEBUG: PlotFigures - Loaded rmsd_matrix.csv. Shape: {data['rmsd_df'].shape}")
            except Exception as e:
                print(f"DEBUG: PlotFigures - Error loading rmsd_matrix.csv: {e}")
        else:
            print("DEBUG: PlotFigures - CRITICAL: RMSD matrix not found in cache or CSV. Many plots might fail.")

    # Load MSA and Distance tables (from grn_assignment.pkl or individual CSVs)
    grn_cache_fname = f'grn_assignment_{chain_id}.pkl'
    cache_path_grn = cache_dir / grn_cache_fname
    grn_data_loaded_from_cache = False
    if cache_path_grn.exists():
        try:
            with open(cache_path_grn, 'rb') as f:
                grn_cache_data = pickle.load(f)
                if isinstance(grn_cache_data, dict):
                    for key in ["residue_table", "distance_table", "ca_residue_table", "ca_distance_table",
                                "msa_table"]:  # msa_table might be an alias
                        if key in grn_cache_data and isinstance(grn_cache_data[key], pd.DataFrame):
                            data[key] = grn_cache_data[key]
                            print(
                                f"  DEBUG: PlotFigures - Loaded '{key}' from {grn_cache_fname}. Shape: {data[key].shape}")
                    grn_data_loaded_from_cache = True
        except Exception as e:
            print(f"DEBUG: PlotFigures - Error loading {grn_cache_fname}: {e}")

    # Fallback to individual CSVs for MSA/Distance tables if not fully loaded from GRN cache
    csv_data_map = {
        'opsin_grn_tables/ca_distance_table_grn.csv': 'ca_distance_table',
        'opsin_grn_tables/distance_table_grn.csv': 'distance_table',
        'opsin_grn_tables/ca_msa_table_grn.csv': 'msa_table',
        'opsin_grn_tables/residue_table_grn.csv': 'residue_table',
    }
    for csv_file, data_key in csv_data_map.items():
        if data_key not in data or (isinstance(data[data_key], pd.DataFrame) and data[data_key].empty):
            csv_path = Path(output_dir_ref) / csv_file
            if csv_path.exists():
                try:
                    data[data_key] = pd.read_csv(csv_path, index_col=0)
                    print(f"DEBUG: PlotFigures - Loaded {csv_file} into '{data_key}'. Shape: {data[data_key].shape}")
                except Exception as e:
                    print(f"DEBUG: PlotFigures - Error loading {csv_file}: {e}")
            else:
                print(f"DEBUG: PlotFigures - CSV file for '{data_key}' not found: {csv_path}")

    # Load structure_errors for error values
    errors_cache_fname = f'structure_errors_{chain_id}.pkl'
    cache_path_errors = cache_dir / errors_cache_fname
    if cache_path_errors.exists():
        try:
            with open(cache_path_errors, 'rb') as f:
                errors_data_cache = pickle.load(f)
                if isinstance(errors_data_cache,
                              dict) and 'structure_errors' in errors_data_cache:  # if pkl saves a dict {'structure_errors': errors_dict}
                    data['structure_errors'] = errors_data_cache['structure_errors']
                elif isinstance(errors_data_cache, dict):  # if pkl directly saves the errors_dict
                    data['structure_errors'] = errors_data_cache
            print(
                f"DEBUG: PlotFigures - Loaded {errors_cache_fname} ({len(data.get('structure_errors', {}))} entries).")
        except Exception as e:
            print(f"DEBUG: PlotFigures - Error loading {errors_cache_fname}: {e}")

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

    for sid in sids_in_analysis:
        mf, domain, avg_error, is_experimental, display_name = "Unknown", "Unknown", None, False, sid
        struct_info = processed_structures_map.get(sid, {})

        if isinstance(struct_info, dict):
            props = struct_info.get('properties', {})
            if isinstance(props, dict):
                mf_val = props.get('molecular_function')
                if pd.notna(mf_val) and str(mf_val).strip() and str(mf_val).lower() != "unknown": mf = str(mf_val)

                domain_val = props.get('domain')
                if pd.notna(domain_val) and str(domain_val).strip() and str(
                    domain_val).lower() != "unknown": domain = str(domain_val)

                is_experimental = props.get('is_experimental', False)
                display_name_prop = props.get('display_name')
                if pd.notna(display_name_prop) and str(display_name_prop).strip(): display_name = str(display_name_prop)

            # Get average_error
            if 'average_error' in struct_info and pd.notna(struct_info['average_error']):
                avg_error = float(struct_info['average_error'])
            elif 'structure_errors' in data and isinstance(data['structure_errors'], dict) and sid in data[
                'structure_errors']:
                s_error_data = data['structure_errors'][sid]
                if isinstance(s_error_data, dict) and 'average_error' in s_error_data and pd.notna(
                        s_error_data['average_error']):
                    avg_error = float(s_error_data['average_error'])

        final_group_dict[sid] = mf
        final_domain_dict_for_plots[sid] = {'domain': domain, 'average_error': avg_error}
        overview_data_list_for_df.append({
            'short_name': display_name, 'molecular_function_normalized': mf,
            'domain': domain, 'experimentally_determined': is_experimental
        })

    data['group_dict'] = final_group_dict
    data['domain_dict_for_plots'] = final_domain_dict_for_plots
    data['overview_df'] = pd.DataFrame(overview_data_list_for_df)

    unknown_mf_final = sum(1 for v in data['group_dict'].values() if v == "Unknown")
    unknown_dom_final = sum(1 for v_dict in data['domain_dict_for_plots'].values() if v_dict['domain'] == "Unknown")
    print(f"DEBUG: PlotFigures - Final group_dict (MF): {len(data['group_dict'])} entries, {unknown_mf_final} Unknown.")
    print(
        f"DEBUG: PlotFigures - Final domain_dict_for_plots (Domain): {len(data['domain_dict_for_plots'])} entries, {unknown_dom_final} Unknown.")
    if not data['overview_df'].empty:
        print(f"DEBUG: PlotFigures - overview_df created with {len(data['overview_df'])} entries.")
    else:
        print("DEBUG: PlotFigures - overview_df is empty.")

    # Fallbacks for msa_table and residue_table
    if 'msa_table' not in data or (isinstance(data['msa_table'], pd.DataFrame) and data['msa_table'].empty):
        if data.get('ca_msa_table_grn') is not None and not data['ca_msa_table_grn'].empty:
            data['msa_table'] = data['ca_msa_table_grn']
            print("  DEBUG: PlotFigures - Used 'ca_msa_table_grn' as 'msa_table'.")
        elif data.get('residue_table_grn') is not None and not data['residue_table_grn'].empty:  # less ideal fallback
            data['msa_table'] = data['residue_table_grn']
            print("  DEBUG: PlotFigures - Used 'residue_table_grn' as 'msa_table'.")

    if 'residue_table' not in data or (isinstance(data['residue_table'], pd.DataFrame) and data['residue_table'].empty):
        if data.get('ca_residue_table_grn') is not None and not data['ca_residue_table_grn'].empty:
            data['residue_table'] = data['ca_residue_table_grn']
            print("  DEBUG: PlotFigures - Used 'ca_residue_table_grn' as 'residue_table'.")
        elif data.get('residue_table_grn') is not None and not data['residue_table_grn'].empty:
            data['residue_table'] = data['residue_table_grn']
            print("  DEBUG: PlotFigures - Used 'residue_table_grn' as 'residue_table'.")
        elif data.get('msa_table') is not None and not data['msa_table'].empty:  # Fallback to msa_table
            data['residue_table'] = data['msa_table']
            print("  DEBUG: PlotFigures - Used 'msa_table' as 'residue_table'.")

    print("DEBUG: PlotFigures - Data loading and property consolidation complete.")
    return data


def generate_summary_csv(data, output_dir):
    summary_data = []
    if 'rmsd_df' not in data or not isinstance(data['rmsd_df'], pd.DataFrame) or data['rmsd_df'].empty:
        print("DEBUG: Skipping summary CSV: 'rmsd_df' not available or empty.")
        return None

    rmsd_df = data['rmsd_df']
    # Use the consolidated property dicts
    domain_dict_for_plots = data.get('domain_dict_for_plots', {})
    group_dict_for_summary = data.get('group_dict', {})
    overview_df_for_names = data.get('overview_df', pd.DataFrame(
        columns=['short_name', 'molecular_function_normalized']))  # For display names

    print("DEBUG: --- Generating Summary CSV ---")
    sids_in_rmsd = rmsd_df.index.tolist()

    for protein_id in sids_in_rmsd:
        rmsd_values = rmsd_df.loc[protein_id].values
        valid_rmsd_values = rmsd_values[~np.isnan(rmsd_values) & (np.abs(rmsd_values) > 1e-6)]
        avg_rmsd_val = np.mean(valid_rmsd_values) if len(valid_rmsd_values) > 0 else np.nan

        # Get display name from overview_df if possible
        display_name = protein_id  # Default
        if not overview_df_for_names.empty and 'short_name' in overview_df_for_names.columns:
            # This assumes sids_in_rmsd can be found by matching overview_df 'short_name' if it was originally the PDB ID
            # or by matching a PDB ID column if overview_df had one.
            # A more robust way is if overview_df was indexed by PDB ID.
            # For now, let's try matching, but this might need adjustment based on overview_df's actual index/columns.
            # A common case: overview_df short_name is display_name, we need mapping from PDB ID to its row.
            # Assuming overview_df's 'short_name' column could be display_name OR original PDB ID for lookup:
            name_row = overview_df_for_names[
                overview_df_for_names['short_name'] == protein_id]  # if short_name is PDB ID
            if not name_row.empty: display_name = name_row['short_name'].iloc[0]
            # If 'processed_structures' was available with 'display_name', that's better.
            # This part needs to ensure `display_name` is correctly fetched.
            # The overview_df generation in load_data attempts to set 'short_name' as display_name or sid.

        domain_plot_info = domain_dict_for_plots.get(protein_id, {'domain': "Unknown", 'average_error': np.nan})
        domain = domain_plot_info['domain']
        error_val = domain_plot_info['average_error']
        function = group_dict_for_summary.get(protein_id, "Unknown")

        summary_data.append([
            display_name,
            round(avg_rmsd_val, 2) if pd.notna(avg_rmsd_val) else 'N/A',
            domain, function,
            round(error_val, 2) if pd.notna(error_val) else 'N/A'
        ])

    summary_df = pd.DataFrame(summary_data,
                              columns=['Protein', 'Average RMSD', 'Domain', 'Molecular Function', 'Error'])
    summary_df = summary_df.sort_values(by=['Molecular Function', 'Domain', 'Average RMSD'], na_position='last')

    csv_path = Path(output_dir) / 'protein_summary.csv'
    summary_df.to_csv(csv_path, index=False)
    print(f"Saved protein summary to {csv_path}")

    all_functions = sorted(
        list(set(s for s in summary_df['Molecular Function'].unique() if s != "Unknown" and pd.notna(s))))
    all_domains = sorted(list(set(s for s in summary_df['Domain'].unique() if s != "Unknown" and pd.notna(s))))

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


def generate_plots(data, output_dir_figures, error_threshold=3.0):
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
            fig, _ = create_opsin_overview_plot(data['overview_df'])  # Adjust if it returns just fig
            if fig:
                fig.savefig(Path(output_dir_figures) / 'A1_opsin_overview.png', dpi=300, bbox_inches='tight')
                plt.close(fig)
                print("  DEBUG: Saved A1_opsin_overview.png")
        except Exception as e:
            print(f"  DEBUG: Error generating opsin overview plot: {e}")
    else:
        print("  DEBUG: Skipped Opsins Overview: 'overview_df' not valid.")

    print("\nDEBUG: 2. Generating RMSD Color Scale Reference...")
    try:
        fig, _ = create_rmsd_color_scale_figure()  # Adjust if it returns just fig
        if fig:
            fig.savefig(Path(output_dir_figures) / 'A2_rmsd_color_scale.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("  DEBUG: Saved A2_rmsd_color_scale.png")
    except Exception as e:
        print(f"  DEBUG: Error generating RMSD color scale: {e}")

    if 'rmsd_df' not in data or not isinstance(data['rmsd_df'], pd.DataFrame) or data['rmsd_df'].empty:
        print("\nDEBUG: CRITICAL: 'rmsd_df' is missing or empty. Skipping RMSD-based plots (B, C).")
    else:
        rmsd_df_orig = data['rmsd_df']
        # --- C. RMSD & Structural Similarity (Simple/Unfiltered) ---
        print("\nDEBUG: 3. Generating Simple Unfiltered Structural Similarity Tree...")
        try:
            fig_simple_tree, _ = plot_similarity_tree(rmsd_df_orig,
                                                      title="Unfiltered Structural Similarity Tree")  # Adjust if returns tuple
            if fig_simple_tree:
                fig_simple_tree.savefig(Path(output_dir_figures) / 'C3_similarity_tree_unfiltered.png', dpi=300,
                                        bbox_inches='tight')
                plt.close(fig_simple_tree)
                print("  DEBUG: Saved C3_similarity_tree_unfiltered.png")
        except Exception as e:
            print(f"  DEBUG: Error generating simple unfiltered similarity tree: {e}")

        # --- B. RMSD & Structural Similarity (Linked & Filtered) ---
        print(f"\nDEBUG: --- Preparing for Linked & Filtered RMSD Plots (Threshold: {error_threshold} Å) ---")
        pdb_list_orig = rmsd_df_orig.index.tolist()
        # Use the consolidated dictionaries from load_data
        group_dict_for_plots = data.get('group_dict', {})
        domain_dict_for_plots_with_error = data.get('domain_dict_for_plots', {})

        kept_ids = []
        # ... (your filtering logic for kept_ids using domain_dict_for_plots_with_error[pdb_id]['average_error']) ...
        # This needs to be careful if domain_dict_for_plots_with_error is not fully populated for all pdb_list_orig
        avg_rmsds_from_matrix = rmsd_df_orig.mean(axis=1)  # Fallback error

        for pdb_id in pdb_list_orig:
            error_val = None
            if pdb_id in domain_dict_for_plots_with_error and \
                    domain_dict_for_plots_with_error[pdb_id].get('average_error') is not None and \
                    pd.notna(domain_dict_for_plots_with_error[pdb_id]['average_error']):
                error_val = domain_dict_for_plots_with_error[pdb_id]['average_error']
            elif pdb_id in avg_rmsds_from_matrix and pd.notna(avg_rmsds_from_matrix[pdb_id]):
                error_val = avg_rmsds_from_matrix[pdb_id]

            if error_val is None or error_val <= error_threshold:  # Keep if no error or below threshold
                kept_ids.append(pdb_id)

        print(f"  DEBUG: Filtered to {len(kept_ids)} structures for linked plots.")

        if len(kept_ids) >= 2:
            filtered_rmsd_df = rmsd_df_orig.loc[kept_ids, kept_ids].copy()
            filtered_group_dict = {k: v for k, v in group_dict_for_plots.items() if k in kept_ids}
            filtered_domain_dict_for_plots_final = {k: v for k, v in domain_dict_for_plots_with_error.items() if
                                                    k in kept_ids}

            Z_linkage = None
            try:
                # ... (Z_linkage calculation from filtered_rmsd_df - your existing robust logic) ...
                matrix_for_linkage = filtered_rmsd_df.fillna(0).values  # Simple fill for now
                np.fill_diagonal(matrix_for_linkage, 0.0)
                condensed_matrix = squareform(matrix_for_linkage, checks=False)  # Turn off checks if issues
                Z_linkage = linkage(condensed_matrix, method='average')
            except Exception as e_link:
                print(f"    DEBUG: ERROR calculating Z_linkage: {e_link}")

            if Z_linkage is not None:
                print("\nDEBUG: 4. Generating Filtered Structural Similarity Tree (Linked)...")
                try:
                    tree_fig, ordered_tree_ids = create_and_visualize_similarity_tree(
                        rmsd_data=filtered_rmsd_df, linkage_matrix=Z_linkage,
                        group_dict=filtered_group_dict, domain_dict=filtered_domain_dict_for_plots_final
                    )
                    if tree_fig:
                        tree_fig.savefig(Path(output_dir_figures) / 'B4_similarity_tree_filtered_linked.png', dpi=300,
                                         bbox_inches='tight')
                        plt.close(tree_fig)
                        print("  DEBUG: Saved B4_similarity_tree_filtered_linked.png")
                except Exception as e_tree:
                    print(f"  DEBUG: Error generating filtered similarity tree: {e_tree}")

                print("\nDEBUG: 5. Generating Filtered RMSD Heatmap (Linked & Clustered)...")
                try:
                    heatmap_clustermap = visualize_rmsd_matrix_improved(
                        rmsd_df=filtered_rmsd_df, linkage_matrix=Z_linkage,
                        group_dict=filtered_group_dict, domain_dict=filtered_domain_dict_for_plots_final,
                        output_file=Path(output_dir_figures) / 'B5_rmsd_heatmap_filtered_linked.png'
                    )
                    if heatmap_clustermap: plt.close(heatmap_clustermap.fig)
                    print("  DEBUG: Saved B5_rmsd_heatmap_filtered_linked.png (by visualize_rmsd_matrix_improved)")
                except Exception as e_heatmap:
                    print(f"  DEBUG: Error generating filtered RMSD heatmap: {e_heatmap}")
        else:
            print(
                "  DEBUG: Skipped Linked Tree and Heatmap (not enough structures after filtering or Z_linkage failed).")

    # --- D. Distance-to-Retinal Plots ---
    # ... (Your D plots logic - ensure data keys like 'ca_distance_table' are checked before use) ...
    # Example for one D plot:
    print("\nDEBUG: Generating CA Distance to Retinal with STD...")
    if 'ca_distance_table' in data and isinstance(data['ca_distance_table'], pd.DataFrame) and not data[
        'ca_distance_table'].empty:
        try:
            fig, _ = plot_distances_with_std(data['ca_distance_table'], title="CA Distance to Retinal by Position",
                                             use_ca=True)  # Adjust if tuple
            if fig:
                fig.savefig(Path(output_dir_figures) / 'D6_ca_distance_retinal_std.png', dpi=300, bbox_inches='tight')
                plt.close(fig)
                print("  DEBUG: Saved D6_ca_distance_retinal_std.png")
        except Exception as e:
            print(f"  DEBUG: Error generating CA Distance STD plot: {e}")
    else:
        print("  DEBUG: Skipped CA Distance STD plot: 'ca_distance_table' invalid.")

    # --- E. Sequence Conservation & Composition ---
    # ... (Your E plots logic - ensure data keys like 'msa_table' are checked) ...
    # Example for one E plot:
    print("\nDEBUG: Generating Residue Conservation Bar Plot...")
    if 'msa_table' in data and isinstance(data['msa_table'], pd.DataFrame) and not data['msa_table'].empty:
        try:
            fig, _ = create_residue_conservation_plot(data['msa_table'], helix_highlighting=True)  # Adjust if tuple
            if fig:
                fig.savefig(Path(output_dir_figures) / 'E12_residue_conservation.png', dpi=300, bbox_inches='tight')
                plt.close(fig)
                print("  DEBUG: Saved E12_residue_conservation.png")
        except Exception as e:
            print(f"  DEBUG: Error generating residue conservation plot: {e}")
    else:
        print("  DEBUG: Skipped Residue Conservation plot: 'msa_table' invalid.")

    # --- F. Combined Plots (Optional) ---
    # ... (Your F plots logic) ...
    print(f"\n--- All plot generation attempts complete. Check {output_dir_figures}. ---")


def verify_filtering_in_summary(summary_csv_path, error_threshold=3.0):
    # ... (your verification logic - no changes needed here for this problem) ...
    pass


def main():
    parser = argparse.ArgumentParser(description='Generate plots for opsin analysis.')
    parser.add_argument('--input-dir', '-i', type=str, default='opsin_output/', help='Directory for input data.')
    parser.add_argument('--output-dir', '-o', type=str, default='opsin_plots_output_debug/',
                        help='Directory for output plots.')  # Changed default
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

    print(f"Opsin Plot Generation Script - DEBUG ENABLED")  # Indicate debug mode
    print(f"Input Dir: {Path(args.input_dir).resolve()}")
    print(f"Output Plots Dir: {output_dir_figures.resolve()}")
    print(f"Output Summaries Dir: {output_dir_summaries.resolve()}")

    dpi_map = {'low': 100, 'medium': 200, 'high': 300}
    plt.rcParams['figure.dpi'] = dpi_map[args.quality]
    plt.rcParams['savefig.dpi'] = dpi_map[args.quality]

    loaded_data = load_data(args.input_dir, args.input_dir, args.chain_id)  # output_dir_ref for CSVs is input_dir

    if not loaded_data or not isinstance(loaded_data, dict):
        print("CRITICAL DEBUG: No data loaded or data is not a dictionary. Exiting plot generation.")
        return

    generate_plots(loaded_data, output_dir_figures, error_threshold=args.error_threshold)
    summary_csv_path = generate_summary_csv(loaded_data, output_dir_summaries)
    # if summary_csv_path:
    #     verify_filtering_in_summary(summary_csv_path, error_threshold=args.error_threshold)

    print(f"\nAll operations complete. Outputs are in {Path(args.output_dir).resolve()}")


if __name__ == "__main__":
    main()