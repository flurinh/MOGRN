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

# --- Import visualization functions ---
# Assuming visualization_functions.py is in the same directory or accessible via PYTHONPATH
# If it's in a sub-module, adjust the import path e.g., from .visualization_module import ...
try:
    from visualization_functions import (
        # Overview & Reference
        create_opsin_overview_plot,  # Note: This was create_overview_plot in plotfigures.py's imports
        create_rmsd_color_scale_figure,

        # RMSD & Structural Similarity (Linked & Filtered)
        create_and_visualize_similarity_tree,
        visualize_rmsd_matrix_improved,

        # RMSD & Structural Similarity (Simple/Unfiltered)
        plot_similarity_tree, # Simple tree
        plot_rmsd_heatmap,    # Simple heatmap (optional, if desired)

        # Distance-to-Retinal Plots
        plot_distances_with_std,
        plot_average_distances_by_helix,
        plot_distance_heatmap,

        # Sequence Conservation & Composition
        create_residue_conservation_plot,
        plot_helix_logo_plots,
        plot_conservation_around_x50,

        # Combined Plots (Optional)
        create_combined_distance_logo_plot,

        # Interactive Plots (Optional)
        visualize_binding_pocket # For Plotly
    )
except ImportError:
    print("Error: visualization_functions.py not found. Make sure it's in the correct path.")
    print("Attempting to import from current directory as a fallback for development.")
    # This is a common structure if plot.py and visualization_functions.py are siblings
    from visualization_functions import (
        create_opsin_overview_plot, create_rmsd_color_scale_figure,
        create_and_visualize_similarity_tree, visualize_rmsd_matrix_improved,
        plot_similarity_tree, plot_rmsd_heatmap,
        plot_distances_with_std, plot_average_distances_by_helix, plot_distance_heatmap,
        create_residue_conservation_plot, plot_helix_logo_plots, plot_conservation_around_x50,
        create_combined_distance_logo_plot, visualize_binding_pocket
    )


# Import the color scheme tools
try:
    from opsin_color_scheme import get_group_colors
except ImportError:
    print("Error: opsin_color_scheme.py not found.")
    # Fallback dummy function if not found, to allow script to run partially
    def get_group_colors(items):
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
    cache_dir = Path(input_dir) / 'cache' # Use Path for better path handling

    print(f"Attempting to load data. Cache directory: {cache_dir}")

    # 1. Try to load from cache files first
    if cache_dir.exists():
        print(f"Found cache directory: {cache_dir}")
        files_to_load_cache = [
            'raw_structures.pkl',
            f'processed_structures_{chain_id}.pkl',
            f'helix_annotations_{chain_id}.pkl',
            f'structure_comparison_{chain_id}.pkl', # Key for rmsd_matrix, group_info
            f'structure_errors_{chain_id}.pkl',    # Key for error filtering
            f'grn_assignment_{chain_id}.pkl'
        ]
        for file_name in files_to_load_cache:
            cache_path = cache_dir / file_name
            if cache_path.exists():
                try:
                    print(f"Loading from cache: {file_name}...")
                    with open(cache_path, 'rb') as f:
                        file_data = pickle.load(f)
                        # Special handling for structure_comparison to extract rmsd_matrix
                        if file_name == f'structure_comparison_{chain_id}.pkl':
                            print(f"  Keys in structure_comparison: {list(file_data.keys())}")
                            if 'rmsd_matrix' in file_data and isinstance(file_data['rmsd_matrix'], pd.DataFrame):
                                data['rmsd_df'] = file_data['rmsd_matrix']
                                print(f"    Loaded 'rmsd_df' from structure_comparison. Shape: {data['rmsd_df'].shape}")
                                # Load other potential keys from this file too
                                for key, value in file_data.items():
                                    if key != 'rmsd_matrix':
                                        data[key] = value
                                continue # Skip generic update for this file
                        data.update(file_data)
                    print(f"  Successfully loaded {file_name}")
                except Exception as e:
                    print(f"  Error loading {file_name}: {e}")
            else:
                print(f"  Warning: Cache file not found: {cache_path}")
    else:
        print(f"Cache directory not found: {cache_dir}. Will rely on CSVs.")

    # 2. Load CSV files (either as primary source or fallback)
    # output_dir_ref is where CSVs are expected (e.g., args.input_dir or a specific 'opsin_output')
    csv_files_to_load = [
        ('rmsd_matrix.csv', 'rmsd_df'), # Load if not already from cache
        ('molecular_functions.csv', 'molecular_functions_df'), # For group_dict
        ('ca_distance_table_grn.csv', 'ca_distance_table'),
        ('distance_table_grn.csv', 'distance_table'), # Sidechain distances
        ('ca_msa_table_grn.csv', 'msa_table'), # For conservation, preferred for CA-based alignment
        ('ca_residue_table_grn.csv', 'residue_table'), # For helix logos
        # Add other potential CSVs if needed by plots
        ('mo_exp_errors.csv', 'mo_exp_errors_df'), # For error filtering fallback
        ('hideaki_errors.csv', 'hideaki_errors_df') # For error filtering fallback
    ]

    for file_name, key_name in csv_files_to_load:
        # Skip loading rmsd_df if already loaded from cache
        if key_name == 'rmsd_df' and 'rmsd_df' in data:
            print(f"  'rmsd_df' already loaded from cache. Skipping {file_name}.")
            continue

        csv_path = Path(output_dir_ref) / file_name
        if csv_path.exists():
            print(f"Loading from CSV: {file_name}...")
            try:
                # Use index_col=0 for tables where first column is the index
                if key_name in ['rmsd_df', 'ca_distance_table', 'distance_table', 'msa_table', 'residue_table']:
                    df = pd.read_csv(csv_path, index_col=0)
                else:
                    df = pd.read_csv(csv_path)
                data[key_name] = df
                print(f"  Successfully loaded {file_name} into '{key_name}'. Shape: {df.shape}")
            except Exception as e:
                print(f"  Error loading {file_name}: {e}")
        else:
            print(f"  Warning: CSV file not found: {csv_path}")

    # 3. Post-processing and data structuring (critical for plotting functions)

    # Ensure 'rmsd_df' exists
    if 'rmsd_df' not in data:
        print("CRITICAL WARNING: 'rmsd_df' (RMSD matrix) not found in cache or CSVs. Many plots will fail.")
    elif not isinstance(data['rmsd_df'], pd.DataFrame):
        print(f"CRITICAL WARNING: 'rmsd_df' is not a DataFrame (type: {type(data['rmsd_df'])}). Many plots will fail.")
    else:
        print(f"RMSD matrix 'rmsd_df' loaded with {len(data['rmsd_df'])} structures.")


    # Populate 'group_dict' (structure_id -> molecular_function)
    # Priority: processed_structures -> molecular_functions_df -> derive from rmsd_df names
    if 'processed_structures' in data and isinstance(data['processed_structures'], dict):
        data['group_dict'] = {}
        for sid, s_data in data['processed_structures'].items():
            if isinstance(s_data, dict) and 'properties' in s_data and \
               isinstance(s_data['properties'], dict) and 'molecular_function' in s_data['properties']:
                data['group_dict'][sid] = s_data['properties']['molecular_function']
            elif 'rmsd_df' in data and sid in data['rmsd_df'].index: # Ensure sid is in RMSD matrix
                 data['group_dict'][sid] = "Unknown" # Fallback if properties missing
        print(f"  Populated 'group_dict' from 'processed_structures' for {len(data['group_dict'])} entries.")
    elif 'molecular_functions_df' in data:
        df_mol_func = data['molecular_functions_df']
        if 'structure_id' in df_mol_func.columns and 'molecular_function' in df_mol_func.columns:
            data['group_dict'] = dict(zip(df_mol_func['structure_id'], df_mol_func['molecular_function']))
            print(f"  Populated 'group_dict' from 'molecular_functions_df' for {len(data['group_dict'])} entries.")
    elif 'rmsd_df' in data and isinstance(data['rmsd_df'], pd.DataFrame):
        print("  Deriving 'group_dict' from structure names in 'rmsd_df' as a fallback.")
        data['group_dict'] = {}
        for sid in data['rmsd_df'].index:
            if 'ChR' in sid or 'channel' in sid.lower(): data['group_dict'][sid] = "Cation channel"
            elif 'HR' in sid or 'pump' in sid.lower() or 'PR' in sid: data['group_dict'][sid] = "Proton pump"
            elif 'ACR' in sid or 'chloride' in sid.lower(): data['group_dict'][sid] = "Chloride pump"
            else: data['group_dict'][sid] = "Unknown"

    # Populate 'domain_dict_for_plots' (structure_id -> {'domain': 'X', 'average_error': Y})
    # This is specifically for the linked tree/heatmap functions.
    # Priority: structure_errors_CHAIN.pkl -> processed_structures_CHAIN.pkl -> mo_exp_errors/hideaki_errors
    data['domain_dict_for_plots'] = {}
    all_sids_for_domain_dict = []
    if 'rmsd_df' in data and isinstance(data['rmsd_df'], pd.DataFrame):
        all_sids_for_domain_dict = data['rmsd_df'].index.tolist()

    for sid in all_sids_for_domain_dict:
        domain_info = {'domain': "Unknown", 'average_error': None}
        # From structure_errors (highest priority for error)
        if f'structure_errors' in data and sid in data[f'structure_errors']:
            s_error_data = data[f'structure_errors'][sid]
            if isinstance(s_error_data, dict) and 'average_error' in s_error_data:
                if isinstance(s_error_data['average_error'], (int, float)):
                    domain_info['average_error'] = float(s_error_data['average_error'])

        # From processed_structures (for domain and potentially other props)
        if f'processed_structures' in data and sid in data[f'processed_structures']:
            s_proc_data = data[f'processed_structures'][sid]
            if isinstance(s_proc_data, dict) and 'properties' in s_proc_data:
                props = s_proc_data['properties']
                if isinstance(props, dict):
                    domain_info['domain'] = str(props.get('domain', domain_info['domain'])) # Ensure string

        # Fallback error data from CSVs if not found in structure_errors
        if domain_info['average_error'] is None:
            for err_df_key in ['mo_exp_errors_df', 'hideaki_errors_df']:
                if err_df_key in data:
                    err_df = data[err_df_key]
                    if 'structure_id' in err_df.columns and sid in err_df['structure_id'].values:
                        row = err_df[err_df['structure_id'] == sid].iloc[0]
                        error_cols = ['backbone_rmsd', 'pocket_rmsd', 'retinal_rmsd']
                        valid_cols = [col for col in error_cols if col in row and pd.notna(row[col])]
                        if valid_cols:
                            avg_err = row[valid_cols].mean()
                            if pd.notna(avg_err):
                                domain_info['average_error'] = float(avg_err)
                                domain_info['domain'] = str(row.get('domain', domain_info['domain'])) # update domain too
                                break # Found error

        data['domain_dict_for_plots'][sid] = domain_info

    if not data['domain_dict_for_plots'] and 'rmsd_df' in data and isinstance(data['rmsd_df'], pd.DataFrame):
         print("  Warning: 'domain_dict_for_plots' is empty. Filling with Unknown domains and no errors.")
         data['domain_dict_for_plots'] = {
             sid: {'domain': "Unknown", 'average_error': None} for sid in data['rmsd_df'].index
         }
    print(f"  Populated 'domain_dict_for_plots' for {len(data['domain_dict_for_plots'])} structures.")


    # Ensure 'msa_table' and 'residue_table' are available for sequence plots
    if 'msa_table' not in data and 'ca_msa_table_grn' in data: # Fallback if specific msa_table missing
        data['msa_table'] = data['ca_msa_table_grn']
        print("  Used 'ca_msa_table_grn' as 'msa_table'.")
    if 'residue_table' not in data and 'msa_table' in data: # Fallback for residue table
        data['residue_table'] = data['msa_table'] # Or ca_residue_table_grn if that's more appropriate
        print("  Used 'msa_table' as 'residue_table'.")

    # For opsin overview plot, structure the input data
    if 'processed_structures' in data and 'rmsd_df' in data and isinstance(data['rmsd_df'], pd.DataFrame):
        overview_data_list = []
        for sid in data['rmsd_df'].index: # Iterate over SIDs in RMSD matrix
            s_data = data['processed_structures'].get(sid, {})
            props = s_data.get('properties', {})
            overview_data_list.append({
                'short_name': s_data.get('display_name', sid),
                'molecular_function_normalized': props.get('molecular_function', 'Unknown'),
                'domain': props.get('domain', 'Unknown'),
                'experimentally_determined': props.get('is_experimental', False) # Assuming this key exists
            })
        data['overview_df'] = pd.DataFrame(overview_data_list)
        print(f"  Prepared 'overview_df' for opsin overview plot with {len(data['overview_df'])} entries.")


    print("Data loading and initial structuring complete.")
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

    for protein_id in rmsd_df.index:
        # Calculate average RMSD to all other proteins
        rmsd_values = rmsd_df.loc[protein_id].values
        # Exclude self-comparison (0) and NaNs if any
        valid_rmsd_values = rmsd_values[~np.isnan(rmsd_values) & (rmsd_values > 1e-6)]
        avg_rmsd_val = np.mean(valid_rmsd_values) if len(valid_rmsd_values) > 0 else np.nan

        s_data = processed_structures.get(protein_id, {})
        display_name = s_data.get('display_name', protein_id)

        domain_plot_info = domain_dict_for_plots.get(protein_id, {'domain': "Unknown", 'average_error': np.nan})
        domain = domain_plot_info['domain']
        error_val = domain_plot_info['average_error']

        # Get molecular function from group_dict
        function = data.get('group_dict', {}).get(protein_id, "Unknown")

        summary_data.append([
            display_name,
            round(avg_rmsd_val, 2) if pd.notna(avg_rmsd_val) else 'N/A',
            domain,
            function,
            round(error_val, 2) if pd.notna(error_val) else 'N/A'
        ])

    summary_df = pd.DataFrame(summary_data, columns=['Protein', 'Average RMSD', 'Domain', 'Molecular Function', 'Error'])
    summary_df = summary_df.sort_values(by='Average RMSD', na_position='last')

    csv_path = Path(output_dir) / 'protein_summary.csv'
    summary_df.to_csv(csv_path, index=False)
    print(f"Saved protein summary to {csv_path}")

    # Save color assignments
    all_functions = sorted(list(set(summary_df['Molecular Function'].unique()) - {'Unknown'}))
    all_domains = sorted(list(set(summary_df['Domain'].unique()) - {'Unknown'}))

    function_colors = get_group_colors(all_functions)
    domain_colors = get_group_colors(all_domains)

    color_data = []
    for func, color in function_colors.items(): color_data.append(['Function', func, color])
    for dom, color in domain_colors.items(): color_data.append(['Domain', dom, color])

    color_df = pd.DataFrame(color_data, columns=['Type', 'Value', 'Color'])
    color_path = Path(output_dir) / 'color_assignments.csv'
    color_df.to_csv(color_path, index=False)
    print(f"Saved color assignments to {color_path}")

    return csv_path


def generate_plots(data, output_dir_figures, error_threshold=3.0):
    """
    Generate all curated plots and save them to the output directory.
    """
    print(f"\n--- Generating plots in {output_dir_figures} ---")
    Path(output_dir_figures).mkdir(parents=True, exist_ok=True)

    # --- A. Overview & Reference ---
    print("\n1. Generating Opsins Overview Plot...")
    if 'overview_df' in data and not data['overview_df'].empty:
        try:
            fig = create_opsin_overview_plot(data['overview_df'])
            fig.savefig(Path(output_dir_figures) / 'A1_opsin_overview.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("  Saved A1_opsin_overview.png")
        except Exception as e:
            print(f"  Error generating opsin overview plot: {e}")
    else:
        print("  Skipped: 'overview_df' not available or empty.")

    print("\n2. Generating RMSD Color Scale Reference...")
    try:
        fig = create_rmsd_color_scale_figure()
        fig.savefig(Path(output_dir_figures) / 'A2_rmsd_color_scale.png', dpi=300, bbox_inches='tight')
        plt.close(fig)
        print("  Saved A2_rmsd_color_scale.png")
    except Exception as e:
        print(f"  Error generating RMSD color scale: {e}")


    # --- Data Availability Check for RMSD-based plots ---
    if 'rmsd_df' not in data or not isinstance(data['rmsd_df'], pd.DataFrame) or data['rmsd_df'].empty:
        print("\nCRITICAL: 'rmsd_df' is missing or empty. Skipping all RMSD-based plots (B, C).")
    else:
        rmsd_df_orig = data['rmsd_df']

        # --- C. RMSD & Structural Similarity (Simple/Unfiltered) ---
        print("\n3. Generating Simple Unfiltered Structural Similarity Tree...")
        try:
            fig_simple_tree = plot_similarity_tree(
                rmsd_df_orig,
                title="Unfiltered Structural Similarity Tree"
            )
            fig_simple_tree.savefig(Path(output_dir_figures) / 'C3_similarity_tree_unfiltered.png', dpi=300, bbox_inches='tight')
            plt.close(fig_simple_tree)
            print("  Saved C3_similarity_tree_unfiltered.png")
        except Exception as e:
            print(f"  Error generating simple unfiltered similarity tree: {e}")

        # Optional: Simple Unfiltered RMSD Heatmap
        # print("\nGenerating Simple Unfiltered RMSD Heatmap...")
        # try:
        #     fig_simple_heatmap = plot_rmsd_heatmap(rmsd_df_orig, title="Unfiltered RMSD Heatmap")
        #     fig_simple_heatmap.savefig(Path(output_dir_figures) / 'C_rmsd_heatmap_unfiltered.png', dpi=300, bbox_inches='tight')
        #     plt.close(fig_simple_heatmap)
        #     print("  Saved C_rmsd_heatmap_unfiltered.png")
        # except Exception as e:
        #     print(f"  Error generating simple unfiltered RMSD heatmap: {e}")


        # --- B. RMSD & Structural Similarity (Linked & Filtered) ---
        print(f"\n--- Preparing for Linked & Filtered RMSD Plots (Threshold: {error_threshold} Å) ---")
        # 1. Centralized Filtering (adapted from plotfigures.py)
        pdb_list_orig = rmsd_df_orig.index.tolist()
        group_dict_orig = data.get('group_dict', {sid: "Unknown" for sid in pdb_list_orig})
        # domain_dict_for_plots has {'domain': NAME, 'average_error': VALUE}
        domain_dict_for_plots = data.get('domain_dict_for_plots', {sid: {'domain': "Unknown", 'average_error': None} for sid in pdb_list_orig})

        # Calculate average RMSDs from the matrix itself as a fallback for error
        avg_rmsds_from_matrix = rmsd_df_orig.mean()

        kept_ids = []
        filtered_out_details = []
        print("  Filtering structures...")
        for pdb_id in pdb_list_orig:
            error_val = None
            source = "N/A"
            if pdb_id in domain_dict_for_plots and domain_dict_for_plots[pdb_id].get('average_error') is not None:
                error_val = domain_dict_for_plots[pdb_id]['average_error']
                source = "explicit_error"
            elif pdb_id in avg_rmsds_from_matrix:
                error_val = avg_rmsds_from_matrix[pdb_id]
                source = "avg_rmsd_from_matrix"

            if error_val is not None and error_val <= error_threshold:
                kept_ids.append(pdb_id)
            elif error_val is not None: # Filtered out
                filtered_out_details.append(f"    - {pdb_id}: value={error_val:.2f} Å (source: {source})")
            else: # No error data, keep by default or filter based on policy
                # Current policy: if no error data, keep it unless strict filtering is on.
                # For now, let's keep them if no error data. Or you can filter them:
                # filtered_out_details.append(f"    - {pdb_id}: No error data, filtered out.")
                kept_ids.append(pdb_id) # Keep if no error data
                print(f"    - {pdb_id}: No error data, kept by default.")


        if filtered_out_details:
            print(f"  Filtered out {len(filtered_out_details)} structures (value > {error_threshold}):")
            for detail in filtered_out_details[:5]: print(detail)
            if len(filtered_out_details) > 5: print(f"    ... and {len(filtered_out_details) - 5} more.")
        else:
            print("  No structures filtered based on error/RMSD threshold.")

        if len(kept_ids) < 2:
            print("  Skipped Linked Plots: Less than 2 structures remaining after filtering.")
        else:
            print(f"  Proceeding with {len(kept_ids)} structures for linked visualizations.")
            filtered_rmsd_df = rmsd_df_orig.loc[kept_ids, kept_ids].copy() # Ensure it's a copy

            # Filter group_dict and domain_dict to only include kept IDs for the plots
            filtered_group_dict = {k: v for k, v in group_dict_orig.items() if k in kept_ids}
            filtered_domain_dict_for_plots = {k: v for k, v in domain_dict_for_plots.items() if k in kept_ids}

            # 2. Calculate Linkage Matrix (Z)
            Z_linkage = None
            try:
                print("  Calculating linkage matrix Z for filtered data...")
                # Clean matrix for linkage (handle NaNs, ensure diagonal is 0, check symmetry)
                matrix_for_linkage = filtered_rmsd_df.values.copy() # Work on a copy
                if np.any(np.isnan(matrix_for_linkage)) or np.any(np.isinf(matrix_for_linkage)):
                    mean_finite = np.nanmean(matrix_for_linkage[np.isfinite(matrix_for_linkage)])
                    fill_val = mean_finite if pd.notna(mean_finite) else 10.0 # Fallback
                    matrix_for_linkage = np.nan_to_num(matrix_for_linkage, nan=fill_val, posinf=fill_val*1.1, neginf=0.0)
                np.fill_diagonal(matrix_for_linkage, 0.0)
                if not np.allclose(matrix_for_linkage, matrix_for_linkage.T, atol=1e-5): # Check symmetry
                    print("    Warning: RMSD matrix for linkage is not perfectly symmetric. Forcing symmetry.")
                    matrix_for_linkage = (matrix_for_linkage + matrix_for_linkage.T) / 2.0
                    np.fill_diagonal(matrix_for_linkage, 0.0)

                condensed_matrix = squareform(matrix_for_linkage, checks=True) # Enable checks
                Z_linkage = linkage(condensed_matrix, method='average')
                print("    Successfully calculated linkage matrix Z.")
            except Exception as e:
                print(f"    ERROR calculating linkage matrix: {e}. Linked plots might fail or be inaccurate.")
                import traceback
                traceback.print_exc()

            if Z_linkage is not None:
                print("\n4. Generating Filtered Structural Similarity Tree (Linked)...")
                try:
                    tree_fig, ordered_tree_ids = create_and_visualize_similarity_tree(
                        rmsd_data=filtered_rmsd_df,
                        linkage_matrix=Z_linkage,
                        group_dict=filtered_group_dict,
                        domain_dict=filtered_domain_dict_for_plots # This has the {'domain': X, 'error': Y} structure
                    )
                    # Saving components if the function doesn't do it.
                    # The function create_and_visualize_similarity_tree should ideally handle saving its components.
                    # For now, just save the main figure.
                    tree_fig.savefig(Path(output_dir_figures) / 'B4_similarity_tree_filtered_linked.png', dpi=300, bbox_inches='tight')
                    # If create_and_visualize_similarity_tree returns paths to components, use them.
                    # Else, you might need to adapt plotfigures.py's legend saving logic here.
                    plt.close(tree_fig)
                    print("  Saved B4_similarity_tree_filtered_linked.png")
                except Exception as e:
                    print(f"  Error generating filtered similarity tree: {e}")
                    import traceback
                    traceback.print_exc()

                print("\n5. Generating Filtered RMSD Heatmap (Linked & Clustered)...")
                try:
                    heatmap_clustermap = visualize_rmsd_matrix_improved(
                        rmsd_df=filtered_rmsd_df,
                        linkage_matrix=Z_linkage, # Crucial
                        group_dict=filtered_group_dict,
                        domain_dict=filtered_domain_dict_for_plots,
                        output_file=Path(output_dir_figures) / 'B5_rmsd_heatmap_filtered_linked.png'
                    )
                    if heatmap_clustermap: # visualize_rmsd_matrix_improved saves internally and returns clustermap
                        plt.close(heatmap_clustermap.fig) # Close the figure associated with clustermap
                        print("  Saved B5_rmsd_heatmap_filtered_linked.png")
                    else:
                        print("  Clustermap generation failed or returned None.")
                except Exception as e:
                    print(f"  Error generating filtered RMSD heatmap: {e}")
                    import traceback
                    traceback.print_exc()
            else:
                print("  Skipped Linked Tree and Heatmap due to Z_linkage calculation failure.")

    # --- D. Distance-to-Retinal Plots ---
    print("\n--- Generating Distance-to-Retinal Plots ---")
    plot_configs_dist = [
        ('ca_distance_table', 'D6_ca_distance_retinal_std.png', "CA Distance to Retinal by Position", True),
        ('distance_table', 'D7_sidechain_distance_retinal_std.png', "Sidechain Distance to Retinal by Position", False),
    ]
    for key, fname, title, use_ca_flag in plot_configs_dist:
        print(f"\n{fname.split('_')[1].upper()} {title}...")
        if key in data and not data[key].empty:
            try:
                fig = plot_distances_with_std(data[key], title=title, use_ca=use_ca_flag)
                fig.savefig(Path(output_dir_figures) / fname, dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"  Saved {fname}")
            except Exception as e:
                print(f"  Error generating {title}: {e}")
        else:
            print(f"  Skipped: Data key '{key}' not available or empty.")

    plot_configs_avg_helix = [
        ('ca_distance_table', 'D8_ca_avg_distance_by_helix.png', "CA Avg Distances by Helix", True),
        ('distance_table', 'D9_sidechain_avg_distance_by_helix.png', "Sidechain Avg Distances by Helix", False),
    ]
    for key, fname, title, use_ca_flag in plot_configs_avg_helix:
        print(f"\n{title}...")
        if key in data and not data[key].empty:
            try:
                fig = plot_average_distances_by_helix(data[key], use_ca=use_ca_flag)
                fig.savefig(Path(output_dir_figures) / fname, dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"  Saved {fname}")
            except Exception as e:
                print(f"  Error generating {title}: {e}")
        else:
            print(f"  Skipped: Data key '{key}' not available or empty.")

    plot_configs_dist_heatmap = [
        ('ca_distance_table', 'D10_ca_distance_heatmap.png', "CA Distance Heatmap"),
        ('distance_table', 'D11_sidechain_distance_heatmap.png', "Sidechain Distance Heatmap"),
    ]
    for key, fname, title in plot_configs_dist_heatmap:
        print(f"\n{title}...")
        if key in data and not data[key].empty:
            try:
                fig = plot_distance_heatmap(data[key]) # title is handled inside
                fig.savefig(Path(output_dir_figures) / fname, dpi=300, bbox_inches='tight')
                plt.close(fig)
                print(f"  Saved {fname}")
            except Exception as e:
                print(f"  Error generating {title}: {e}")
        else:
            print(f"  Skipped: Data key '{key}' not available or empty.")

    # --- E. Sequence Conservation & Composition ---
    print("\n--- Generating Sequence Conservation & Composition Plots ---")
    print("\n12. Residue Conservation Bar Plot...")
    if 'msa_table' in data and not data['msa_table'].empty:
        try:
            fig = create_residue_conservation_plot(data['msa_table'], helix_highlighting=True)
            fig.savefig(Path(output_dir_figures) / 'E12_residue_conservation.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("  Saved E12_residue_conservation.png")
        except Exception as e:
            print(f"  Error generating residue conservation plot: {e}")
    else:
        print("  Skipped: 'msa_table' not available or empty.")

    print("\n13. Helix Sequence Logos (X.50 +/- 4)...")
    if 'residue_table' in data and not data['residue_table'].empty:
        try:
            fig = plot_helix_logo_plots(data['residue_table'])
            fig.savefig(Path(output_dir_figures) / 'E13_helix_logos_x50.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("  Saved E13_helix_logos_x50.png")
        except Exception as e:
            print(f"  Error generating helix logo plots: {e}")
    else:
        print("  Skipped: 'residue_table' not available or empty.")

    print("\n14. Conservation-Weighted Sequence Logos (X.50 +/- 4)...")
    if 'residue_table' in data and not data['residue_table'].empty:
        try:
            fig = plot_conservation_around_x50(data['residue_table'])
            fig.savefig(Path(output_dir_figures) / 'E14_conservation_weighted_logos_x50.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("  Saved E14_conservation_weighted_logos_x50.png")
        except Exception as e:
            print(f"  Error generating conservation-weighted logos: {e}")
    else:
        print("  Skipped: 'residue_table' not available or empty.")

    # --- F. Combined Plots (Optional) ---
    print("\n--- Generating Combined Plots (Optional) ---")
    print("\n15. Combined Distance Line Plot & Sequence Logo...")
    if 'distance_table' in data and not data['distance_table'].empty and \
       'msa_table' in data and not data['msa_table'].empty:
        try:
            # Assuming distance_table is for sidechains; choose ca_distance_table if preferred
            fig = create_combined_distance_logo_plot(data['distance_table'], data['msa_table'])
            fig.savefig(Path(output_dir_figures) / 'F15_combined_distance_logo.png', dpi=300, bbox_inches='tight')
            plt.close(fig)
            print("  Saved F15_combined_distance_logo.png")
        except Exception as e:
            print(f"  Error generating combined distance logo plot: {e}")
            import traceback
            traceback.print_exc()
    else:
        print("  Skipped: 'distance_table' or 'msa_table' not available or empty.")

    # --- G. Interactive Plots (Example - save as HTML) ---
    # print("\n--- Generating Interactive Plots (Optional) ---")
    # print("\n16. Binding Pocket 3D Visualization (Example)...")
    # if 'processed_structures' in data and 'SOME_STRUCTURE_ID' in data['processed_structures']:
    #     try:
    #         # This requires specific data for a single structure
    #         # Example: structure_df = data['processed_structures']['SOME_STRUCTURE_ID']['atom_df']
    #         #          retinal_df = data['processed_structures']['SOME_STRUCTURE_ID']['retinal_df']
    #         #          residue_ids_pocket = [...]
    #         # fig = visualize_binding_pocket(structure_df, residue_ids_pocket, retinal_df)
    #         # fig.write_html(Path(output_dir_figures) / 'G16_binding_pocket_interactive.html')
    #         # print("  Saved G16_binding_pocket_interactive.html (Example - needs specific data setup)")
    #         print("  Skipped: Example interactive plot, requires specific data setup for a single structure.")
    #     except Exception as e:
    #         print(f"  Error generating binding pocket plot: {e}")
    # else:
    #     print("  Skipped: Data for example binding pocket plot not available.")

    print(f"\n--- All plot generation attempts complete. Check {output_dir_figures}. ---")


def verify_filtering_in_summary(summary_csv_path, error_threshold=3.0):
    """
    Verify that the protein_summary.csv doesn't contain structures that should have been filtered.
    """
    if not summary_csv_path or not Path(summary_csv_path).exists():
        print("\nVerification: Protein summary CSV not found.")
        return

    print(f"\n--- Verifying Filtering in Protein Summary ({summary_csv_path}) ---")
    df = pd.read_csv(summary_csv_path)
    violations = 0

    # Check Average RMSD column
    if 'Average RMSD' in df.columns:
        # Convert to numeric, coercing errors to NaN
        df['Average RMSD Numeric'] = pd.to_numeric(df['Average RMSD'], errors='coerce')
        high_rmsd_entries = df[df['Average RMSD Numeric'] > error_threshold]
        if not high_rmsd_entries.empty:
            print(f"  WARNING: Found {len(high_rmsd_entries)} entries with Average RMSD > {error_threshold}:")
            for _, row in high_rmsd_entries.iterrows():
                print(f"    - {row['Protein']}: Avg RMSD = {row['Average RMSD']}")
            violations += len(high_rmsd_entries)

    # Check Error column
    if 'Error' in df.columns:
        df['Error Numeric'] = pd.to_numeric(df['Error'], errors='coerce')
        high_error_entries = df[df['Error Numeric'] > error_threshold]
        if not high_error_entries.empty:
            print(f"  WARNING: Found {len(high_error_entries)} entries with Error > {error_threshold}:")
            for _, row in high_error_entries.iterrows():
                print(f"    - {row['Protein']}: Error = {row['Error']}")
            violations += len(high_error_entries)
            # Note: This might double count if a structure has high avg RMSD AND high explicit error.

    if violations == 0:
        print(f"  Verification PASSED: No entries found in summary exceeding threshold {error_threshold} Å (based on available numeric values).")
    else:
        print(f"  Verification FAILED: {violations} potential filtering violations found.")
    print("  Note: This verification checks the summary CSV. Linked plots use their own live filtering.")


def main():
    parser = argparse.ArgumentParser(description='Generate plots for opsin analysis.')
    parser.add_argument('--input-dir', '-i', type=str, default='.',
                        help='Directory containing input data files (cache subdir, CSVs). Default: current directory.')
    parser.add_argument('--output-dir', '-o', type=str, default='opsin_plots_output',
                        help='Directory to save output plots and summaries. Default: opsin_plots_output')
    parser.add_argument('--chain-id', '-c', type=str, default='A',
                        help='Chain ID used in the analysis (default: A)')
    parser.add_argument('--quality', '-q', type=str, choices=['low', 'medium', 'high'], default='high',
                        help='Figure quality (affects DPI). Default: high')
    parser.add_argument('--error-threshold', '-et', type=float, default=3.0,
                        help='RMSD/Error threshold in Angstrom for filtering linked plots. Default: 3.0')
    args = parser.parse_args()

    # Ensure output directory exists for plots and summaries
    output_dir_figures = Path(args.output_dir) / 'figures'
    output_dir_summaries = Path(args.output_dir) / 'summaries'
    output_dir_figures.mkdir(parents=True, exist_ok=True)
    output_dir_summaries.mkdir(parents=True, exist_ok=True)

    print(f"Opsin Plot Generation Script")
    print(f"Input Data Directory: {Path(args.input_dir).resolve()}")
    print(f"Output Directory (Plots): {output_dir_figures.resolve()}")
    print(f"Output Directory (Summaries): {output_dir_summaries.resolve()}")
    print(f"Chain ID: {args.chain_id}")
    print(f"Figure Quality: {args.quality}")
    print(f"Error/RMSD Filtering Threshold for Linked Plots: {args.error_threshold} Å")

    # Set Matplotlib figure quality parameters
    dpi_map = {'low': 100, 'medium': 200, 'high': 300}
    plt.rcParams['figure.dpi'] = dpi_map[args.quality]
    plt.rcParams['savefig.dpi'] = dpi_map[args.quality]

    # Load data
    # The load_data function will look for cache in 'input_dir/cache'
    # and CSVs in 'input_dir' (or a specified path if you modify load_data)
    loaded_data = load_data(args.input_dir, args.input_dir, args.chain_id)

    if not loaded_data:
        print("No data loaded. Exiting.")
        return

    # Generate plots
    generate_plots(loaded_data, output_dir_figures, error_threshold=args.error_threshold)

    # Generate protein summary CSV
    summary_csv_path = generate_summary_csv(loaded_data, output_dir_summaries)

    # Verify filtering in the summary (optional sanity check)
    if summary_csv_path:
        verify_filtering_in_summary(summary_csv_path, error_threshold=args.error_threshold)

    print(f"\nAll operations complete. Outputs are in {Path(args.output_dir).resolve()}")

if __name__ == "__main__":
    main()