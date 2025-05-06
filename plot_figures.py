#!/usr/bin/env python
"""
Script to load precomputed data and generate plots for opsin analysis.
This script focuses only on visualization without any additional processing.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import pickle
import argparse
from pathlib import Path

# Import visualization functions
from projects.opsin_analysis.visualization_functions import (
    # RMSD Visualizations
    visualize_rmsd_heatmap,
    create_and_visualize_similarity_tree,
    visualize_rmsd_matrix_improved,
    plot_rmsd_heatmap,
    plot_similarity_tree,
    create_rmsd_color_scale_figure,
    
    # Distance and structure visualizations
    plot_distances_with_std,
    plot_helix_logo_plots,
    plot_average_distances_by_helix,
    plot_distance_heatmap,
    
    # Residue composition and sequence plots
    create_residue_conservation_plot,
    print_residue_composition,
    
    # New conservation-based visualization
    plot_conservation_around_x50
)

# Import the color scheme tools
from projects.opsin_analysis.opsin_color_scheme import get_group_colors

def load_data(input_dir, output_dir=None, chain_id='A'):
    """
    Load precomputed data from workflow cache or CSV files.
    
    Args:
        input_dir: Directory containing the data files
        output_dir: Directory to save output files (defaults to input_dir)
        chain_id: Chain ID to use (default: 'A')
        
    Returns:
        Dictionary with loaded data
    """
    if output_dir is None:
        output_dir = input_dir
    
    # Make output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize data dictionary
    data = {}
    
    # Check for cache directory and files
    cache_dir = os.path.join(input_dir, 'cache')
    
    # First try to load from our new cache files
    if os.path.exists(cache_dir):
        print(f"Found cache directory: {cache_dir}")
        
        # Cache files to check in order of dependency - same as in run_new_visualizations.py
        files_to_load = [
            'raw_structures.pkl',
            f'processed_structures_{chain_id}.pkl',
            f'helix_annotations_{chain_id}.pkl',
            f'structure_comparison_{chain_id}.pkl',
            f'structure_errors_{chain_id}.pkl',
            f'grn_assignment_{chain_id}.pkl'
        ]
        
        # Try to load each cache file
        for file_name in files_to_load:
            cache_path = os.path.join(cache_dir, file_name)
            if os.path.exists(cache_path):
                try:
                    print(f"Loading {file_name}...")
                    with open(cache_path, 'rb') as f:
                        file_data = pickle.load(f)
                        
                        # Check if this is the structure_comparison file and examine keys
                        if file_name == f'structure_comparison_{chain_id}.pkl':
                            print(f"structure_comparison_{chain_id}.pkl keys: {list(file_data.keys())}")
                            # Check if the file contains RMSD data directly
                            rmsd_keys = [k for k in file_data.keys() if 'rmsd' in k.lower()]
                            if rmsd_keys:
                                print(f"Found RMSD-related keys in structure_comparison: {rmsd_keys}")
                                # If 'rmsd_matrix' is present, use it directly
                                if 'rmsd_matrix' in file_data:
                                    print("Using rmsd_matrix from structure_comparison directly")
                                    data['rmsd_df'] = file_data['rmsd_matrix']
                                    # Still update other keys from this file
                                    for key, value in file_data.items():
                                        if key != 'rmsd_matrix':
                                            data[key] = value
                                    continue  # Skip the general update below
                        
                        data.update(file_data)
                    print(f"Successfully loaded {file_name}")
                except Exception as e:
                    print(f"Error loading {file_name}: {e}")
            else:
                print(f"Warning: Cache file not found: {file_name}")
        
        # Skip loading rmsd_matrix.csv if we already have rmsd_df from the cache file
        if 'rmsd_df' in data:
            print("Using RMSD data from cache file")
            # Other CSV files that might contain necessary data
            csv_files = [
                ('ca_distance_table_grn.csv', 'ca_distance_table'),
                ('sidechain_distance_table_grn.csv', 'distance_table'),
                ('distance_table_grn.csv', 'distance_table_grn'),
                ('msa_table_grn.csv', 'msa_table'),
                ('ca_msa_table_grn.csv', 'ca_msa_table_grn'),
            ]
        else:
            print("No RMSD data found in cache files, checking CSV file")
            # Include rmsd_matrix.csv in files to load
            csv_files = [
                ('rmsd_matrix.csv', 'rmsd_df'),
                ('ca_distance_table_grn.csv', 'ca_distance_table'),
                ('sidechain_distance_table_grn.csv', 'distance_table'),
                ('distance_table_grn.csv', 'distance_table_grn'),
                ('msa_table_grn.csv', 'msa_table'),
                ('ca_msa_table_grn.csv', 'ca_msa_table_grn'),
            ]
        
        for file_name, key_name in csv_files:
            csv_path = os.path.join(output_dir, file_name)
            if os.path.exists(csv_path):
                print(f"Loading {file_name}...")
                try:
                    df = pd.read_csv(csv_path, index_col=0)
                    data[key_name] = df
                    print(f"Successfully loaded {file_name} into {key_name}")
                except Exception as e:
                    print(f"Error loading {file_name}: {e}")
        
        # Check if we loaded the key data we need
        if 'rmsd_df' in data:
            print(f"Successfully loaded RMSD data with {len(data['rmsd_df'])} structures")
            
            # Debug output
            if 'structure_errors' in data:
                print("\nDEBUG: Structure error information:")
                error_count = sum(1 for sid, err in data['structure_errors'].items() 
                                 if 'average_error' in err and err['average_error'] is not None)
                print(f"Found error data for {error_count} structures")
                
                # Check which structures have errors above threshold
                high_error_count = sum(1 for sid, err in data['structure_errors'].items() 
                                     if 'average_error' in err and 
                                     isinstance(err['average_error'], (int, float)) and 
                                     err['average_error'] > 3.0)
                print(f"Found {high_error_count} structures with error > 3.0 that should be filtered out")
                
                if high_error_count > 0:
                    print("\nStructures that should be filtered out:")
                    filtered_structs = []
                    for sid, err in data['structure_errors'].items():
                        if ('average_error' in err and 
                            isinstance(err['average_error'], (int, float)) and 
                            err['average_error'] > 3.0):
                            filtered_structs.append((sid, err['average_error']))
                    
                    # Show top structures by error
                    for sid, error in sorted(filtered_structs, key=lambda x: x[1], reverse=True)[:5]:
                        print(f"  {sid}: error={error:.2f}")
                
                # Show a few sample structures
                if error_count > 0:
                    print("\nSample structure errors:")
                    sample_size = min(3, error_count)
                    sample_structures = list(data['structure_errors'].keys())[:sample_size]
                    
                    for struct_id in sample_structures:
                        if 'average_error' in data['structure_errors'][struct_id]:
                            error = data['structure_errors'][struct_id]['average_error']
                            print(f"  {struct_id}: error={error}")
            
            # Return the loaded data
            return data, output_dir
        else:
            print("RMSD data not found in cache, falling back to CSV files")
    
    # Try to load error data from the file system if available
    mo_exp_errors_path = os.path.join(input_dir, 'mo_exp_errors.csv')
    hideaki_errors_path = os.path.join(input_dir, 'hideaki_errors.csv')
    errors_data = {}
    
    if os.path.exists(mo_exp_errors_path):
        try:
            print(f"Loading MO experiment errors from {mo_exp_errors_path}")
            mo_errors_df = pd.read_csv(mo_exp_errors_path)
            for _, row in mo_errors_df.iterrows():
                structure_id = row.get('structure_id')
                if structure_id:
                    # Calculate average error for this structure
                    error_cols = ['backbone_rmsd', 'pocket_rmsd', 'retinal_rmsd']
                    avg_error = row[error_cols].mean()
                    errors_data[structure_id] = {
                        'domain': row.get('domain', 'Unknown'),
                        'molecular_function': row.get('molecular_function', 'Unknown'),
                        'average_error': avg_error,
                        'backbone_rmsd': row.get('backbone_rmsd', 0),
                        'pocket_rmsd': row.get('pocket_rmsd', 0),
                        'retinal_rmsd': row.get('retinal_rmsd', 0)
                    }
            print(f"Loaded errors for {len(errors_data)} structures")
        except Exception as e:
            print(f"Error loading MO experiment errors: {e}")
    
    if os.path.exists(hideaki_errors_path):
        try:
            print(f"Loading Hideaki errors from {hideaki_errors_path}")
            hideaki_errors_df = pd.read_csv(hideaki_errors_path)
            for _, row in hideaki_errors_df.iterrows():
                structure_id = row.get('structure_id')
                if structure_id:
                    # Calculate average error for this structure
                    error_cols = ['backbone_rmsd', 'pocket_rmsd', 'retinal_rmsd']
                    avg_error = row[error_cols].mean()
                    errors_data[structure_id] = {
                        'domain': row.get('domain', 'Unknown'),
                        'molecular_function': row.get('molecular_function', 'Unknown'),
                        'average_error': avg_error,
                        'backbone_rmsd': row.get('backbone_rmsd', 0),
                        'pocket_rmsd': row.get('pocket_rmsd', 0),
                        'retinal_rmsd': row.get('retinal_rmsd', 0)
                    }
        except Exception as e:
            print(f"Error loading Hideaki errors: {e}")
    
    # Store error data in the main data dictionary
    if errors_data:
        data['properties'] = errors_data
    
    # If we couldn't load from cache, try the traditional CSV approach
    print("Looking for CSV files in the input directory...")
    
    # Load RMSD matrix
    rmsd_path = os.path.join(input_dir, 'rmsd_matrix.csv')
    if os.path.exists(rmsd_path):
        try:
            data['rmsd_df'] = pd.read_csv(rmsd_path, index_col=0)
            print(f"Loaded RMSD matrix from {rmsd_path}")
        except Exception as e:
            print(f"Error loading RMSD matrix: {e}")
    else:
        print(f"RMSD matrix file not found at {rmsd_path}")

    # Load molecular function data for coloring
    func_path = os.path.join(input_dir, 'molecular_functions.csv')
    if os.path.exists(func_path):
        try:
            func_df = pd.read_csv(func_path)
            data['group_dict'] = dict(zip(func_df['structure_id'], func_df['molecular_function']))
            print(f"Loaded molecular function data from {func_path}")
        except Exception as e:
            print(f"Error loading molecular function data: {e}")
            # Create default group dictionary based on structure names
            if 'rmsd_df' in data:
                data['group_dict'] = {}
                for sid in data['rmsd_df'].index:
                    # Try to extract function from name (e.g. "channel" in "CrChR2_channel")
                    if '_' in sid:
                        data['group_dict'][sid] = sid.split('_')[-1]
                    else:
                        data['group_dict'][sid] = "Unknown"
    else:
        print(f"Molecular function file not found at {func_path}")
        # Create default group dictionary
        if 'rmsd_df' in data:
            data['group_dict'] = {}
            for sid in data['rmsd_df'].index:
                # Try to extract function from name (common naming patterns)
                if 'ChR' in sid or 'channel' in sid.lower():
                    data['group_dict'][sid] = "Cation channel"
                elif 'HR' in sid or 'pump' in sid.lower() or 'PR' in sid:
                    data['group_dict'][sid] = "Proton pump"
                elif 'ACR' in sid or 'chloride' in sid.lower():
                    data['group_dict'][sid] = "Chloride pump"
                else:
                    data['group_dict'][sid] = "Unknown"

    # Load CA distance table
    distance_path = os.path.join(input_dir, 'ca_distance_table_grn.csv')
    if os.path.exists(distance_path):
        try:
            data['ca_distance_table'] = pd.read_csv(distance_path, index_col=0)
            print(f"Loaded CA distance table from {distance_path}")
        except Exception as e:
            print(f"Error loading CA distance table: {e}")
    else:
        print(f"CA distance table file not found at {distance_path}")
        # Look for pickle files as an alternative
        distance_pkl = os.path.join(input_dir, 'ca_distance_table_grn.pkl')
        if os.path.exists(distance_pkl):
            try:
                with open(distance_pkl, 'rb') as f:
                    data['ca_distance_table'] = pickle.load(f)
                print(f"Loaded CA distance table from {distance_pkl}")
            except Exception as e:
                print(f"Error loading CA distance table from pickle: {e}")
                
    # Load sidechain distance table
    sc_distance_path = os.path.join(input_dir, 'distance_table_grn.csv')
    if os.path.exists(sc_distance_path):
        try:
            data['distance_table'] = pd.read_csv(sc_distance_path, index_col=0)
            print(f"Loaded sidechain distance table from {sc_distance_path}")
        except Exception as e:
            print(f"Error loading sidechain distance table: {e}")
    else:
        print(f"Sidechain distance table file not found at {sc_distance_path}")
        # Look for pickle files as an alternative
        sc_distance_pkl = os.path.join(input_dir, 'distance_table_grn.pkl')
        if os.path.exists(sc_distance_pkl):
            try:
                with open(sc_distance_pkl, 'rb') as f:
                    data['distance_table'] = pickle.load(f)
                print(f"Loaded sidechain distance table from {sc_distance_pkl}")
            except Exception as e:
                print(f"Error loading sidechain distance table from pickle: {e}")

    # Load MSA table
    msa_path = os.path.join(input_dir, 'ca_msa_table_grn.csv')
    if os.path.exists(msa_path):
        try:
            data['msa_table'] = pd.read_csv(msa_path, index_col=0)
            # For compatibility with both formats
            data['msa_df'] = data['msa_table']
            print(f"Loaded MSA table from {msa_path}")
        except Exception as e:
            print(f"Error loading MSA table: {e}")
    else:
        print(f"MSA table file not found at {msa_path}")
        # Look for pickle files as an alternative
        msa_pkl = os.path.join(input_dir, 'ca_msa_table_grn.pkl')
        if os.path.exists(msa_pkl):
            try:
                with open(msa_pkl, 'rb') as f:
                    data['msa_table'] = pickle.load(f)
                    # For compatibility with both formats
                    data['msa_df'] = data['msa_table']
                print(f"Loaded MSA table from {msa_pkl}")
            except Exception as e:
                print(f"Error loading MSA table from pickle: {e}")
    
    # Load residue table (needed for residue logos)
    residue_path = os.path.join(input_dir, 'ca_residue_table_grn.csv')
    if os.path.exists(residue_path):
        try:
            data['residue_table'] = pd.read_csv(residue_path, index_col=0)
            print(f"Loaded residue table from {residue_path}")
        except Exception as e:
            print(f"Error loading residue table: {e}")
    else:
        print(f"Residue table file not found at {residue_path}")
        # Look for pickle files as an alternative
        residue_pkl = os.path.join(input_dir, 'ca_residue_table_grn.pkl')
        if os.path.exists(residue_pkl):
            try:
                with open(residue_pkl, 'rb') as f:
                    data['residue_table'] = pickle.load(f)
                print(f"Loaded residue table from {residue_pkl}")
            except Exception as e:
                print(f"Error loading residue table from pickle: {e}")
        # Fall back to using MSA table if available
        elif 'msa_table' in data:
            data['residue_table'] = data['msa_table'].copy()
            print("Using MSA table as fallback for residue table")
    
    return data, output_dir

def generate_summary_csv(data, output_dir):
    """
    Generate a 4-column CSV file with protein info:
    1. Protein name/display name
    2. Average RMSD
    3. Domain
    4. Molecular function
    
    Args:
        data: Dictionary with loaded data
        output_dir: Directory to save the CSV file
        
    Returns:
        Path to the generated CSV file
    """
    import csv
    
    # Create a list to hold the data
    csv_data = []
    
    # Extract protein info
    if 'rmsd_df' in data and 'processed_structures' in data:
        rmsd_df = data['rmsd_df']
        
        # Calculate average RMSD for each protein
        avg_rmsd = {}
        for protein_id in rmsd_df.index:
            # Calculate average RMSD to all other proteins
            rmsd_values = rmsd_df.loc[protein_id].values
            avg_rmsd[protein_id] = np.mean(rmsd_values[rmsd_values > 0])  # Exclude self-comparison (0)
        
        # First, identify high RMSD structures that should be excluded
        high_rmsd_structures = []
        for protein_id, avg in avg_rmsd.items():
            if avg > 3.0:
                high_rmsd_structures.append(protein_id)
                
        if high_rmsd_structures:
            print(f"Found {len(high_rmsd_structures)} structures with average RMSD > 3.0 that should be excluded from summary:")
            for pid in high_rmsd_structures[:5]:  # Show first 5
                print(f"  {pid}: average RMSD = {avg_rmsd[pid]:.2f}")
            if len(high_rmsd_structures) > 5:
                print(f"  ... and {len(high_rmsd_structures) - 5} more")
        
        # Get domain and function info for all proteins except high RMSD ones
        processed_proteins = set(rmsd_df.index) - set(high_rmsd_structures)
        for protein_id in processed_proteins:
            display_name = protein_id  # Default to ID if no display name
            
            # Try to get display name from processed_structures
            if protein_id in data['processed_structures']:
                if 'display_name' in data['processed_structures'][protein_id]:
                    display_name = data['processed_structures'][protein_id]['display_name']
            
            # Get domain and molecular function
            domain = "Unknown"
            function = "Unknown"
            
            if protein_id in data['processed_structures']:
                if 'properties' in data['processed_structures'][protein_id]:
                    properties = data['processed_structures'][protein_id]['properties']
                    domain = properties.get('domain', "Unknown")
                    function = properties.get('molecular_function', "Unknown")
            
            # Also get error if available
            error = "N/A"
            if 'structure_errors' in data and protein_id in data['structure_errors']:
                if 'average_error' in data['structure_errors'][protein_id]:
                    error_val = data['structure_errors'][protein_id]['average_error']
                    if error_val is not None:
                        error = round(error_val, 2)
            
            # Add to CSV data
            csv_data.append([
                display_name,
                round(avg_rmsd.get(protein_id, 0), 2),
                domain,
                function, 
                error
            ])
        
        print(f"Adding {len(csv_data)} structures to protein summary CSV (excluding high RMSD structures)")
    
    # Sort by average RMSD
    csv_data.sort(key=lambda x: x[1] if isinstance(x[1], (int, float)) else float('inf'))
    
    # Write to CSV file
    csv_path = os.path.join(output_dir, 'protein_summary.csv')
    with open(csv_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Write header
        csv_writer.writerow(['Protein', 'Average RMSD', 'Domain', 'Molecular Function', 'Error'])
        # Write data
        csv_writer.writerows(csv_data)
    
    # Save color assignments
    # Extract unique functions and domains
    functions = set()
    domains = set()
    for row in csv_data:
        functions.add(row[3])  # Index 3 is molecular function
        domains.add(row[2])    # Index 2 is domain
    
    # Get color assignments
    from projects.opsin_analysis.opsin_color_scheme import get_group_colors
    function_colors = get_group_colors(sorted(functions))
    domain_colors = get_group_colors(sorted(domains))
    
    # Write color assignments to CSV
    color_path = os.path.join(output_dir, 'color_assignments.csv')
    with open(color_path, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        # Write header
        csv_writer.writerow(['Type', 'Value', 'Color'])
        # Write function colors
        for func, color in function_colors.items():
            csv_writer.writerow(['Function', func, color])
        # Write domain colors
        for domain, color in domain_colors.items():
            csv_writer.writerow(['Domain', domain, color])
    
    print(f"Saved protein summary CSV to {csv_path}")
    print(f"Saved color assignments to {color_path}")
    return csv_path

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram

def generate_plots(data, output_dir):
    """
    Generate all plots and save them to the output directory.
    Ensures consistent ordering between RMSD heatmap and similarity tree.

    Args:
        data: Dictionary with loaded data
        output_dir: Directory to save output files
    """
    # --- Standard Plots (Unaffected by RMSD/Tree linking) ---

    # Generate Overview plot
    print("Generating overview plot...")
    try:
        from projects.opsin_analysis.visualization_plots import create_overview_plot
        create_overview_plot(data, os.path.join(output_dir, 'opsin_overview.png'))
    except Exception as e:
        print(f"Error generating overview plot: {e}")
        import traceback
        traceback.print_exc()

    # Generate RMSD color scale figure first
    print("Generating RMSD color scale reference...")
    try:
        color_scale_fig = create_rmsd_color_scale_figure()
        color_scale_path = os.path.join(output_dir, 'rmsd_color_scale.png')
        color_scale_fig.savefig(color_scale_path, dpi=300, bbox_inches='tight')
        plt.close(color_scale_fig)
        print(f"Saved RMSD color scale to {color_scale_path}")
    except Exception as e:
        print(f"Error generating RMSD color scale: {e}")

    # Generate simple similarity tree (unlinked to the main heatmap)
    if 'rmsd_df' in data:
        print("Generating simple similarity tree...")
        try:
            # Use the original, unfiltered data for the simple tree
            fig_simple = plot_similarity_tree(
                data['rmsd_df'],
                title="Simple Structural Similarity Tree (Unfiltered)"
            )
            output_path_simple = os.path.join(output_dir, 'simple_similarity_tree.png')
            fig_simple.savefig(output_path_simple, dpi=300, bbox_inches='tight')
            plt.close(fig_simple)
            print(f"Saved simple similarity tree to {output_path_simple}")
        except Exception as e:
            print(f"Error generating simple similarity tree: {e}")
            import traceback
            traceback.print_exc()

    # Generate CA distance plot
    if 'ca_distance_table' in data:
        print("Generating CA distance plot...")
        try:
            fig = plot_distances_with_std(
                data['ca_distance_table'],
                title="CA Distance to Retinal by Position",
                figsize=(18, 10),
                use_ca=True
            )
            output_path = os.path.join(output_dir, 'ca_distance_plot.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved CA distance plot to {output_path}")
        except Exception as e:
            print(f"Error generating CA distance plot: {e}")

    # Generate sidechain distance plot
    if 'distance_table' in data:
        print("Generating sidechain distance plot...")
        try:
            fig = plot_distances_with_std(
                data['distance_table'],
                title="Sidechain Distance to Retinal by Position",
                figsize=(18, 10),
                use_ca=False
            )
            output_path = os.path.join(output_dir, 'sidechain_distance_plot.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved sidechain distance plot to {output_path}")
        except Exception as e:
            print(f"Error generating sidechain distance plot: {e}")

    # Generate helix logo plots
    if 'residue_table' in data:
        print("Generating helix logo plots...")
        try:
            fig = plot_helix_logo_plots(
                data['residue_table'],
                figsize=(20, 5)
            )
            output_path = os.path.join(output_dir, 'helix_logo_plots.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved helix logo plots to {output_path}")
        except Exception as e:
            print(f"Error generating helix logo plots: {e}")

    # Generate average CA distances by helix
    if 'ca_distance_table' in data:
        print("Generating average CA distances by helix...")
        try:
            fig = plot_average_distances_by_helix(data['ca_distance_table'], use_ca=True)
            output_path = os.path.join(output_dir, 'ca_average_distances_by_helix.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved average CA distances by helix to {output_path}")
        except Exception as e:
            print(f"Error generating average CA distances by helix: {e}")

    # Generate average sidechain distances by helix
    if 'distance_table' in data:
        print("Generating average sidechain distances by helix...")
        try:
            fig = plot_average_distances_by_helix(data['distance_table'], use_ca=False)
            output_path = os.path.join(output_dir, 'sidechain_average_distances_by_helix.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved average sidechain distances by helix to {output_path}")
        except Exception as e:
            print(f"Error generating average sidechain distances by helix: {e}")

    # Generate CA distance heatmap
    if 'ca_distance_table' in data:
        print("Generating CA distance heatmap...")
        try:
            fig = plot_distance_heatmap(data['ca_distance_table'])
            output_path = os.path.join(output_dir, 'ca_distance_heatmap.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved CA distance heatmap to {output_path}")
        except Exception as e:
            print(f"Error generating CA distance heatmap: {e}")

    # Generate sidechain distance heatmap
    if 'distance_table' in data:
        print("Generating sidechain distance heatmap...")
        try:
            fig = plot_distance_heatmap(data['distance_table'])
            output_path = os.path.join(output_dir, 'sidechain_distance_heatmap.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved sidechain distance heatmap to {output_path}")
        except Exception as e:
            print(f"Error generating sidechain distance heatmap: {e}")

    # Generate residue conservation plot
    if 'msa_table' in data:
        print("Generating residue conservation plot...")
        try:
            fig = create_residue_conservation_plot(
                data['msa_table'],
                helix_highlighting=True,
                figsize=(14, 8)
            )
            output_path = os.path.join(output_dir, 'residue_conservation.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved residue conservation plot to {output_path}")
        except Exception as e:
            print(f"Error generating residue conservation plot: {e}")

    # Generate conservation-based logo plots around X.50 positions
    if 'residue_table' in data:
        print("Generating conservation logo plots around X.50 positions...")
        try:
            fig = plot_conservation_around_x50(
                data['residue_table'],
                figsize=(20, 10)
            )
            output_path = os.path.join(output_dir, 'conservation_around_x50.png')
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
            plt.close(fig)
            print(f"Saved conservation logo plots around X.50 to {output_path}")
        except Exception as e:
            print(f"Error generating conservation logo plots around X.50: {e}")


    # --- Linked RMSD Matrix and Similarity Tree Section ---
    if 'rmsd_df' in data:
        print("\nProcessing RMSD matrix and similarity tree with consistent ordering...")

        # --- 1. Centralized Filtering ---
        rmsd_df = data['rmsd_df']
        original_pdb_list = rmsd_df.index.tolist()

        # Prepare domain_dict and group_dict, ensuring they cover all original IDs
        domain_dict = {} # Structure: {'domain': 'X', 'average_error': Y} or {'domain': 'X'}
        group_dict = {}  # Structure: 'molecular_function' (string)

        # Populate domain_dict and group_dict from processed_structures and structure_errors
        if 'processed_structures' in data:
            print("DEBUG: Populating group/domain info from processed_structures...")
            count_props = 0
            count_errors = 0
            for struct_id, struct_data in data['processed_structures'].items():
                if struct_id in original_pdb_list: # Only consider structures in the RMSD matrix
                    properties = struct_data.get('properties', {})
                    domain = properties.get('domain', 'Unknown')
                    function = properties.get('molecular_function', 'Unknown')
                    count_props += 1

                    # Initialize domain entry as a dictionary
                    domain_entry = {'domain': str(domain)} # Ensure string

                    # Add error data if available from structure_errors cache
                    if 'structure_errors' in data and struct_id in data['structure_errors']:
                        error_val = data['structure_errors'][struct_id].get('average_error')
                        if error_val is not None and isinstance(error_val, (int, float)):
                            domain_entry['average_error'] = float(error_val)
                            count_errors += 1

                    domain_dict[struct_id] = domain_entry
                    group_dict[struct_id] = str(function) # Ensure string
            print(f"DEBUG: Found properties for {count_props} structures, error values for {count_errors} structures.")

        # Fallback for missing group/domain info
        print("DEBUG: Applying fallbacks for missing group/domain info...")
        unknown_group_count = 0
        unknown_domain_count = 0
        for sid in original_pdb_list:
            if sid not in group_dict:
                group_dict[sid] = "Unknown"
                unknown_group_count += 1
            if sid not in domain_dict:
                 # Ensure domain_dict entries are always dictionaries
                domain_dict[sid] = {'domain': "Unknown"}
                unknown_domain_count +=1
            elif not isinstance(domain_dict[sid], dict):
                 # Convert simple string domain entries to dict format if they somehow occurred
                 print(f"DEBUG: Converting non-dict domain entry for {sid}")
                 domain_dict[sid] = {'domain': str(domain_dict[sid])} # Ensure string domain
        print(f"DEBUG: Assigned 'Unknown' group to {unknown_group_count}, 'Unknown' domain to {unknown_domain_count}.")


        # Define error threshold
        error_threshold = 3.0
        print(f"Using error/RMSD threshold for filtering: {error_threshold} Å")

        # Filter based on average RMSD calculated from the matrix itself OR explicit error value
        try:
            average_rmsds = rmsd_df.mean()
        except Exception as avg_e:
             print(f"ERROR calculating average RMSDs: {avg_e}. Cannot perform filtering based on matrix average.")
             average_rmsds = pd.Series(dtype=float) # Empty series


        kept_ids = []
        filtered_out_details = []

        print("DEBUG: Starting filtering loop...")
        filter_value_sources = {'error_dict': 0, 'avg_rmsd': 0, 'skipped': 0}

        for pdb_id in original_pdb_list:
            # Check explicit error from domain_dict first
            error_val = None
            if isinstance(domain_dict.get(pdb_id), dict):
                 error_val = domain_dict[pdb_id].get('average_error') # Returns None if key missing

            filter_value = None
            source = 'skipped'

            if error_val is not None:
                filter_value = error_val
                source = 'error_dict'
                filter_value_sources['error_dict'] += 1
            elif pdb_id in average_rmsds.index:
                filter_value = average_rmsds.loc[pdb_id]
                source = 'avg_rmsd'
                filter_value_sources['avg_rmsd'] += 1
            else:
                 filter_value_sources['skipped'] += 1
                 print(f"DEBUG: Skipping {pdb_id}, no error value and not in average_rmsds.")
                 # Decide whether to keep or filter skipped ones, let's filter them
                 filter_value = error_threshold + 1 # Ensure it gets filtered

            # Perform the check
            if filter_value is not None and filter_value <= error_threshold:
                kept_ids.append(pdb_id)
                # Optional: print kept structures for debug
                # print(f"DEBUG: Keeping {pdb_id} (value={filter_value:.2f}, source={source})")
            elif filter_value is not None:
                filtered_out_details.append(f"  - {pdb_id}: value={filter_value:.2f} Å (source: {source})")
                # Optional: print filtered structures for debug
                # print(f"DEBUG: Filtering {pdb_id} (value={filter_value:.2f}, source={source})")


        print(f"DEBUG: Filtering sources: {filter_value_sources}")

        if filtered_out_details:
             print(f"Filtering out {len(filtered_out_details)} structures (value > {error_threshold}):")
             # Sort details for clearer output
             filtered_out_details.sort(key=lambda x: float(x.split('=')[1].split(' ')[0]), reverse=True)
             print("\n".join(filtered_out_details[:15])) # Show up to 15
             if len(filtered_out_details) > 15: print(f"  ... and {len(filtered_out_details) - 15} more")
        else:
             print("No structures filtered based on threshold.")

        # --- Proceed only if enough structures remain ---
        if len(kept_ids) < 2:
            print("\nNot enough structures remain after filtering (< 2). Skipping linked RMSD heatmap and similarity tree.")
        else:
            print(f"\nProceeding with {len(kept_ids)} structures for linked visualizations.")
            filtered_rmsd_df = rmsd_df.loc[kept_ids, kept_ids]

            # Filter group_dict and domain_dict to only include kept IDs
            filtered_group_dict = {k: v for k, v in group_dict.items() if k in kept_ids}
            filtered_domain_dict = {k: v for k, v in domain_dict.items() if k in kept_ids}
            print(f"DEBUG: Filtered group_dict size: {len(filtered_group_dict)}")
            print(f"DEBUG: Filtered domain_dict size: {len(filtered_domain_dict)}")


            # --- 2. Calculate Linkage Once ---
            Z_linkage = None # Initialize
            try:
                print("Calculating linkage matrix Z...")
                # Ensure matrix is clean for linkage calculation
                # Use fillna with a high value (e.g., mean of existing values + std dev or just a fixed high number)
                # This is generally safer than fillna(0) for distance matrices
                matrix_for_linkage = filtered_rmsd_df.values
                if np.any(~np.isfinite(matrix_for_linkage)):
                    mean_finite = np.nanmean(matrix_for_linkage[np.isfinite(matrix_for_linkage)])
                    std_finite = np.nanstd(matrix_for_linkage[np.isfinite(matrix_for_linkage)])
                    fill_val = mean_finite + 3*std_finite if not np.isnan(mean_finite) else 10.0 # Fallback fill value
                    print(f"DEBUG: Filling NaNs/Infs with {fill_val:.2f}")
                    matrix_for_linkage = np.nan_to_num(matrix_for_linkage, nan=fill_val, posinf=fill_val*1.1, neginf=0.0)
                else:
                    matrix_for_linkage = matrix_for_linkage.copy() # Ensure it's a copy

                np.fill_diagonal(matrix_for_linkage, 0.0) # Essential

                # Check for negative values after cleaning
                if np.any(matrix_for_linkage < 0):
                     print("WARNING: Negative values found in distance matrix after cleaning. Clipping to zero.")
                     matrix_for_linkage[matrix_for_linkage < 0] = 0.0

                # Check symmetry
                if not np.allclose(matrix_for_linkage, matrix_for_linkage.T):
                    print("WARNING: Matrix is not symmetric after cleaning. Forcing symmetry.")
                    matrix_for_linkage = (matrix_for_linkage + matrix_for_linkage.T) / 2.0
                    np.fill_diagonal(matrix_for_linkage, 0.0) # Re-ensure diagonal is zero


                condensed_matrix = squareform(matrix_for_linkage, checks=False) # Disable checks for robustness
                # Use the same method as the tree function ('average')
                Z_linkage = linkage(condensed_matrix, method='average')
                print("Successfully calculated linkage matrix Z.")

            except Exception as e:
                 print(f"ERROR calculating linkage matrix: {e}")
                 import traceback
                 traceback.print_exc()
                 Z_linkage = None # Ensure it's None if calculation failed

            # --- Proceed only if Linkage Matrix was calculated ---
            if Z_linkage is not None:
                # --- 3. Generate Similarity Tree (using Z) ---
                print("\nGenerating similarity tree (using calculated Z)...")
                try:
                    tree_output_path = os.path.join(output_dir, 'similarity_tree_linked.png')
                    tree_fig, ordered_tree_ids = create_and_visualize_similarity_tree(
                        rmsd_data=filtered_rmsd_df, # Pass filtered df
                        linkage_matrix=Z_linkage,    # Pass pre-calculated Z
                        group_dict=filtered_group_dict,
                        domain_dict=filtered_domain_dict,
                        # No error_threshold needed here
                    )
                    # Add saving logic for tree components if needed
                    tree_fig.savefig(tree_output_path, dpi=300, bbox_inches='tight')

                    # Save separate components for the tree
                    # 1. Save content-only version (no legends)
                    content_path_tree = tree_output_path.replace('.png', '_content_only.png')
                    # Temporarily remove legends for saving
                    original_legends_tree = list(tree_fig.legends) # Store legends
                    for legend in original_legends_tree:
                        legend.remove()
                    tree_fig.savefig(content_path_tree, dpi=300, bbox_inches='tight')
                    # Restore legends if needed later, or just close the figure
                    # for legend in original_legends_tree: tree_fig.legends.append(legend) # Restore (optional)

                    # 2. Create and save legends-only version for the tree
                    legend_path_tree = tree_output_path.replace('.png', '_legends_only.png')
                    legend_fig_tree = plt.figure(figsize=(6, 4)) # Adjust size as needed
                    legend_ax_tree = legend_fig_tree.add_subplot(111)
                    legend_ax_tree.axis('off')

                    # Recreate function legend (using filtered data)
                    if filtered_group_dict:
                         unique_funcs_filt = sorted(set(filtered_group_dict.values()))
                         func_colors_filt = get_group_colors(unique_funcs_filt)
                         func_handles = [plt.Rectangle((0, 0), 1, 1, facecolor=func_colors_filt.get(f, '#CCCCCC'), edgecolor='none', label=f) for f in unique_funcs_filt]
                         leg1_tree = legend_ax_tree.legend(handles=func_handles, title="Molecular Function", loc='upper center', fontsize=10)
                         legend_ax_tree.add_artist(leg1_tree)

                    # Recreate domain legend (using filtered data)
                    if filtered_domain_dict:
                        domain_values_filt = [str(d.get('domain', 'Unknown') if isinstance(d, dict) else d) for d in filtered_domain_dict.values()]
                        unique_doms_filt = sorted(set(domain_values_filt))
                        dom_colors_filt = get_group_colors(unique_doms_filt)
                        dom_handles = [plt.Rectangle((0, 0), 1, 1, facecolor=dom_colors_filt.get(d, '#CCCCCC'), edgecolor='none', label=d) for d in unique_doms_filt]
                        # Place domain legend below function legend
                        leg2_tree = legend_ax_tree.legend(handles=dom_handles, title="Domain", loc='lower center', fontsize=10)
                        legend_ax_tree.add_artist(leg2_tree) # Ensure second legend is added correctly

                    legend_fig_tree.tight_layout()
                    legend_fig_tree.savefig(legend_path_tree, dpi=300, bbox_inches='tight')
                    plt.close(legend_fig_tree)

                    # Close the main tree figure
                    plt.close(tree_fig)
                    print(f"Saved linked similarity tree (and components) to {tree_output_path}")
                    print(f"  Ordered IDs from tree: {len(ordered_tree_ids)}")

                    # --- 4. Generate Improved RMSD Matrix (using Z) ---
                    print("\nGenerating improved RMSD matrix (using calculated Z)...")
                    try:
                        output_path_heatmap = os.path.join(output_dir, 'rmsd_matrix_improved_linked.png')
                        heatmap_clustermap = visualize_rmsd_matrix_improved(
                            rmsd_df=filtered_rmsd_df,    # Pass filtered df
                            linkage_matrix=Z_linkage,     # Pass pre-calculated Z REQUIRED
                            group_dict=filtered_group_dict,
                            domain_dict=filtered_domain_dict,
                            output_file=output_path_heatmap, # Function saves internally now
                            # No error_threshold needed here
                            figsize=(14, 12)
                        )
                        if heatmap_clustermap:
                             # Figure is saved within the function, just need to close it
                             plt.close(heatmap_clustermap.fig)
                             print(f"Saved linked improved RMSD matrix (and components) to {output_path_heatmap}")
                        else:
                             print("Skipped saving RMSD matrix (clustermap generation failed or returned None).")

                    except Exception as e:
                        print(f"Error generating improved RMSD matrix: {e}")
                        import traceback
                        traceback.print_exc()

                except Exception as e:
                    print(f"Error generating similarity tree: {e}")
                    import traceback
                    traceback.print_exc()
            else:
                 print("Skipping linked RMSD/Tree plots because Linkage Matrix calculation failed.")

    # --- End Linked RMSD Matrix and Similarity Tree Section ---

    print(f"\nAll requested plots generated and saved to {output_dir}")


def main():
    parser = argparse.ArgumentParser(description='Generate plots for opsin analysis')
    parser.add_argument('--input-dir', '-i', type=str, default='projects/opsin_analysis/opsin_output',
                        help='Directory containing input data files')
    parser.add_argument('--output-dir', '-o', type=str, default='projects/opsin_analysis/opsin_output',
                        help='Directory to save output plots (defaults to input-dir)')
    parser.add_argument('--chain-id', '-c', type=str, default='A',
                        help='Chain ID used in the analysis (default: A)')
    parser.add_argument('--quality', '-q', type=str, choices=['low', 'medium', 'high'], default='high',
                        help='Figure quality (affects DPI and size)')
    args = parser.parse_args()
    
    # Load data
    print(f"Loading data from {args.input_dir} for chain {args.chain_id}...")
    data, output_dir = load_data(args.input_dir, args.output_dir, args.chain_id)
    
    # Create figures directory
    figures_dir = os.path.join(output_dir, 'figures')
    os.makedirs(figures_dir, exist_ok=True)
    
    # Set figure quality parameters
    if args.quality == 'high':
        plt.rcParams['figure.dpi'] = 300
        plt.rcParams['savefig.dpi'] = 300
    elif args.quality == 'medium':
        plt.rcParams['figure.dpi'] = 200
        plt.rcParams['savefig.dpi'] = 200
    else:  # low
        plt.rcParams['figure.dpi'] = 100
        plt.rcParams['savefig.dpi'] = 100
    
    print(f"Using {args.quality} quality for figures")
        
    # Generate plots
    if data:
        print("Generating plots...")
        generate_plots(data, figures_dir)
        
        # Generate protein summary CSV
        print("Generating protein summary CSV...")
        generate_summary_csv(data, figures_dir)
        
        print(f"All plots saved to {figures_dir}")
    else:
        print("No data loaded, cannot generate plots.")
        print("Make sure you've run the workflow first with:")
        print(f"  python -m projects.opsin_analysis.opsin_analysis_workflow --output-dir {args.input_dir} --chain-id {args.chain_id}")

def verify_filtering(output_dir):
    """
    Verify that the RMSD visualizations don't contain structures with high RMSD values.
    This is a sanity check to run after generating visualizations.
    
    Args:
        output_dir: Directory containing the output files
    """
    import os
    import pandas as pd
    
    # Check the protein summary CSV to verify filtering worked
    summary_path = os.path.join(output_dir, 'protein_summary.csv')
    if os.path.exists(summary_path):
        print("\nVerifying filtering in output files:")
        df = pd.read_csv(summary_path)
        
        # Look for any proteins with average RMSD > 3.0
        if 'Average RMSD' in df.columns:
            high_rmsd = df[df['Average RMSD'] > 3.0]
            if high_rmsd.empty:
                print("  ✓ No structures with average RMSD > 3.0 found in protein summary")
            else:
                print(f"  ✗ Found {len(high_rmsd)} structures with average RMSD > 3.0 in protein summary:")
                for _, row in high_rmsd.iterrows():
                    print(f"    {row['Protein']}: RMSD = {row['Average RMSD']}")
    
    print("\nFiltering verification complete. If any high-RMSD structures were found,")
    print("check that the visualizations correctly show the filtered data.")
    
if __name__ == "__main__":
    main()
    # After generating all visualizations, verify filtering worked
    print("\n==== VERIFICATION STEP ====")
    verify_filtering(os.path.join('projects', 'opsin_analysis', 'opsin_output', 'figures'))