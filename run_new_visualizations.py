#!/usr/bin/env python
"""
Script to run the new visualization_plots.py functionality.
This replaces the deprecated plot_figures.py script.
"""
import os
import pickle
import pandas as pd
from visualization_plots import *

def load_analysis_data(cache_dir, output_dir):
    """Load analysis data from cache files and CSV files"""
    print(f"Loading analysis data from cache: {cache_dir}")
    
    # Gather cache files
    files_to_load = [
        'raw_structures.pkl',
        'processed_structures_A.pkl',
        'helix_annotations_A.pkl',
        'structure_comparison_A.pkl',
        'structure_errors_A.pkl',
        'grn_assignment_A.pkl'
    ]
    
    # Load data from cache files
    data = {}
    for file_name in files_to_load:
        cache_path = os.path.join(cache_dir, file_name)
        if os.path.exists(cache_path):
            print(f"Loading {file_name}...")
            try:
                with open(cache_path, 'rb') as f:
                    file_data = pickle.load(f)
                    data.update(file_data)
                print(f"Successfully loaded {file_name}")
            except Exception as e:
                print(f"Error loading {file_name}: {e}")
        else:
            print(f"Warning: Cache file not found: {file_name}")
    
    # Check for and load any CSV files that might contain necessary data
    csv_files = [
        ('rmsd_matrix.csv', 'rmsd_df'),
        ('ca_distance_table_grn.csv', 'ca_distance_table_grn'),
        ('sidechain_distance_table_grn.csv', 'sidechain_distance_table_grn'),
        ('distance_table_grn.csv', 'distance_table_grn'),
        ('msa_table_grn.csv', 'msa_table_grn'),  # Add MSA table for helix logo plots
        ('ca_msa_table_grn.csv', 'ca_msa_table_grn'),  # Add CA MSA table as fallback
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
    
    # If we have a raw RMSD matrix but no rmsd_df in the data
    if 'rmsd_matrix.csv' in [os.path.basename(f) for f in os.listdir(output_dir) if os.path.isfile(os.path.join(output_dir, f))] and 'rmsd_df' not in data:
        rmsd_path = os.path.join(output_dir, 'rmsd_matrix.csv')
        try:
            print(f"Loading RMSD matrix from {rmsd_path}")
            rmsd_df = pd.read_csv(rmsd_path, index_col=0)
            data['rmsd_df'] = rmsd_df
            print(f"Successfully loaded RMSD matrix")
            
            # Debug output
            print(f"RMSD matrix shape: {rmsd_df.shape}")
            print(f"RMSD matrix has {len(rmsd_df.index)} structures")
            print(f"RMSD matrix error threshold analysis:")

            # If we have structure_errors
            if 'structure_errors' in data:
                errors_by_structure = {}
                for struct_id in rmsd_df.index:
                    # Check if structure is in errors dict
                    if struct_id in data['structure_errors']:
                        error = data['structure_errors'][struct_id].get('average_error', float('inf'))
                        errors_by_structure[struct_id] = error
                
                if errors_by_structure:
                    error_values = list(errors_by_structure.values())
                    below_threshold = sum(1 for e in error_values if e <= 3.0)
                    print(f"  {below_threshold} of {len(error_values)} structures have error <= 3.0")
                    print(f"  Min error: {min(error_values)}, Max error: {max(error_values)}")
                    
                    # If no structures are below threshold, adjust it
                    if below_threshold == 0:
                        print("WARNING: No structures meet the error threshold of 3.0.")
                        print("  Suggesting a higher threshold based on data distribution.")
                        # Find a good threshold that would include ~20% of structures
                        sorted_errors = sorted(error_values)
                        idx = min(len(sorted_errors) // 5, len(sorted_errors) - 1)
                        suggested_threshold = sorted_errors[idx]
                        print(f"  Suggested threshold: {suggested_threshold}")
                        # Store the suggested threshold in data
                        data['suggested_error_threshold'] = suggested_threshold
        except Exception as e:
            print(f"Error processing RMSD matrix: {e}")
    
    return data

def main():
    # Set paths
    base_dir = os.path.dirname(__file__)
    output_dir = os.path.join(base_dir, 'opsin_output')
    cache_dir = os.path.join(output_dir, 'cache')
    figures_dir = os.path.join(output_dir, 'figures')
    
    # Create directories if needed
    os.makedirs(figures_dir, exist_ok=True)
    
    # Load analysis data
    data = load_analysis_data(cache_dir, output_dir)
    
    if not data:
        print("No data loaded from cache files. Cannot proceed.")
        return
    
    # Set RMSD filtering threshold
    error_threshold = 4.0  # Use 4.0 as specified in your comment
    print(f"Using RMSD filtering threshold: {error_threshold}")
    
    # Generate visualizations
    print(f"Generating visualizations in: {figures_dir}")
    
    # Debug the structure error values
    if 'structure_errors' in data:
        print("\nDEBUG: Examining structure errors:")
        error_values = []
        for struct_id, error_data in data['structure_errors'].items():
            if 'average_error' in error_data and error_data['average_error'] is not None:
                error_values.append((struct_id, error_data['average_error']))
        
        if error_values:
            sorted_errors = sorted(error_values, key=lambda x: x[1])
            print(f"  Total structures with error data: {len(sorted_errors)}")
            print(f"  Min error: {sorted_errors[0][1]} (ID: {sorted_errors[0][0]})")
            print(f"  Max error: {sorted_errors[-1][1]} (ID: {sorted_errors[-1][0]})")
            print(f"  Structures with error ≤ 3.0: {sum(1 for _, err in sorted_errors if err <= 3.0)}")
            print(f"  Structures with error ≤ 8.0: {sum(1 for _, err in sorted_errors if err <= 8.0)}")
            print(f"  Structures with error ≤ 15.0: {sum(1 for _, err in sorted_errors if err <= 15.0)}")
            
            # Find a good threshold that includes at least some structures
            for threshold in [5.0, 8.0, 10.0, 15.0, 20.0]:
                count = sum(1 for _, err in sorted_errors if err <= threshold)
                if count >= 10:  # At least 10 structureplot_helix_logo_plotss
                    print(f"  Recommended threshold: {threshold} (includes {count} structures)")
                    error_threshold = threshold
                    break
    
    # Generate standard visualizations
    figures = generate_visualizations(
        data, 
        output_dir=figures_dir,
        visualize_ca=True,
        visualize_sidechain=True,
        error_threshold=error_threshold
    )
    
    # Generate helix logo plots with 15% frequency threshold
    logo_plot_path = os.path.join(figures_dir, 'helix_logo_plots_2.png')
    print(f"Generating logo plots with 15% minimum frequency threshold: {logo_plot_path}")
    if 'msa_table_grn' in data:
        logo_fig = plot_helix_logo_plots(
            data['msa_table_grn'], 
            output_path=logo_plot_path,
            min_frequency=0.15
        )
        print(f"Successfully generated filtered logo plots")
        # Check if figures is a list before appending
        if isinstance(figures, list):
            figures.append(logo_plot_path)
        else:
            print(f"Warning: figures is a {type(figures).__name__}, not a list")
    
    print(f"Successfully generated {len(figures)} visualizations")
    print(f"Visualization files saved to: {figures_dir}")

if __name__ == "__main__":
    main()

