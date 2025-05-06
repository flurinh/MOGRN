"""
New visualization script for opsin analysis that creates standardized plots.
This replaces the deprecated plot_figures.py script.
"""
import os
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, dendrogram
import logomaker
from matplotlib.colors import BoundaryNorm
from matplotlib.patches import Patch
from collections import Counter
import matplotlib.patches as mpatches

# Import global color scheme and visualization utilities
from projects.opsin_analysis.opsin_color_scheme import (
    OPSIN_COLORS, RMSD_CMAP, DIVERGING_CMAP, RMSD_DISCRETE_CMAP, RMSD_BOUNDS,
    HELIX_COLORS, HELIX_COLORS_STR, HELIX_COLORS_LIST, get_group_colors)

def _save_separate_content_and_legends(fig, output_path, unique_funcs, function_colors_dict, 
                                     function_map, unique_domains, domain_color_map, 
                                     title="Visualization"):
    """Helper function to save separate content and legend files for a figure"""
    # Create and save a version without legends (content only)
    content_path = output_path.replace('.png', '_content_only.png')
    
    # Hide all legends from the figure
    for legend in fig.legends:
        legend.remove()
    
    # Hide any annotation text
    for text in fig.texts:
        text.set_visible(False)
        
    # Save content-only version
    fig.savefig(content_path, dpi=300, bbox_inches='tight')
    print(f"Saved content-only version to {os.path.basename(content_path)}")
    
    # Create and save legends-only figure
    legends_path = output_path.replace('.png', '_legends_only.png')
    legend_fig = plt.figure(figsize=(10, 6))
    
    # Recreate all legends on a clean figure
    legend_ax = legend_fig.add_subplot(111)
    legend_ax.axis('off')
    
    # Function legend
    function_handles = [
        plt.Line2D([0],[0], color=function_colors_dict.get(f, OPSIN_COLORS['gray_light']),
                  marker='o', label=function_map.get(f, f.replace('_', ' ')), 
                  markersize=10, linestyle='None')
        for f in unique_funcs
    ]
    
    # Domain legend
    domain_handles = [
        plt.Line2D([0],[0], color=domain_color_map.get(d, OPSIN_COLORS['gray_light']),
                  marker='o', label=d, markersize=10, linestyle='None')
        for d in unique_domains
    ]
    
    # Add title
    legend_fig.suptitle(f"Legend for {title}", fontsize=14, y=0.95)
    
    # Create separate legends with better spacing
    if function_handles:
        legend1 = legend_ax.legend(handles=function_handles, 
                                 loc='upper left', 
                                 title="Molecular Function",
                                 fontsize=12)
        legend_ax.add_artist(legend1)  # Add first legend
    
    if domain_handles:
        legend2 = legend_ax.legend(handles=domain_handles, 
                                 loc='upper right', 
                                 title="Domain",
                                 fontsize=12)
        legend_ax.add_artist(legend2)  # Add second legend
    
    # Save legends-only figure
    legend_fig.tight_layout()
    legend_fig.savefig(legends_path, dpi=300, bbox_inches='tight')
    plt.close(legend_fig)
    print(f"Saved legends-only version to {os.path.basename(legends_path)}")


def create_overview_plot(data, output_path=None, figsize=(16, 14)):
    """
    Create a circular overview visualization of the dataset with 4 rings:
    1. Inner ring: molecular_function (color coded)
    2. Second ring: domain (color coded)
    3. Third ring: experimental and predicted structure availability (dots)
    4. Outer ring: display name of each structure
    
    Args:
        data: Dictionary containing analysis results
        output_path: Path to save the figure (optional)
        figsize: Figure size tuple (width, height)
        
    Returns:
        Matplotlib figure
    """
    # Debug: Understand the structure_mapping format
    print("Debug: Structure mapping analysis")
    print(f"Structure mapping type: {type(data.get('structure_mapping', {}))}")
    if 'structure_mapping' in data and data['structure_mapping']:
        sample_key = next(iter(data['structure_mapping'].keys()))
        sample_value = data['structure_mapping'][sample_key]
        print(f"Sample key: {sample_key}, type: {type(sample_key)}")
        print(f"Sample value: {sample_value}, type: {type(sample_value)}")
        
    # Extract the structure data and mapping from cache
    processed_structures = data.get('processed_structures', {})
    structure_mapping = data.get('structure_mapping', {})
    
    # Create a DataFrame for plotting
    plot_data = []
    
    for struct_id, struct_data in processed_structures.items():
        properties = struct_data.get('properties', {})
        
        # Check if this structure has an experimental counterpart
        has_exp_structure = False
        
        # Based on the debug output, the structure mapping is: {exp_id: pred_id}
        # We need to check if this structure is a PREDICTED structure (appears as a value)
        if struct_id in structure_mapping.values():
            has_exp_structure = True
            
        # Get short display name
        short_name = struct_data.get('display_name', struct_id)
        if len(short_name) > 15:
            short_name = short_name[:12] + '...'
            
        # Get molecular function and domain
        molecular_function = properties.get('molecular_function', 'Unknown')
        domain = properties.get('domain', 'Unknown')
        
        # Normalize molecular function and domain strings
        if not molecular_function or pd.isna(molecular_function):
            molecular_function = 'Unknown'
        if not domain or pd.isna(domain):
            domain = 'Unknown'
            
        # Add to data list
        plot_data.append({
            'structure_id': struct_id,
            'short_name': short_name,
            'molecular_function': molecular_function,
            'molecular_function_normalized': molecular_function.replace(' ', '_'),
            'domain': domain,
            'experimentally_determined': has_exp_structure
        })
    
    # Convert to DataFrame
    df = pd.DataFrame(plot_data)
    
    # Sort by molecular function
    df_sorted = df.sort_values('molecular_function_normalized', ignore_index=True)
    N = len(df_sorted)
    
    # Normalize function names for consistency with color mapping
    # Map normalized function names to display format
    function_map = {
        'Proton_pump': 'Proton pump',
        'proton_pump': 'Proton pump',
        'Chloride_pump': 'Chloride pump',
        'chloride_pump': 'Chloride pump',
        'Cation_channel': 'Cation channel',
        'cation_channel': 'Cation channel',
        'Anion_channel': 'Anion channel',
        'anion_channel': 'Anion channel',
        'Photosensor': 'Photosensor',
        'photosensor': 'Photosensor',
        'Phototaxis': 'Phototaxis',
        'phototaxis': 'Phototaxis',
        'unknown': 'Unknown',
        'Unknown': 'Unknown'
    }
    
