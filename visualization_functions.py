"""Updated visualization functions for opsin analysis that use the global color scheme."""
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd
from scipy.spatial.distance import squareform, cdist
from scipy.cluster.hierarchy import linkage, dendrogram
import plotly.graph_objects as go
import logomaker
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

import os

# Import our color scheme
from opsin_color_scheme import (
    OPSIN_COLORS, RMSD_CMAP, RMSD_COMPACT_CMAP, DIVERGING_CMAP, RMSD_DISCRETE_CMAP, RMSD_BOUNDS,
    HELIX_COLORS, HELIX_COLORS_STR, HELIX_COLORS_LIST, get_group_colors)

def create_rmsd_color_scale_figure(output_path=None):
    """
    Create a standalone figure showing the RMSD color scale.
    
    Args:
        output_path: Optional path to save the figure
        
    Returns:
        Matplotlib figure
    """
    # Create figure
    fig, ax = plt.subplots(figsize=(6, 2))
    
    # Create a gradient from 0 to 3.0 and above
    gradient = np.linspace(0, 3.5, 100).reshape(1, -1)
    
    # Show the gradient with our colormap
    im = ax.imshow(gradient, aspect='auto', cmap=RMSD_COMPACT_CMAP, 
                  extent=[0, 3.5, 0, 1], vmin=0, vmax=3.0)
    
    # Add ticks at key points
    ax.set_xticks([0, 1.5, 3.0, 3.5])
    ax.set_xticklabels(['0.0', '1.5', '3.0', '>3.0'])
    ax.set_yticks([])
    
    # Add annotations for the color regions
    ax.annotate('Dark Blue', xy=(0.75, 0.5), ha='center', va='center', color='white')
    ax.annotate('Light Blue', xy=(2.25, 0.5), ha='center', va='center', color='black')
    ax.annotate('Yellow', xy=(3.25, 0.5), ha='center', va='center', color='black')
    
    # Add title and labels
    ax.set_title('RMSD Color Scale')
    ax.set_xlabel('RMSD (Å)')
    
    # Save if requested
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    
    return fig

def plot_rmsd_heatmap(rmsd_df, title="All-vs-All RMSD Heatmap"):
    """Visualize the RMSD matrix as a heatmap using the global color scheme."""
    plt.figure(figsize=(10, 8))
    # Use our custom compact colormap for RMSD values
    vmax = 3.0  # Max value for color scale (anything above will be yellow)
    sns.heatmap(rmsd_df, annot=True, fmt=".2f", cmap=RMSD_COMPACT_CMAP, vmin=0, vmax=vmax)
    plt.title(title)
    plt.xlabel("Structure")
    plt.ylabel("Structure")
    plt.tight_layout()
    return plt.gcf()

def plot_similarity_tree(rmsd_df, title="Structural Similarity Tree"):
    """Generate a dendrogram based on RMSD values with consistent colors."""
    matrix = rmsd_df.fillna(rmsd_df.max().max()).values
    condensed = squareform(matrix)
    Z = linkage(condensed, method='average')
    plt.figure(figsize=(12, 8))
    
    # Use a consistent function for branch colors
    def link_color_func(k):
        return OPSIN_COLORS['gray_dark']
    
    dendrogram(
        Z,
        labels=rmsd_df.index,
        leaf_rotation=45,
        link_color_func=link_color_func
    )
    plt.title(title)
    plt.xlabel("Structure")
    plt.ylabel("RMSD (Å)")
    plt.tight_layout()
    return plt.gcf()

def create_residue_conservation_plot(msa_df, helix_highlighting=True, figsize=(12, 8)):
    """
    Create a heatmap showing conservation levels at each position in the MSA.
    Handles both dot notation (1.50) and x notation (1x50) GRN formats.
    
    Args:
        msa_df: MSA DataFrame with GRN column names
        helix_highlighting: Whether to color-code positions by helix
        figsize: Figure size tuple
    Returns:
        Matplotlib figure
    """
    # Extract residue types from MSA (first character of each cell)
    residue_df = msa_df.applymap(lambda x: x[0] if isinstance(x, str) and len(x) > 0 and x != '-' else '-')
    
    # Calculate conservation percentage
    conservation = {}
    for col in residue_df.columns:
        # Count residue types at this position
        counts = residue_df[col].value_counts()
        total = (residue_df[col] != '-').sum()  # Don't count gaps
        
        if total > 0:
            # Calculate percentage of most common residue
            top_res = counts.index[0] if not counts.empty else '-'
            if top_res != '-':
                top_pct = (counts[top_res] / total) * 100
            else:
                # If top residue is a gap, use the next one
                top_res = counts.index[1] if len(counts) > 1 else '-'
                top_pct = (counts[top_res] / total) * 100 if top_res != '-' else 0
            
            conservation[col] = {
                'residue': top_res,
                'percentage': top_pct
            }
        else:
            conservation[col] = {
                'residue': '-',
                'percentage': 0
            }
    
    # Create a DataFrame with position, helix, and conservation data
    plot_data = []
    for pos, data in conservation.items():
        # Handle different GRN formats
        if 'x' in str(pos):
            # Handle 1x50 format (transmembrane helices)
            try:
                helix, num = pos.split('x')
                if helix.isdigit():
                    plot_data.append({
                        'Position': pos,
                        'Helix': helix,
                        'Position_Num': float(num),
                        'Residue': data['residue'],
                        'Conservation': data['percentage']
                    })
                else:
                    # Non-numeric helix
                    plot_data.append({
                        'Position': pos,
                        'Helix': 'Other',
                        'Position_Num': float(num) if num.isdigit() else 0,
                        'Residue': data['residue'],
                        'Conservation': data['percentage']
                    })
            except:
                # Fallback for malformed positions
                plot_data.append({
                    'Position': pos,
                    'Helix': 'Other',
                    'Position_Num': 0,
                    'Residue': data['residue'],
                    'Conservation': data['percentage']
                })
        elif '.' in str(pos):
            # Handle multiple dot notation formats
            parts = str(pos).split('.')
            
            if pos.startswith('L.'):
                # Loop region with L. prefix
                plot_data.append({
                    'Position': pos,
                    'Helix': 'L',
                    'Position_Num': float(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0,
                    'Residue': data['residue'],
                    'Conservation': data['percentage']
                })
            elif pos.startswith('n.'):
                # N-terminal region
                plot_data.append({
                    'Position': pos,
                    'Helix': 'N',
                    'Position_Num': float(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0,
                    'Residue': data['residue'],
                    'Conservation': data['percentage']
                })
            elif pos.startswith('c.'):
                # C-terminal region
                plot_data.append({
                    'Position': pos,
                    'Helix': 'C',
                    'Position_Num': float(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0,
                    'Residue': data['residue'],
                    'Conservation': data['percentage']
                })
            elif len(parts) == 2 and parts[0].isdigit():
                # Standard helix.position format (1.50)
                plot_data.append({
                    'Position': pos,
                    'Helix': parts[0],
                    'Position_Num': float(parts[1]) if parts[1].isdigit() else 0,
                    'Residue': data['residue'],
                    'Conservation': data['percentage']
                })
            elif len(parts) >= 2 and len(parts[0]) == 2 and parts[0].isdigit():
                # Loop region in AB.CCC format (12.003)
                plot_data.append({
                    'Position': pos,
                    'Helix': 'L',
                    'Position_Num': float('0.' + parts[1]) if parts[1].isdigit() else 0,
                    'Residue': data['residue'],
                    'Conservation': data['percentage']
                })
            else:
                # Other dot notation formats
                plot_data.append({
                    'Position': pos,
                    'Helix': 'Other',
                    'Position_Num': 0,
                    'Residue': data['residue'],
                    'Conservation': data['percentage']
                })
        else:
            # Handle numeric or other position labels
            plot_data.append({
                'Position': pos,
                'Helix': 'Other',
                'Position_Num': float(pos) if str(pos).isdigit() else 0,
                'Residue': data['residue'],
                'Conservation': data['percentage']
            })
    
    # Convert to DataFrame and sort
    conservation_df = pd.DataFrame(plot_data)
    
    # Sort by helix and position
    if 'Helix' in conservation_df.columns and not conservation_df['Helix'].isna().all():
        def helix_sort_key(h):
            if h.isdigit():
                return int(h)  # TM helices first, sorted by number
            elif h == 'N':
                return -1  # N-terminal before helices
            elif h == 'L':
                return 10  # Loop regions after TM helices
            elif h == 'C':
                return 11  # C-terminal after loops
            else:
                return 12  # Other regions last
        
        conservation_df['Helix_Sort'] = conservation_df['Helix'].apply(helix_sort_key)
        conservation_df = conservation_df.sort_values(['Helix_Sort', 'Position_Num'])
    else:
        # Fall back to sorting by position string
        conservation_df = conservation_df.sort_values('Position')
    
    # Create figure
    fig, ax = plt.subplots(figsize=figsize)
    
    # Set up positions for bars
    x = np.arange(len(conservation_df))
    
    # Create bars
    bars = ax.bar(x, conservation_df['Conservation'], width=0.8)
    
    # Color bars by helix if requested
    if helix_highlighting and 'Helix' in conservation_df.columns:
        for i, (_, row) in enumerate(conservation_df.iterrows()):
            helix = row['Helix']
            if helix.isdigit() and int(helix) in HELIX_COLORS:
                bars[i].set_color(HELIX_COLORS[int(helix)])
            elif helix == 'L':
                bars[i].set_color(OPSIN_COLORS['gray_light'])
            elif helix == 'N':
                bars[i].set_color(OPSIN_COLORS['blue_light'])
            elif helix == 'C':
                bars[i].set_color(OPSIN_COLORS['green'])
            else:
                bars[i].set_color(OPSIN_COLORS['gray_dark'])
    
    # Annotate bars with residue letters
    for i, (_, row) in enumerate(conservation_df.iterrows()):
        if row['Conservation'] > 5:  # Only annotate if there's enough space
            ax.text(i, row['Conservation'] + 2, row['Residue'],
                  ha='center', va='bottom', fontweight='bold')
    
    # Add labels and title
    ax.set_ylabel('Conservation (%)', fontsize=12)
    ax.set_xlabel('Position', fontsize=12)
    ax.set_title('Residue Conservation by Position', fontsize=14)
    
    # Set x-tick labels to position names
    ax.set_xticks(x)
    ax.set_xticklabels(conservation_df['Position'], rotation=90)
    
    # Add a reference line at 50% conservation
    ax.axhline(y=50, color=OPSIN_COLORS['gray_light'], linestyle='--', alpha=0.7)
    
    # Set y-axis limits
    ax.set_ylim(0, 105)
    
    # Add grid
    ax.grid(axis='y', alpha=0.3)
    
    # Add legend for helix colors
    if helix_highlighting:
        legend_elements = []
        for i in range(1, 8):
            if i in HELIX_COLORS:
                legend_elements.append(
                    Patch(facecolor=HELIX_COLORS[i], label=f'Helix {i}')
                )
        legend_elements.extend([
            Patch(facecolor=OPSIN_COLORS['blue_light'], label='N-terminal'),
            Patch(facecolor=OPSIN_COLORS['gray_light'], label='Loop'),
            Patch(facecolor=OPSIN_COLORS['green'], label='C-terminal')
        ])
        ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    return fig

def visualize_rmsd_heatmap(rmsd_df, structure_ids, group_dict=None, domain_dict=None, name_dict=None,
                         annot=False, font_scale=1.0, group_by='molecular_function', error_threshold=3.0):
    """
    Visualize RMSD matrix as a heatmap with our consistent color scheme.
    Args:
        rmsd_df: DataFrame with RMSD values
        structure_ids: List of structure IDs
        group_dict: Dictionary mapping structure IDs to molecular function groups
        domain_dict: Dictionary mapping structure IDs to domain groups
        name_dict: Dictionary mapping structure IDs to display names
        annot: Whether to annotate the heatmap with values (default: False)
        font_scale: Scale factor for font sizes
        group_by: Property to use for primary grouping (default: 'molecular_function')
        error_threshold: Remove entries with average error above this threshold
    Returns:
        Matplotlib figure
    """
    sns.set_context("notebook", font_scale=font_scale)
    plt.figure(figsize=(16, 14), dpi=100)
    
    # Filter out structures with high error if domain_dict contains average error values
    filtered_structure_indices = list(range(len(structure_ids)))
    if domain_dict is not None and any('average_error' in domain_dict.get(sid, {}) for sid in structure_ids):
        # Keep only structures with average error below threshold
        filtered_structure_indices = [
            i for i, sid in enumerate(structure_ids) 
            if 'average_error' not in domain_dict.get(sid, {}) or 
               domain_dict.get(sid, {}).get('average_error', 0) <= error_threshold
        ]
        
        if len(filtered_structure_indices) < len(structure_ids):
            print(f"Filtered out {len(structure_ids) - len(filtered_structure_indices)} structures with average error > {error_threshold}")
            # Update structure_ids to only include filtered ones
            structure_ids = [structure_ids[i] for i in filtered_structure_indices]
    
    # Create domain dictionary if only string values are provided
    domain_dict_local = {}
    if domain_dict is not None:
        for sid in structure_ids:
            if sid in domain_dict:
                # Get domain directly if it's a string, or from the domain key if it's a dict
                if isinstance(domain_dict[sid], dict) and 'domain' in domain_dict[sid]:
                    domain_dict_local[sid] = domain_dict[sid]['domain']
                elif isinstance(domain_dict[sid], str):
                    domain_dict_local[sid] = domain_dict[sid]
                else:
                    domain_dict_local[sid] = "Unknown"
            else:
                domain_dict_local[sid] = "Unknown"
    
    # If a group_dict is provided, sort the indices by group so that rows/columns become contiguous.
    if group_dict is not None:
        # Sort first by molecular function, then by domain for structures with the same function
        if domain_dict_local:
            # Create a composite key for sorting
            sort_key = lambda i: (
                group_dict.get(structure_ids[i], "Unknown"), 
                domain_dict_local.get(structure_ids[i], "Unknown")
            )
        else:
            sort_key = lambda i: group_dict.get(structure_ids[i], "Unknown")
        
        sorted_indices = sorted(range(len(structure_ids)), key=sort_key)
        structure_ids_sorted = [structure_ids[i] for i in sorted_indices]
        rmsd_matrix = rmsd_df.loc[structure_ids_sorted, structure_ids_sorted].values
    else:
        structure_ids_sorted = structure_ids
        rmsd_matrix = rmsd_df.loc[structure_ids_sorted, structure_ids_sorted].values
    
    # Replace NaN and infinite values.
    max_rmsd = np.nanmax(rmsd_matrix[np.isfinite(rmsd_matrix)])
    rmsd_matrix = np.nan_to_num(rmsd_matrix, nan=0.0, posinf=max_rmsd, neginf=0.0)
    
    # Define color boundaries and use our discrete colormap
    rmsd_bounds = RMSD_BOUNDS.copy()
    if max_rmsd > rmsd_bounds[-1]:
        rmsd_bounds.append(max_rmsd)
    norm = BoundaryNorm(rmsd_bounds, RMSD_DISCRETE_CMAP.N)
    
    # Plot the heatmap
    ax = plt.gca()
    im = ax.imshow(rmsd_matrix, cmap=RMSD_DISCRETE_CMAP, norm=norm)
    
    # Add colorbar on the left
    try:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("left", size="5%", pad=0.5)
        cbar = plt.colorbar(im, cax=cax, ticks=rmsd_bounds)
        cbar.set_label('RMSD (Å)', rotation=90, fontsize=14, labelpad=15)
        cax.yaxis.set_ticks_position('left')
        cax.yaxis.set_label_position('left')
    except (ImportError, FileNotFoundError) as e:
        print(f"Warning: Could not create colorbar with axes_grid1: {e}")
        # Create a simpler colorbar without axes_grid1
        try:
            cbar = plt.colorbar(im, ax=ax, orientation='vertical', pad=0.05)
            cbar.set_label('RMSD (Å)', rotation=90, fontsize=14, labelpad=15)
        except Exception as e2:
            print(f"Warning: Could not create alternate colorbar: {e2}")
    
    # Hide axis ticks and labels
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.set_title('Structural Similarity Matrix (RMSD)', pad=20, fontsize=18, fontweight='bold')
    
    # If group_dict is provided, compute group boundaries and add visual indicators
    if group_dict is not None:
        # Since structure_ids_sorted are now sorted by group, create a list of group names
        molecular_functions = [group_dict.get(sid, "Unknown") for sid in structure_ids_sorted]
        domains = [domain_dict_local.get(sid, "Unknown") for sid in structure_ids_sorted]
        
        # Create unique groups and boundaries for molecular functions
        unique_functions = []
        function_sizes = []
        current_function = None
        count = 0
        
        for func in molecular_functions:
            if func != current_function:
                if current_function is not None:
                    unique_functions.append(current_function)
                    function_sizes.append(count)
                current_function = func
                count = 1
            else:
                count += 1
        
        if current_function is not None:
            unique_functions.append(current_function)
            function_sizes.append(count)
        
        function_boundaries = np.cumsum(function_sizes)
        
        # Draw white separator lines between groups for molecular functions
        for boundary in function_boundaries[:-1]:
            ax.axhline(y=boundary - 0.5, color=OPSIN_COLORS['white'], linewidth=2.5)
            ax.axvline(x=boundary - 0.5, color=OPSIN_COLORS['white'], linewidth=2.5)
        
        # Now handle domains - create boundaries within molecular function groups
        if domain_dict_local:
            domain_subgroups = []  # Stores (domain, start_idx, size) tuples
            
            start_idx = 0
            for func_size, func_name in zip(function_sizes, unique_functions):
                # Get all domains within this function group
                func_range = list(range(start_idx, start_idx + func_size))
                func_domains = [domains[i] for i in func_range]
                
                # Count domains within this function
                domain_counts = {}
                last_domain = None
                for i, domain in enumerate(func_domains):
                    if domain != last_domain:
                        domain_counts[domain] = domain_counts.get(domain, 0) + 1
                        last_domain = domain
                        domain_subgroups.append((domain, start_idx + i, 1))
                    else:
                        # Increment the size of the last domain subgroup
                        _, subgroup_start, subgroup_size = domain_subgroups[-1]
                        domain_subgroups[-1] = (domain, subgroup_start, subgroup_size + 1)
                
                start_idx += func_size
            
            # Draw lighter separator lines between domains within function groups
            for _, subgroup_start, subgroup_size in domain_subgroups:
                end_idx = subgroup_start + subgroup_size
                if end_idx < len(structure_ids_sorted) and end_idx not in function_boundaries:
                    ax.axhline(y=end_idx - 0.5, color=OPSIN_COLORS['white'], linewidth=1.0, alpha=0.5)
                    ax.axvline(x=end_idx - 0.5, color=OPSIN_COLORS['white'], linewidth=1.0, alpha=0.5)
        
        # Get consistent colors for functions and domains
        function_color_map = get_group_colors(sorted(set(molecular_functions)))
        domain_color_map = get_group_colors(sorted(set(domains)))
        
        # Add colored bars for functions along the top edge
        pos = 0
        for func, size in zip(unique_functions, function_sizes):
            # Top edge color bar for function
            ax.add_patch(plt.Rectangle(
                (pos, -0.05*rmsd_matrix.shape[0]), size, 0.02*rmsd_matrix.shape[0],
                facecolor=function_color_map[func], edgecolor='none', alpha=0.8,
                transform=ax.transData, clip_on=False
            ))
            pos += size
        
        # Add colored bars for domains along the right edge
        if domain_dict_local:
            for domain, start_idx, size in domain_subgroups:
                # Right edge color bar for domain
                ax.add_patch(plt.Rectangle(
                    (rmsd_matrix.shape[1], start_idx), 0.02*rmsd_matrix.shape[1], size,
                    facecolor=domain_color_map[domain], edgecolor='none', alpha=0.8,
                    transform=ax.transData, clip_on=False
                ))
        
        # Create legend elements for both functions and domains
        function_legend_elements = [
            plt.Rectangle((0, 0), 1, 1, facecolor=function_color_map[func], 
                         edgecolor='none', label=func)
            for func in sorted(set(molecular_functions))
        ]
        
        domain_legend_elements = [
            plt.Rectangle((0, 0), 1, 1, facecolor=domain_color_map[domain], 
                         edgecolor='none', label=domain)
            for domain in sorted(set(domains))
        ]
        
        # Create a new figure for the legends to avoid reusing artists
        legend_fig = plt.figure(figsize=(3, 5))
        legend_ax = legend_fig.add_subplot(111)
        legend_ax.axis('off')
        
        # Add a legend for molecular functions
        leg1 = legend_ax.legend(
            handles=function_legend_elements,
            title="Molecular Function",
            loc='upper center',
            fontsize=12,
            title_fontsize=14
        )
        leg1.get_frame().set_alpha(0.8)
        
        # Add a legend for domains if available
        if domain_dict_local:
            leg2 = legend_ax.legend(
                handles=domain_legend_elements,
                title="Domain",
                loc='lower center',
                fontsize=12,
                title_fontsize=14
            )
            leg2.get_frame().set_alpha(0.8)
            
        # Save the legend figure
        legend_path = 'rmsd_heatmap_legends.png'
        legend_fig.savefig(legend_path, dpi=300, bbox_inches='tight')
        plt.close(legend_fig)
        
        # Add text to indicate legends are in a separate file
        ax.text(0.95, 0.05, "Legends saved separately", 
                transform=ax.transAxes, ha='right', va='bottom',
                fontsize=10, color='black', bbox=dict(facecolor='white', alpha=0.7))
    
    plt.tight_layout()
    return plt.gcf()

def create_and_visualize_similarity_tree(rmsd_data, group_dict=None, domain_dict=None,
                                         name_dict=None, font_scale=1.0, leaf_font_size=12,
                                         group_by='molecular_function',
                                         linkage_matrix=None):  # Added linkage_matrix param

    # Removed error_threshold
    """
    Create a similarity tree (dendrogram) from the RMSD matrix with consistent colors.
    Args:
        rmsd_data: RMSD DataFrame (should be pre-filtered)
        group_dict: Dictionary mapping structure IDs to molecular function groups
        domain_dict: Dictionary mapping structure IDs to domain groups (dict or string)
        name_dict: Dictionary mapping structure IDs to display names
        font_scale: Scale factor for font sizes
        leaf_font_size: Font size for leaf labels
        group_by: Property to use for primary coloring (default: 'molecular_function')
        linkage_matrix: Pre-calculated linkage matrix (Z) (optional)
    Returns:
        Tuple: (Matplotlib figure, list of ordered structure IDs)
    """
    sns.set_context("notebook", font_scale=font_scale)
    fig = plt.figure(figsize=(16, 12), dpi=100)  # Create figure instance
    plt.grid(False)

    # Extract structure IDs and RMSD matrix from input (assuming rmsd_data is DataFrame)
    structure_ids = rmsd_data.index.tolist()
    rmsd_matrix = rmsd_data.values
    N = len(structure_ids)

    # Check dimensions
    if rmsd_matrix.shape != (N, N):
        raise ValueError(f"RMSD matrix dimensions ({rmsd_matrix.shape}) must match the number of structures ({N}).")

    # Use provided linkage matrix Z or calculate it if not provided
    if linkage_matrix is None:
        print("DEBUG (Tree): Calculating linkage matrix internally.")
        # Clean the matrix: replace NaN/Inf and force zeros on the diagonal.
        if np.any(~np.isfinite(rmsd_matrix)):
            mask = np.isfinite(rmsd_matrix)
            if not np.any(mask):
                max_rmsd = 1.0
            else:
                max_rmsd = np.max(rmsd_matrix[mask])
            clean_matrix = np.nan_to_num(rmsd_matrix, nan=0.0, posinf=max_rmsd, neginf=0.0)
        else:
            clean_matrix = rmsd_matrix.copy()

        np.fill_diagonal(clean_matrix, 0.0)

        # Compute the condensed distance matrix and linkage.
        condensed_matrix = squareform(clean_matrix, checks=False)
        Z = linkage(condensed_matrix, method='average')
    else:
        print("DEBUG (Tree): Using provided linkage matrix.")
        Z = linkage_matrix  # Use the pre-calculated Z


    # Force branch colors to black
    def link_color_func(k):
        return OPSIN_COLORS['gray_dark']


    # Create dendrogram without labels first (we'll add them later with colors)
    ax = plt.gca()  # Get current axes associated with the figure
    dendro = dendrogram(
        Z,
        ax=ax,  # Explicitly pass axes
        labels=None,  # No labels initially
        orientation='right',
        link_color_func=link_color_func,
        above_threshold_color=OPSIN_COLORS['gray_dark'],
        show_leaf_counts=False,
        no_labels=True  # Hide labels
    )

    ax.tick_params(axis='x', labelsize=12)
    ax.tick_params(axis='y', labelcolor='none')  # Hide y-tick labels

    # Get the molecular function and domain for each structure
    # ... (logic for function_dict, domain_dict_local, function_color_map, domain_color_map remains the same) ...
    # Make sure this uses the structure_ids from the passed (filtered) rmsd_data
    function_dict = {}
    domain_dict_local = {}

    if group_dict is not None:
        for sid in structure_ids:
            function_dict[sid] = str(group_dict.get(sid, "Unknown"))  # Ensure string
    else:
        function_dict = {sid: "Unknown" for sid in structure_ids}

    if domain_dict is not None:
        for sid in structure_ids:
            domain_info = domain_dict.get(sid, "Unknown")
            if isinstance(domain_info, dict):
                domain_dict_local[sid] = str(domain_info.get('domain', "Unknown"))  # Ensure string
            else:  # Assume it's already a string or needs conversion
                domain_dict_local[sid] = str(domain_info)
    else:
        domain_dict_local = {sid: "Unknown" for sid in structure_ids}

    # Get unique function and domain categories
    unique_functions = sorted(set(function_dict.values()))
    unique_domains = sorted(set(domain_dict_local.values()))

    # Get consistent colors for functions and domains
    function_color_map = get_group_colors(unique_functions)
    domain_color_map = get_group_colors(unique_domains)

    # The dendrogram returns an ordering of the original indices
    leaf_indices = dendro['leaves']
    ordered_structure_ids = [structure_ids[i] for i in leaf_indices]  # Calculate the order

    # Add colored function and domain labels on the right side
    # ... (patch drawing logic remains the same, using ax) ...
    for i, leaf_idx in enumerate(leaf_indices):
        if leaf_idx < len(structure_ids):
            sid = structure_ids[leaf_idx]
            function = function_dict.get(sid, "Unknown")
            domain = domain_dict_local.get(sid, "Unknown")

            # Add colored rectangular patch for molecular function
            y_pos = i * 10  # Dendrogram y-coordinates are distances; use index for positioning patches
            rect_height = 7  # Adjust as needed based on dendrogram scaling

            # Use ax.transData for positioning relative to data coordinates
            # The y-coordinates from dendrogram need careful handling for patch placement.
            # A simpler approach is to get the y-coords directly from the dendrogram output if possible,
            # or iterate through the plotted lines. Let's use the index `i` for simplicity, assuming leaves are evenly spaced visually.

            # We need the y-coordinates of the leaves. dendrogram plots them at 10, 20, 30...
            y_coord = (i * 10) + 5  # Center of the leaf position visually

            # Molecular function rectangle (first from the right)
            # Adjust x position based on data limits
            x_limit_max = ax.get_xlim()[1]
            function_rect_width = 0.07 * x_limit_max
            function_rect = plt.Rectangle(
                (x_limit_max, y_coord - rect_height / 2),  # Use calculated y_coord
                function_rect_width, rect_height,
                facecolor=function_color_map.get(function, OPSIN_COLORS['gray_light']),  # Use .get for safety
                edgecolor='none', alpha=0.8,
                clip_on=False  # Allow drawing outside axes
            )
            ax.add_patch(function_rect)

            # Domain rectangle (second from the right)
            domain_rect_width = 0.04 * x_limit_max
            domain_rect = plt.Rectangle(
                (x_limit_max + function_rect_width + 0.01 * x_limit_max, y_coord - rect_height / 2),
                # Use calculated y_coord
                domain_rect_width, rect_height,
                facecolor=domain_color_map.get(domain, OPSIN_COLORS['gray_light']),  # Use .get for safety
                edgecolor='none', alpha=0.8,
                clip_on=False  # Allow drawing outside axes
            )
            ax.add_patch(domain_rect)

    # Create legends
    # ... (legend creation logic remains the same, using ax) ...
    # Create a legend for molecular functions
    function_legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor=function_color_map.get(func, OPSIN_COLORS['gray_light']),
                      edgecolor='none', label=func)
        for func in unique_functions
    ]

    # Create a legend for domains
    domain_legend_elements = [
        plt.Rectangle((0, 0), 1, 1, facecolor=domain_color_map.get(domain, OPSIN_COLORS['gray_light']),
                      edgecolor='none', label=domain)
        for domain in unique_domains
    ]

    # Add legends on the right side
    # First legend for molecular function
    leg1 = ax.legend(  # Use ax.legend
        handles=function_legend_elements,
        title="Molecular Function",
        bbox_to_anchor=(1.35, 0.7),
        loc='center left',
        fontsize=leaf_font_size,
        title_fontsize=leaf_font_size + 2
    )
    leg1.get_frame().set_alpha(0.8)
    ax.add_artist(leg1)  # Add legend back to the axes

    # Second legend for domain
    if domain_dict_local:  # Only add if domains are present
        leg2 = ax.legend(  # Use ax.legend
            handles=domain_legend_elements,
            title="Domain",
            bbox_to_anchor=(1.35, 0.3),
            loc='center left',
            fontsize=leaf_font_size,
            title_fontsize=leaf_font_size + 2
        )
        leg2.get_frame().set_alpha(0.8)
        # Add the second legend manually AFTER the first one has been added
        ax.add_artist(leg2)

    # Add colorbar for RMSD distances on the left
    # ... (colorbar logic remains the same, using ax) ...
    # Use fixed maximum value of 3.0 for consistent color scale
    # Calculate max RMSD from the original matrix used for linkage
    max_dist_for_norm = np.max(Z[:, 2]) if Z is not None and len(Z) > 0 else 3.0
    norm = plt.Normalize(0, min(3.0, max_dist_for_norm))
    cmap = plt.cm.get_cmap(RMSD_COMPACT_CMAP)

    # Create colorbar with proper placement
    try:
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("left", size="5%", pad=0.5)
        cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax)
        cb.set_label('RMSD (Å)', rotation=90, fontsize=12, labelpad=15)
        cax.yaxis.set_ticks_position('left')
        cax.yaxis.set_label_position('left')
    except Exception as e:  # Catch more general exceptions
        print(f"Warning: Could not create colorbar for similarity tree: {e}")
        # Create a simpler colorbar without axes_grid1
        try:
            cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),  # Use fig.colorbar
                              ax=ax, orientation='vertical', pad=0.05, aspect=40)  # Adjust aspect
            cb.set_label('RMSD (Å)', rotation=90, fontsize=12, labelpad=15)
        except Exception as e2:
            print(f"Warning: Could not create alternate colorbar: {e2}")

    # Set title and labels
    ax.set_title('Structural Similarity Tree', fontsize=18, pad=20, fontweight='bold')
    ax.set_xlabel('Distance (RMSD, Å)', fontsize=14)

    # Remove the top and right spines for a cleaner look
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    fig.tight_layout(rect=[0.05, 0, 0.85, 1])  # Adjust layout to make space for legends

    return fig, ordered_structure_ids  # Return the figure and the calculated order

def visualize_rmsd_matrix_improved(rmsd_df, group_dict=None, name_dict=None,
                                       output_file=None,
                                       # cluster_method='average', # Method defined externally now
                                       figsize=(14, 12), domain_dict=None,
                                       linkage_matrix=None):  # Added linkage_matrix param
    # Removed error_threshold
    """
    Improved visualization of RMSD matrix using a pre-calculated clustering.
    Uses the global color scheme for consistent visualization.
    Args:
        rmsd_df: DataFrame with RMSD values (should be pre-filtered)
        group_dict: Dictionary mapping structure ID to group/type (for y-axis)
        name_dict: Dictionary mapping structure ID to display name
        output_file: Optional file path to save the figure
        figsize: Figure size
        domain_dict: Dictionary mapping structure ID to domain (for x-axis, dict or string)
        linkage_matrix: Pre-calculated linkage matrix (Z) REQUIRED
    Returns:
        Seaborn clustermap object or None if error
    """
    # Set up plotting - DO NOT create a new figure here, clustermap does it.
    # plt.figure(figsize=figsize) # REMOVE THIS

    # Get structure IDs from the filtered DataFrame
    pdb_list = rmsd_df.index.tolist()

    # --- Filtering is done BEFORE calling this function ---
    # --- Linkage matrix Z is calculated BEFORE calling this function ---

    if linkage_matrix is None:
        print("ERROR (Heatmap): Linkage matrix (Z) is required for visualize_rmsd_matrix_improved.")
        # Optionally, calculate Z internally as a fallback, but it might not match the tree
        # print("WARNING (Heatmap): Calculating linkage matrix internally as fallback.")
        # clean_matrix = rmsd_df.fillna(0).values
        # np.fill_diagonal(clean_matrix, 0)
        # if np.any(~np.isfinite(clean_matrix)):
        #     max_val = np.nanmax(clean_matrix[np.isfinite(clean_matrix)]) if np.any(np.isfinite(clean_matrix)) else 1.0
        #     clean_matrix = np.nan_to_num(clean_matrix, nan=max_val, posinf=max_val, neginf=0.0)
        # condensed_dist = squareform(clean_matrix, checks=False)
        # linkage_matrix = linkage(condensed_dist, method='average') # Or use default passed method
        return None  # Recommended to fail if Z is not provided

    # Check if we have enough structures to continue
    if len(pdb_list) < 2:
        print("Not enough structures to create RMSD matrix visualization.")
        return None

    # clustermap will derive the order from the provided linkage matrix Z
    # We need to prepare the row and column colors based on the *final* order
    # that clustermap will determine from Z.

    # Perform a dummy dendrogram calculation just to get the order derived from Z
    try:
        dendrogram_out = dendrogram(linkage_matrix, no_plot=True)
        ordered_indices = dendrogram_out['leaves']
        ordered_ids = [pdb_list[i] for i in ordered_indices]
        print(f"DEBUG (Heatmap): Order derived from provided Z: {len(ordered_ids)} structures.")
    except Exception as e:
        print(f"ERROR (Heatmap): Could not derive order from provided linkage matrix: {e}")
        ordered_ids = pdb_list  # Fallback to original order if dendrogram fails

    # Group colors for y-axis (molecular function) based on the derived order
    row_colors_data = None
    row_color_map = {}
    if group_dict:
        unique_groups = sorted(set(str(v) for k, v in group_dict.items() if k in ordered_ids))  # Use filtered dict keys
        if 'Unknown' not in unique_groups: unique_groups.append('Unknown')
        row_color_map = get_group_colors(unique_groups)
        row_colors_list = [row_color_map.get(str(group_dict.get(id, 'Unknown')), row_color_map['Unknown']) for id in
                           ordered_ids]
        row_colors_data = pd.Series(row_colors_list, index=ordered_ids, name='Function')  # Use Series for clustermap

    # Domain colors for x-axis (domain) based on the derived order
    col_colors_data = None
    col_color_map = {}
    if domain_dict:
        domain_values = []
        for struct_id in ordered_ids:
            domain_info = domain_dict.get(struct_id, "Unknown")
            if isinstance(domain_info, dict):
                domain_values.append(str(domain_info.get('domain', 'Unknown')))
            else:
                domain_values.append(str(domain_info))

        unique_domains = sorted(set(domain_values))
        if 'Unknown' not in unique_domains: unique_domains.append('Unknown')
        col_color_map = get_group_colors(unique_domains)
        col_colors_list = [col_color_map.get(domain_values[i], col_color_map['Unknown']) for i in range(len(ordered_ids))]
        col_colors_data = pd.Series(col_colors_list, index=ordered_ids, name='Domain')  # Use Series for clustermap

    # Use the same colors for columns if domain dict not provided but group dict is
    elif group_dict:
        col_colors_data = row_colors_data.rename('Group')  # Use same colors but different label
        col_color_map = row_color_map

    # --- Plot using clustermap with pre-calculated Z ---
    print(f"DEBUG (Heatmap): Plotting clustermap with {rmsd_df.shape} matrix and provided Z.")
    try:
        g = sns.clustermap(
            rmsd_df,  # Pass the original filtered (but unordered) DataFrame
            row_linkage=linkage_matrix,  # Provide Z for rows
            col_linkage=linkage_matrix,  # Provide Z for columns
            row_colors=row_colors_data,  # Provide Series if available
            col_colors=col_colors_data,  # Provide Series if available
            xticklabels=False, yticklabels=False,  # Hide labels to declutter the plot
            figsize=figsize, cmap=RMSD_COMPACT_CMAP,
            vmin=0, vmax=3.0,  # Fixed range for color scale (0-3.0)
            cbar_kws={'label': 'RMSD (Å)'},
            dendrogram_ratio=(.2, .2)  # Adjust dendrogram size if needed
        )

        # --- Add Legends ---
        # This part is tricky with clustermap as it uses its own figure structure.
        # We need to add legends to the main figure (g.fig) or specific axes.

        handles = []
        labels = []
        # Molecular Function Legend (using row colors)
        if row_colors_data is not None and row_color_map:
            handles.extend([plt.Rectangle((0, 0), 1, 1, color=color) for group, color in sorted(row_color_map.items())])
            labels.extend([group for group, color in sorted(row_color_map.items())])

        # Domain Legend (using column colors)
        if col_colors_data is not None and col_color_map and col_colors_data.name != row_colors_data.name:  # Add only if different from row colors
            # Add a separator visually if needed
            handles.append(plt.Rectangle((0, 0), 0, 0, color='white'))  # Invisible separator
            labels.append("")
            handles.extend([plt.Rectangle((0, 0), 1, 1, color=color) for domain, color in sorted(col_color_map.items())])
            labels.extend([domain for domain, color in sorted(col_color_map.items())])

        if handles:
            try:
                # Add legend to the figure - position might need tweaking
                legend_title = "Function / Domain"  # Generic title
                if row_colors_data is not None and col_colors_data is not None and col_colors_data.name != row_colors_data.name:
                    legend_title = f"{row_colors_data.name} / {col_colors_data.name}"
                elif row_colors_data is not None:
                    legend_title = row_colors_data.name
                elif col_colors_data is not None:
                    legend_title = col_colors_data.name

                g.fig.legend(handles, labels, title=legend_title,
                             loc='center left', bbox_to_anchor=(0.9, 0.5),  # Adjust position
                             fontsize=8, title_fontsize=10)
                # Adjust layout to prevent overlap
                g.fig.subplots_adjust(right=0.85)  # Make space for legend

            except Exception as leg_e:
                print(f"Warning (Heatmap): Could not add legend: {leg_e}")

        # Add overall title
        g.fig.suptitle("Structure Similarity Matrix (RMSD)", fontsize=16, y=0.98)  # Use g.fig.suptitle

        # Save the figure if requested
        if output_file:
            print(f"DEBUG (Heatmap): Saving clustermap figure to {output_file}")
            # Save full visualization
            # Use bbox_inches='tight' carefully with clustermap legends
            g.fig.savefig(output_file, dpi=300, bbox_inches='tight')

            # --- Saving components separately (more complex with clustermap) ---
            # It might be easier to just save the main figure generated by clustermap.
            # If separate components are strictly needed, consider the sns.heatmap approach when order is fixed.
            print(f"Saved heatmap figure to {os.path.basename(output_file)}")
            # Optional: try saving components if needed, might require accessing g.ax_heatmap, g.ax_row/col_dendrogram etc.

        return g  # Return the clustermap object

    except Exception as e:
        print(f"ERROR (Heatmap): Failed to generate clustermap: {e}")
        import traceback

        traceback.print_exc()
        return None

def visualize_binding_pocket(structure, residue_ids, retinal=None, highlight_residues=None,
                           distance_cutoff=4.0, figsize=(800, 600)):
    """
    Create an interactive 3D visualization of a binding pocket with retinal.
    Args:
        structure: DataFrame with structure coordinates
        residue_ids: List of residue IDs to include in visualization
        retinal: DataFrame with retinal coordinates (optional)
        highlight_residues: Dict mapping residue IDs to highlight colors
        distance_cutoff: Distance cutoff for showing atom contacts
        figsize: Figure size tuple (width, height)
    Returns:
        Plotly figure object
    """
    # Filter structure to selected residues
    pocket_df = structure[structure['auth_seq_id'].isin(residue_ids)]
    
    # Set up traces for visualization
    traces = []
    
    # Add residues in pocket
    for res_id, res_df in pocket_df.groupby('auth_seq_id'):
        # Determine color based on highlighting or helix assignment
        if highlight_residues and res_id in highlight_residues:
            color = highlight_residues[res_id]
        elif 'helix_num' in res_df.columns and not res_df['helix_num'].isna().all():
            helix_num = res_df['helix_num'].iloc[0]
            if helix_num in HELIX_COLORS:
                color = HELIX_COLORS[helix_num]
            else:
                color = OPSIN_COLORS['gray_dark']
        else:
            color = OPSIN_COLORS['gray_dark']
        
        # Get residue type
        if 'res_name3l' in res_df.columns:
            res_type = res_df['res_name3l'].iloc[0]
        elif 'res_name1l' in res_df.columns:
            res_type = res_df['res_name1l'].iloc[0]
        else:
            res_type = '?'
        
        # Add trace for this residue
        traces.append(go.Scatter3d(
            x=res_df['x'],
            y=res_df['y'],
            z=res_df['z'],
            mode='markers',
            marker=dict(
                size=5,
                color=color,
                opacity=0.8
            ),
            text=[f"{res_type}{res_id} - {atom}" for atom in res_df['res_atom_name']],
            hoverinfo='text',
            name=f"{res_type}{res_id}"
        ))
    
    # Add retinal if provided
    if retinal is not None and not retinal.empty:
        traces.append(go.Scatter3d(
            x=retinal['x'],
            y=retinal['y'],
            z=retinal['z'],
            mode='markers',
            marker=dict(
                size=6,
                color=HELIX_COLORS['retinal'],
                symbol='diamond',
                opacity=0.9
            ),
            text=[f"RET - {atom}" for atom in retinal['res_atom_name']],
            hoverinfo='text',
            name="Retinal"
        ))
        
        # Add contact lines between retinal and nearby residues
        if distance_cutoff > 0:
            # Calculate distances between retinal and pocket atoms
            ret_coords = retinal[['x', 'y', 'z']].values
            pocket_coords = pocket_df[['x', 'y', 'z']].values
            
            # IDs for mapping back
            pocket_indices = pocket_df.index.tolist()
            
            # Calculate distances
            distances = cdist(ret_coords, pocket_coords)
            
            # Find contacts below cutoff
            contacts = np.where(distances < distance_cutoff)
            
            # Create line traces for contacts
            if len(contacts[0]) > 0:
                contact_x, contact_y, contact_z = [], [], []
                for ret_idx, pocket_idx in zip(contacts[0], contacts[1]):
                    # Get coordinates
                    rx, ry, rz = ret_coords[ret_idx]
                    px, py, pz = pocket_coords[pocket_idx]
                    
                    # Add line segment (with nan to separate segments)
                    contact_x.extend([rx, px, None])
                    contact_y.extend([ry, py, None])
                    contact_z.extend([rz, pz, None])
                
                # Add contact lines trace
                traces.append(go.Scatter3d(
                    x=contact_x,
                    y=contact_y,
                    z=contact_z,
                    mode='lines',
                    line=dict(
                        color=OPSIN_COLORS['yellow'],
                        width=2,
                        dash='dot'
                    ),
                    opacity=0.6,
                    hoverinfo='none',
                    name="Contacts",
                    showlegend=True
                ))
    
    # Create figure
    fig = go.Figure(data=traces)
    
    # Update layout
    fig.update_layout(
        title="Binding Pocket Visualization",
        scene=dict(
            xaxis_title='X',
            yaxis_title='Y',
            zaxis_title='Z',
            aspectmode='data',
            xaxis=dict(showbackground=True, backgroundcolor=OPSIN_COLORS['white']),
            yaxis=dict(showbackground=True, backgroundcolor=OPSIN_COLORS['white']),
            zaxis=dict(showbackground=True, backgroundcolor=OPSIN_COLORS['white'])
        ),
        width=figsize[0],
        height=figsize[1],
        margin=dict(l=0, r=0, b=0, t=40),
        legend=dict(
            x=0.01,
            y=0.99,
            bgcolor=OPSIN_COLORS['white'],
            bordercolor=OPSIN_COLORS['gray_light']
        )
    )
    
    return fig

def plot_distances_with_std(distance_table, title="Distance to Retinal by Position", figsize=(18, 10), use_ca=False):
    """
    Plot the average distances with standard deviation error bars for TM residues only (helices 1-7).

    Args:
        distance_table: DataFrame with distances
        title: Title for the plot
        figsize: Figure size tuple
        use_ca: Whether to label as CA distances or all-atom distances

    Returns:
        Matplotlib figure
    """
    # Handle empty table case
    if distance_table.empty:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No data available",
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.set_title(title, fontsize=16)
        return fig

    # Filter for TM residues only (positions with format N.XX where N is 1-7)
    tm_positions = []
    for pos in distance_table.columns:
        pos_str = str(pos)
        if '.' in pos_str:
            helix_num = pos_str.split('.')[0]
            if helix_num.isdigit() and 1 <= int(helix_num) <= 7:
                tm_positions.append(pos)

    # If no TM positions found, show warning
    if not tm_positions:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No TM residues found in data",
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        ax.set_title(title, fontsize=16)
        return fig

    # Create filtered table with only TM positions
    tm_distance_table = distance_table[tm_positions]

    # Calculate means and standard deviations
    mean_distances = tm_distance_table.mean(skipna=True)
    std_distances = tm_distance_table.std(skipna=True)

    # Create figure
    fig, ax = plt.subplots(figsize=figsize)

    # Simple sorting by helix number and position
    def custom_sort_key(pos):
        pos_str = str(pos)
        if '.' in pos_str:
            parts = pos_str.split('.')
            if parts[0].isdigit():
                return (int(parts[0]), float(parts[1]))
        return (999, 999)

    # Sort positions by helix number and position
    positions = sorted(mean_distances.index, key=custom_sort_key)

    # Prepare data for plotting
    x_vals = range(len(positions))
    y_vals = [mean_distances[pos] for pos in positions]
    err_vals = [std_distances[pos] for pos in positions]

    # Create point colors based on helix
    point_colors = []
    for pos in positions:
        helix_num = int(str(pos).split('.')[0])
        point_colors.append(HELIX_COLORS.get(helix_num, OPSIN_COLORS['gray_dark']))

    # Plot points with helix-specific colors
    for i, (x, y, color) in enumerate(zip(x_vals, y_vals, point_colors)):
        ax.plot(x, y, 'o', color=color, markersize=8)  # Slightly larger markers

    # Add error bars with increased visibility
    for i, (x, y, err, color) in enumerate(zip(x_vals, y_vals, err_vals, point_colors)):
        ax.errorbar(x, y, yerr=err, fmt='none', ecolor=color, alpha=0.8, capsize=5,
                   elinewidth=1.5, capthick=1.5)  # More visible error bars

    # Connect points within same helix with lines
    current_helix = None
    helix_x_vals = []
    helix_y_vals = []
    helix_color = None

    for i, pos in enumerate(positions):
        pos_str = str(pos)
        helix_num = int(pos_str.split('.')[0])

        if current_helix is None:
            current_helix = helix_num
            helix_color = HELIX_COLORS.get(helix_num, OPSIN_COLORS['gray_dark'])

        if helix_num == current_helix:
            helix_x_vals.append(x_vals[i])
            helix_y_vals.append(y_vals[i])
        else:
            # Plot connected line for previous helix
            ax.plot(helix_x_vals, helix_y_vals, '-', color=helix_color, linewidth=2.0, alpha=0.7)

            # Start new helix
            current_helix = helix_num
            helix_color = HELIX_COLORS.get(helix_num, OPSIN_COLORS['gray_dark'])
            helix_x_vals = [x_vals[i]]
            helix_y_vals = [y_vals[i]]

    # Plot the last helix line
    if helix_x_vals:
        ax.plot(helix_x_vals, helix_y_vals, '-', color=helix_color, linewidth=2.0, alpha=0.7)

    # Set custom x-axis labels
    ax.set_xticks(x_vals)
    ax.set_xticklabels(positions, rotation=90, fontsize=8)

    # Add reference lines for key positions (X.50)
    for i, pos in enumerate(positions):
        pos_str = str(pos)
        if pos_str.endswith('.50') and pos_str[0].isdigit():
            helix = pos_str.split('.')[0]
            if int(helix) in HELIX_COLORS:
                color = HELIX_COLORS[int(helix)]

                # Add vertical line
                ax.axvline(x=i, color=color, linestyle='--', alpha=0.6, linewidth=1.5)

                # Add helix label
                ax.text(i, max(y_vals) * 1.05, f"H{helix}",
                       ha='center', va='bottom', fontsize=12, color=color, fontweight='bold')

    # Add grid and labels
    ax.grid(True, alpha=0.3)
    ax.set_xlabel('Position (GRN Format)', fontsize=14)
    ax.set_ylabel('Distance to Retinal (Å)', fontsize=14)

    # Set appropriate y-axis limits
    if y_vals:
        max_y = max(y_vals) + max(err_vals) if err_vals else max(y_vals)
        ax.set_ylim(0, max_y * 1.1)
    else:
        ax.set_ylim(0, 10)  # Default if no data

    # Add a title
    subtitle = "CA Atoms" if use_ca else "All Atoms"
    ax.set_title(f"{title} ({subtitle})", fontsize=18)

    # Add legend for helix colors
    import matplotlib.patches as mpatches
    legend_handles = []
    for helix_num in range(1, 8):
        color = HELIX_COLORS.get(helix_num)
        legend_handles.append(mpatches.Patch(color=color, label=f"Helix {helix_num}"))

    ax.legend(handles=legend_handles, loc='upper right', ncol=2, fontsize=10, frameon=True, framealpha=0.8)

    plt.tight_layout()
    return fig

def plot_helix_logo_plots(residue_table, figsize=(20, 5)):
    """
    Create sequence logo plots for positions around X.50 for each of the 7 transmembrane helices.
    Shows amino acid composition variability for each position.
    
    Args:
        residue_table: DataFrame with residue information (ca_residue_table_grn)
        figsize: Figure size tuple (width, height)
        
    Returns:
        Matplotlib figure with 7 logo subplots
    """
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import logomaker
    from collections import Counter
    
    # Check if table is provided
    if residue_table is None or residue_table.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.text(0.5, 0.5, "No data available", 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        return fig
    
    # Define the window size (positions before and after X.50)
    window_size = 4
    
    # Create a figure with 7 subplots (one for each helix)
    fig, axes = plt.subplots(1, 7, figsize=figsize, sharey=True)
    
    # Define amino acid properties for coloring
    aa_properties = {
        'A': 'hydrophobic', 'I': 'hydrophobic', 'L': 'hydrophobic', 'M': 'hydrophobic',
        'F': 'hydrophobic', 'W': 'hydrophobic', 'V': 'hydrophobic', 'P': 'hydrophobic',
        'S': 'polar', 'T': 'polar', 'C': 'polar', 'N': 'polar', 'Q': 'polar', 'Y': 'polar',
        'D': 'negative', 'E': 'negative',
        'K': 'positive', 'R': 'positive', 'H': 'positive',
        'G': 'special'
    }
    
    aa_colors = {
        'hydrophobic': OPSIN_COLORS['blue_dark'],
        'polar': OPSIN_COLORS['green'],
        'negative': OPSIN_COLORS['red'],
        'positive': OPSIN_COLORS['orange'],
        'special': OPSIN_COLORS['purple']
    }
    
    # Filter residue table to only include helices 1-7 with proper GRN format
    valid_columns = []
    for col in residue_table.columns:
        col_str = str(col)
        # Only include columns with format N.XX where N is 1-7
        if '.' in col_str:
            try:
                helix_num = int(col_str.split('.')[0])
                if 1 <= helix_num <= 7:
                    valid_columns.append(col)
            except ValueError:
                # Skip columns with non-numeric helix values
                continue
                
    # Create filtered table with only valid TM helix positions
    filtered_table = residue_table[valid_columns]
    
    # Process each helix
    for helix_num in range(1, 8):
        ax = axes[helix_num - 1]
        
        # Get the X.50 position for this helix
        pivot_pos = f"{helix_num}.50"
        
        # Find all positions for this helix
        helix_positions = [col for col in filtered_table.columns 
                          if str(col).startswith(f"{helix_num}.")]
        
        # Sort positions by residue number
        def position_sort_key(pos):
            return float(str(pos).split('.')[1])
        
        helix_positions = sorted(helix_positions, key=position_sort_key)
        
        # Find the index of the pivot position in sorted positions
        try:
            pivot_idx = helix_positions.index(pivot_pos)
        except ValueError:
            # If pivot position not found, estimate it based on closest position
            closest_positions = sorted([(abs(position_sort_key(pos) - 50.0), pos) 
                                      for pos in helix_positions])
            if closest_positions:
                pivot_pos = closest_positions[0][1]
                pivot_idx = helix_positions.index(pivot_pos)
                print(f"[WARNING] Position {helix_num}.50 not found, using {pivot_pos} as pivot")
            else:
                # Skip this helix if no positions found
                ax.text(0.5, 0.5, f"No data for Helix {helix_num}", 
                        ha='center', va='center', transform=ax.transAxes)
                continue
        
        # Select positions within window around pivot
        start_idx = max(0, pivot_idx - window_size)
        end_idx = min(len(helix_positions), pivot_idx + window_size + 1)
        window_positions = helix_positions[start_idx:end_idx]
        
        # Create a position weight matrix with float values for frequencies
        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # 20 standard amino acids
        logo_df = pd.DataFrame(0.0, index=range(len(window_positions)), columns=list(amino_acids))
        
        # Store position mapping for reference
        position_mapping = {i: pos for i, pos in enumerate(window_positions)}
        
        # Count amino acid frequencies at each position
        for i, pos in enumerate(window_positions):
            # Process all cells at this position to extract amino acids
            aa_values = []
            for cell in filtered_table[pos].values:
                if pd.isna(cell) or cell == '-':
                    continue
                    
                # Extract just the first character for amino acid type
                if isinstance(cell, str) and len(cell) > 0:
                    aa = cell[0]
                    if aa in amino_acids:  # Only consider standard amino acids
                        aa_values.append(aa)
            
            # Count frequencies
            if aa_values:
                aa_counter = Counter(aa_values)
                total = sum(aa_counter.values())
                
                # Calculate frequencies
                if total > 0:
                    freq_dict = {aa: count / total for aa, count in aa_counter.items()}
                    
                    # Ensure frequencies sum to 1.0 (handling any floating point issues)
                    sum_freqs = sum(freq_dict.values())
                    if sum_freqs > 0 and abs(sum_freqs - 1.0) > 1e-6:
                        # Normalize to ensure sum is exactly 1.0
                        freq_dict = {aa: freq / sum_freqs for aa, freq in freq_dict.items()}
                    
                    # Fill the dataframe with normalized frequencies
                    for aa, freq in freq_dict.items():
                        logo_df.at[i, aa] = freq
        
        # Create position labels - distance from X.50
        position_labels = []
        pivot_value = float(pivot_pos.split('.')[1])
        
        for pos in window_positions:
            pos_value = float(str(pos).split('.')[1])
            offset = int(pos_value - pivot_value)
            position_labels.append(offset)
        
        # Create the logo plot
        logo = logomaker.Logo(logo_df, ax=ax, color_scheme='chemistry')
        
        # Configure subplot
        ax.set_title(f"Helix {helix_num}", color=HELIX_COLORS[helix_num], fontweight='bold')
        ax.set_xlabel("Offset from X.50", fontsize=10)
        ax.set_xticks(range(len(window_positions)))
        ax.set_xticklabels(position_labels, fontsize=9)
        
        # Add vertical line at X.50 position
        try:
            x50_idx = window_positions.index(pivot_pos)
            ax.axvline(x=x50_idx, color=HELIX_COLORS[helix_num], linestyle='--', alpha=0.5)
        except ValueError:
            pass
        
        # Highlight conserved residues (>80% conservation)
        for i, pos in enumerate(window_positions):
            max_freq = logo_df.loc[i].max()
            if max_freq > 0.8:  # Highly conserved
                max_aa = logo_df.loc[i].idxmax()
                ax.text(i, -0.1, max_aa, ha='center', fontweight='bold', fontsize=12)
    
    # Set consistent y-axis limits
    for ax in axes:
        ax.set_ylim(0, 2)
        
    # Add y-axis label only to the leftmost subplot
    axes[0].set_ylabel("Information Content", fontsize=12)
    
    # Add figure title
    plt.suptitle("Conserved Amino Acids around Position X.50 in Transmembrane Helices",
                fontsize=14, y=0.98)
    
    plt.tight_layout()
    return fig


def plot_average_distances_by_helix(distance_table, use_ca=True):
    """
    Plot the average distances per position, grouped by helix with consistent helix colors.
    Handles both dot notation (1.50) and x notation (1x50) GRN formats.
    
    Args:
        distance_table: DataFrame with distances
        use_ca: Boolean indicating whether the distances are from CA atoms (True) or sidechains (False)
    Returns:
        Matplotlib figure
    """
    # Calculate column means
    mean_distances = distance_table.mean(skipna=True)
    
    # Filter out loop regions (L.xx) and terminal regions (n.xx, c.xx)
    tm_distances = mean_distances[~mean_distances.index.str.startswith('L.') & 
                                  ~mean_distances.index.str.startswith('n.') &
                                  ~mean_distances.index.str.startswith('c.')]

    # Create a DataFrame for plotting
    plot_data = []
    for position, distance in tm_distances.items():
        # Handle both dot notation and x notation GRN formats
        if 'x' in position:
            # Handle 1x50 format
            helix, pos = position.split('x')
            plot_data.append({
                'Helix': helix,
                'Position': float(pos),
                'Distance (Å)': distance
            })
        elif '.' in position and not position.startswith('L.'):
            # Check if this is in helix.position format (like 1.50)
            parts = position.split('.')
            if len(parts) == 2 and parts[0].isdigit():
                helix, pos = parts
                plot_data.append({
                    'Helix': helix,
                    'Position': float(pos),
                    'Distance (Å)': distance
                })
    
    if not plot_data:
        # If no valid data was parsed, create a default plot
        fig, ax = plt.subplots(figsize=(14, 8))
        ax.text(0.5, 0.5, "No valid TM helix positions found in the data",
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=14)
        plt.tight_layout()
        return fig
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Plot each helix as a separate line with consistent helix colors
    for helix, group in plot_df.groupby('Helix'):
        sorted_group = group.sort_values('Position')
        
        # Get the color from our helix color map
        if helix in HELIX_COLORS_STR:
            color = HELIX_COLORS_STR[helix]
        else:
            # Fallback for non-standard helix numbers
            color = OPSIN_COLORS['gray_dark']
        
        ax.plot(sorted_group['Position'], sorted_group['Distance (Å)'],
              marker='o', linestyle='-', label=f'Helix {helix}',
              color=color, linewidth=2)
    
    # Add labels and title
    ax.set_xlabel('Position Number', fontsize=12)
    ax.set_ylabel('Average Distance to Retinal (Å)', fontsize=12)
    
    # Set title based on distance type
    distance_type = "Cα Atoms" if use_ca else "Sidechains"
    ax.set_title(f'Average Distance to Retinal by Residue Position ({distance_type})', fontsize=14)
    
    # Add vertical lines at the X.50 position for each helix
    for helix in sorted(plot_df['Helix'].unique()):
        if helix.isdigit() and int(helix) in HELIX_COLORS:
            ax.axvline(x=50, color=HELIX_COLORS[int(helix)], linestyle='--', alpha=0.3)
    
    # Add legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10)
    
    # Set grid for better readability
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # Add text annotation for X.50 positions
    ax.text(50, ax.get_ylim()[1] * 0.95, "X.50", 
            ha='center', va='top', fontsize=10, 
            bbox=dict(facecolor='white', alpha=0.8, edgecolor='gray'))
    
    # Set some plot styling
    plt.tight_layout()
    
    return fig

def plot_distance_heatmap(distance_table):
    """
    Create a heatmap of average distances ordered by helix and position using the consistent color scheme.
    Handles both dot notation (1.50) and x notation (1x50) GRN formats.
    
    Args:
        distance_table: DataFrame with distances
    Returns:
        Matplotlib figure
    """
    # Calculate column means
    mean_distances = distance_table.mean(skipna=True)
    
    # Filter out loop regions (L.xx) and terminal regions (n.xx, c.xx)
    tm_distances = mean_distances[~mean_distances.index.str.startswith('L.') & 
                                  ~mean_distances.index.str.startswith('n.') &
                                  ~mean_distances.index.str.startswith('c.')]
    
    # Create a DataFrame for plotting
    plot_data = []
    for position, distance in tm_distances.items():
        # Handle both dot notation and x notation GRN formats
        if 'x' in position:
            # Handle 1x50 format
            helix, pos = position.split('x')
            try:
                helix_num = int(helix)
                plot_data.append({
                    'Helix': helix_num,
                    'Position': float(pos),
                    'Distance (Å)': distance
                })
            except ValueError:
                # Skip non-numeric helix values
                continue
        elif '.' in position:
            # Check if this is in helix.position format (like 1.50)
            parts = position.split('.')
            if len(parts) == 2 and parts[0].isdigit():
                helix, pos = parts
                try:
                    helix_num = int(helix)
                    plot_data.append({
                        'Helix': helix_num,
                        'Position': float(pos),
                        'Distance (Å)': distance
                    })
                except ValueError:
                    # Skip non-numeric helix values
                    continue
    
    if not plot_data:
        # If no valid data was parsed, create a default plot
        fig, ax = plt.subplots(figsize=(16, 8))
        ax.text(0.5, 0.5, "No valid TM helix positions found in the data",
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=14)
        plt.tight_layout()
        return fig
    
    plot_df = pd.DataFrame(plot_data)
    
    # Create pivot table for heatmap
    pivot_df = plot_df.pivot(index='Helix', columns='Position', values='Distance (Å)')
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 8))
    
    # Create heatmap with our custom colormap - reversed for distances
    # where darker colors represent closer distances
    reversed_cmap = plt.cm.get_cmap(RMSD_CMAP).reversed()
    
    # Adjust vmin and vmax for better color contrast
    vmin = min(3.0, pivot_df.min().min())  # Minimum 3Å or actual minimum
    vmax = min(20.0, pivot_df.max().max())  # Maximum 20Å or actual maximum
    
    sns.heatmap(pivot_df, cmap=reversed_cmap, ax=ax,
              cbar_kws={'label': 'Distance to Retinal (Å)'},
              vmin=vmin, vmax=vmax)
    
    # Add title and labels
    ax.set_title('Distance to Retinal by Position', fontsize=14)
    ax.set_xlabel('Position Number', fontsize=12)
    ax.set_ylabel('Helix', fontsize=12)
    
    # Highlight the X.50 position
    if 50 in pivot_df.columns:
        ax.axvline(x=pivot_df.columns.get_loc(50) + 0.5,
                 color=OPSIN_COLORS['white'], linestyle='--', alpha=0.7)
    
    # Customize y-tick labels with helix colors
    plt.yticks(
        np.arange(len(pivot_df.index)) + 0.5,
        pivot_df.index,
        fontweight='bold'
    )
    
    # Apply helix-specific colors to y-tick labels
    for i, tick in enumerate(ax.get_yticklabels()):
        helix_num = pivot_df.index[i]
        if helix_num in HELIX_COLORS:
            tick.set_color(HELIX_COLORS[helix_num])
    
    return fig

def print_residue_composition(composition_dict, highlight_threshold=20.0):
    """
    Pretty-print the residue composition results with highlighting of significant values.
    Args:
        composition_dict: Dictionary from analyze_residue_composition
        highlight_threshold: Percentage threshold above which to highlight residues
    """
    print("\nResidue Composition at Key Positions:")
    print("=====================================")
    
    for position, composition in composition_dict.items():
        print(f"\nPosition {position}:")
        print("-" * (len(str(position)) + 10))
        
        # Skip positions that contain error messages
        if isinstance(composition, dict) and 'error' in composition:
            print(f"  {composition['error']}")
            continue
            
        # Handle case when frequencies are available
        if isinstance(composition, dict) and 'frequencies' in composition:
            # Use the frequencies directly
            for aa, freq in composition['frequencies'].items():
                percentage = freq * 100
                if percentage >= highlight_threshold:
                    print(f"  \033[1m{aa}: {percentage:.1f}%\033[0m")  # Bold in terminal
                else:
                    print(f"  {aa}: {percentage:.1f}%")
        # Handle case when sorted frequencies are available
        elif isinstance(composition, dict) and 'sorted' in composition:
            # Use the sorted list of (aa, freq) pairs
            for aa, freq in composition['sorted']:
                percentage = freq * 100
                if percentage >= highlight_threshold:
                    print(f"  \033[1m{aa}: {percentage:.1f}%\033[0m")  # Bold in terminal
                else:
                    print(f"  {aa}: {percentage:.1f}%")
        # Handle direct percentage values (older format)
        elif isinstance(composition, dict):
            for residue, value in composition.items():
                # Skip non-percentage values or sub-dictionaries
                if not isinstance(value, (int, float)):
                    continue
                    
                # Format with highlighting for high percentages
                if value >= highlight_threshold:
                    print(f"  \033[1m{residue}: {value:.1f}%\033[0m")  # Bold in terminal
                else:
                    print(f"  {residue}: {value:.1f}%")
        else:
            print("  [Error: Unexpected composition format]")

def calculate_helix_distances(distance_table):
    """
    Calculate the mean and standard deviation of distances to retinal for each position,
    grouped by helix, using the GRN column names (N.YY format).
    
    Args:
        distance_table: DataFrame with distances (columns are GRN positions)
        
    Returns:
        dict: Dictionary with helix statistics
    """
    # Calculate column means and std deviations
    mean_distances = distance_table.mean(skipna=True)
    std_distances = distance_table.std(skipna=True)
    
    # Group by helix
    helix_stats = {}
    
    for position in mean_distances.index:
        pos_str = str(position)
        
        # Process TM helix positions in N.YY format
        if '.' in pos_str and len(pos_str.split('.')) == 2:
            prefix, suffix = pos_str.split('.')
            
            # Check if this is a TM helix position
            if prefix.isdigit() and int(prefix) >= 1 and int(prefix) <= 7:
                helix = prefix
                pos_num = float(suffix)
                
                if helix not in helix_stats:
                    helix_stats[helix] = {'positions': [], 'means': [], 'stds': []}
                
                helix_stats[helix]['positions'].append(position)
                helix_stats[helix]['means'].append(mean_distances[position])
                helix_stats[helix]['stds'].append(std_distances[position])
                
    # Sort positions within each helix and find closest position
    for helix, data in helix_stats.items():
        # Sort by position number
        sorted_idx = sorted(range(len(data['positions'])), 
                          key=lambda i: float(str(data['positions'][i]).split('.')[1]))
        
        data['positions'] = [data['positions'][i] for i in sorted_idx]
        data['means'] = [data['means'][i] for i in sorted_idx]
        data['stds'] = [data['stds'][i] for i in sorted_idx]
        
        # Find the position with minimum mean distance
        if data['means']:
            min_idx = data['means'].index(min(data['means']))
            data['closest_position'] = data['positions'][min_idx]
            data['closest_mean'] = data['means'][min_idx]
        else:
            data['closest_position'] = None
            data['closest_mean'] = None
    
    return helix_stats

def visualize_msa_distances(msa_results, output_dir="./", prefix=""):
    """
    Creates visualization plots for MSA distance data with standard deviation.
    
    Args:
        msa_results: Results dictionary from generate_grn_msa_tables function
        output_dir: Directory to save plots
        prefix: Prefix for output filenames
        
    Returns:
        Dictionary with paths to generated plots
    """
    import os
    import matplotlib.pyplot as plt
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Access distance tables
    ca_distance_table = msa_results["ca_distance_table"]
    distance_table = msa_results["distance_table"]
    
    # Calculate helix stats
    helix_stats = calculate_helix_distances(distance_table)
    ca_helix_stats = calculate_helix_distances(ca_distance_table)
    
    # Print information about the closest residues
    print("Closest residues to retinal (sidechain atoms):")
    for helix in sorted(helix_stats.keys(), key=int):
        if helix_stats[helix]['closest_position']:
            print(f"Helix {helix}: Position {helix_stats[helix]['closest_position']} (distance: {helix_stats[helix]['closest_mean']:.2f}Å)")
    
    print("\nClosest residues to retinal (CA atoms):")
    for helix in sorted(ca_helix_stats.keys(), key=int):
        if ca_helix_stats[helix]['closest_position']:
            print(f"Helix {helix}: Position {ca_helix_stats[helix]['closest_position']} (distance: {ca_helix_stats[helix]['closest_mean']:.2f}Å)")
    
    # Create visualizations
    output_files = {}
    
    # Create plots
    # CA atoms plot
    ca_plot = plot_distances_with_std(
        ca_distance_table, 
        title="CA Atom Distances to Retinal", 
        use_ca=True
    )

    ca_path = os.path.join(output_dir, f"{prefix}ca_distances.png")
    ca_plot.savefig(ca_path, dpi=300, bbox_inches='tight')
    output_files["ca_distances"] = ca_path
    plt.close(ca_plot)
        
    # All atoms plot
    all_plot = plot_distances_with_std(
        distance_table, 
        title="All-Atom Distances to Retinal", 
        use_ca=False
    )
    all_path = os.path.join(output_dir, f"{prefix}all_distances.png")
    all_plot.savefig(all_path, dpi=300, bbox_inches='tight')
    output_files["all_distances"] = all_path
    plt.close(all_plot)
    
    # Create sequence logo plots for residues around position X.50 in each helix
    try:
        # Try to install logomaker if not available
        try:
            import logomaker
        except ImportError:
            print("[INFO] Installing logomaker package for sequence logo plots")
            import subprocess
            subprocess.check_call(["pip", "install", "logomaker"])
            import logomaker
            
        # Create logo plots
        print("[INFO] Creating sequence logo plots for positions around X.50")
        logo_plot = plot_helix_logo_plots(
            residue_table=msa_results["ca_residue_table"],
            figsize=(20, 5)
        )
        
        logo_path = os.path.join(output_dir, f"{prefix}sequence_logos.png")
        logo_plot.savefig(logo_path, dpi=300, bbox_inches='tight')
        output_files["sequence_logos"] = logo_path
        plt.close(logo_plot)
        print(f"[INFO] Sequence logo plots saved to {logo_path}")
    except Exception as e:
        print(f"[WARNING] Could not create sequence logo plots: {e}")
        
    print(f"Distance visualizations saved to {output_dir}")
    return output_files

def create_combined_distance_logo_plot(distance_table, msa_df):
    """
    Create a combined figure with a distance line plot on top and a sequence logo below,
    with precisely aligned x-axes for both dot notation (1.50) and x notation (1x50) GRN formats.
    
    Args:
        distance_table: DataFrame with distances for plotting (columns as GRN positions)
        msa_df: Multiple sequence alignment DataFrame (columns as GRN positions)
    Returns:
        Matplotlib figure
    """
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import numpy as np
    import logomaker
    
    # Create figure with gridspec for layout control
    fig = plt.figure(figsize=(15, 10))
    gs = gridspec.GridSpec(2, 1, height_ratios=[1, 0.6], hspace=0.05)
    
    # Create subplots
    ax_top = fig.add_subplot(gs[0])
    ax_bottom = fig.add_subplot(gs[1], sharex=ax_top)
    
    # ---- Process data for distance plot ----
    # Calculate mean distances across structures
    mean_distances = distance_table.mean()
    
    # Define a custom sorting function for GRN positions
    def grn_sort_key(pos):
        # Primary sorting by helix number, secondary by position within helix
        if isinstance(pos, str):
            # Handle N-terminal regions
            if pos.startswith('n.'):
                return (-1, float(pos.split('.')[1]) if len(pos.split('.')) > 1 else 0)
                
            # Handle TM helices with x notation (1x50)
            if 'x' in pos:
                parts = pos.split('x')
                if parts[0].isdigit() and parts[1].isdigit():
                    return (int(parts[0]), int(parts[1]))
            
            # Handle TM helices with dot notation (1.50)
            if '.' in pos and not pos.startswith('L.') and not pos.startswith('c.'):
                parts = pos.split('.')
                if len(parts) == 2 and parts[0].isdigit() and parts[1].isdigit():
                    return (int(parts[0]), int(parts[1]))
            
            # Handle loop regions (AB.CCC) - sort after TM helices
            if '.' in pos and len(pos) >= 6 and pos[0:2].isdigit():
                parts = pos.split('.')
                if parts[0].isdigit() and len(parts[0]) == 2:
                    return (int(parts[0][0]) * 10 + int(parts[0][1]), float('0.' + parts[1]))
            
            # Handle loop regions with L. prefix
            if pos.startswith('L.'):
                return (80, float(pos.split('.')[1]) if len(pos.split('.')) > 1 else 0)
                
            # Handle C-terminal regions
            if pos.startswith('c.'):
                return (90, float(pos.split('.')[1]) if len(pos.split('.')) > 1 else 0)
        
        # Default for unknown formats
        return (999, 999)
    
    # Sort GRN positions
    try:
        # Try to use the GRN utilities for sorting if available
        from protos.processing.grn.grn_utils import sort_grns_str
        sorted_grns = sort_grns_str(mean_distances.index)
    except ImportError:
        # Fall back to custom sorting logic
        sorted_grns = sorted(mean_distances.index, key=grn_sort_key)
    
    # Create plotting data
    x_values = np.arange(len(sorted_grns))
    y_values = [mean_distances[grn] for grn in sorted_grns]
    
    # Plot distance data as a single line
    ax_top.plot(x_values, y_values, 'o-', color='blue', linewidth=2)
    
    # Add labels and styling to top plot
    ax_top.set_title('Distance to Retinal and Sequence Conservation', fontsize=14)
    ax_top.set_ylabel('Average Distance to Retinal (Å)', fontsize=12)
    ax_top.grid(True, alpha=0.3)
    
    # Add reference lines for key positions (X.50 or Xx50)
    for i, grn in enumerate(sorted_grns):
        if ('x50' in grn) or ('.50' in grn and grn[0].isdigit()):
            ax_top.axvline(x=i, color='gray', linestyle='--', alpha=0.7)
            ax_top.text(i, ax_top.get_ylim()[1] * 0.95, grn,
                        rotation=90, ha='center', va='top', fontsize=8)
    
    # Hide x-tick labels on top plot
    plt.setp(ax_top.get_xticklabels(), visible=False)
    
    # ---- Process data for sequence logo ----
    # Extract common GRN positions between distance data and MSA
    common_grns = [grn for grn in sorted_grns if grn in msa_df.columns]
    
    # Create filtered MSA with only the positions we have distance data for
    filtered_msa = msa_df[common_grns]
    
    # Create position weight matrix for the logo
    # First extract just the first character from each sequence at each position
    aa_data = {}
    for grn in common_grns:
        position_idx = sorted_grns.index(grn)
        aa_counts = {}
        total = 0
        
        # Count amino acids at this position
        for aa_str in filtered_msa[grn]:
            # Skip empty, None, NaN, or '-' values
            if not isinstance(aa_str, str) or not aa_str or aa_str[0] == '-':
                continue
            
            aa = aa_str[0]
            aa_counts[aa] = aa_counts.get(aa, 0) + 1
            total += 1
        
        # Store the counts for this position
        aa_data[position_idx] = (aa_counts, total)
    
    # Create a proper position weight matrix for Logomaker
    # with positions as rows and amino acids as columns
    # Include all standard amino acids but exclude '-' for the PWM
    amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
    
    import pandas as pd
    pwm_data = {}
    
    # Initialize with zeros
    for aa in amino_acids:
        pwm_data[aa] = [0.0] * len(aa_data)
    
    # Fill in frequencies
    for idx, (aa_counts, total) in aa_data.items():
        if total > 0:
            for aa in amino_acids:
                # Skip the gap character '-' entirely
                if aa == '-':
                    pwm_data[aa][idx] = 0.0
                elif aa in aa_counts:
                    pwm_data[aa][idx] = aa_counts[aa] / total
    
    # Create the PWM DataFrame with numeric indices
    pwm = pd.DataFrame(pwm_data, index=range(len(aa_data)))
    
    # Create sequence logo in bottom plot
    if not pwm.empty and pwm.sum().sum() > 0:
        # Use OPSIN color scheme for amino acids
        aa_colors = {
            'A': OPSIN_COLORS['blue_light'],
            'C': OPSIN_COLORS['yellow'],
            'D': OPSIN_COLORS['red'],
            'E': OPSIN_COLORS['red'],
            'F': OPSIN_COLORS['blue_dark'],
            'G': OPSIN_COLORS['orange'],
            'H': OPSIN_COLORS['teal'],
            'I': OPSIN_COLORS['blue_light'],
            'K': OPSIN_COLORS['blue_dark'],
            'L': OPSIN_COLORS['blue_light'],
            'M': OPSIN_COLORS['blue_light'],
            'N': OPSIN_COLORS['green'],
            'P': OPSIN_COLORS['yellow'],
            'Q': OPSIN_COLORS['green'],
            'R': OPSIN_COLORS['blue_dark'],
            'S': OPSIN_COLORS['green'],
            'T': OPSIN_COLORS['green'],
            'V': OPSIN_COLORS['blue_light'],
            'W': OPSIN_COLORS['blue_dark'],
            'Y': OPSIN_COLORS['blue_dark'],
            '-': OPSIN_COLORS['gray_light']
        }
        
        try:
            # Create the logo with logomaker
            logo = logomaker.Logo(pwm, ax=ax_bottom, color_scheme=aa_colors)
            
            # Style the logo
            logo.style_spines(visible=False)
            logo.style_xticks(rotation=90, fmt='%d')
            
            # Map indices back to GRN positions for x-axis labels
            idx_to_grn = {sorted_grns.index(grn): grn for grn in common_grns}
            ax_bottom.set_xticks(list(idx_to_grn.keys()))
            ax_bottom.set_xticklabels([idx_to_grn[i] for i in sorted(idx_to_grn.keys())], rotation=90)
        except Exception as e:
            # Fallback if logomaker fails
            print(f"Error creating logo: {e}")
            
            # Create a placeholder bar chart instead
            ax_bottom.bar(range(len(common_grns)), [0.5] * len(common_grns), alpha=0.2)
            ax_bottom.set_xticks(range(len(common_grns)))
            ax_bottom.set_xticklabels(common_grns, rotation=90)
    else:
        # Fallback if we have no valid PWM data
        ax_bottom.text(0.5, 0.5, "Insufficient data for sequence logo",
                       ha='center', va='center', transform=ax_bottom.transAxes)
        ax_bottom.set_xticks(x_values)
        ax_bottom.set_xticklabels(sorted_grns, rotation=90)
    
    # Add labels to bottom plot
    ax_bottom.set_xlabel('Position (GRN)', fontsize=12)
    ax_bottom.set_ylabel('Information Content', fontsize=12)
    
    # Add reference lines to match top plot
    for i, grn in enumerate(sorted_grns):
        if (('x50' in grn) or ('.50' in grn and grn[0].isdigit())) and i in x_values:
            ax_bottom.axvline(x=i, color='gray', linestyle='--', alpha=0.7)
    
    # Set the same x-limits for both plots
    xlim = (-0.5, len(sorted_grns) - 0.5)
    ax_top.set_xlim(xlim)
    ax_bottom.set_xlim(xlim)
    
    # Adjust layout
    plt.tight_layout()
    
    return fig


def create_opsin_overview_plot(df, output_path=None, figsize=(16, 16)):
    """
    Create a circular overview plot of opsin structures with multiple rings showing 
    molecular function, domain, and experimental/predicted status.
    
    Args:
        df: DataFrame with opsin information, must include columns:
            - molecular_function_normalized: Function category
            - domain: Taxonomic domain
            - experimentally_determined: Boolean flag for experimental structures
            - short_name: Short identifier for each structure
        output_path: Path to save the figure (optional)
        figsize: Figure size tuple (width, height)
        
    Returns:
        Matplotlib figure with circular plot
    """
    # Ensure DataFrame is properly formatted
    df = df.copy()
    df['molecular_function_normalized'] = df['molecular_function_normalized'].fillna('Unknown')
    df['domain'] = df['domain'].fillna('Unknown')
    df_sorted = df.sort_values('molecular_function_normalized', ignore_index=True)
    N = len(df_sorted)
    
    # Define consistent colors for molecular functions using our global color scheme
    function_colors_dict = get_group_colors(df_sorted['molecular_function_normalized'].unique())
    
    # Define consistent colors for domains
    domain_colors_dict = {
        'Bacteria':   OPSIN_COLORS['red'],
        'Eukaryota':  OPSIN_COLORS['green'],
        'Archaea':    OPSIN_COLORS['yellow'],
        'Virus':      OPSIN_COLORS['blue_dark'],
        'Synthetic':  OPSIN_COLORS['orange'],
        'Unknown':    OPSIN_COLORS['gray_dark']
    }
    
    # Prepare colors for each ring
    ring1_colors = [function_colors_dict.get(f, OPSIN_COLORS['gray_light']) for f in df_sorted['molecular_function_normalized']]
    ring2_colors = [domain_colors_dict.get(d, OPSIN_COLORS['gray_light']) for d in df_sorted['domain']]
    ring3_labels = ['•' if x else '' for x in df_sorted['experimentally_determined']]
    ring4_labels = ['•'] * N  # All predicted
    outer_labels = df_sorted['short_name'].tolist()
    
    # Create figure & axes
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes([0.05, 0.1, 2, 0.65])
    
    # Define radii for each ring
    R1_RADIUS, R2_RADIUS, R3_RADIUS, R4_RADIUS, R5_RADIUS = 0.65, 0.85, 0.93, 0.98, 1.06
    SLICE_WIDTH = 0.18
    DOT_WIDTH   = 0.03
    
    # Plot the rings
    wedges1, _ = ax.pie(
        [1]*N, colors=ring1_colors, radius=R1_RADIUS,
        wedgeprops=dict(width=SLICE_WIDTH, edgecolor='white'), startangle=0
    )
    wedges2, _ = ax.pie(
        [1]*N, colors=ring2_colors, radius=R2_RADIUS,
        wedgeprops=dict(width=SLICE_WIDTH, edgecolor='white'), startangle=0
    )
    wedges3, texts3 = ax.pie(
        [1]*N, labels=ring3_labels, radius=R3_RADIUS, labeldistance=1.0,
        wedgeprops=dict(width=DOT_WIDTH, edgecolor='white', color='white'), startangle=0
    )
    for t in texts3:
        t.set_size(16)
        t.set_color('black')
    
    wedges4, texts4 = ax.pie(
        [1]*N, labels=ring4_labels, radius=R4_RADIUS, labeldistance=1.0,
        wedgeprops=dict(width=DOT_WIDTH, edgecolor='white', color='white'), startangle=0
    )
    for t in texts4:
        t.set_size(16)
        t.set_color('green')
    
    wedges5, texts5 = ax.pie(
        [1]*N, labels=outer_labels, radius=R5_RADIUS, labeldistance=1.1,
        wedgeprops=dict(width=DOT_WIDTH, edgecolor='white', color='white'), startangle=0
    )
    for i, text in enumerate(texts5):
        text.set_size(8)
        theta = np.deg2rad(wedges5[i].theta1 + (wedges5[i].theta2 - wedges5[i].theta1)/2)
        rotation = np.rad2deg(theta)
        if 90 < rotation <= 270:
            rotation -= 180
        text.set_rotation(rotation)
        text.set_color('black')
    
    ax.set_title("Opsin Structures", pad=20, fontsize=18)
    
    # Add legends
    unique_funcs = df_sorted['molecular_function_normalized'].unique()
    unique_domains = df_sorted['domain'].unique()
    """
    # Function legend
    func_legend = ax.legend(
        [
            plt.Line2D([0],[0], color=function_colors_dict.get(f, OPSIN_COLORS['gray_light']),
                       marker='o', label=f, markersize=8, linestyle='None')
            for f in unique_funcs
        ],
        unique_funcs,
        title="Molecular Function",
        title_fontsize=10,
        fontsize=8,
        loc='center left',
        bbox_to_anchor=(1.05, 0.95),
        ncol=2,
        labelspacing=0.4,
        columnspacing=0.8,
        handletextpad=0.5,
    )
    ax.add_artist(func_legend)
    
    # Domain legend
    dom_legend = ax.legend(
        [
            plt.Line2D([0],[0], color=domain_colors_dict.get(d, OPSIN_COLORS['gray_light']),
                       marker='o', label=d, markersize=8, linestyle='None')
            for d in unique_domains
        ],
        unique_domains,
        title="Domain",
        title_fontsize=10,
        fontsize=8,
        loc='center left',
        bbox_to_anchor=(1.05, 0.8),
        ncol=1,
        labelspacing=0.4,
        handletextpad=0.5
    )
    ax.add_artist(dom_legend)
    
    # Structure type legend
    struct_legend = ax.legend(
        [
            plt.Line2D([0],[0], color='black', marker='o',
                       label='Experimental (• if True)', markersize=8, linestyle='None'),
            plt.Line2D([0],[0], color='green', marker='o',
                       label='Predicted (• for all)', markersize=8, linestyle='None')
        ],
        ['Experimental','Predicted'],
        title="Structure",
        title_fontsize=10,
        fontsize=8,
        loc='center left',
        bbox_to_anchor=(1.05, 0.7),
        ncol=1,
        labelspacing=0.4,
        handletextpad=0.5
    )
    ax.add_artist(struct_legend)
    """
    print("hello from opsin overview")
    # Save the figure if output path is provided
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
    
    return fig

def plot_conservation_around_x50(residue_table, conservation_scores=None, figsize=(20, 10)):
    """
    Create sequence logo plots for 4 residues around X.50 positions in each helix (1-7),
    based on residue conservation.
    
    Args:
        residue_table: DataFrame with residue information (ca_residue_table_grn)
        conservation_scores: Dictionary of conservation scores by position (optional)
        figsize: Figure size tuple (width, height)
        
    Returns:
        Matplotlib figure with 7 logo subplots showing conservation
    """
    import pandas as pd
    import numpy as np
    import matplotlib.pyplot as plt
    import logomaker
    from collections import Counter
    
    # Check if table is provided
    if residue_table is None or residue_table.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.text(0.5, 0.5, "No data available", 
                ha='center', va='center', transform=ax.transAxes, fontsize=14)
        return fig
    
    # Define the window size (positions before and after X.50)
    window_size = 4
    
    # Create a figure with 7 subplots (one for each helix)
    fig, axes = plt.subplots(1, 7, figsize=figsize, sharey=True)
    
    # Define amino acid properties for coloring
    aa_properties = {
        'A': 'hydrophobic', 'I': 'hydrophobic', 'L': 'hydrophobic', 'M': 'hydrophobic',
        'F': 'hydrophobic', 'W': 'hydrophobic', 'V': 'hydrophobic', 'P': 'hydrophobic',
        'S': 'polar', 'T': 'polar', 'C': 'polar', 'N': 'polar', 'Q': 'polar', 'Y': 'polar',
        'D': 'negative', 'E': 'negative',
        'K': 'positive', 'R': 'positive', 'H': 'positive',
        'G': 'special'
    }
    
    aa_colors = {
        'hydrophobic': OPSIN_COLORS['blue_dark'],
        'polar': OPSIN_COLORS['green'],
        'negative': OPSIN_COLORS['red'],
        'positive': OPSIN_COLORS['orange'],
        'special': OPSIN_COLORS['purple']
    }
    
    # Filter residue table to only include helices 1-7 with proper GRN format
    valid_columns = []
    for col in residue_table.columns:
        col_str = str(col)
        # Only include columns with format N.XX where N is 1-7
        if '.' in col_str:
            try:
                helix_num = int(col_str.split('.')[0])
                if 1 <= helix_num <= 7:
                    valid_columns.append(col)
            except ValueError:
                # Skip columns with non-numeric helix values
                continue
                
    # Create filtered table with only valid TM helix positions
    filtered_table = residue_table[valid_columns]
    
    # Calculate conservation scores if not provided
    if conservation_scores is None:
        conservation_scores = {}
        for pos in filtered_table.columns:
            # Process all cells at this position to extract amino acids
            aa_values = []
            for cell in filtered_table[pos].values:
                if pd.isna(cell) or cell == '-':
                    continue
                    
                # Extract just the first character for amino acid type
                if isinstance(cell, str) and len(cell) > 0:
                    aa = cell[0]
                    if aa in 'ACDEFGHIKLMNPQRSTVWY':  # Only consider standard amino acids
                        aa_values.append(aa)
            
            # Count frequencies and get conservation score
            if aa_values:
                aa_counter = Counter(aa_values)
                total = sum(aa_counter.values())
                max_freq = max(aa_counter.values()) / total
                conservation_scores[pos] = max_freq
            else:
                # Handle positions with no valid amino acids
                conservation_scores[pos] = 0.0
    
    # Process each helix
    for helix_num in range(1, 8):
        ax = axes[helix_num - 1]
        
        # Get the X.50 position for this helix
        pivot_pos = f"{helix_num}.50"
        
        # Find all positions for this helix
        helix_positions = [col for col in filtered_table.columns 
                          if str(col).startswith(f"{helix_num}.")]
        
        # Sort positions by residue number
        def position_sort_key(pos):
            return float(str(pos).split('.')[1])
        
        helix_positions = sorted(helix_positions, key=position_sort_key)
        
        # Find the index of the pivot position in sorted positions
        try:
            pivot_idx = helix_positions.index(pivot_pos)
        except ValueError:
            # If pivot position not found, estimate it based on closest position
            closest_positions = sorted([(abs(position_sort_key(pos) - 50.0), pos) 
                                      for pos in helix_positions])
            if closest_positions:
                pivot_pos = closest_positions[0][1]
                pivot_idx = helix_positions.index(pivot_pos)
                print(f"[WARNING] Position {helix_num}.50 not found, using {pivot_pos} as pivot")
            else:
                # Skip this helix if no positions found
                ax.text(0.5, 0.5, f"No data for Helix {helix_num}", 
                        ha='center', va='center', transform=ax.transAxes)
                continue
        
        # Select positions within window around pivot
        start_idx = max(0, pivot_idx - window_size)
        end_idx = min(len(helix_positions), pivot_idx + window_size + 1)
        window_positions = helix_positions[start_idx:end_idx]
        
        # Create a position weight matrix with float values for frequencies
        amino_acids = 'ACDEFGHIKLMNPQRSTVWY'  # 20 standard amino acids
        logo_df = pd.DataFrame(0.0, index=range(len(window_positions)), columns=list(amino_acids))
        
        # Store position mapping for reference
        position_mapping = {i: pos for i, pos in enumerate(window_positions)}
        
        # For each position in the window, determine amino acid frequencies
        for i, pos in enumerate(window_positions):
            # Process all cells at this position to extract amino acids
            aa_values = []
            for cell in filtered_table[pos].values:
                if pd.isna(cell) or cell == '-':
                    continue
                    
                # Extract just the first character for amino acid type
                if isinstance(cell, str) and len(cell) > 0:
                    aa = cell[0]
                    if aa in amino_acids:  # Only consider standard amino acids
                        aa_values.append(aa)
            
            # Count frequencies and normalize
            if aa_values:
                aa_counter = Counter(aa_values)
                total = sum(aa_counter.values())
                
                # Calculate raw frequencies
                raw_freqs = {aa: count / total for aa, count in aa_counter.items()}
                
                # Get conservation score
                cons_score = conservation_scores.get(pos, 0)
                
                # Weight the frequencies by conservation
                # Higher conservation = larger letter height
                weighted_freqs = {aa: freq * cons_score * 2 for aa, freq in raw_freqs.items()}
                
                # Normalize the weighted frequencies to sum to 1.0
                if weighted_freqs:
                    weight_sum = sum(weighted_freqs.values())
                    if weight_sum > 0:
                        normalized_freqs = {aa: freq / weight_sum for aa, freq in weighted_freqs.items()}
                        
                        # Fill the dataframe with normalized frequencies
                        for aa, freq in normalized_freqs.items():
                            logo_df.at[i, aa] = freq
        
        # Create position labels - distance from X.50
        position_labels = []
        pivot_value = float(pivot_pos.split('.')[1])
        
        for pos in window_positions:
            pos_value = float(str(pos).split('.')[1])
            offset = int(pos_value - pivot_value)
            position_labels.append(offset)
        
        # Create the logo plot
        logo = logomaker.Logo(logo_df, ax=ax, color_scheme='chemistry')
        
        # Configure subplot
        ax.set_title(f"Helix {helix_num}", color=HELIX_COLORS[helix_num], fontweight='bold')
        ax.set_xlabel("Offset from X.50", fontsize=10)
        ax.set_xticks(range(len(window_positions)))
        ax.set_xticklabels(position_labels, fontsize=9)
        
        # Add vertical line at X.50 position
        try:
            x50_idx = window_positions.index(pivot_pos)
            ax.axvline(x=x50_idx, color=HELIX_COLORS[helix_num], linestyle='--', alpha=0.5)
        except ValueError:
            pass
        
        # Add conservation score labels below the plot
        for i, pos in enumerate(window_positions):
            cons_score = conservation_scores.get(pos, 0)
            if cons_score > 0:
                # Show percent conservation
                ax.text(i, -0.1, f"{int(cons_score*100)}%", 
                        ha='center', fontsize=8, rotation=90, va='top')
                # Add the most frequent amino acid
                if cons_score > 0.4:  # Only show for reasonably conserved positions
                    max_aa = logo_df.loc[i].idxmax()
                    ax.text(i, -0.5, max_aa, ha='center', fontweight='bold', fontsize=10)
    
    # Set consistent y-axis limits
    for ax in axes:
        ax.set_ylim(0, 2)
        
    # Add y-axis label only to the leftmost subplot
    axes[0].set_ylabel("Conservation-Weighted\nInformation Content", fontsize=12)
    
    # Add figure title
    plt.suptitle("Amino Acid Conservation around X.50 Positions (TM Helices 1-7)",
                fontsize=14, y=0.98)
    
    plt.tight_layout()
    return fig
