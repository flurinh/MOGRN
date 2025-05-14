"""
Opsin structure visualization color scheme.
This module defines consistent color palettes for all opsin structure visualizations.
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import seaborn as sns
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

# =============================================================================
# MAIN COLOR PALETTES
# =============================================================================

# Core color palette - used as the basis for many visualizations
# Based on the AlphaFold confidence color scheme with some modifications
OPSIN_COLORS = {
    'blue_dark': '#0053D6',  # Deep blue - most reliable/closest
    'blue_light': '#65CBF3',  # Light blue - reliable/close
    'green': '#00D45A',  # Green - moderately reliable/medium distance
    'yellow': '#FFD662',  # Yellow - less reliable/distant
    'orange': '#FF7D45',  # Orange - uncertain/very distant
    'red': '#FF4C4C',  # Red - least reliable/furthest
    'purple': '#9D6CFF',  # Purple - Used for special highlights
    'teal': '#00B5AD',  # Teal - Alternative highlight
    'gray_light': '#E0E0E0',  # Light gray - background elements
    'gray_dark': '#404040',  # Dark gray - text and lines
    'white': '#FFFFFF',  # White - background
    'black': '#000000'  # Black - text and outlines
}

# =============================================================================
# DERIVED COLOR MAPS
# =============================================================================

# Sequential colormap for RMSD and distance values (blue → yellow → red)
RMSD_COLORS = [
    OPSIN_COLORS['blue_dark'],
    OPSIN_COLORS['blue_light'],
    OPSIN_COLORS['green'],
    OPSIN_COLORS['yellow'],
    OPSIN_COLORS['orange'],
    OPSIN_COLORS['red']
]

# Create a proper colormap for continuous values
RMSD_CMAP = LinearSegmentedColormap.from_list('opsin_rmsd', RMSD_COLORS)

# Custom RMSD colormap with limited variation (dark blue -> light blue -> yellow)
# Dark blue up to 1.5, light blue up to 3.0, yellow above
RMSD_COMPACT_COLORS = [
    OPSIN_COLORS['blue_dark'],   # 0.0
    OPSIN_COLORS['blue_dark'],   # 0.75
    OPSIN_COLORS['blue_light'],  # 1.5
    OPSIN_COLORS['blue_light'],  # 2.25
    OPSIN_COLORS['yellow'],      # 3.0+
]
RMSD_COMPACT_CMAP = LinearSegmentedColormap.from_list('opsin_rmsd_compact', RMSD_COMPACT_COLORS)

# For diverging color schemes (e.g., deviation from average)
DIVERGING_COLORS = [
    OPSIN_COLORS['blue_dark'],
    OPSIN_COLORS['blue_light'],
    OPSIN_COLORS['white'],
    OPSIN_COLORS['orange'],
    OPSIN_COLORS['red']
]

DIVERGING_CMAP = LinearSegmentedColormap.from_list('opsin_diverging', DIVERGING_COLORS)

# Discrete colormap for specific RMSD ranges
RMSD_BOUNDS = [0, 1.5, 3.0, 5.0]  # Updated bounds to match new color scheme
RMSD_DISCRETE_CMAP = ListedColormap([OPSIN_COLORS['blue_dark'], OPSIN_COLORS['blue_light'], OPSIN_COLORS['yellow']])  # Simplified color scheme

# =============================================================================
# HELIX COLORS
# =============================================================================

# Consistent colors for the 7 transmembrane helices
HELIX_COLORS = {
    1: OPSIN_COLORS['blue_dark'],  # Helix 1
    2: OPSIN_COLORS['blue_light'],  # Helix 2
    3: OPSIN_COLORS['green'],  # Helix 3
    4: OPSIN_COLORS['teal'],  # Helix 4
    5: OPSIN_COLORS['yellow'],  # Helix 5
    6: OPSIN_COLORS['orange'],  # Helix 6
    7: OPSIN_COLORS['red'],  # Helix 7
    'retinal': OPSIN_COLORS['purple']  # Retinal
}

# String versions for when helix numbers are strings
HELIX_COLORS_STR = {str(k): v for k, v in HELIX_COLORS.items() if isinstance(k, int)}
HELIX_COLORS_STR['retinal'] = HELIX_COLORS['retinal']

# List version for sequential access
HELIX_COLORS_LIST = [HELIX_COLORS[i] for i in range(1, 8)]

# =============================================================================
# GROUP COLORS
# =============================================================================

# Pre-defined colors for up to 10 different protein groups
GROUP_COLORS = {
    'Rhodopsin': OPSIN_COLORS['blue_dark'],
    'Cone_Opsin': OPSIN_COLORS['blue_light'],
    'Visual_Opsin': OPSIN_COLORS['teal'],
    'Non_visual_Opsin': OPSIN_COLORS['green'],
    'Melanopsin': OPSIN_COLORS['yellow'],
    'Photoisomerase': OPSIN_COLORS['orange'],
    'Cnidopsin': OPSIN_COLORS['red'],
    'Bistable': OPSIN_COLORS['purple'],
    'Other': OPSIN_COLORS['gray_dark'],
    'Unknown': OPSIN_COLORS['gray_light']
}


# Function to dynamically assign colors to arbitrary group names
def get_group_colors(group_names):
    """
    Dynamically assign colors to group names.
    Uses predefined colors if available, otherwise creates new colors.

    Args:
        group_names: List of group names

    Returns:
        Dictionary mapping group names to colors
    """
    # Start with predefined colors
    group_colors = {}

    # Get unique groups and sort them for consistency
    unique_groups = sorted(set(group_names))

    # Assign colors
    for group in unique_groups:
        if group in GROUP_COLORS:
            group_colors[group] = GROUP_COLORS[group]
        else:
            # If not predefined, we'll assign dynamically
            pass

    # For groups without predefined colors, use a colormap
    remaining_groups = [g for g in unique_groups if g not in group_colors]
    if remaining_groups:
        n_remaining = len(remaining_groups)
        if n_remaining <= 6:
            # Use main colors
            remaining_colors = [OPSIN_COLORS['blue_dark'], OPSIN_COLORS['green'],
                                OPSIN_COLORS['yellow'], OPSIN_COLORS['orange'],
                                OPSIN_COLORS['red'], OPSIN_COLORS['purple']][:n_remaining]
        else:
            # For more groups, generate from colormap
            remaining_colors = [RMSD_CMAP(i / (n_remaining - 1)) for i in range(n_remaining)]

        # Assign remaining colors
        for i, group in enumerate(remaining_groups):
            group_colors[group] = remaining_colors[i]

    return group_colors


# =============================================================================
# UTILITIES
# =============================================================================



def get_helix_cmap():
    """Get a discrete colormap for the 7 transmembrane helices"""
    return ListedColormap(HELIX_COLORS_LIST)


def register_colormaps():
    """
    Register custom colormaps with matplotlib.
    This function is compatible with different matplotlib versions.
    """
    cmap_dict = {
        'opsin_rmsd': RMSD_CMAP,
        'opsin_diverging': DIVERGING_CMAP,
        'opsin_discrete': RMSD_DISCRETE_CMAP
    }

    try:
        # Method for newer matplotlib versions
        for name, cmap in cmap_dict.items():
            plt.colormaps.register(name=name, cmap=cmap)
    except AttributeError:
        try:
            # Method for older matplotlib versions
            from matplotlib.cm import register_cmap
            for name, cmap in cmap_dict.items():
                register_cmap(name=name, cmap=cmap)
        except (ImportError, AttributeError):
            # Fallback method
            print("[WARNING] Could not register custom colormaps. Visualizations will use default colors.")


# Register the colormaps when the module is imported
register_colormaps()

# Modified application of the style to prevent errors
def apply_style():
    """Apply consistent styling to matplotlib"""
    plt.style.use('default')  # Reset to default first

    style_params = {
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10,
        'axes.labelsize': 12,
        'axes.titlesize': 14,
        'xtick.labelsize': 10,
        'ytick.labelsize': 10,
        'legend.fontsize': 10,
        'figure.titlesize': 16,
        'axes.spines.top': False,
        'axes.spines.right': False,
        'axes.grid': True,
        'grid.alpha': 0.3,
        'lines.linewidth': 2,
        'axes.labelpad': 8,
        'axes.titlepad': 10,
        'legend.frameon': False,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight'
    }

    try:
        plt.rcParams.update(style_params)
        # Use seaborn for better aesthetics but don't change our color schemes
        sns.set_style('ticks')
    except Exception as e:
        print(f"[WARNING] Could not fully apply style settings: {e}")

    # Return the style dict in case it's needed
    return plt.rcParams.copy()


apply_style()


# =============================================================================
# EXAMPLE USAGE
# =============================================================================

if __name__ == "__main__":
    # Demo of the color schemes
    import matplotlib.pyplot as plt

    # Example plot to show the colormaps
    fig, axs = plt.subplots(3, 1, figsize=(10, 8))

    # Show continuous RMSD colormap
    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    axs[0].imshow(gradient, aspect='auto', cmap=RMSD_CMAP)
    axs[0].set_title('RMSD Colormap')
    axs[0].set_yticks([])

    # Show diverging colormap
    axs[1].imshow(gradient, aspect='auto', cmap=DIVERGING_CMAP)
    axs[1].set_title('Diverging Colormap')
    axs[1].set_yticks([])

    # Show helix colors
    helix_colors = np.array([[i] for i in range(7)]).reshape(1, -1)
    axs[2].imshow(helix_colors, aspect='auto', cmap=get_helix_cmap())
    axs[2].set_title('Helix Colors')
    axs[2].set_yticks([])
    axs[2].set_xticks(range(7))
    axs[2].set_xticklabels([f'Helix {i + 1}' for i in range(7)])

    plt.tight_layout()
    plt.show()

    # Example of group colors
    example_groups = ['Rhodopsin', 'Cone_Opsin', 'NewGroup1', 'NewGroup2']
    group_colors = get_group_colors(example_groups)

    plt.figure(figsize=(8, 4))
    for i, (group, color) in enumerate(group_colors.items()):
        plt.bar(i, 1, color=color, label=group)

    plt.legend()
    plt.title('Group Color Assignment Example')
    plt.xticks([])
    plt.yticks([])
    plt.show()