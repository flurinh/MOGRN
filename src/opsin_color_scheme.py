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
# MAIN COLOR PALETTES (Core Hex Codes)
# =============================================================================

OPSIN_COLORS = {
    # Domain Palette: Blue -> Cyan -> Green
    'domain_blue_dark': '#08306B',    # Deep Blue
    'domain_blue_medium': '#08519C',
    'domain_blue_light': '#2171B5',
    'domain_cyan_medium': '#41B6C4',   # Cyan
    'domain_cyan_light': '#7FCDBB',
    'domain_green_medium': '#4CAF50',   # Green
    'domain_green_light': '#A1D99B',

    # Molecular Function Palette: Yellow -> Orange -> Red -> Dark Purple
    'mol_func_yellow_light': '#FFFFD4', # Lightest Yellow
    'mol_func_yellow_medium': '#FED976',
    'mol_func_orange_light': '#FEB24C',
    'mol_func_orange_medium': '#FD8D3C',
    'mol_func_red_medium': '#FC4E2A',
    'mol_func_red_dark': '#E31A1C',
    'mol_func_purple_medium': '#800080', # Purple
    'mol_func_purple_dark': '#4B0082',   # Indigo/Dark Purple

    # Grayscale: White to Black (Unchanged)
    'gray_1': '#FFFFFF', 'gray_2': '#F0F0F0', 'gray_3': '#D9D9D9',
    'gray_4': '#BDBDBD', 'gray_5': '#969696', 'gray_6': '#737373',
    'gray_7': '#525252', 'gray_8': '#252525', 'gray_9': '#000000',

    # Utility Colors (Can be adjusted or expanded)
    'teal': '#4DB6AC',
    'pink': '#F48FB1',
    'white': '#FFFFFF',
    'black': '#000000'
}

# Backwards compatibility & aliases (adjust as needed)
OPSIN_COLORS['blue_dark'] = OPSIN_COLORS['domain_blue_dark']
OPSIN_COLORS['blue_light'] = OPSIN_COLORS['domain_blue_light']
OPSIN_COLORS['yellow'] = OPSIN_COLORS['mol_func_yellow_medium']
OPSIN_COLORS['orange'] = OPSIN_COLORS['mol_func_orange_medium']
OPSIN_COLORS['red'] = OPSIN_COLORS['mol_func_red_dark']
OPSIN_COLORS['purple'] = OPSIN_COLORS['mol_func_purple_medium']
OPSIN_COLORS['green'] = OPSIN_COLORS['domain_green_medium']
OPSIN_COLORS['gray_light'] = OPSIN_COLORS['gray_3']
OPSIN_COLORS['gray_dark'] = OPSIN_COLORS['gray_7']


# =============================================================================
# PROPERTY COLOR MAPS (NEW DEFINITIONS)
# =============================================================================

DOMAIN_PROPERTY_COLORS = [
    OPSIN_COLORS['mol_func_yellow_light'],
    OPSIN_COLORS['mol_func_yellow_medium'],
    OPSIN_COLORS['mol_func_orange_light'],
    OPSIN_COLORS['mol_func_orange_medium'],
    OPSIN_COLORS['mol_func_red_medium'],
    OPSIN_COLORS['mol_func_red_dark'],
    OPSIN_COLORS['mol_func_purple_medium'],
    OPSIN_COLORS['mol_func_purple_dark'],
]
DOMAIN_PROPERTY_CMAP = LinearSegmentedColormap.from_list('domain_property_map', DOMAIN_PROPERTY_COLORS)

MOL_FUNC_PROPERTY_COLORS = [
    OPSIN_COLORS['domain_blue_dark'],
    OPSIN_COLORS['domain_blue_medium'],
    OPSIN_COLORS['domain_blue_light'],
    OPSIN_COLORS['domain_cyan_medium'],
    OPSIN_COLORS['domain_cyan_light'],
    OPSIN_COLORS['domain_green_medium'],
    OPSIN_COLORS['domain_green_light'],
]
MOL_FUNC_PROPERTY_CMAP = LinearSegmentedColormap.from_list('mol_func_property_map', MOL_FUNC_PROPERTY_COLORS)

# =============================================================================
# RMSD COLOR MAPS (USING SPECTRAL OR GRAYSCALE AS BEFORE)
# =============================================================================
# Spectral RMSD: (Cold -> Warm, e.g., Blue -> Yellow/Orange)
# This can remain independent of the new property colors, or you can adapt it.
# For now, keeping the previous spectral definition:
_rmsd_spectral_temp = [
    OPSIN_COLORS['domain_blue_dark'], # Low RMSD
    OPSIN_COLORS['domain_blue_light'],
    OPSIN_COLORS['white'],
    OPSIN_COLORS['mol_func_yellow_medium'],
    OPSIN_COLORS['mol_func_orange_medium'] # High RMSD
]
RMSD_CMAP = LinearSegmentedColormap.from_list('rmsd_spectral', _rmsd_spectral_temp)

RMSD_COMPACT_COLORS_SPECTRAL = [
    OPSIN_COLORS['domain_blue_dark'], # e.g., 0 - 1.5 Å
    OPSIN_COLORS['domain_blue_light'],# e.g., 1.5 - 3.0 Å
    OPSIN_COLORS['mol_func_orange_medium']   # e.g., > 3.0 Å
]
RMSD_BOUNDS = [0, 1.5, 3.0, 5.0]
RMSD_COMPACT_CMAP = ListedColormap(RMSD_COMPACT_COLORS_SPECTRAL)
RMSD_DISCRETE_CMAP = RMSD_COMPACT_CMAP

# Grayscale RMSD (Unchanged)
RMSD_GRAYSCALE_COLORS = [OPSIN_COLORS[f'gray_{i}'] for i in range(1,10) if f'gray_{i}' in OPSIN_COLORS]
RMSD_GRAYSCALE_CMAP = LinearSegmentedColormap.from_list('rmsd_grayscale', RMSD_GRAYSCALE_COLORS)
RMSD_GRAYSCALE_COMPACT_COLORS = [OPSIN_COLORS['gray_1'], OPSIN_COLORS['gray_5'], OPSIN_COLORS['gray_9']]
RMSD_GRAYSCALE_COMPACT_CMAP = ListedColormap(RMSD_GRAYSCALE_COMPACT_COLORS)

# =============================================================================
# DISTANCE COLOR MAPS (ADAPTED)
# =============================================================================
# Using reversed Molecular Function palette (Dark Purple for far, Yellow for close)
DISTANCE_CMAP_MOL_FUNC_REV = LinearSegmentedColormap.from_list(
    'distance_mol_func_rev', MOL_FUNC_PROPERTY_COLORS[::-1]
)
# Using reversed Domain palette (Green for far, Dark Blue for close)
DISTANCE_CMAP_DOMAIN_REV = LinearSegmentedColormap.from_list(
    'distance_domain_rev', DOMAIN_PROPERTY_COLORS[::-1]
)
# Default distance map - choose one, e.g., mol_func_rev
DISTANCE_CMAP = DISTANCE_CMAP_MOL_FUNC_REV

# =============================================================================
# DIVERGING COLOR MAP (ADAPTED)
# =============================================================================
DIVERGING_COLORS = [
    OPSIN_COLORS['domain_blue_dark'], # Strong negative
    OPSIN_COLORS['domain_blue_light'],# Weak negative
    OPSIN_COLORS['white'],            # Neutral
    OPSIN_COLORS['mol_func_orange_light'],# Weak positive
    OPSIN_COLORS['mol_func_red_dark']  # Strong positive
]
DIVERGING_CMAP = LinearSegmentedColormap.from_list('opsin_diverging_new', DIVERGING_COLORS)

# =============================================================================
# HELIX COLORS (NEW DEFINITION: Blue -> Cyan -> Green -> Yellow -> Orange -> Red -> Purple)
# =============================================================================
HELIX_COLORS = {
    1: OPSIN_COLORS['domain_blue_dark'],     # Start with Domain Palette
    2: OPSIN_COLORS['domain_cyan_medium'],
    3: OPSIN_COLORS['domain_green_medium'],
    4: OPSIN_COLORS['mol_func_yellow_medium'],# Transition to Mol Func Palette
    5: OPSIN_COLORS['mol_func_orange_medium'],
    6: OPSIN_COLORS['mol_func_red_dark'],
    7: OPSIN_COLORS['mol_func_purple_dark'],  # End with Mol Func Palette
    'retinal': OPSIN_COLORS['pink'] # A distinct color for retinal
}
HELIX_COLORS_STR = {str(k): v for k, v in HELIX_COLORS.items() if isinstance(k, int)}
HELIX_COLORS_STR['retinal'] = HELIX_COLORS['retinal']
HELIX_COLORS_LIST = [HELIX_COLORS[i] for i in range(1, 8)]

def get_helix_cmap():
    return ListedColormap(HELIX_COLORS_LIST)

# =============================================================================
# CATEGORICAL COLORS (GROUPS, DOMAINS - PREDEFINED)
# =============================================================================
# For Molecular Function (uses new Yellow -> Purple palette)
GROUP_COLORS_PREDEFINED = {
    # Key groups mapped to distinct points in the new Mol Func palette
    'Rhodopsin': OPSIN_COLORS['mol_func_red_dark'],
    'Cone_Opsin': OPSIN_COLORS['mol_func_orange_medium'],
    'Visual_Opsin': OPSIN_COLORS['mol_func_yellow_medium'],
    'Non_visual_Opsin': OPSIN_COLORS['mol_func_yellow_light'],
    'Melanopsin': OPSIN_COLORS['mol_func_purple_medium'], # Distinct end
    'Photoisomerase': OPSIN_COLORS['domain_cyan_medium'], # Cross-palette for distinction
    'Cnidopsin': OPSIN_COLORS['domain_green_medium'],   # Cross-palette for distinction
    'Bistable': OPSIN_COLORS['teal'],                  # Utility color
    'Other': OPSIN_COLORS['gray_5'],
    'Unknown': OPSIN_COLORS['gray_3']
}

# For Domains (uses new Blue -> Green palette)
DOMAIN_COLORS_PREDEFINED = {
    'Bacteria': OPSIN_COLORS['domain_blue_dark'],
    'Eukaryota': OPSIN_COLORS['domain_blue_medium'],
    'Archaea': OPSIN_COLORS['domain_cyan_medium'],
    'Virus': OPSIN_COLORS['domain_green_light'],
    'Synthetic': OPSIN_COLORS['mol_func_orange_light'], # Cross-palette for distinction
    'Unknown': OPSIN_COLORS['gray_3']
}

# Function to dynamically assign colors (Updated logic for fallback)
def get_group_colors(group_names, palette_type='default', custom_predefined=None):
    """
    Dynamically assigns colors to a list of group names.
    Args:
        palette_type (str): 'mol_func', 'domain', 'gray', 'helix', or 'default' (acts like 'mol_func').
    """
    if isinstance(group_names, dict):
        unique_groups = sorted(list(set(group_names.keys())))
    else:
        unique_groups = sorted(list(set(group_names)))

    assigned_colors = {}

    # Determine the primary set of predefined colors and the fallback palette
    if custom_predefined is not None:
        predefined = custom_predefined
        # Fallback logic for custom_predefined needs to be context-aware or generic
        if palette_type == 'mol_func' or palette_type == 'default':
            fallback_palette_colors = MOL_FUNC_PROPERTY_COLORS
        elif palette_type == 'domain':
            fallback_palette_colors = DOMAIN_PROPERTY_COLORS
        elif palette_type == 'gray':
            fallback_palette_colors = RMSD_GRAYSCALE_COLORS
        elif palette_type == 'helix':
            fallback_palette_colors = HELIX_COLORS_LIST + [OPSIN_COLORS['pink']] # Added retinal color
        else: # Generic fallback for custom_predefined
            fallback_palette_colors = MOL_FUNC_PROPERTY_COLORS + DOMAIN_PROPERTY_COLORS # Mix
    elif palette_type == 'mol_func' or palette_type == 'default':
        predefined = GROUP_COLORS_PREDEFINED
        fallback_palette_colors = MOL_FUNC_PROPERTY_COLORS
    elif palette_type == 'domain':
        predefined = DOMAIN_COLORS_PREDEFINED
        fallback_palette_colors = DOMAIN_PROPERTY_COLORS
    elif palette_type == 'gray':
        predefined = {}
        fallback_palette_colors = RMSD_GRAYSCALE_COLORS
    elif palette_type == 'helix':
        predefined = {} # Usually HELIX_COLORS is used directly for the 7TMs
        fallback_palette_colors = HELIX_COLORS_LIST + [OPSIN_COLORS['pink']]
    else: # Unknown palette_type string
        predefined = {}
        fallback_palette_colors = MOL_FUNC_PROPERTY_COLORS + DOMAIN_PROPERTY_COLORS # Mix

    # Assign predefined colors
    for group in unique_groups:
        if group in predefined:
            assigned_colors[group] = predefined[group]

    # Assign from fallback_palette for remaining groups
    remaining_groups = [g for g in unique_groups if g not in assigned_colors]
    if not remaining_groups:
        return assigned_colors

    num_remaining = len(remaining_groups)
    colors_to_assign = []
    if num_remaining <= len(fallback_palette_colors):
        colors_to_assign = fallback_palette_colors[:num_remaining]
    else:
        for i in range(num_remaining):
            colors_to_assign.append(fallback_palette_colors[i % len(fallback_palette_colors)])

    for i, group in enumerate(remaining_groups):
        assigned_colors[group] = colors_to_assign[i]
    return assigned_colors


# =============================================================================
# AMINO ACID COLORS FOR LOGO PLOTS (ADAPTED)
# =============================================================================
AMINO_ACID_COLORS_COMBINED = {
    # Hydrophobic (Mol Func - Yellow/Orange)
    'A': OPSIN_COLORS['mol_func_yellow_light'], 'V': OPSIN_COLORS['mol_func_yellow_medium'],
    'I': OPSIN_COLORS['mol_func_orange_light'], 'L': OPSIN_COLORS['mol_func_orange_medium'],
    'M': OPSIN_COLORS['mol_func_red_medium'], 'F': OPSIN_COLORS['mol_func_red_dark'],
    'Y': OPSIN_COLORS['mol_func_purple_medium'], 'W': OPSIN_COLORS['mol_func_purple_dark'],
    # Polar Neutral (Domain - Blue/Cyan/Green)
    'S': OPSIN_COLORS['domain_green_light'], 'T': OPSIN_COLORS['domain_green_medium'],
    'N': OPSIN_COLORS['domain_cyan_light'], 'Q': OPSIN_COLORS['domain_cyan_medium'],
    # Special (Neutrals or distinct utility)
    'G': OPSIN_COLORS['gray_3'], 'P': OPSIN_COLORS['teal'],
    'C': OPSIN_COLORS['pink'], # Was purple_medium, now pink for more distinction
    # Charged Basic (+) (Domain - Blues)
    'K': OPSIN_COLORS['domain_blue_light'], 'R': OPSIN_COLORS['domain_blue_medium'],
    'H': OPSIN_COLORS['domain_blue_dark'],
    # Charged Acidic (-) (Mol Func - Reds, distinct from hydrophobic oranges)
    'D': OPSIN_COLORS['mol_func_red_dark'], 'E': OPSIN_COLORS['mol_func_red_medium'], # Using more distinct reds
    '-': OPSIN_COLORS['gray_1']
}
_standard_aas = "ACDEFGHIKLMNPQRSTVWY-"
for aa in _standard_aas:
    if aa not in AMINO_ACID_COLORS_COMBINED:
        AMINO_ACID_COLORS_COMBINED[aa] = OPSIN_COLORS['gray_5']

# =============================================================================
# UTILITIES & REGISTRATION
# =============================================================================
def register_colormaps():
    """Register custom colormaps with matplotlib."""
    cmap_dict = {
        'opsin_domain_prop': DOMAIN_PROPERTY_CMAP,
        'opsin_mol_func_prop': MOL_FUNC_PROPERTY_CMAP,
        'opsin_rmsd_spectral': RMSD_CMAP,
        'opsin_rmsd_spectral_compact': RMSD_COMPACT_CMAP,
        'opsin_rmsd_grayscale': RMSD_GRAYSCALE_CMAP,
        'opsin_distance_default': DISTANCE_CMAP,
        'opsin_distance_mol_func_rev': DISTANCE_CMAP_MOL_FUNC_REV,
        'opsin_distance_domain_rev': DISTANCE_CMAP_DOMAIN_REV,
        'opsin_diverging': DIVERGING_CMAP,
        'opsin_helices': get_helix_cmap()
    }
    for name, cmap in cmap_dict.items():
        try: plt.colormaps.register(name=name, cmap=cmap)
        except AttributeError: # Older matplotlib
            from matplotlib.cm import register_cmap
            register_cmap(name=name, cmap=cmap)
        except Exception as e: print(f"[WARNING] Could not register colormap '{name}': {e}")

register_colormaps()

def apply_style(): # (Keep your existing apply_style function)
    plt.style.use('default')
    style_params = {
        'font.family': 'sans-serif', 'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
        'font.size': 10, 'axes.labelsize': 12, 'axes.titlesize': 14,
        'xtick.labelsize': 10, 'ytick.labelsize': 10, 'legend.fontsize': 10,
        'figure.titlesize': 16, 'axes.spines.top': False, 'axes.spines.right': False,
        'axes.grid': True, 'grid.alpha': 0.3, 'grid.linestyle': '--',
        'lines.linewidth': 1.5, 'lines.markersize': 6,
        'axes.labelpad': 8, 'axes.titlepad': 12, 'legend.frameon': False,
        'savefig.dpi': 300, 'savefig.bbox': 'tight', 'figure.facecolor': 'white'
    }
    try:
        plt.rcParams.update(style_params)
        sns.set_style("ticks", rc=style_params)
    except Exception as e:
        print(f"[WARNING] Could not fully apply style settings: {e}")
    return plt.rcParams.copy()

apply_style()

# =============================================================================
# EXAMPLE USAGE
# =============================================================================
if __name__ == "__main__":
    fig, axs = plt.subplots(8, 1, figsize=(10, 16)) # Increased to 8 for new maps
    gradient = np.linspace(0, 1, 256).reshape(1, -1)

    axs[0].imshow(gradient, aspect='auto', cmap=DOMAIN_PROPERTY_CMAP)
    axs[0].set_title('Domain Property (Blue -> Cyan -> Green)')
    axs[1].imshow(gradient, aspect='auto', cmap=MOL_FUNC_PROPERTY_CMAP)
    axs[1].set_title('Molecular Function Property (Yellow -> Red -> Purple)')
    axs[2].imshow(gradient, aspect='auto', cmap=RMSD_CMAP)
    axs[2].set_title('RMSD (Spectral)')
    axs[3].imshow(gradient, aspect='auto', cmap=DISTANCE_CMAP)
    axs[3].set_title('Distance (Default - Mol Func Reversed)')
    axs[4].imshow(gradient, aspect='auto', cmap=DIVERGING_CMAP)
    axs[4].set_title('Diverging Colormap')

    helix_indices = np.array([[i] for i in range(7)]).reshape(1, -1)
    axs[5].imshow(helix_indices, aspect='auto', cmap=get_helix_cmap())
    axs[5].set_title('Helix Colors (New Combined Gradient)'); axs[5].set_xticks(range(7)); axs[5].set_xticklabels([f'H{i+1}' for i in range(7)])

    aa_labels = list(AMINO_ACID_COLORS_COMBINED.keys())
    aa_colors_list = [AMINO_ACID_COLORS_COMBINED[aa] for aa in aa_labels]
    aa_indices = np.array([[i] for i in range(len(aa_labels))]).reshape(1, -1)
    axs[6].imshow(aa_indices, aspect='auto', cmap=ListedColormap(aa_colors_list))
    axs[6].set_title('Amino Acid Combined Colors (New Palettes)')
    axs[6].set_xticks(range(len(aa_labels))); axs[6].set_xticklabels(aa_labels, rotation=45, ha="right")

    # Demo the grayscale RMSD
    axs[7].imshow(gradient, aspect='auto', cmap=RMSD_GRAYSCALE_CMAP)
    axs[7].set_title('RMSD (Grayscale)')


    for ax in axs: ax.set_yticks([])
    plt.tight_layout()
    plt.show()

    # Example of group colors
    example_mol_funcs = ['Rhodopsin', 'Cone_Opsin', 'NewFuncA', 'NewFuncB', 'Unknown']
    example_domains = ['Bacteria', 'Eukaryota', 'NewDomainX', 'Unknown']

    plt.figure(figsize=(12, 6))
    plt.subplot(1, 2, 1)
    colors1 = get_group_colors(example_mol_funcs, palette_type='mol_func')
    for i, (group, color) in enumerate(colors1.items()): plt.bar(i, 1, color=color, label=group)
    plt.title('Molecular Function Group Colors'); plt.legend(); plt.xticks([])

    plt.subplot(1, 2, 2)
    colors2 = get_group_colors(example_domains, palette_type='domain')
    for i, (group, color) in enumerate(colors2.items()): plt.bar(i, 1, color=color, label=group)
    plt.title('Domain Group Colors'); plt.legend(); plt.xticks([])
    plt.tight_layout()
    plt.show()