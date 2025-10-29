"""
Opsin structure visualization color scheme.
This module defines consistent color palettes and color assignment logic
for all opsin structure visualizations.
"""

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np # Used in example usage
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

# =============================================================================
# I. CORE HEX CODES (Primary Palette Definitions)
# =============================================================================

OPSIN_COLORS = {
    # A. WARM Palette (e.g., for Molecular Function, Hydrophobic AAs)
    # Yellow -> Orange -> Red -> Dark Purple
    'warm_yellow_lightest': '#FFFFD4',
    'warm_yellow_light': '#FFFFB2',
    'warm_yellow_medium': '#FED976',
    'warm_orange_light': '#FEB24C',
    'warm_orange_medium': '#FD8D3C',
    'warm_red_medium': '#FC4E2A',
    'warm_red_dark': '#E31A1C',
    'warm_purple_medium': '#BD0026',
    'warm_purple_dark': '#800026',
    'warm_magenta_dark': '#5C001F',

    # B. COLD Palette (e.g., for Domain, Polar/Charged AAs)
    # Dark Blue -> Medium Blue -> Light Blue -> Cyan -> Light Green -> Medium Green
    'cold_blue_darkest': '#08306B',
    'cold_blue_dark': '#08519C',
    'cold_blue_medium': '#2171B5',
    'cold_blue_light': '#4292C6',
    'cold_cyan_light': '#7FCDBB',
    'cold_cyan_medium': '#41B6C4',
    'cold_green_light': '#A1D99B',
    'cold_green_medium': '#4CAF50',
    'cold_green_dark': '#228B22',

    # C. GRAYSCALE Palette (White to Black)
    'gray_1_white': '#FFFFFF',
    'gray_2_lightest': '#F0F0F0',
    'gray_3_light': '#D9D9D9',
    'gray_4_light_mid': '#BDBDBD',
    'gray_5_mid': '#969696',
    'gray_6_dark_mid': '#737373',
    'gray_7_dark': '#525252',
    'gray_8_darkest': '#252525',
    'gray_9_black': '#000000',

    # D. UTILITY & ACCENT Colors (Distinct from Warm/Cold/Gray)
    'utility_teal': '#4DB6AC',
    'utility_pink': '#F48FB1',
    'utility_lime': '#AFE1AF',
    'utility_lavender': '#E6E6FA',
    'white': '#FFFFFF', # Alias for gray_1_white
    'black': '#000000', # Alias for gray_9_black
    # Muted comparison palette entries
    'muted_blue_gray': '#89A8B2',
    'muted_coral': '#E4A5A5',
}
OPSIN_COLORS['white'] = OPSIN_COLORS['gray_1_white'] # Explicit alias
OPSIN_COLORS['black'] = OPSIN_COLORS['gray_9_black'] # Explicit alias


# =============================================================================
# II. SEQUENTIAL COLOR LISTS (Derived from Core Hex Codes)
# =============================================================================

WARM_SEQUENTIAL_COLORS = [
    OPSIN_COLORS['warm_yellow_lightest'], OPSIN_COLORS['warm_yellow_medium'],
    OPSIN_COLORS['warm_orange_light'], OPSIN_COLORS['warm_orange_medium'],
    OPSIN_COLORS['warm_red_medium'], OPSIN_COLORS['warm_red_dark'],
    OPSIN_COLORS['warm_purple_medium'], OPSIN_COLORS['warm_purple_dark']
]

COLD_SEQUENTIAL_COLORS = [
    OPSIN_COLORS['cold_blue_darkest'], OPSIN_COLORS['cold_blue_dark'],
    OPSIN_COLORS['cold_blue_medium'], OPSIN_COLORS['cold_blue_light'],
    OPSIN_COLORS['cold_cyan_medium'], OPSIN_COLORS['cold_cyan_light'],
    OPSIN_COLORS['cold_green_medium'], OPSIN_COLORS['cold_green_dark']
]

GRAYSCALE_SEQUENTIAL_COLORS = [
    OPSIN_COLORS['gray_1_white'], OPSIN_COLORS['gray_2_lightest'], OPSIN_COLORS['gray_3_light'],
    OPSIN_COLORS['gray_4_light_mid'], OPSIN_COLORS['gray_5_mid'], OPSIN_COLORS['gray_6_dark_mid'],
    OPSIN_COLORS['gray_7_dark'], OPSIN_COLORS['gray_8_darkest'], OPSIN_COLORS['gray_9_black']
]

# =============================================================================
# III. CATEGORICAL COLOR DEFINITIONS (Predefined Mappings)
# =============================================================================

# A. For Property 1 (Molecular Function - uses WARM palette)
# Property 1 values: ['Sensor / Regulatory', 'Proton Pump', 'Anion Channel',
#                     'Chloride Pump', 'Unknown', 'Cation Channel', 'Sodium Pump']
PROPERTY1_COLORS_PREDEFINED = {
    'Sensor / Regulatory': OPSIN_COLORS['warm_yellow_medium'],
    'Proton Pump':         OPSIN_COLORS['warm_orange_medium'],
    'Anion Channel':       OPSIN_COLORS['warm_red_medium'],
    'Chloride Pump':       OPSIN_COLORS['warm_red_dark'],      # Darker shade for related pump
    'Cation Channel':      OPSIN_COLORS['warm_purple_medium'],
    'Sodium Pump':         OPSIN_COLORS['warm_purple_dark'],   # Darker shade for related pump
    'Unknown':             OPSIN_COLORS['gray_3_light'],     # Specific key for 'Unknown'
    'Other_Property1':     OPSIN_COLORS['gray_5_mid']        # For a generic 'Other' if it appears
}
# Ensure 'Unknown_Property1' can be used as an alias if needed by some parts of code
PROPERTY1_COLORS_PREDEFINED['Unknown_Property1'] = PROPERTY1_COLORS_PREDEFINED['Unknown']


# B. For Property 2 (Domain - uses COLD palette)
# Property 2 values: ['Eukaryota', 'Bacteria', 'Synthetic', 'Archaea', 'Unknown', 'Virus']
PROPERTY2_COLORS_PREDEFINED = {
    'Eukaryota':           OPSIN_COLORS['cold_blue_medium'],
    'Bacteria':            OPSIN_COLORS['cold_blue_darkest'], # Often many, give it a strong base color
    'Archaea':             OPSIN_COLORS['cold_cyan_medium'],
    'Virus':               OPSIN_COLORS['cold_green_light'],
    'Synthetic':           OPSIN_COLORS['cold_green_dark'],   # Distinct cold color
    'Unknown':             OPSIN_COLORS['gray_3_light'],      # Specific key for 'Unknown'
    'Other_Property2':     OPSIN_COLORS['gray_5_mid']         # For a generic 'Other' if it appears
}
# Ensure 'Unknown_Property2' can be used as an alias
PROPERTY2_COLORS_PREDEFINED['Unknown_Property2'] = PROPERTY2_COLORS_PREDEFINED['Unknown']


# C. For Helix Colors (combines COLD to WARM spectrum)
HELIX_NUMBER_COLORS = {
    1: OPSIN_COLORS['cold_blue_darkest'],   # Helix 1 (Cold start)
    2: OPSIN_COLORS['cold_blue_medium'],
    3: OPSIN_COLORS['cold_cyan_medium'],
    4: OPSIN_COLORS['warm_yellow_medium'],  # Mid-point transition to Warm
    5: OPSIN_COLORS['warm_orange_medium'],
    6: OPSIN_COLORS['warm_red_dark'],
    7: OPSIN_COLORS['warm_purple_dark'],    # Helix 7 (Warm end)
    'retinal': OPSIN_COLORS['utility_pink'] # Special ligand color
}
# String version for convenience if helix numbers are strings in data
HELIX_STRING_COLORS = {str(k): v for k, v in HELIX_NUMBER_COLORS.items()}

# D. For Amino Acid Logo (categorized by property, using WARM, COLD, and UTILITY)
AMINO_ACID_LOGO_COLORS = {
    # Hydrophobic (WARM palette primarily)
    'A': OPSIN_COLORS['warm_yellow_lightest'], 'V': OPSIN_COLORS['warm_yellow_medium'],
    'I': OPSIN_COLORS['warm_orange_light'], 'L': OPSIN_COLORS['warm_orange_medium'],
    'M': OPSIN_COLORS['warm_red_medium'],
    'F': OPSIN_COLORS['warm_red_dark'], 'Y': OPSIN_COLORS['warm_purple_medium'], 'W': OPSIN_COLORS['warm_purple_dark'],
    # Polar Neutral (COLD palette primarily)
    'S': OPSIN_COLORS['cold_green_light'], 'T': OPSIN_COLORS['cold_green_medium'],
    'N': OPSIN_COLORS['cold_cyan_light'], 'Q': OPSIN_COLORS['cold_cyan_medium'],
    # Special (UTILITY or distinct GRAYSCALE)
    'G': OPSIN_COLORS['gray_4_light_mid'], 'P': OPSIN_COLORS['utility_teal'],
    'C': OPSIN_COLORS['utility_pink'],
    # Charged Basic (+) (COLD palette - blues)
    'K': OPSIN_COLORS['cold_blue_light'], 'R': OPSIN_COLORS['cold_blue_medium'], 'H': OPSIN_COLORS['cold_blue_dark'],
    # Charged Acidic (-) (WARM palette - distinct reds/magentas)
    'D': OPSIN_COLORS['warm_magenta_dark'], 'E': OPSIN_COLORS['warm_purple_dark'],
    '-': OPSIN_COLORS['gray_1_white'] # Gap color
}
# Ensure all standard AAs have a color (fallback to a mid-gray)
_STANDARD_AAS = "ACDEFGHIKLMNPQRSTVWY-"
for aa_code in _STANDARD_AAS:
    if aa_code not in AMINO_ACID_LOGO_COLORS:
        AMINO_ACID_LOGO_COLORS[aa_code] = OPSIN_COLORS['gray_5_mid']

# E. Specific Status Colors
STATUS_EXPERIMENTAL_COLOR = OPSIN_COLORS['black']
STATUS_PREDICTED_COLOR = OPSIN_COLORS['utility_lime']


# =============================================================================
# IV. COLORMAPS (for continuous data or specific visualizations)
# =============================================================================

# A. Property Colormaps (derived from sequential lists)
WARM_PROPERTY_CMAP = LinearSegmentedColormap.from_list('opsin_warm_property_map', WARM_SEQUENTIAL_COLORS)
COLD_PROPERTY_CMAP = LinearSegmentedColormap.from_list('opsin_cold_property_map', COLD_SEQUENTIAL_COLORS)

# B. RMSD Colormaps
RMSD_GRAYSCALE_CMAP = LinearSegmentedColormap.from_list('opsin_rmsd_grayscale', GRAYSCALE_SEQUENTIAL_COLORS)
RMSD_GRAYSCALE_COMPACT_COLORS = [OPSIN_COLORS['gray_1_white'], OPSIN_COLORS['gray_5_mid'], OPSIN_COLORS['gray_9_black']]
RMSD_GRAYSCALE_COMPACT_CMAP = ListedColormap(RMSD_GRAYSCALE_COMPACT_COLORS)
_rmsd_spectral_colors = [
    OPSIN_COLORS['cold_blue_darkest'], OPSIN_COLORS['cold_blue_light'], OPSIN_COLORS['white'],
    OPSIN_COLORS['warm_yellow_medium'], OPSIN_COLORS['warm_orange_medium']
]

RMSD_CUSTOM_GRAY_COLORS = [
    OPSIN_COLORS['gray_2_lightest'], # Lightest gray (almost white)
    OPSIN_COLORS['gray_5_mid'],      # Medium gray
    OPSIN_COLORS['gray_8_darkest']      # Dark gray (not black)
]
RMSD_CUSTOM_GRAY_CMAP = ListedColormap(RMSD_CUSTOM_GRAY_COLORS)

RMSD_SPECTRAL_CMAP = LinearSegmentedColormap.from_list('opsin_rmsd_spectral', _rmsd_spectral_colors)
RMSD_SPECTRAL_COMPACT_CMAP = ListedColormap([
    OPSIN_COLORS['cold_blue_darkest'], OPSIN_COLORS['cold_blue_light'], OPSIN_COLORS['warm_orange_medium']
])

RMSD_WHITE_TO_DARKGRAY_COLORS = [
    OPSIN_COLORS['gray_1_white'],      # For min RMSD value
    OPSIN_COLORS['gray_8_darkest']       # For max RMSD value
]   
RMSD_WHITE_TO_DARKGRAY_CMAP = LinearSegmentedColormap.from_list(
    'opsin_rmsd_white_to_darkgray',
    RMSD_WHITE_TO_DARKGRAY_COLORS
)

RMSD_BOUNDS = [0.0, 1.0, 2.5, 5.0]

# C. Distance Colormaps
DISTANCE_WARM_REVERSED_CMAP = LinearSegmentedColormap.from_list(
    'opsin_distance_warm_rev', WARM_SEQUENTIAL_COLORS[::-1]
)
DISTANCE_COLD_REVERSED_CMAP = LinearSegmentedColormap.from_list(
    'opsin_distance_cold_rev', COLD_SEQUENTIAL_COLORS[::-1]
)
DEFAULT_DISTANCE_CMAP = DISTANCE_WARM_REVERSED_CMAP

# D. Diverging Colormap
DIVERGING_COLORS_LIST = [
    OPSIN_COLORS['cold_blue_darkest'], OPSIN_COLORS['cold_blue_light'], OPSIN_COLORS['white'],
    OPSIN_COLORS['warm_orange_light'], OPSIN_COLORS['warm_red_dark']
]
DIVERGING_CMAP = LinearSegmentedColormap.from_list('opsin_diverging', DIVERGING_COLORS_LIST)

# E. Helix Colormap
HELIX_CMAP_LIST = [HELIX_NUMBER_COLORS[i] for i in range(1, 8)]
HELIX_CMAP = ListedColormap(HELIX_CMAP_LIST)


# =============================================================================
# V. DYNAMIC COLOR ASSIGNMENT FUNCTION (for categorical data)
# =============================================================================

def get_categorical_colors(item_names, property_type='property1', custom_predefined=None):
    """
    Assigns colors to a list of item names based on a specified property type.
    Uses predefined colors first, then falls back to a sequential palette.

    Args:
        item_names (list or set): Unique names of items to color.
        property_type (str):
            'property1' (e.g., Molecular Function - uses WARM palette).
            'property2' (e.g., Domain - uses COLD palette).
            'grayscale' (uses GRAYSCALE palette).
            'helix' (uses HELIX_STRING_COLORS directly if item_names are helix numbers/strings).
            'custom' (requires `custom_predefined` to be set).
        custom_predefined (dict, optional): A custom dictionary of item_name:color pairs.

    Returns:
        dict: A dictionary mapping item_name to its assigned hex color string.
    """
    if isinstance(item_names, dict):
        unique_items = sorted(list(set(item_names.keys())))
    else:
        unique_items = sorted(list(set(str(item) for item in item_names))) # Ensure strings for dict keys

    assigned_colors = {}

    if property_type == 'property1':
        predefined_map = PROPERTY1_COLORS_PREDEFINED
        fallback_palette = WARM_SEQUENTIAL_COLORS
    elif property_type == 'property2':
        predefined_map = PROPERTY2_COLORS_PREDEFINED
        fallback_palette = COLD_SEQUENTIAL_COLORS
    elif property_type == 'grayscale':
        predefined_map = {}
        fallback_palette = GRAYSCALE_SEQUENTIAL_COLORS
    elif property_type == 'helix':
        predefined_map = HELIX_STRING_COLORS
        fallback_palette = HELIX_CMAP_LIST
    elif property_type == 'custom' and custom_predefined is not None:
        predefined_map = custom_predefined
        fallback_palette = WARM_SEQUENTIAL_COLORS + COLD_SEQUENTIAL_COLORS
    else:
        print(f"Warning: Unknown property_type '{property_type}'. Using Property1/Warm palette as default.")
        predefined_map = PROPERTY1_COLORS_PREDEFINED
        fallback_palette = WARM_SEQUENTIAL_COLORS

    for item in unique_items:
        if item in predefined_map:
            assigned_colors[item] = predefined_map[item]
        elif str(item) in predefined_map: # Check for string version too
             assigned_colors[item] = predefined_map[str(item)]


    remaining_items = [item for item in unique_items if item not in assigned_colors]
    if not remaining_items:
        return assigned_colors

    num_remaining = len(remaining_items)
    colors_to_assign_from_fallback = []
    if not fallback_palette: # Handle empty fallback_palette case
        print(f"Warning: Fallback palette for property_type '{property_type}' is empty. Remaining items will get gray.")
        for _ in range(num_remaining):
            colors_to_assign_from_fallback.append(OPSIN_COLORS['gray_5_mid'])
    elif num_remaining <= len(fallback_palette):
        colors_to_assign_from_fallback = fallback_palette[:num_remaining]
    else:
        for i in range(num_remaining):
            colors_to_assign_from_fallback.append(fallback_palette[i % len(fallback_palette)])

    for i, item in enumerate(remaining_items):
        assigned_colors[item] = colors_to_assign_from_fallback[i]

    return assigned_colors


# =============================================================================
# VI. MATPLOTLIB STYLING AND REGISTRATION
# =============================================================================

def register_opsin_colormaps():
    """Registers all custom colormaps with Matplotlib."""
    colormaps_to_register = {
        'opsin_warm_property': WARM_PROPERTY_CMAP,
        'opsin_cold_property': COLD_PROPERTY_CMAP,
        'opsin_rmsd_gray': RMSD_GRAYSCALE_CMAP,
        'opsin_rmsd_custom_gray': RMSD_CUSTOM_GRAY_CMAP,
        'opsin_rmsd_gray_compact': RMSD_GRAYSCALE_COMPACT_CMAP,
        'opsin_rmsd_spectral': RMSD_SPECTRAL_CMAP,
        'opsin_rmsd_spectral_compact': RMSD_SPECTRAL_COMPACT_CMAP,
        'opsin_distance_warm_rev': DISTANCE_WARM_REVERSED_CMAP,
        'opsin_distance_cold_rev': DISTANCE_COLD_REVERSED_CMAP,
        'opsin_distance_default': DEFAULT_DISTANCE_CMAP,
        'opsin_diverging': DIVERGING_CMAP,
        'opsin_helices': HELIX_CMAP
    }
    for name, cmap in colormaps_to_register.items():
        try:
            plt.colormaps.register(name=name, cmap=cmap)
        except ValueError:
            pass
        except Exception as e:
            print(f"[opsin_color_scheme] Warning: Could not register colormap '{name}': {e}")
    plt.colormaps.register(name='opsin_rmsd_white_to_darkgray', cmap=RMSD_WHITE_TO_DARKGRAY_CMAP)

def apply_opsin_style():
    """Applies a consistent Matplotlib style for opsin visualizations."""
    plt.style.use('default')
    style_params = {
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans', 'Bitstream Vera Sans'],
        'font.size': 10, 'axes.labelsize': 11, 'axes.titlesize': 12,
        'xtick.labelsize': 9, 'ytick.labelsize': 9, 'legend.fontsize': 9,
        'legend.title_fontsize': 10, 'figure.titlesize': 14,
        'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True,
        'grid.alpha': 0.4, 'grid.linestyle': ':', 'lines.linewidth': 1.5,
        'lines.markersize': 5, 'axes.labelpad': 6, 'axes.titlepad': 10,
        'legend.frameon': False, 'savefig.dpi': 300, 'savefig.bbox': 'tight',
        'figure.facecolor': OPSIN_COLORS['white']
    }
    try:
        plt.rcParams.update(style_params)
    except Exception as e:
        print(f"[opsin_color_scheme] Warning: Could not fully apply style settings: {e}")
    return plt.rcParams.copy()

register_opsin_colormaps()
_ = apply_opsin_style()


# =============================================================================
# VII. EXAMPLE USAGE (for testing and demonstration)
# =============================================================================
if __name__ == "__main__":
    print("--- Opsin Color Scheme Demo ---")

    print("\nTesting get_categorical_colors:")
    prop1_actual_values = ['Sensor / Regulatory', 'Proton Pump', 'Anion Channel',
                           'Chloride Pump', 'Unknown', 'Cation Channel', 'Sodium Pump', 'NewFunctionX']
    prop1_item_colors = get_categorical_colors(prop1_actual_values, property_type='property1')
    print(f"Property 1 (Function) Colors for {prop1_actual_values}:")
    for item, color in prop1_item_colors.items(): print(f"  {item}: {color}")

    prop2_actual_values = ['Eukaryota', 'Bacteria', 'Synthetic', 'Archaea', 'Unknown', 'Virus', 'NewDomainY']
    prop2_item_colors = get_categorical_colors(prop2_actual_values, property_type='property2')
    print(f"\nProperty 2 (Domain) Colors for {prop2_actual_values}:")
    for item, color in prop2_item_colors.items(): print(f"  {item}: {color}")

    helix_items = ['1', '3', '7', 'retinal', 'loop12']
    helix_item_colors = get_categorical_colors(helix_items, property_type='helix')
    print("\nHelix Item Colors:", helix_item_colors)

    cmaps_to_show = {
        'Warm Property': WARM_PROPERTY_CMAP, 'Cold Property': COLD_PROPERTY_CMAP,
        'RMSD Grayscale': RMSD_GRAYSCALE_CMAP, 'RMSD Grayscale Compact': RMSD_GRAYSCALE_COMPACT_CMAP,
        'Distance (Warm Rev)': DISTANCE_WARM_REVERSED_CMAP, 'Diverging': DIVERGING_CMAP,
        'Helices (1-7)': HELIX_CMAP
    }
    fig, axs = plt.subplots(len(cmaps_to_show), 1, figsize=(8, 2 * len(cmaps_to_show)))
    if len(cmaps_to_show) == 1: axs = [axs]

    gradient = np.linspace(0, 1, 256).reshape(1, -1)
    for i, (title, cmap) in enumerate(cmaps_to_show.items()):
        axs[i].imshow(gradient, aspect='auto', cmap=cmap)
        axs[i].set_title(title); axs[i].set_yticks([])
    plt.tight_layout(); plt.suptitle("Opsin Colormap Overview", y=1.02, fontsize=16); plt.show()

    plt.figure(figsize=(12, 3))
    aa_labels_sorted = sorted(AMINO_ACID_LOGO_COLORS.keys())
    aa_colors_list_sorted = [AMINO_ACID_LOGO_COLORS[aa] for aa in aa_labels_sorted]
    bar_positions = np.arange(len(aa_labels_sorted))
    plt.bar(bar_positions, np.ones(len(aa_labels_sorted)), color=aa_colors_list_sorted, tick_label=aa_labels_sorted)
    plt.title('Amino Acid Logo Colors'); plt.xticks(rotation=45, ha="right"); plt.tight_layout(); plt.show()

    print("\n--- End of Demo ---")
