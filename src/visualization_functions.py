"""Updated visualization functions for opsin analysis that use the global color scheme."""

import os
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from collections import Counter
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec # For create_combined_distance_logo_plot

import pandas as pd
import numpy as np
from collections import Counter
import logomaker


import os
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform, cdist
import logomaker

# Import from our updated color scheme
from src.opsin_color_scheme import (
    OPSIN_COLORS,
    get_categorical_colors,
    # Specific categorical maps (can be used directly or via get_categorical_colors)
    # PROPERTY1_COLORS_PREDEFINED,
    # PROPERTY2_COLORS_PREDEFINED,
    HELIX_NUMBER_COLORS,
    HELIX_STRING_COLORS, # If helix IDs are often strings
    AMINO_ACID_LOGO_COLORS,
    STATUS_EXPERIMENTAL_COLOR,
    STATUS_PREDICTED_COLOR,
    # Colormaps (can be used by string name if registered, or direct object)
    RMSD_GRAYSCALE_COMPACT_CMAP,
    RMSD_GRAYSCALE_CMAP,
    DEFAULT_DISTANCE_CMAP, # This is DISTANCE_WARM_REVERSED_CMAP
    # Other cmaps if needed by name later:
    # 'opsin_rmsd_gray_compact', 'opsin_default_distance', etc.
    RMSD_BOUNDS # For BoundaryNorm with compact RMSD maps
)


import numpy as np
import pandas as pd
import json
from pathlib import Path
from scipy.cluster.hierarchy import fcluster
from itertools import combinations

def _pairwise_values_from_square(df: pd.DataFrame):
    m = df.values.astype(float)
    n = m.shape[0]
    # exclude self (diag) and take upper triangle
    triu = m[np.triu_indices(n, k=1)]
    return triu[~np.isnan(triu)]


def _cluster_masks(labels: np.ndarray):
    """
    Given an array of cluster labels (e.g. from fcluster),
    return a dict mapping cluster_id -> boolean mask.
    """
    labs = np.asarray(labels)
    return {c: (labs == c) for c in np.unique(labs)}

def compute_rmsd_metrics(
    rmsd_df: pd.DataFrame,
    linkage_matrix,
    thresholds=(2.0, 2.5, 3.0),
    n_clusters=2,
    outdir: Path = Path("metrics")
):
    outdir.mkdir(parents=True, exist_ok=True)
    names = rmsd_df.index.to_list()
    M = rmsd_df.values.astype(float)

    # -------- Global distribution --------
    all_pairs = _pairwise_values_from_square(rmsd_df)
    def pct_below(x):
        return {f"pct_pairs_lt_{t:.1f}A": float(np.mean(x < t) * 100.0) for t in thresholds}
    global_stats = {
        "n_structures": int(M.shape[0]),
        "n_pairs": int(len(all_pairs)),
        "median_RMSD": float(np.nanmedian(all_pairs)),
        "IQR_RMSD": [float(np.nanpercentile(all_pairs, 25)), float(np.nanpercentile(all_pairs, 75))],
        **pct_below(all_pairs),
    }

    # -------- Clusters (hierarchical cut) --------
    labels = fcluster(linkage_matrix, n_clusters, criterion="maxclust")
    cluster_masks = _cluster_masks(labels)
    cluster_stats = {}
    for cid, mask in cluster_masks.items():
        idx = np.where(mask)[0]
        sub = M[np.ix_(idx, idx)]
        sub_pairs = _pairwise_values_from_square(pd.DataFrame(sub))
        base_stats = {
            "size": int(mask.sum()),
            "median_intra": float(np.nanmedian(sub_pairs)) if len(sub_pairs) else np.nan,
            "IQR_intra": [float(np.nanpercentile(sub_pairs,25)) if len(sub_pairs) else np.nan,
                          float(np.nanpercentile(sub_pairs,75)) if len(sub_pairs) else np.nan],
        }
        cluster_stats[int(cid)] = base_stats | (pct_below(sub_pairs) if len(sub_pairs) else {})

    # Inter-cluster (all pairs where labels differ)
    inter_vals = []
    for (c1, m1), (c2, m2) in combinations(cluster_masks.items(), 2):
        sub = M[np.ix_(np.where(m1)[0], np.where(m2)[0])]
        inter_vals.append(sub.reshape(-1))
    inter_vals = np.concatenate(inter_vals) if len(inter_vals) else np.array([])
    inter_stats = {
        "median_inter": float(np.nanmedian(inter_vals)) if inter_vals.size else np.nan,
        "IQR_inter": [float(np.nanpercentile(inter_vals,25)) if inter_vals.size else np.nan,
                      float(np.nanpercentile(inter_vals,75)) if inter_vals.size else np.nan],
        **(pct_below(inter_vals) if inter_vals.size else {})
    }

    # -------- Reference-like & outliers --------
    # Reference-like score: mean RMSD to the OTHER cluster(s)
    ref_scores = {}
    for i in range(M.shape[0]):
        other = M[i, labels != labels[i]]
        other = other[~np.isnan(other)]
        ref_scores[names[i]] = float(np.nanmean(other)) if other.size else np.nan
    # Outlier score: mean RMSD to ALL others
    out_scores = {}
    for i in range(M.shape[0]):
        row = M[i, :]
        row = row[~np.isnan(row)]
        out_scores[names[i]] = float(np.nanmean(row)) if row.size else np.nan

    # rank
    ref_like = sorted(ref_scores.items(), key=lambda kv: kv[1])[:10]  # lower is better
    outliers = sorted(out_scores.items(), key=lambda kv: kv[1], reverse=True)[:10]  # higher is worse

    # -------- Save to disk --------
    (outdir / "rmsd_global_stats.json").write_text(json.dumps(global_stats, indent=2))
    (outdir / "rmsd_cluster_stats.json").write_text(json.dumps({
        "n_clusters": int(n_clusters),
        "clusters": cluster_stats,
        "inter": inter_stats
    }, indent=2))
    pd.DataFrame({"name": list(ref_scores.keys()), "ref_like_mean_crosscluster_rmsd": list(ref_scores.values())})\
      .sort_values("ref_like_mean_crosscluster_rmsd").to_csv(outdir / "reference_like_top.csv", index=False)
    pd.DataFrame({"name": list(out_scores.keys()), "outlier_mean_global_rmsd": list(out_scores.values())})\
      .sort_values("outlier_mean_global_rmsd", ascending=False).to_csv(outdir / "outliers_top.csv", index=False)

    # Compose short overlay text for the figure
    overlay_lines = [
        f"Global median RMSD = {global_stats['median_RMSD']:.2f} Å (IQR {global_stats['IQR_RMSD'][0]:.2f}–{global_stats['IQR_RMSD'][1]:.2f})",
    ] + [f"{k}: {v:.1f}%" for k, v in global_stats.items() if k.startswith("pct_pairs_lt_")]
    overlay_lines.append(f"Inter-cluster median = {inter_stats['median_inter']:.2f} Å")
    for cid, st in cluster_stats.items():
        overlay_lines.append(f"Cluster {cid} (n={st['size']}): median intra = {st['median_intra']:.2f} Å; "
                             f"IQR {st['IQR_intra'][0]:.2f}–{st['IQR_intra'][1]:.2f}")
    overlay_text = "\n".join(overlay_lines)

    return {
        "global": global_stats,
        "clusters": cluster_stats,
        "inter": inter_stats,
        "labels": labels.tolist(),
        "overlay_text": overlay_text,
        "ref_like_top10": ref_like,
        "outliers_top10": outliers,
    }


def _annotate_metrics_on_clustergrid(clustergrid, text: str, loc=(0.01, 0.99)):
    """
    clustergrid: seaborn ClusterGrid or object exposing .ax_heatmap
    loc: axes fraction coordinates (x,y) for the top-left of the text box
    """
    ax = getattr(clustergrid, "ax_heatmap", None) or clustergrid  # support Figure fallback
    ax.text(
        loc[0], loc[1], text,
        transform=ax.transAxes, ha="left", va="top",
        fontsize=9, family="monospace",
        bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="black", lw=0.8, alpha=0.85)
    )

def plot_rmsd_heatmap(rmsd_df, title="All-vs-All RMSD Heatmap"):
    """Visualize the RMSD matrix as a heatmap using the grayscale color scheme."""
    plt.figure(figsize=(10, 8))
    vmax = 3.0  # Max value for color scale for RMSD_GRAYSCALE_COMPACT_CMAP if using its default bounds
                # If using RMSD_GRAYSCALE_CMAP (continuous), this vmax is still relevant.
    sns.heatmap(rmsd_df, annot=True, fmt=".2f", cmap=RMSD_GRAYSCALE_COMPACT_CMAP,
                vmin=0, vmax=vmax) # Or use 'opsin_rmsd_gray_compact'
    plt.title(title)
    plt.xlabel("Structure")
    plt.ylabel("Structure")
    plt.tight_layout()
    return plt.gcf()


def plot_similarity_tree(rmsd_df, title="Structural Similarity Tree"):
    """Generate a dendrogram based on RMSD values with consistent colors."""
    matrix = rmsd_df.fillna(rmsd_df.max().max()).values # Fill NaNs for linkage
    condensed = squareform(matrix)
    Z = linkage(condensed, method='weighted')

    plt.figure(figsize=(12, 8))

    def link_color_func(k):
        return OPSIN_COLORS['gray_7_dark'] # Slightly darker gray for good visibility

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
    Create a bar plot showing conservation levels at each position in the MSA.
    Uses helix colors for highlighting TM regions.
    """
    residue_df = msa_df.map(lambda x: x[0] if isinstance(x, str) and len(x) > 0 and x != '-' else '-')
    conservation = {}
    for col in residue_df.columns:
        counts = residue_df[col].value_counts()
        total = (residue_df[col] != '-').sum()
        if total > 0:
            top_res = counts.index[0] if not counts.empty else '-'
            if top_res != '-': top_pct = (counts[top_res] / total) * 100
            else:
                top_res = counts.index[1] if len(counts) > 1 else '-'
                top_pct = (counts[top_res] / total) * 100 if top_res != '-' else 0
            conservation[col] = {'residue': top_res, 'percentage': top_pct}
        else:
            conservation[col] = {'residue': '-', 'percentage': 0}

    plot_data = []
    for pos, data in conservation.items():
        pos_str = str(pos)
        helix_label = 'Other'
        pos_num = 0.0
        if 'x' in pos_str:
            try:
                helix_part, num_part = pos_str.split('x')
                if helix_part.isdigit(): helix_label = helix_part
                if num_part.replace('.', '', 1).isdigit(): pos_num = float(num_part)
            except ValueError: pass
        elif '.' in pos_str:
            parts = pos_str.split('.')
            num_after_dot = float(parts[1]) if len(parts) > 1 and parts[1].replace('.', '', 1).isdigit() else 0
            if pos_str.startswith('L.'): helix_label, pos_num = 'L', num_after_dot
            elif pos_str.startswith('n.'): helix_label, pos_num = 'N', num_after_dot
            elif pos_str.startswith('c.'): helix_label, pos_num = 'C', num_after_dot
            elif len(parts) == 2 and parts[0].isdigit(): helix_label, pos_num = parts[0], num_after_dot
            elif len(parts) >= 2 and len(parts[0]) == 2 and parts[0].isdigit(): # Loop AB.CCC
                helix_label, pos_num = 'L', float('0.' + parts[1]) if parts[1].isdigit() else 0
        elif pos_str.replace('.', '', 1).isdigit():
            pos_num = float(pos_str) # Helix remains 'Other'
        plot_data.append({'Position': pos, 'Helix': helix_label, 'Position_Num': pos_num,
                          'Residue': data['residue'], 'Conservation': data['percentage']})

    conservation_df = pd.DataFrame(plot_data)
    if 'Helix' in conservation_df.columns and not conservation_df['Helix'].isna().all():
        def helix_sort_key(h_val):
            h_str = str(h_val)
            if h_str.isdigit(): return int(h_str)
            elif h_str == 'N': return -1
            elif h_str == 'L': return 8 # Loops after TM7
            elif h_str == 'C': return 9 # C-term after loops
            else: return 10 # Other last
        conservation_df['Helix_Sort'] = conservation_df['Helix'].apply(helix_sort_key)
        conservation_df = conservation_df.sort_values(['Helix_Sort', 'Position_Num'])
    else:
        conservation_df = conservation_df.sort_values('Position')

    fig, ax = plt.subplots(figsize=figsize)
    x_coords = np.arange(len(conservation_df)) # Renamed x to x_coords
    bars = ax.bar(x_coords, conservation_df['Conservation'], width=0.8)

    if helix_highlighting and 'Helix' in conservation_df.columns:
        for i, (_, row) in enumerate(conservation_df.iterrows()):
            helix_str = str(row['Helix'])
            bar_color = OPSIN_COLORS['gray_6_dark_mid'] # Default for 'Other'
            if helix_str.isdigit() and int(helix_str) in HELIX_NUMBER_COLORS:
                bar_color = HELIX_NUMBER_COLORS[int(helix_str)]
            elif helix_str == 'L': bar_color = OPSIN_COLORS['gray_4_light_mid'] # Loop color
            elif helix_str == 'N': bar_color = OPSIN_COLORS['cold_blue_light']   # N-term color
            elif helix_str == 'C': bar_color = OPSIN_COLORS['cold_green_medium'] # C-term color
            bars[i].set_color(bar_color)

    for i, (_, row) in enumerate(conservation_df.iterrows()):
        if row['Conservation'] > 5:
            ax.text(i, row['Conservation'] + 2, row['Residue'], ha='center', va='bottom', fontweight='bold')

    ax.set_ylabel('Conservation (%)', fontsize=12)
    ax.set_xlabel('Position', fontsize=12)
    ax.set_title('Residue Conservation by Position', fontsize=14)
    ax.set_xticks(x_coords)
    ax.set_xticklabels(conservation_df['Position'], rotation=90)
    ax.axhline(y=50, color=OPSIN_COLORS['gray_3_light'], linestyle='--', alpha=0.7)
    ax.set_ylim(0, 105)
    ax.grid(axis='y', alpha=0.3)

    if helix_highlighting:
        legend_elements = []
        for i in range(1, 8):
            if i in HELIX_NUMBER_COLORS:
                legend_elements.append(Patch(facecolor=HELIX_NUMBER_COLORS[i], label=f'Helix {i}'))
        legend_elements.extend([
            Patch(facecolor=OPSIN_COLORS['cold_blue_light'], label='N-terminal'),
            Patch(facecolor=OPSIN_COLORS['gray_4_light_mid'], label='Loop'),
            Patch(facecolor=OPSIN_COLORS['cold_green_medium'], label='C-terminal'),
            Patch(facecolor=OPSIN_COLORS['gray_6_dark_mid'], label='Other')
        ])
        ax.legend(handles=legend_elements, loc='upper right')
    plt.tight_layout()
    return fig


def visualize_rmsd_heatmap(rmsd_df, structure_ids, group_dict=None, domain_dict=None, name_dict=None,
                           annot=False, font_scale=1.0, group_by='molecular_function', error_threshold=3.0):
    """Visualize RMSD matrix with side color annotations using the global color scheme."""
    sns.set_context("notebook", font_scale=font_scale)
    fig = plt.figure(figsize=(16, 14), dpi=100)

    # Filter by error_threshold (assuming domain_dict might contain error info)
    if domain_dict is not None and any('average_error' in domain_dict.get(sid, {}) for sid in structure_ids):
        structure_ids = [
            sid for sid in structure_ids
            if 'average_error' not in domain_dict.get(sid, {}) or
               domain_dict.get(sid, {}).get('average_error', 0) <= error_threshold
        ]
        if not structure_ids: # All filtered out
            ax = plt.gca()
            ax.text(0.5, 0.5, "All structures filtered by error threshold.", ha='center', va='center')
            return fig
        rmsd_df = rmsd_df.loc[structure_ids, structure_ids] # Re-filter rmsd_df

    # Prepare domain_dict_local (Property2)
    domain_data_local = {} # Renamed to avoid conflict
    if domain_dict is not None:
        for sid in structure_ids:
            info = domain_dict.get(sid, "Unknown")
            domain_data_local[sid] = str(info.get('domain', info)) if isinstance(info, dict) else str(info)
    else:
        domain_data_local = {sid: "Unknown" for sid in structure_ids}

    # Prepare group_dict_local (Property1)
    group_data_local = {}
    if group_dict is not None:
        for sid in structure_ids: group_data_local[sid] = str(group_dict.get(sid, "Unknown"))
    else:
        group_data_local = {sid: "Unknown" for sid in structure_ids}


    # Sorting
    if group_dict is not None:
        sort_key = lambda sid_val: (group_data_local.get(sid_val, "Unknown"), domain_data_local.get(sid_val, "Unknown"))
        structure_ids_sorted = sorted(structure_ids, key=sort_key)
        rmsd_matrix_values = rmsd_df.loc[structure_ids_sorted, structure_ids_sorted].values
    else:
        structure_ids_sorted = structure_ids
        rmsd_matrix_values = rmsd_df.loc[structure_ids_sorted, structure_ids_sorted].values

    max_rmsd_val = np.nanmax(rmsd_matrix_values[np.isfinite(rmsd_matrix_values)]) if np.any(np.isfinite(rmsd_matrix_values)) else 3.0
    rmsd_matrix_values = np.nan_to_num(rmsd_matrix_values, nan=0.0, posinf=max_rmsd_val, neginf=0.0)

    # Use BoundaryNorm with RMSD_GRAYSCALE_COMPACT_CMAP
    # RMSD_BOUNDS should define the edges of the bins for RMSD_GRAYSCALE_COMPACT_CMAP
    norm = BoundaryNorm(RMSD_BOUNDS, RMSD_GRAYSCALE_COMPACT_CMAP.N)

    ax = plt.gca()
    im = ax.imshow(rmsd_matrix_values, cmap=RMSD_GRAYSCALE_COMPACT_CMAP, norm=norm) # Or 'opsin_rmsd_gray_compact'

    try:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("left", size="5%", pad=0.5)
        cbar = plt.colorbar(im, cax=cax, ticks=RMSD_BOUNDS) # Use RMSD_BOUNDS for ticks
        cbar.set_label('RMSD (Å)', rotation=90, fontsize=14, labelpad=15)
        cax.yaxis.set_ticks_position('left'); cax.yaxis.set_label_position('left')
    except Exception as e:
        print(f"Warning: Could not create colorbar with axes_grid1: {e}")
        # Fallback colorbar (might not align perfectly with BoundaryNorm)
        cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=RMSD_GRAYSCALE_COMPACT_CMAP),
                            ax=ax, orientation='vertical', pad=0.05, ticks=RMSD_BOUNDS)
        cbar.set_label('RMSD (Å)', rotation=90, fontsize=14, labelpad=15)

    ax.set_xticks([]); ax.set_yticks([])
    ax.set_title('Structural Similarity Matrix (RMSD)', pad=20, fontsize=18, fontweight='bold')

    if group_dict is not None:
        prop1_values = [group_data_local.get(sid, "Unknown") for sid in structure_ids_sorted]
        prop2_values = [domain_data_local.get(sid, "Unknown") for sid in structure_ids_sorted]

        unique_prop1 = []; prop1_sizes = []; current_prop1 = None; count = 0
        for p1 in prop1_values:
            if p1 != current_prop1:
                if current_prop1 is not None: unique_prop1.append(current_prop1); prop1_sizes.append(count)
                current_prop1 = p1; count = 1
            else: count += 1
        if current_prop1 is not None: unique_prop1.append(current_prop1); prop1_sizes.append(count)
        prop1_boundaries = np.cumsum(prop1_sizes)

        for boundary in prop1_boundaries[:-1]:
            ax.axhline(y=boundary - 0.5, color=OPSIN_COLORS['white'], linewidth=2.5)
            ax.axvline(x=boundary - 0.5, color=OPSIN_COLORS['white'], linewidth=2.5)

        prop1_color_map = get_categorical_colors(unique_prop1, property_type='property1')
        prop2_color_map = get_categorical_colors(list(set(prop2_values)), property_type='property2')

        pos_val = 0
        for p1, size in zip(unique_prop1, prop1_sizes): # Add Property1 color strip at bottom
            ax.add_patch(plt.Rectangle(
                (pos_val, -0.05 * rmsd_matrix_values.shape[0]), size, 0.02 * rmsd_matrix_values.shape[0],
                facecolor=prop1_color_map.get(p1, OPSIN_COLORS['gray_3_light']),
                edgecolor='none', alpha=0.8, transform=ax.transData, clip_on=False))
            pos_val += size

        # Add Property2 color strip at right (iterate through sorted structures for correct coloring)
        # This requires careful indexing if there are subgroups within Property1 for Property2
        # For simplicity, let's color based on the sorted order directly for the side strip:
        current_y = 0
        for sid in structure_ids_sorted: # Iterate through the sorted IDs to get their Prop2 color
            p2_val = domain_data_local.get(sid, "Unknown")
            ax.add_patch(plt.Rectangle(
                (rmsd_matrix_values.shape[1], current_y), 0.02 * rmsd_matrix_values.shape[1], 1, # Height is 1 unit (per structure)
                facecolor=prop2_color_map.get(p2_val, OPSIN_COLORS['gray_3_light']),
                edgecolor='none', alpha=0.8, transform=ax.transData, clip_on=False))
            current_y += 1


        # Create separate legend figure
        legend_elements_prop1 = [Patch(facecolor=color, label=p1) for p1, color in sorted(prop1_color_map.items())]
        legend_elements_prop2 = [Patch(facecolor=color, label=p2) for p2, color in sorted(prop2_color_map.items())]

        legend_fig = plt.figure(figsize=(4, max(len(legend_elements_prop1), len(legend_elements_prop2)) * 0.5 + 2 ), dpi=100)
        gs_legend = legend_fig.add_gridspec(2,1, hspace=0.4)
        ax_leg1 = legend_fig.add_subplot(gs_legend[0])
        ax_leg1.legend(handles=legend_elements_prop1, title="Property 1 (e.g., Function)", loc='center', frameon=False)
        ax_leg1.axis('off')
        ax_leg2 = legend_fig.add_subplot(gs_legend[1])
        ax_leg2.legend(handles=legend_elements_prop2, title="Property 2 (e.g., Domain)", loc='center', frameon=False)
        ax_leg2.axis('off')
        legend_path = 'rmsd_heatmap_legends.png'
        legend_fig.savefig(legend_path, bbox_inches='tight')
        plt.close(legend_fig)
        ax.text(0.98, 0.02, "Legends: rmsd_heatmap_legends.png", transform=ax.transAxes, ha='right', va='bottom',
                fontsize=8, bbox=dict(facecolor=OPSIN_COLORS['white'], alpha=0.7, pad=0.2))

    plt.tight_layout(rect=[0.05, 0.05, 0.95, 0.95]) # Adjust to make space if colorbars are moved
    return fig


def create_and_visualize_similarity_tree(rmsd_data, group_dict=None, domain_dict=None,
                                         name_dict=None, font_scale=1.0, leaf_font_size=10,
                                         linkage_matrix=None):
    """Create a similarity tree (dendrogram) with side color bars for properties."""
    sns.set_context("notebook", font_scale=font_scale)
    fig = plt.figure(figsize=(16, max(12, len(rmsd_data) * 0.3) ), dpi=100) # Adjust height dynamically
    plt.grid(False)

    structure_ids = rmsd_data.index.tolist()
    N = len(structure_ids)

    if linkage_matrix is None:
        clean_matrix = np.nan_to_num(rmsd_data.values, nan=0.0,
                                     posinf=np.nanmax(rmsd_data.values[np.isfinite(rmsd_data.values)]) if np.any(np.isfinite(rmsd_data.values)) else 1.0,
                                     neginf=0.0)
        np.fill_diagonal(clean_matrix, 0.0)
        condensed_matrix = squareform(clean_matrix, checks=False)
        Z_linkage = linkage(condensed_matrix, method='weighted') # Renamed Z
    else:
        Z_linkage = linkage_matrix

    ax = plt.gca()
    dendro = dendrogram(Z_linkage, ax=ax, orientation='right',
                        link_color_func=lambda k: OPSIN_COLORS['gray_7_dark'],
                        above_threshold_color=OPSIN_COLORS['gray_7_dark'],
                        no_labels=True) # Labels will be custom

    ax.tick_params(axis='x', labelsize=10)
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False); ax.spines['left'].set_visible(False)
    ax.tick_params(left=False) # Remove y-ticks

    # Prepare Property1 and Property2 data
    prop1_data = {sid: str(group_dict.get(sid, "Unknown")) for sid in structure_ids} if group_dict else {sid: "Unknown" for sid in structure_ids}
    prop2_data = {}
    if domain_dict:
        for sid in structure_ids:
            info = domain_dict.get(sid, "Unknown")
            prop2_data[sid] = str(info.get('domain', info)) if isinstance(info, dict) else str(info)
    else:
        prop2_data = {sid: "Unknown" for sid in structure_ids}

    unique_prop1 = sorted(list(set(prop1_data.values())))
    unique_prop2 = sorted(list(set(prop2_data.values())))
    prop1_color_map = get_categorical_colors(unique_prop1, property_type='property1')
    prop2_color_map = get_categorical_colors(unique_prop2, property_type='property2')

    leaf_indices_ordered = dendro['leaves'] # Renamed leaf_indices
    ordered_sids = [structure_ids[i] for i in leaf_indices_ordered]

    # Add custom labels and color bars
    y_coords = np.arange(N) * 10 + 5 # y-positions for labels/bars
    label_x_pos = ax.get_xlim()[0] - 0.02 * (ax.get_xlim()[1] - ax.get_xlim()[0]) # Position labels to the left
    bar_start_x = ax.get_xlim()[1] # Start bars at the right edge of dendrogram
    bar_width_prop1 = 0.05 * (ax.get_xlim()[1] - ax.get_xlim()[0])
    bar_width_prop2 = 0.05 * (ax.get_xlim()[1] - ax.get_xlim()[0])
    bar_gap = 0.01 * (ax.get_xlim()[1] - ax.get_xlim()[0])

    for i, sid_val in enumerate(ordered_sids): # Renamed sid
        p1_val = prop1_data.get(sid_val, "Unknown") # Renamed p1
        p2_val = prop2_data.get(sid_val, "Unknown") # Renamed p2
        display_name = name_dict.get(sid_val, sid_val) if name_dict else sid_val

        # Leaf label (structure name)
        ax.text(label_x_pos, y_coords[i], display_name, ha='right', va='center', fontsize=leaf_font_size)

        # Property 1 color bar
        ax.add_patch(plt.Rectangle((bar_start_x, y_coords[i] - 4), bar_width_prop1, 8,
                                   facecolor=prop1_color_map.get(p1_val, OPSIN_COLORS['gray_3_light']),
                                   edgecolor=OPSIN_COLORS['white'], clip_on=False))
        # Property 2 color bar
        ax.add_patch(plt.Rectangle((bar_start_x + bar_width_prop1 + bar_gap, y_coords[i] - 4), bar_width_prop2, 8,
                                   facecolor=prop2_color_map.get(p2_val, OPSIN_COLORS['gray_3_light']),
                                   edgecolor=OPSIN_COLORS['white'], clip_on=False))

    ax.set_ylim(-5, N * 10 + 5) # Adjust y-limits for labels/bars

    # Legends
    legend_handles_p1 = [Patch(facecolor=color, label=p1) for p1, color in sorted(prop1_color_map.items()) if p1 != "Unknown"]
    legend_handles_p2 = [Patch(facecolor=color, label=p2) for p2, color in sorted(prop2_color_map.items()) if p2 != "Unknown"]

    # Position legends to the right of the color bars
    legend_x_anchor = bar_start_x + bar_width_prop1 + bar_gap + bar_width_prop2 + 2 * bar_gap
    legend_x_anchor_norm = (legend_x_anchor - ax.get_xlim()[0]) / (ax.get_xlim()[1] - ax.get_xlim()[0]) # Normalize for bbox_to_anchor

    leg1 = ax.legend(handles=legend_handles_p1, title="Property 1", loc='upper right',
                     bbox_to_anchor=(legend_x_anchor_norm + 0.2, 0.95), # Adjust bbox to move outside plot area
                     fontsize=leaf_font_size * 0.8, title_fontsize=leaf_font_size * 0.9, frameon=False)
    if legend_handles_p2:
        leg2 = ax.legend(handles=legend_handles_p2, title="Property 2", loc='center right',
                         bbox_to_anchor=(legend_x_anchor_norm + 0.2, 0.5),
                         fontsize=leaf_font_size * 0.8, title_fontsize=leaf_font_size * 0.9, frameon=False)
        ax.add_artist(leg1) # Re-add first legend

    # Colorbar for RMSD scale
    max_dist_norm = np.max(Z_linkage[:, 2]) if Z_linkage is not None and Z_linkage.shape[0] > 0 else 3.0
    norm_cb = plt.Normalize(0, min(3.0, max_dist_norm))
    try:
        from mpl_toolkits.axes_grid1 import make_axes_locatable
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("bottom", size="5%", pad=0.3) # Colorbar at bottom
        cb = plt.colorbar(plt.cm.ScalarMappable(norm=norm_cb, cmap=RMSD_GRAYSCALE_CMAP), # Use continuous grayscale
                          cax=cax, orientation='horizontal')
        cb.set_label('RMSD (Å)', fontsize=10)
        cax.tick_params(labelsize=8)
    except Exception as e: print(f"Warning: Could not create colorbar for similarity tree: {e}")

    ax.set_title('Structural Similarity Tree', fontsize=16, pad=15)
    ax.set_xlabel('RMSD (Å)', fontsize=12)
    fig.tight_layout(rect=[0.15, 0.1, 0.80, 0.95]) # Adjust for labels and legends
    return fig, ordered_sids


def visualize_rmsd_matrix_improved(rmsd_df, group_dict=None, name_dict=None,
                                   output_file=None, figsize=(16, 12),  # Adjusted for right legend panel
                                   domain_dict=None, linkage_matrix=None,
                                   rmsd_vmin=0.0, rmsd_vmax=None,
                                   color_mode='continuous', step_cutoffs=None,
                                   step_colors=None):
    """
    Clustermap with continuous white-to-dark_gray RMSD.
    Categorical property legends are at the top-right, with RMSD colorbar directly below them.
    Names next to property color bars are removed.
    """
    pdb_ids_list = rmsd_df.index.tolist()
    if linkage_matrix is None:
        print("ERROR (Clustermap): Linkage matrix (Z) is required.")
        return None
    if len(pdb_ids_list) < 2:
        print("Not enough structures for Clustermap visualization.")
        return None

    try:
        dendrogram_out = dendrogram(linkage_matrix, no_plot=True)
        ordered_ids_for_legend_cats = [pdb_ids_list[i] for i in dendrogram_out['leaves']]
    except Exception as e:
        print(f"Warning (Clustermap): Could not derive order from linkage for legend: {e}.")
        ordered_ids_for_legend_cats = pdb_ids_list

    # --- Prepare Row Colors (Property 1) ---
    row_colors_series, row_legend_map = None, {}
    prop1_name = 'Molecular Function'  # Legend title
    if group_dict:
        unique_prop1 = sorted(
            list(set(str(group_dict.get(id_val, 'Unknown')) for id_val in ordered_ids_for_legend_cats)))
        row_legend_map = get_categorical_colors(unique_prop1, property_type='property1')
        colors_row = [row_legend_map.get(str(group_dict.get(id_val, 'Unknown')), OPSIN_COLORS['gray_3_light'])
                      for id_val in rmsd_df.index]
        row_colors_series = pd.Series(colors_row, index=rmsd_df.index)
        row_colors_series.name = None

    # --- Prepare Column Colors (Property 2) ---
    col_colors_series, col_legend_map = None, {}
    prop2_name = 'Source'  # Legend title
    if domain_dict:
        prop2_vals = [str(domain_dict.get(id_val, {}).get('domain', domain_dict.get(id_val, "Unknown"))) for id_val in
                      ordered_ids_for_legend_cats]
        unique_prop2 = sorted(list(set(prop2_vals)))
        col_legend_map = get_categorical_colors(unique_prop2, property_type='property2')
        colors_col = []
        for id_val in rmsd_df.index:
            info = domain_dict.get(id_val, "Unknown")
            val_str = str(info.get('domain', info)) if isinstance(info, dict) else str(info)
            colors_col.append(col_legend_map.get(val_str, OPSIN_COLORS['gray_3_light']))
        col_colors_series = pd.Series(colors_col, index=rmsd_df.index)
        col_colors_series.name = None
    elif group_dict and row_colors_series is not None:  # Fallback
        col_colors_series = row_colors_series.copy()
        col_legend_map = row_legend_map

    # --- Determine vmin and vmax for RMSD heatmap ---
    if rmsd_vmax is None:
        finite_rmsd_values = rmsd_df.values[np.isfinite(rmsd_df.values)]
        if finite_rmsd_values.size > 0:
            rmsd_vmax_calc = np.max(finite_rmsd_values);
            rmsd_vmax = min(rmsd_vmax_calc, 6.0)
        else:
            rmsd_vmax = 3.0
    if rmsd_vmin is None: rmsd_vmin = 0.0

    # --- Configure colormap / normalization ---
    cmap = 'opsin_rmsd_white_to_darkgray'
    norm = None
    step_tick_positions = []
    step_tick_labels = []

    if color_mode not in {'continuous', 'step'}:
        raise ValueError("color_mode must be 'continuous' or 'step'")

    if color_mode == 'step':
        step_cutoffs = step_cutoffs or [0.5, 1.5, 2.5]
        if sorted(step_cutoffs) != step_cutoffs:
            raise ValueError("step_cutoffs must be in ascending order")
        if step_colors is None:
            step_colors = [
                OPSIN_COLORS['gray_1_white'],
                OPSIN_COLORS['gray_3_light'],
                OPSIN_COLORS['gray_6_dark_mid'],
                OPSIN_COLORS['gray_9_black']
            ]
        num_bins = len(step_colors)
        if num_bins != len(step_cutoffs) + 1:
            raise ValueError("step_colors must have exactly one more entry than step_cutoffs")
        upper_bound = max(rmsd_vmax, step_cutoffs[-1] + 0.5)
        boundaries = [rmsd_vmin] + step_cutoffs + [upper_bound]
        norm = BoundaryNorm(boundaries, num_bins, clip=True)
        cmap = ListedColormap(step_colors)
        rmsd_vmax = upper_bound
        cbar_tick_centers = [(low + high) / 2 for low, high in zip(boundaries[:-1], boundaries[1:])]
        step_tick_positions = cbar_tick_centers
        label_segments = []
        for idx, center in enumerate(cbar_tick_centers):
            if idx == 0:
                label_segments.append(f"≤{step_cutoffs[0]:.1f}")
            elif idx == len(cbar_tick_centers) - 1:
                label_segments.append(f">{step_cutoffs[-1]:.1f}")
            else:
                label_segments.append(f"{step_cutoffs[idx-1]:.1f}–{step_cutoffs[idx]:.1f}")
        step_tick_labels = label_segments
        cbar_kws = {'label': 'RMSD (Å)', 'ticks': step_tick_positions}
    else:
        cbar_kws = {'label': 'RMSD (Å)', 'format': '%.1f'}

    # --- Generate Clustermap ---
    try:
        g = sns.clustermap(
            rmsd_df.fillna(rmsd_vmax + 1),
            row_linkage=linkage_matrix, col_linkage=linkage_matrix,
            row_colors=row_colors_series, col_colors=col_colors_series,
            xticklabels=False, yticklabels=False,
            figsize=figsize,
            cmap=cmap,
            standard_scale=None,
            vmin=rmsd_vmin, vmax=rmsd_vmax,
            norm=norm,
            cbar_kws=cbar_kws,
            dendrogram_ratio=(0.12, 0.12)
        )

        g.fig.suptitle("Structure Similarity Clustermap (RMSD)", fontsize=16)  # y default is fine

        # --- Introduce subtle gap between property bars/dendrograms and heatmap ---
        gap_main = 0.008
        gap_vertical = 0.008

        if hasattr(g, 'ax_heatmap') and g.ax_heatmap is not None:
            heatmap_pos = list(g.ax_heatmap.get_position().bounds)
            gap_x = min(gap_main, heatmap_pos[2] * 0.2)
            gap_y = min(gap_vertical, heatmap_pos[3] * 0.2)
            heatmap_pos[0] += gap_x
            heatmap_pos[2] -= gap_x
            heatmap_pos[3] -= gap_y
            g.ax_heatmap.set_position(heatmap_pos)

        if hasattr(g, 'ax_row_dendrogram') and g.ax_row_dendrogram is not None:
            rd_pos = list(g.ax_row_dendrogram.get_position().bounds)
            rd_gap = min(gap_main * 0.6, rd_pos[2] * 0.5)
            rd_pos[2] -= rd_gap
            g.ax_row_dendrogram.set_position(rd_pos)

        if hasattr(g, 'ax_row_colors') and g.ax_row_colors is not None:
            rc_pos = list(g.ax_row_colors.get_position().bounds)
            rc_gap = min(gap_main * 0.6, rc_pos[2] * 0.5)
            rc_pos[2] -= rc_gap
            g.ax_row_colors.set_position(rc_pos)

        if hasattr(g, 'ax_col_dendrogram') and g.ax_col_dendrogram is not None:
            cd_pos = list(g.ax_col_dendrogram.get_position().bounds)
            cd_gap = min(gap_vertical * 0.6, cd_pos[3] * 0.5)
            cd_pos[3] -= cd_gap
            g.ax_col_dendrogram.set_position(cd_pos)

        if hasattr(g, 'ax_col_colors') and g.ax_col_colors is not None:
            cc_pos = list(g.ax_col_colors.get_position().bounds)
            cc_gap = min(gap_vertical * 0.6, cc_pos[3] * 0.5)
            cc_pos[3] -= cc_gap
            g.ax_col_colors.set_position(cc_pos)

        # --- Remove Names from Property Color Bars ---
        if hasattr(g, 'ax_row_colors') and g.ax_row_colors is not None:
            g.ax_row_colors.set_ylabel('')  # Remove y-label (name) if it exists
            g.ax_row_colors.tick_params(labelleft=False, left=False)  # Remove ticks and labels
        if hasattr(g, 'ax_col_colors') and g.ax_col_colors is not None:
            g.ax_col_colors.set_xlabel('')  # Remove x-label (name) if it exists
            g.ax_col_colors.tick_params(labelbottom=False, bottom=False)  # Remove ticks and labels

        # --- Legend Panel Definition & Positioning ---
        # Adjust main plot area to make space on the right for ALL legends
        # [left, bottom, right, top] in figure coordinates
        # Reduce 'right' to make space. Reduce 'top' a bit for suptitle.
        g.fig.subplots_adjust(left=0.05, bottom=0.05, right=0.75, top=0.93)  # Tune 'right' (e.g. 0.75)

        # Define the full legend panel area on the right
        legend_panel_left = g.fig.subplotpars.right + 0.03  # Start after (adjusted) clustermap area
        legend_panel_bottom = 0.15  # Bottom of legend panel (adjust as needed)
        legend_panel_width = 0.97 - legend_panel_left  # Width of legend panel (up to fig edge)
        legend_panel_height = 0.75  # Height of legend panel (adjust as needed, relative to suptitle and bottom)

        # --- Build separate legends for Molecular Function and Source ---
        handles_prop1 = []
        if row_colors_series is not None and row_legend_map:
            for name, color in sorted(row_legend_map.items()):
                if name == "Unknown":
                    continue
                handles_prop1.append(Patch(facecolor=color, edgecolor=OPSIN_COLORS['gray_7_dark'], label=name))
            if "Unknown" in row_legend_map:
                handles_prop1.append(Patch(facecolor=row_legend_map["Unknown"],
                                           edgecolor=OPSIN_COLORS['gray_7_dark'],
                                           label="Others"))

        handles_prop2 = []
        if col_colors_series is not None and col_legend_map:
            for name, color in sorted(col_legend_map.items()):
                if name == "Unknown":
                    continue
                handles_prop2.append(Patch(facecolor=color, edgecolor=OPSIN_COLORS['gray_7_dark'], label=name))
            if "Unknown" in col_legend_map:
                handles_prop2.append(Patch(facecolor=col_legend_map["Unknown"],
                                           edgecolor=OPSIN_COLORS['gray_7_dark'],
                                           label="Unknown Source"))

        legend_current_top = legend_panel_bottom + legend_panel_height
        legend_gap = 0.015
        cbar_min_height = 0.1
        legend_axes = []

        def _add_legend(handles, title):
            nonlocal legend_current_top
            if not handles:
                return None
            max_height = legend_current_top - (legend_panel_bottom + cbar_min_height)
            if max_height <= 0:
                max_height = legend_current_top - legend_panel_bottom
            if max_height <= 0:
                return None
            legend_height = min(0.035 * len(handles) + 0.06, max_height)
            legend_height = max(min(max_height, legend_height), 0.05)
            legend_current_top -= legend_height
            ax_leg = g.fig.add_axes([
                legend_panel_left,
                legend_current_top,
                legend_panel_width,
                legend_height
            ])
            ax_leg.axis('off')
            legend_obj = ax_leg.legend(handles=handles, title=title, loc='upper left',
                                       fontsize=9, title_fontsize=10, frameon=False)
            legend_axes.append(legend_obj)
            if legend_current_top - legend_gap > legend_panel_bottom + cbar_min_height:
                legend_current_top -= legend_gap
            else:
                legend_current_top = legend_panel_bottom + cbar_min_height
            return legend_obj

        _add_legend(handles_prop1, prop1_name)
        _add_legend(handles_prop2, prop2_name)

        available_cbar_height = legend_current_top - legend_panel_bottom
        if available_cbar_height < cbar_min_height:
            available_cbar_height = cbar_min_height
        if available_cbar_height > legend_panel_height:
            available_cbar_height = legend_panel_height

        # --- Reposition RMSD Colorbar beneath legends ---
        if g.cax is not None:
            g.cax.set_position([
                legend_panel_left + legend_panel_width * 0.1,
                legend_panel_bottom,
                0.035,
                available_cbar_height
            ])
            g.cax.tick_params(labelsize=8)
            g.cax.set_ylabel(g.cax.get_ylabel(), fontsize=9, rotation=90, labelpad=15)
            g.cax.yaxis.set_label_position('left')
            g.cax.yaxis.set_ticks_position('left')

            if color_mode == 'step' and step_tick_positions:
                g.cax.set_yticks(step_tick_positions)
                g.cax.set_yticklabels(step_tick_labels, fontsize=8)

        if output_file:
            g.fig.savefig(output_file, dpi=300)
            print(f"Saved Clustermap figure to {os.path.basename(str(output_file))}")
        return g
    except Exception as e:
        print(f"ERROR (Clustermap): Failed to generate clustermap: {e}")
        import traceback
        traceback.print_exc()
        return None

def visualize_binding_pocket(structure_df, residue_ids_list, retinal_df=None, highlight_residues_dict=None, # Renamed args
                             distance_cutoff=4.0, figsize=(800, 600)):
    """Interactive 3D visualization of a binding pocket."""
    pocket_data = structure_df[structure_df['auth_seq_id'].isin(residue_ids_list)] # Renamed
    traces = []
    for res_id, res_data in pocket_data.groupby('auth_seq_id'): # Renamed
        color = OPSIN_COLORS['gray_6_dark_mid'] # Default
        if highlight_residues_dict and res_id in highlight_residues_dict:
            color = highlight_residues_dict[res_id]
        elif 'helix_num' in res_data.columns and pd.notna(res_data['helix_num'].iloc[0]):
            helix_val = int(res_data['helix_num'].iloc[0]) # Renamed helix_num
            if helix_val in HELIX_NUMBER_COLORS: color = HELIX_NUMBER_COLORS[helix_val]
        res_name_display = str(res_data['res_name3l'].iloc[0] if 'res_name3l' in res_data and pd.notna(res_data['res_name3l'].iloc[0]) else '?')

        traces.append(go.Scatter3d(
            x=res_data['x'], y=res_data['y'], z=res_data['z'], mode='markers',
            marker=dict(size=5, color=color, opacity=0.8),
            text=[f"{res_name_display}{res_id} - {atom}" for atom in res_data['res_atom_name']],
            hoverinfo='text', name=f"{res_name_display}{res_id}"))

    if retinal_df is not None and not retinal_df.empty:
        retinal_color = HELIX_NUMBER_COLORS.get('retinal', OPSIN_COLORS['utility_pink'])
        traces.append(go.Scatter3d(
            x=retinal_df['x'], y=retinal_df['y'], z=retinal_df['z'], mode='markers',
            marker=dict(size=6, color=retinal_color, symbol='diamond', opacity=0.9),
            text=[f"RET - {atom}" for atom in retinal_df['res_atom_name']],
            hoverinfo='text', name="Retinal"))

        if distance_cutoff > 0 and not pocket_data.empty:
            ret_coords = retinal_df[['x', 'y', 'z']].values
            pocket_atom_coords = pocket_data[['x', 'y', 'z']].values # Renamed
            distances = cdist(ret_coords, pocket_atom_coords)
            contacts_indices = np.where(distances < distance_cutoff) # Renamed
            if len(contacts_indices[0]) > 0:
                cx, cy, cz = [], [], [] # Renamed
                for ret_idx, pocket_idx in zip(contacts_indices[0], contacts_indices[1]):
                    rx, ry, rz = ret_coords[ret_idx]
                    px, py, pz = pocket_atom_coords[pocket_idx]
                    cx.extend([rx, px, None]); cy.extend([ry, py, None]); cz.extend([rz, pz, None])
                traces.append(go.Scatter3d(
                    x=cx, y=cy, z=cz, mode='lines',
                    line=dict(color=OPSIN_COLORS['warm_yellow_medium'], width=2, dash='dot'), # Contact line color
                    opacity=0.6, hoverinfo='none', name="Contacts", showlegend=True))

    fig = go.Figure(data=traces)
    fig.update_layout(title="Binding Pocket Visualization",
                      scene=dict(xaxis_title='X', yaxis_title='Y', zaxis_title='Z', aspectmode='data',
                                 xaxis=dict(backgroundcolor=OPSIN_COLORS['gray_1_white']),
                                 yaxis=dict(backgroundcolor=OPSIN_COLORS['gray_1_white']),
                                 zaxis=dict(backgroundcolor=OPSIN_COLORS['gray_1_white'])),
                      width=figsize[0], height=figsize[1], margin=dict(l=0, r=0, b=0, t=40),
                      legend=dict(x=0.01, y=0.99, bgcolor=OPSIN_COLORS['gray_2_lightest'], bordercolor=OPSIN_COLORS['gray_4_light_mid']))
    return fig


def plot_distances_with_std(distance_table, title="Distance to Retinal by Position", figsize=(12, 8), use_ca=False):
    """
    Plot the average distances with standard deviation error bars for TM residues only (helices 1-7).
    Includes horizontal lines at 1.5, 3.9, and 6.0 Angstrom.
    X-axis ticks and labels are removed. Columns with 1 or less non-NaN entries are excluded.
    Text labels are slightly larger for better downscaling.

    Args:
        distance_table: DataFrame with distances
        title: Title for the plot
        figsize: Figure size tuple
        use_ca: Whether to label as CA distances or all-atom distances

    Returns:
        Matplotlib figure
    """
    if distance_table.empty:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No data available (input table empty)",
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        ax.set_title(title, fontsize=18)
        return fig

    # Filter columns with 1 or less non-NaN entries
    original_cols_count = distance_table.shape[1]
    # Use dropna on axis=1 (columns), thresh=2 means keep columns with at least 2 non-NaNs
    distance_table_cleaned = distance_table.dropna(axis=1, thresh=2)
    filtered_cols_count = distance_table_cleaned.shape[1]

    if filtered_cols_count < original_cols_count:
        print(f"Filtered out {original_cols_count - filtered_cols_count} columns with less than 2 non-NaN entries.")

    if distance_table_cleaned.empty:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No data available after filtering columns with <2 non-NaN entries",
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        ax.set_title(title, fontsize=18)
        return fig

    tm_positions_list = []
    for pos_val_str in distance_table_cleaned.columns:  # Iterate over columns of the cleaned table
        s = str(pos_val_str)
        helix_part = ""
        if '.' in s:  # Standard N.YY format
            parts = s.split('.')
            if len(parts) == 2 and parts[0].isdigit() and 1 <= int(parts[0]) <= 7 and parts[1].replace('.', '',
                                                                                                       1).isdigit():
                helix_part = parts[0]
        elif 'x' in s:  # Ballesteros-Weinstein XxYY format
            parts = s.split('x')
            if len(parts) == 2 and parts[0].isdigit() and 1 <= int(parts[0]) <= 7 and parts[1].replace('.', '',
                                                                                                       1).isdigit():
                helix_part = parts[0]

        if helix_part:  # If helix_part was successfully extracted (i.e., it's a TM helix 1-7)
            tm_positions_list.append(pos_val_str)

    if not tm_positions_list:
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, "No TM residues (H1-H7) found in data after filtering",
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        ax.set_title(title, fontsize=18)
        return fig

    tm_dist_table_final = distance_table_cleaned[tm_positions_list]  # Renamed
    mean_dists_final = tm_dist_table_final.mean(skipna=True)  # Renamed
    std_dists_final = tm_dist_table_final.std(skipna=True)  # Renamed

    fig, ax = plt.subplots(figsize=figsize)

    def universal_sort_key(pos_str_val):  # Renamed
        s = str(pos_str_val)
        if '.' in s:
            parts = s.split('.'); return (int(parts[0]), float(parts[1]))
        elif 'x' in s:
            parts = s.split('x'); return (int(parts[0]), float(parts[1]))
        return (99, 0)  # Fallback

    sorted_positions_final = sorted(mean_dists_final.index, key=universal_sort_key)  # Renamed
    x_indices_final = np.arange(len(sorted_positions_final))  # Renamed
    y_means_final = [mean_dists_final[p] for p in sorted_positions_final]  # Renamed
    y_errors_final = [max(0, std_dists_final.get(p, 0)) if pd.notna(std_dists_final.get(p)) else 0 for p in
                      sorted_positions_final]  # Renamed

    # Store points per helix for connecting lines
    helix_data_map = {}

    for i, pos_val_str in enumerate(sorted_positions_final):  # Renamed
        s = str(pos_val_str)
        helix_num_str_part = ""  # Renamed
        if '.' in s:
            helix_num_str_part = s.split('.')[0]
        elif 'x' in s:
            helix_num_str_part = s.split('x')[0]

        helix_idx_val = int(helix_num_str_part) if helix_num_str_part.isdigit() else 0  # Renamed
        point_plot_color = HELIX_NUMBER_COLORS.get(helix_idx_val,
                                                   OPSIN_COLORS['gray_6_dark_mid'])  # Renamed & default color

        ax.errorbar(x_indices_final[i], y_means_final[i], yerr=y_errors_final[i], fmt='o', color=point_plot_color,
                    ecolor=point_plot_color, markersize=8, alpha=0.9, capsize=5,
                    elinewidth=1.5, capthick=1.5)

        if helix_idx_val not in helix_data_map: helix_data_map[helix_idx_val] = []
        helix_data_map[helix_idx_val].append((x_indices_final[i], y_means_final[i]))

    # Plot connecting lines per helix
    for helix_idx_val, points_list in sorted(helix_data_map.items()):  # Renamed
        if points_list:
            line_x_coords, line_y_coords = zip(
                *sorted(points_list, key=lambda pt: pt[0]))  # Ensure points are x-sorted for line
            line_plot_color = HELIX_NUMBER_COLORS.get(helix_idx_val, OPSIN_COLORS['gray_6_dark_mid'])  # Renamed
            ax.plot(line_x_coords, line_y_coords, '-', color=line_plot_color, linewidth=2.0, alpha=0.7)

    # Remove x-axis ticks and labels
    ax.set_xticks([])
    ax.set_xticklabels([])

    max_y_for_text_label = (max(y_means_final) if y_means_final else 10)  # Renamed
    for i, pos_val_str in enumerate(sorted_positions_final):  # Renamed
        s = str(pos_val_str)
        is_x50_type = (s.endswith('.50') and s.split('.')[0].isdigit()) or \
                      (s.endswith('x50') and s.split('x')[0].isdigit())  # BW-like Xx50

        if is_x50_type:
            helix_str_part = s.split('.')[0] if '.' in s else s.split('x')[0]  # Renamed
            if helix_str_part.isdigit():
                helix_idx_val_x50 = int(helix_str_part)  # Renamed
                color_x50 = HELIX_NUMBER_COLORS.get(helix_idx_val_x50, OPSIN_COLORS['gray_5_mid'])  # Renamed
                ax.axvline(x=x_indices_final[i], color=color_x50, linestyle=':', alpha=0.7,
                           linewidth=1.5)  # Changed to ':'
                ax.text(x_indices_final[i], max_y_for_text_label * 1.02, f"H{helix_str_part}",
                        # Reduced y-offset slightly
                        ha='center', va='bottom', fontsize=14, color=color_x50, fontweight='bold')

    # Add 3 horizontal lines
    h_line_style = {'color': OPSIN_COLORS['gray_5_mid'], 'linestyle': '--', 'alpha': 0.7, 'linewidth': 1.2}  # Renamed
    ax.axhline(y=1.5, **h_line_style)
    ax.axhline(y=3.9, **h_line_style)
    ax.axhline(y=6.0, **h_line_style)

    # Text labels for horizontal lines (optional, can be added if needed)
    # current_xlims = ax.get_xlim() # Renamed
    # text_x_pos_lines = current_xlims[1] * 0.98 # Renamed
    # ax.text(text_x_pos_lines, 1.5, '1.5 Å', va='center_baseline', ha='right', color=h_line_style['color'], fontsize=11)
    # ax.text(text_x_pos_lines, 3.9, '3.9 Å', va='center_baseline', ha='right', color=h_line_style['color'], fontsize=11)
    # ax.text(text_x_pos_lines, 6.0, '6.0 Å', va='center_baseline', ha='right', color=h_line_style['color'], fontsize=11)

    ax.grid(True, alpha=0.3, axis='y', linestyle=':')  # Keep y-axis grid, make it dotted
    ax.set_xlabel('GRN Position', fontsize=14, fontweight='bold')  # Cleaner X-label
    ax.set_ylabel('Distance (Å)', fontsize=14, fontweight='bold')  # Cleaner Y-label, same size as H1-H7

    # Make y-axis tick labels bigger (same as H1-H7 labels)
    ax.tick_params(axis='y', labelsize=14)
    ax.tick_params(axis='x', labelsize=14)

    if y_means_final:
        max_y_err_val = max(y_errors_final) if y_errors_final else 0  # Renamed
        data_max_y = max(y_means_final) + max_y_err_val  # Renamed
        plot_upper_y_limit = max(data_max_y, 6.5)  # Ensure 6.0 line is visible # Renamed
        ax.set_ylim(0, plot_upper_y_limit * 1.1 if plot_upper_y_limit > 0 else 10)
    else:
        ax.set_ylim(0, 10)

    # Cleaner title
    atom_type = "Cα" if use_ca else "All-Atom"
    ax.set_title(f"Distance to Retinal ({atom_type})", fontsize=16, fontweight='bold')

    legend_elements_list = [  # Renamed
        Patch(color=HELIX_NUMBER_COLORS[h_idx], label=f"H{h_idx}")
        for h_idx in range(1, 8) if h_idx in HELIX_NUMBER_COLORS
    ]
    if legend_elements_list:
        # Place legend outside the plot area, to the upper right
        ax.legend(handles=legend_elements_list, loc='upper left', bbox_to_anchor=(1.01, 1.0),
                  ncol=1, fontsize=14, title="Helix", title_fontsize=14,
                  frameon=True, framealpha=0.85)

    plt.tight_layout(rect=[0, 0, 0.88, 1])  # Make room for legend on the right
    return fig


def plot_helix_logo_plots(
        residue_table_data,
        figsize=(18, 12),
        frequency_threshold=0.15
):
    """
    Sequence logos around X.50 for TM helices.
    Uses AMINO_ACID_LOGO_COLORS. Amino acids with original frequency below
    frequency_threshold at a position are excluded.
    The heights of the displayed letters represent their original absolute frequencies.
    Y-axis is hardcoded to 0-1 range, and y-axis title AND tick labels are removed.
    X-axis tick labels are the actual GRN column names.
    """
    if residue_table_data is None or residue_table_data.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.text(0.5, 0.5, "No data available for logos", ha='center', va='center')
        return fig

    window_size = 4  # Number of positions to show on either side of X.50
    fig = plt.figure(figsize=figsize)
    gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.3)

    layout_positions = {
        1: (0, 0),
        2: (0, 1),
        3: (0, 2),
        4: (1, 0),
        5: (1, 1),
        6: (2, 0),
        7: (2, 1),
    }

    axes_dict = {}
    shared_axis = None
    for helix_idx in range(1, 8):
        row, col = layout_positions[helix_idx]
        sharey_kwargs = {'sharey': shared_axis} if shared_axis is not None else {}
        current_ax = fig.add_subplot(gs[row, col], **sharey_kwargs)
        if shared_axis is None:
            shared_axis = current_ax
        axes_dict[helix_idx] = current_ax

    for row, col in [(1, 2), (2, 2)]:
        fig.add_subplot(gs[row, col]).axis('off')

    tm_cols = [
        c for c in residue_table_data.columns
        if str(c).count('.') == 1 and
           str(c).split('.')[0].isdigit() and
           1 <= int(str(c).split('.')[0]) <= 7
    ]
    filtered_residue_table = residue_table_data[tm_cols]

    if filtered_residue_table.empty:
        for helix_idx, ax in axes_dict.items():
            ax.text(0.5, 0.5, f"No TM data for Helix {helix_idx}", ha='center', va='center',
                    transform=ax.transAxes)
            ax.set_yticks([])  # Remove y-ticks as well if no data
        fig.tight_layout()
        return fig

    for helix_idx in range(1, 8):
        ax = axes_dict[helix_idx]
        pivot_grn_target = f"{helix_idx}.50"

        helix_grn_cols_sorted = sorted(
            [c for c in filtered_residue_table.columns if str(c).startswith(f"{helix_idx}.")],
            key=lambda p_str: float(str(p_str).split('.')[1])
        )

        actual_pivot_grn = pivot_grn_target
        pivot_col_index_in_helix_list = -1

        if not helix_grn_cols_sorted:
            ax.text(0.5, 0.5, f"No data for H{helix_idx}", ha='center', va='center', transform=ax.transAxes)
            ax.set_ylim(0, 1.0);
            ax.set_ylabel("")
            ax.set_yticks([])  # Remove y-axis ticks and their labels
            continue

        try:
            pivot_col_index_in_helix_list = helix_grn_cols_sorted.index(pivot_grn_target)
        except ValueError:
            closest_grn_match = sorted(
                [(abs(float(str(p_str).split('.')[1]) - 50.0), p_str) for p_str in helix_grn_cols_sorted]
            )
            if closest_grn_match:
                actual_pivot_grn = closest_grn_match[0][1]
                pivot_col_index_in_helix_list = helix_grn_cols_sorted.index(actual_pivot_grn)
            else:
                ax.text(0.5, 0.5, f"No data for H{helix_idx}", ha='center', va='center', transform=ax.transAxes)
                ax.set_ylim(0, 1.0);
                ax.set_ylabel("")
                ax.set_yticks([])
                continue

        start_idx_in_helix_list = max(0, pivot_col_index_in_helix_list - window_size)
        end_idx_in_helix_list = min(len(helix_grn_cols_sorted), pivot_col_index_in_helix_list + window_size + 1)
        window_grn_ids_to_plot = helix_grn_cols_sorted[start_idx_in_helix_list:end_idx_in_helix_list]

        aa_alphabet_for_logo = list(AMINO_ACID_LOGO_COLORS.keys())
        if '-' in aa_alphabet_for_logo:
            aa_alphabet_for_logo.remove('-')

        logo_data_for_plot_df = pd.DataFrame(0.0, index=range(len(window_grn_ids_to_plot)),
                                             columns=aa_alphabet_for_logo)

        for i_window_col, grn_id_in_window in enumerate(window_grn_ids_to_plot):
            aas_at_current_pos = [
                str(cell)[0].upper() for cell in filtered_residue_table[grn_id_in_window].dropna()
                if isinstance(cell, str) and cell and cell[0] != '-' and cell[0].upper() in aa_alphabet_for_logo
            ]
            if aas_at_current_pos:
                aa_counts_map = Counter(aas_at_current_pos)
                total_aa_count_at_pos = sum(aa_counts_map.values())

                if total_aa_count_at_pos > 0:
                    for aa_code_val, num_count in aa_counts_map.items():
                        original_frequency = num_count / total_aa_count_at_pos
                        if original_frequency >= frequency_threshold:
                            logo_data_for_plot_df.loc[i_window_col, aa_code_val] = original_frequency

        # X-axis tick labels are now the GRN IDs themselves
        x_tick_labels_for_plot = window_grn_ids_to_plot

        if not logo_data_for_plot_df.empty and logo_data_for_plot_df.sum().sum() > 1e-6:
            try:
                logomaker.Logo(logo_data_for_plot_df, ax=ax, color_scheme=AMINO_ACID_LOGO_COLORS)

                ax.set_title(f"Helix {helix_idx}", color=HELIX_NUMBER_COLORS.get(helix_idx, OPSIN_COLORS['black']),
                             fontweight='bold')
                # X-axis label can be more generic or removed if GRNs are self-explanatory
                ax.set_xlabel("Position (GRN)", fontsize=10)  # Or ax.set_xlabel("")
                ax.set_xticks(range(len(window_grn_ids_to_plot)))
                ax.set_xticklabels(x_tick_labels_for_plot, fontsize=8, rotation=45, ha="right")  # Smaller font, rotated

                try:
                    pivot_in_window_plot_idx = window_grn_ids_to_plot.index(actual_pivot_grn)
                    ax.axvline(x=pivot_in_window_plot_idx,
                               color=HELIX_NUMBER_COLORS.get(helix_idx, OPSIN_COLORS['gray_5_mid']), linestyle=':',
                               alpha=0.6)
                except ValueError:
                    pass

            except Exception as e_logo:
                ax.text(0.5, 0.5, f"Logo Error H{helix_idx}", ha='center', va='center', transform=ax.transAxes,
                        fontsize=8)
                print(f"Error creating logo for Helix {helix_idx}: {e_logo}")
        else:
            ax.text(0.5, 0.5, f"No AA data (freq >= {frequency_threshold * 100:.0f}%) for H{helix_idx}", ha='center',
                    va='center', transform=ax.transAxes)

        ax.set_ylim(0, 1.0)
        ax.set_ylabel("")  # Remove y-axis title
        ax.set_yticks([])  # Remove y-axis ticks and their labels
        # If sharey=True, ensure ticks removed on shared axis as well.

    for ax_iter in axes_dict.values():  # Redundant but ensures all axes are consistent
        ax_iter.set_yticks([])

    fig.suptitle(f"Prominent Amino Acids (Freq >= {frequency_threshold * 100:.0f}%) around X.50", fontsize=14, y=0.98)
    fig.tight_layout(rect=[0.03, 0.03, 1, 0.95])  # May need to adjust rect bottom if x-labels are rotated
    return fig


def plot_average_distances_by_helix(distance_table_data, use_ca=True): # Renamed
    """Plot average distances per position, grouped by helix, using helix colors."""
    if distance_table_data.empty:
        fig, ax = plt.subplots(figsize=(14,8)); ax.text(0.5,0.5,"No data",ha='center',va='center'); return fig

    mean_dists_all = distance_table_data.mean(skipna=True) # Renamed
    tm_mean_dists = mean_dists_all[ # Renamed
        ~mean_dists_all.index.astype(str).str.contains(r'^[Lnc]\.') # Regex for L. n. c.
    ]
    plot_points = [] # Renamed
    for pos_grn, dist_val in tm_mean_dists.items(): # Renamed
        pos_str = str(pos_grn)
        helix_id_str, pos_num_str = ('', '') # Renamed
        if 'x' in pos_str: helix_id_str, pos_num_str = pos_str.split('x', 1)
        elif '.' in pos_str: helix_id_str, pos_num_str = pos_str.split('.', 1)
        if helix_id_str.isdigit() and pos_num_str.replace('.','',1).isdigit():
            plot_points.append({'Helix': helix_id_str, 'Position': float(pos_num_str), 'Distance': dist_val}) # Renamed

    if not plot_points:
        fig, ax = plt.subplots(figsize=(14,8)); ax.text(0.5,0.5,"No TM data",ha='center',va='center'); return fig

    plot_data_df = pd.DataFrame(plot_points) # Renamed
    fig, ax = plt.subplots(figsize=(14, 8))
    for helix_id, group_df in plot_data_df.groupby('Helix'): # Renamed
        sorted_group_df = group_df.sort_values('Position') # Renamed
        color = HELIX_STRING_COLORS.get(str(helix_id), OPSIN_COLORS['gray_6_dark_mid']) # Use string for HELIX_STRING_COLORS
        ax.plot(sorted_group_df['Position'], sorted_group_df['Distance'], marker='o', linestyle='-', label=f'Helix {helix_id}', color=color, linewidth=1.8)

    ax.set_xlabel('Position Number (within helix)', fontsize=12) # Clarified X-axis
    ax.set_ylabel(f'Avg Distance to Retinal (Å) - {"CA" if use_ca else "All Atom"}', fontsize=12)
    ax.set_title(f'Average Distance to Retinal by Residue ({ "CA" if use_ca else "All Atom"})', fontsize=14)

    # Highlight .50 position for all plotted helices
    plotted_helix_nums = [int(h) for h in plot_data_df['Helix'].unique() if h.isdigit()]
    if 50.0 in plot_data_df['Position'].unique(): # Check if 50 is a position number
        for h_num in plotted_helix_nums:
            ax.axvline(x=50, color=HELIX_NUMBER_COLORS.get(h_num, OPSIN_COLORS['gray_4_light_mid']), linestyle=':', alpha=0.4)
        ax.text(50, ax.get_ylim()[1]*0.95, "X.50", ha='center', va='top', bbox=dict(facecolor=OPSIN_COLORS['white'], alpha=0.7))

    ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left', fontsize=9); ax.grid(True, linestyle=':', alpha=0.5)
    plt.tight_layout(rect=[0,0,0.88,1]); return fig


def plot_distance_heatmap(distance_table_data): # Renamed
    """Heatmap of average distances, ordered by helix/position, using distance colormap."""
    if distance_table_data.empty:
        fig, ax = plt.subplots(figsize=(16,8)); ax.text(0.5,0.5,"No data",ha='center',va='center'); return fig

    mean_dists_all = distance_table_data.mean(skipna=True) # Renamed
    tm_mean_dists = mean_dists_all[~mean_dists_all.index.astype(str).str.contains(r'^[Lnc]\.')]
    plot_points = []
    for pos_grn, dist_val in tm_mean_dists.items():
        pos_str = str(pos_grn); helix_id_str, pos_num_str = ('', '')
        if 'x' in pos_str: helix_id_str, pos_num_str = pos_str.split('x', 1)
        elif '.' in pos_str: helix_id_str, pos_num_str = pos_str.split('.', 1)
        if helix_id_str.isdigit() and pos_num_str.replace('.','',1).isdigit():
            plot_points.append({'Helix': int(helix_id_str), 'Position': float(pos_num_str), 'Distance': dist_val})

    if not plot_points:
        fig, ax = plt.subplots(figsize=(16,8)); ax.text(0.5,0.5,"No TM data for heatmap",ha='center',va='center'); return fig

    plot_data_df = pd.DataFrame(plot_points).sort_values(by=['Helix', 'Position'])
    pivot_data_df = plot_data_df.pivot(index='Helix', columns='Position', values='Distance') # Renamed

    fig, ax = plt.subplots(figsize=(16, 8))
    # Use DEFAULT_DISTANCE_CMAP (which is DISTANCE_WARM_REVERSED_CMAP)
    sns.heatmap(pivot_data_df, cmap=DEFAULT_DISTANCE_CMAP, ax=ax, # Or 'opsin_distance_default'
                cbar_kws={'label': 'Avg Distance to Retinal (Å)'},
                vmin=0, vmax=max(10, pivot_data_df.max().max() if not pivot_data_df.empty else 10), annot=False) # Adjust vmax

    ax.set_title('Distance to Retinal Heatmap (TM Helices)', fontsize=14) # Renamed
    ax.set_xlabel('Position Number (within helix)', fontsize=12); ax.set_ylabel('Helix', fontsize=12)

    if 50.0 in pivot_data_df.columns:
        try:
            col_idx_50 = pivot_data_df.columns.get_loc(50.0)
            ax.axvline(x=col_idx_50 + 0.5, color=OPSIN_COLORS['white'], linestyle=':', alpha=0.8, linewidth=1.2)
        except KeyError: pass

    ax.set_yticks(np.arange(len(pivot_data_df.index)) + 0.5)
    ax.set_yticklabels(pivot_data_df.index, fontweight='bold')
    for i, tick_label_obj in enumerate(ax.get_yticklabels()): # Renamed
        helix_val_int = pivot_data_df.index[i] # Renamed
        if helix_val_int in HELIX_NUMBER_COLORS: tick_label_obj.set_color(HELIX_NUMBER_COLORS[helix_val_int])
    plt.tight_layout(); return fig


def print_residue_composition(composition_results_dict, highlight_thresh=20.0): # Renamed args
    """Pretty-print residue composition with ANSI highlighting."""
    print("\nResidue Composition at Key Positions:\n" + "="*37)
    for pos_key, comp_data in composition_results_dict.items(): # Renamed
        print(f"\nPosition {pos_key}:\n" + "-"*(len(str(pos_key))+10))
        if isinstance(comp_data, dict) and 'error' in comp_data: print(f"  {comp_data['error']}"); continue

        freq_map_to_print = None # Renamed
        if isinstance(comp_data, dict):
            if 'frequencies' in comp_data and isinstance(comp_data['frequencies'], dict): freq_map_to_print = comp_data['frequencies']
            elif 'sorted' in comp_data and isinstance(comp_data['sorted'], list): freq_map_to_print = dict(comp_data['sorted'])
            else: freq_map_to_print = {k:v for k,v in comp_data.items() if isinstance(v, (int,float))} # Fallback for flat dict

        if freq_map_to_print:
            for aa_code, freq_val in sorted(freq_map_to_print.items(), key=lambda item: item[1], reverse=True):
                percent_val = freq_val * 100 # Renamed
                highlight_on = "\033[1m" if percent_val >= highlight_thresh else ""
                highlight_off = "\033[0m" if percent_val >= highlight_thresh else ""
                print(f"  {highlight_on}{aa_code}: {percent_val:.1f}%{highlight_off}")
        else: print(f"  [Info: No parsable frequency data for {pos_key}]")


def calculate_helix_distances(distance_table_data): # Renamed
    """Calculate mean/std distances by helix (TM 1-7)."""
    if distance_table_data.empty: return {}
    mean_dists = distance_table_data.mean(skipna=True); std_dists = distance_table_data.std(skipna=True)
    helix_analysis_stats = {} # Renamed
    for pos_grn in mean_dists.index:
        pos_str = str(pos_grn)
        if pos_str.count('.') == 1:
            helix_id_str, pos_num_part = pos_str.split('.',1) # Renamed
            if helix_id_str.isdigit() and 1 <= int(helix_id_str) <= 7 and pos_num_part.replace('.','',1).isdigit():
                if helix_id_str not in helix_analysis_stats:
                    helix_analysis_stats[helix_id_str] = {'positions':[],'means':[],'stds':[],'grn_pos':[]}
                if pos_grn in mean_dists: # Check if pos_grn exists in mean_dists
                    helix_analysis_stats[helix_id_str]['grn_pos'].append(pos_grn) # Store original GRN
                    helix_analysis_stats[helix_id_str]['positions'].append(float(pos_num_part)) # Store numeric part for sorting
                    helix_analysis_stats[helix_id_str]['means'].append(mean_dists[pos_grn])
                    helix_analysis_stats[helix_id_str]['stds'].append(std_dists.get(pos_grn, np.nan))

    for helix_id_str, data_map in helix_analysis_stats.items(): # Renamed
        if data_map['positions']:
            # Sort by numeric position part, then reconstruct original GRNs in sorted order
            sorted_indices = sorted(range(len(data_map['positions'])), key=lambda i: data_map['positions'][i])
            data_map['grn_pos'] = [data_map['grn_pos'][i] for i in sorted_indices] # Sorted original GRNs
            data_map['positions'] = [data_map['positions'][i] for i in sorted_indices] # Sorted numeric parts
            data_map['means'] = [data_map['means'][i] for i in sorted_indices]
            data_map['stds'] = [data_map['stds'][i] for i in sorted_indices]

            valid_means_list = [m for m in data_map['means'] if pd.notna(m)] # Renamed
            if valid_means_list:
                min_mean = min(valid_means_list)
                min_idx_val = data_map['means'].index(min_mean) # Renamed
                data_map['closest_grn'] = data_map['grn_pos'][min_idx_val] # Use sorted original GRN
                data_map['closest_mean_dist'] = min_mean # Renamed
            else: data_map['closest_grn'] = None; data_map['closest_mean_dist'] = None
    return helix_analysis_stats


def visualize_msa_distances(msa_results_dict, output_dir_path="./", file_prefix=""): # Renamed
    """Orchestrates distance visualizations."""
    os.makedirs(output_dir_path, exist_ok=True)
    ca_dist_table = msa_results_dict.get("ca_distance_table", pd.DataFrame()) # Renamed
    all_atom_dist_table = msa_results_dict.get("distance_table", pd.DataFrame()) # Renamed
    plot_paths_dict = {} # Renamed

    if not all_atom_dist_table.empty:
        stats = calculate_helix_distances(all_atom_dist_table)
        print("\nClosest residues (All-atom):")
        for h_id in sorted(stats.keys(), key=int): # Renamed
            if stats[h_id].get('closest_grn'): print(f"  H{h_id}: {stats[h_id]['closest_grn']} ({stats[h_id]['closest_mean_dist']:.2f}Å)")
        fig_all = plot_distances_with_std(all_atom_dist_table, title="All-Atom Distances to Retinal", use_ca=False)
        path_all = os.path.join(output_dir_path, f"{file_prefix}all_atom_distances_std.png"); fig_all.savefig(path_all); plt.close(fig_all)
        plot_paths_dict["all_atom_distances_std"] = path_all
    else: print("Skipping All-atom distance plot: table empty.")

    if not ca_dist_table.empty:
        stats_ca = calculate_helix_distances(ca_dist_table) # Renamed
        print("\nClosest residues (CA-atom):")
        for h_id in sorted(stats_ca.keys(), key=int):
            if stats_ca[h_id].get('closest_grn'): print(f"  H{h_id}: {stats_ca[h_id]['closest_grn']} ({stats_ca[h_id]['closest_mean_dist']:.2f}Å)")
        fig_ca = plot_distances_with_std(ca_dist_table, title="CA-Atom Distances to Retinal", use_ca=True)
        path_ca = os.path.join(output_dir_path, f"{file_prefix}ca_atom_distances_std.png"); fig_ca.savefig(path_ca); plt.close(fig_ca)
        plot_paths_dict["ca_atom_distances_std"] = path_ca
    else: print("Skipping CA-atom distance plot: table empty.")

    msa_data_for_logo = msa_results_dict.get("residue_table", msa_results_dict.get("msa_table", pd.DataFrame())) # Renamed
    if not msa_data_for_logo.empty:
        try:
            fig_logo = plot_helix_logo_plots(msa_data_for_logo, figsize=(20,5)) # Renamed
            path_logo = os.path.join(output_dir_path, f"{file_prefix}sequence_logos_x50.png"); fig_logo.savefig(path_logo); plt.close(fig_logo)
            plot_paths_dict["sequence_logos_x50"] = path_logo
        except Exception as e: print(f"WARNING: Logo plot failed: {e}"); import traceback; traceback.print_exc()
    else: print("Skipping X.50 logo plots: MSA data empty.")
    print(f"\nDistance visualizations saved to {output_dir_path}"); return plot_paths_dict


def create_combined_distance_logo_plot(distance_table_data, msa_data_df): # Renamed
    """Combined distance line plot and sequence logo."""
    if distance_table_data.empty or msa_data_df.empty:
        fig, ax = plt.subplots(figsize=(15,10)); ax.text(0.5,0.5,"Empty data",ha='center',va='center'); return fig

    fig = plt.figure(figsize=(15,10)); gs = gridspec.GridSpec(2,1,height_ratios=[1,0.6],hspace=0.05)
    ax_top = fig.add_subplot(gs[0]); ax_bottom = fig.add_subplot(gs[1], sharex=ax_top)

    mean_dists_all = distance_table_data.mean(skipna=True) # Renamed
    common_grns = mean_dists_all.index.intersection(msa_data_df.columns) # Renamed
    if common_grns.empty:
        ax_top.text(0.5,0.5,"No common GRNs",ha='center',va='center'); return fig

    def grn_universal_sort_key(grn_str_val): # Renamed
        s = str(grn_str_val)
        if s.startswith('n.'): p = s.split('.'); return (-2, float(p[1]) if len(p)>1 and p[1].replace('.','',1).isdigit() else 0)
        if s.count('.')==1 and s.split('.')[0].isdigit(): p=s.split('.'); return (int(p[0]), float(p[1])) # TM: 1.50
        if 'x' in s: p=s.split('x'); return(int(p[0]),float(p[1])) if p[0].isdigit() and p[1].replace('.','',1).isdigit() else (8,0) # TMx: 1x50
        if s.startswith('L.'): p=s.split('.'); return (80, float(p[1]) if len(p)>1 and p[1].replace('.','',1).isdigit() else 0) # Loop
        if s.startswith('c.'): p=s.split('.'); return (90, float(p[1]) if len(p)>1 and p[1].replace('.','',1).isdigit() else 0) # C-term
        return (999, 0) # Fallback

    sorted_common_grns_list = sorted(common_grns, key=grn_universal_sort_key) # Renamed
    x_indices_plot = np.arange(len(sorted_common_grns_list)) # Renamed
    y_dist_values = [mean_dists_all[g] for g in sorted_common_grns_list] # Renamed

    ax_top.plot(x_indices_plot, y_dist_values, 'o-', color=OPSIN_COLORS['black'], linewidth=1.5) # Neutral color
    ax_top.set_title('Distance Profile & Sequence Conservation', fontsize=14); ax_top.set_ylabel('Avg Distance (Å)', fontsize=12)
    ax_top.grid(True, alpha=0.3, linestyle=':')

    max_y_val_dist = (max(y_dist_values) if y_dist_values else 10) # Renamed
    for i, grn_val in enumerate(sorted_common_grns_list): # Renamed
        s = str(grn_val)
        if (s.endswith('.50') and s.split('.')[0].isdigit()) or ('x50' in s and s.split('x')[0].isdigit()):
            ax_top.axvline(x=i, color=OPSIN_COLORS['gray_5_mid'], linestyle=':', alpha=0.7)
            ax_top.text(i, max_y_val_dist*0.95, s, rotation=90, ha='center', va='top', fontsize=8, color=OPSIN_COLORS['gray_7_dark'])
    plt.setp(ax_top.get_xticklabels(), visible=False)

    # Sequence Logo
    msa_for_logo_plot = msa_data_df[sorted_common_grns_list] # Renamed
    aa_alphabet = list(AMINO_ACID_LOGO_COLORS.keys()); aa_alphabet.remove('-') # Renamed
    pwm_data_df = pd.DataFrame(0.0, index=x_indices_plot, columns=aa_alphabet) # Renamed
    for i, grn_col_id in enumerate(sorted_common_grns_list): # Renamed
        aas_at_pos = [str(cell)[0].upper() for cell in msa_for_logo_plot[grn_col_id].dropna() if isinstance(cell,str) and cell and cell[0]!='-' and cell[0].upper() in aa_alphabet]
        if aas_at_pos:
            counts = Counter(aas_at_pos); total = sum(counts.values())
            for aa_code, num in counts.items(): pwm_data_df.loc[i, aa_code] = num/total
    if not pwm_data_df.empty and pwm_data_df.sum().sum() > 0:
        try:
            logomaker.Logo(pwm_data_df, ax=ax_bottom, color_scheme=AMINO_ACID_LOGO_COLORS) # Use our scheme
            ax_bottom.set_xticks(x_indices_plot); ax_bottom.set_xticklabels(sorted_common_grns_list, rotation=90, fontsize=8)
        except Exception as e: ax_bottom.text(0.5,0.5,"Logo Error",ha='center',va='center'); print(f"Combined Logo Error: {e}")
    else: ax_bottom.text(0.5,0.5,"No data for logo",ha='center',va='center')
    ax_bottom.set_xlabel('Position (GRN)',fontsize=12); ax_bottom.set_ylabel('Information (bits)',fontsize=12)
    ax_bottom.set_ylim(0, np.log2(len(aa_alphabet)))
    for i, grn_val in enumerate(sorted_common_grns_list):
        s = str(grn_val)
        if (s.endswith('.50') and s.split('.')[0].isdigit()) or ('x50' in s and s.split('x')[0].isdigit()):
            ax_bottom.axvline(x=i, color=OPSIN_COLORS['gray_5_mid'], linestyle=':', alpha=0.7)
    if x_indices_plot.size > 0:
        ax_top.set_xlim(-0.5, x_indices_plot[-1]+0.5); # ax_bottom shares xlim
    plt.tight_layout(rect=[0,0.03,1,0.95]); return fig


def create_opsin_overview_plot(opsin_df, output_path_str=None, figsize_tuple=(16,16)): # Renamed
    """Circular overview plot highlighting molecular function and source domain."""
    df_proc = opsin_df.copy() # Renamed
    df_proc['molecular_function'] = df_proc['molecular_function'].fillna('Unknown').astype(str)
    df_proc['domain'] = df_proc['domain'].fillna('Unknown').astype(str)
    df_sorted_plot = df_proc.sort_values(by=['molecular_function','domain','short_name'],ignore_index=True) # Renamed
    N_items = len(df_sorted_plot) # Renamed
    if N_items == 0:
        fig, ax = plt.subplots(figsize=figsize_tuple); ax.text(0.5,0.5,"No data",ha='center',va='center'); return fig

    # Use get_categorical_colors for consistent property coloring
    prop1_colors_map = get_categorical_colors(df_sorted_plot['molecular_function'].unique(), property_type='property1')
    prop2_colors_map = get_categorical_colors(df_sorted_plot['domain'].unique(), property_type='property2')

    ring1_plot_colors = [prop1_colors_map.get(f, OPSIN_COLORS['gray_3_light']) for f in df_sorted_plot['molecular_function']]
    ring2_plot_colors = [prop2_colors_map.get(d, OPSIN_COLORS['gray_3_light']) for d in df_sorted_plot['domain']]
    # Ensure 'experimentally_determined' is boolean
    if 'experimentally_determined' in df_sorted_plot.columns and df_sorted_plot['experimentally_determined'].dtype != bool:
        df_sorted_plot['experimentally_determined'] = df_sorted_plot['experimentally_determined'].astype(str).str.lower().map({'true':True,'1':True,'yes':True}).fillna(False)
    elif 'experimentally_determined' not in df_sorted_plot.columns:
        df_sorted_plot['experimentally_determined'] = False

    ring3_dot_labels = ['•' if x else '' for x in df_sorted_plot['experimentally_determined']]
    # Assuming 'is_predicted' or similar for Ring 4. If not present, all dots will appear.
    # For now, let's assume if 'is_predicted' is not there, we show all as if predicted.
    if 'is_predicted' in df_sorted_plot.columns and df_sorted_plot['is_predicted'].dtype != bool:
        df_sorted_plot['is_predicted'] = df_sorted_plot['is_predicted'].astype(str).str.lower().map({'true':True,'1':True,'yes':True}).fillna(False)
    elif 'is_predicted' not in df_sorted_plot.columns:
         df_sorted_plot['is_predicted'] = True # Fallback: assume predicted if column missing

    # Ring 4 should display a dot for every entry to reflect that all structures are modelled
    ring4_dot_labels = ['•'] * N_items


    fig = plt.figure(figsize=figsize_tuple); ax = fig.add_axes([0.05,0.05,0.7,0.9])
    R1, R2, R3, R4 = 0.65, 0.85, 0.93, 0.98 # Renamed radii
    SLICE_W, DOT_W = 0.18, 0.02 # Renamed widths (smaller dots)
    DOT_FONT_SIZE = 18
    ring_dot_textprops = {'va':'center','ha':'center','fontsize':DOT_FONT_SIZE}

    ax.pie([1]*N_items, colors=ring1_plot_colors, radius=R1, wedgeprops=dict(width=SLICE_W, edgecolor=OPSIN_COLORS['white']), startangle=90, counterclock=False)
    ax.pie([1]*N_items, colors=ring2_plot_colors, radius=R2, wedgeprops=dict(width=SLICE_W, edgecolor=OPSIN_COLORS['white']), startangle=90, counterclock=False)
    ax.pie([1]*N_items, labels=ring3_dot_labels, radius=R3, labeldistance=1.02,
           wedgeprops=dict(width=DOT_W, edgecolor=OPSIN_COLORS['white'], facecolor='none'), startangle=90, counterclock=False,
           textprops=ring_dot_textprops | {'color': STATUS_EXPERIMENTAL_COLOR})
    ax.pie([1]*N_items, labels=ring4_dot_labels, radius=R4, labeldistance=1.02,
           wedgeprops=dict(width=DOT_W, edgecolor=OPSIN_COLORS['white'], facecolor='none'), startangle=90, counterclock=False,
           textprops=ring_dot_textprops | {'color': STATUS_PREDICTED_COLOR})
    ax.set_title("Opsin Structures Overview", pad=25, fontsize=18) # Adjusted pad

    # Legends
    leg_elements_p1 = [Patch(facecolor=c, label=p1) for p1,c in sorted(prop1_colors_map.items()) if p1 != 'Unknown' and p1 in df_sorted_plot['molecular_function'].unique()]
    leg_elements_p2 = [Patch(facecolor=c, label=p2) for p2,c in sorted(prop2_colors_map.items()) if p2 != 'Unknown' and p2 in df_sorted_plot['domain'].unique()]
    # Add Unknown if present
    if 'Unknown' in prop1_colors_map and 'Unknown' in df_sorted_plot['molecular_function'].unique():
        leg_elements_p1.append(Patch(facecolor=prop1_colors_map['Unknown'], label='Unknown Molecular Function'))
    if 'Unknown' in prop2_colors_map and 'Unknown' in df_sorted_plot['domain'].unique():
        leg_elements_p2.append(Patch(facecolor=prop2_colors_map['Unknown'], label='Unknown Source'))


    leg_elements_status_dots = [
        Line2D([0], [0], marker='o', color='none', markerfacecolor=STATUS_EXPERIMENTAL_COLOR,
               markeredgecolor=STATUS_EXPERIMENTAL_COLOR, markersize=4, label='Experimental'),
        Line2D([0], [0], marker='o', color='none', markerfacecolor=STATUS_PREDICTED_COLOR,
               markeredgecolor=STATUS_PREDICTED_COLOR, markersize=4, label='Modelled')
    ]

    leg_ax_plot = fig.add_axes([0.72, 0.15, 0.25, 0.7]); leg_ax_plot.axis('off') # Renamed
    all_created_legends = [] # Renamed
    if leg_elements_p1:
        leg1_plot = leg_ax_plot.legend(handles=leg_elements_p1, title="Molecular Function", loc="upper left", fontsize=9, title_fontsize=10, frameon=False)
        all_created_legends.append(leg1_plot)
    current_y = 0.75 # Adjusted for potentially more items
    if leg_elements_p2:
        leg2_plot = leg_ax_plot.legend(handles=leg_elements_p2, title="Source", loc="upper left", bbox_to_anchor=(0,current_y), fontsize=9, title_fontsize=10, frameon=False)
        all_created_legends.append(leg2_plot)
        current_y -= (len(leg_elements_p2) * 0.035 + 0.1) # Dynamic offset attempt
    if leg_elements_status_dots:
        leg3_plot = leg_ax_plot.legend(handles=leg_elements_status_dots, loc="upper left", bbox_to_anchor=(0,current_y), fontsize=9, frameon=False)
        all_created_legends.append(leg3_plot)

    for i_leg in range(len(all_created_legends) -1 ): # Re-add all but the last
        if all_created_legends[i_leg]: leg_ax_plot.add_artist(all_created_legends[i_leg])

    if output_path_str: plt.savefig(output_path_str, dpi=300, bbox_inches='tight')
    return fig


def plot_error_violin(set_a_errors, set_b_errors,
                      labels=("Set A (Training)", "Set B (Validation)"),
                      error_metric_label="Retinal RMSD (Å)",
                      title="Prediction Error Distributions",
                      output_path=None):
    """Render side-by-side violin plots comparing two error distributions."""
    set_a_series = pd.Series(set_a_errors, dtype=float).dropna()
    set_b_series = pd.Series(set_b_errors, dtype=float).dropna()

    if set_a_series.empty and set_b_series.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No error data available", ha='center', va='center', fontsize=12)
        ax.axis('off')
        if output_path:
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
        return fig

    data_frames = []
    if not set_a_series.empty:
        data_frames.append(pd.DataFrame({
            'Dataset': labels[0],
            'Error': set_a_series
        }))
    if not set_b_series.empty:
        data_frames.append(pd.DataFrame({
            'Dataset': labels[1],
            'Error': set_b_series
        }))

    combined_df = pd.concat(data_frames, ignore_index=True)

    fig, ax = plt.subplots(figsize=(7, 6))
    palette = [OPSIN_COLORS['warm_orange_medium'], OPSIN_COLORS['cold_blue_medium']]
    sns.violinplot(data=combined_df, x='Dataset', y='Error', palette=palette[:len(data_frames)],
                   cut=0, inner='quartile', linewidth=1.1, saturation=0.85, ax=ax)
    sns.stripplot(data=combined_df, x='Dataset', y='Error', color=OPSIN_COLORS['gray_7_dark'],
                  size=3, alpha=0.45, jitter=True, ax=ax)

    ax.set_xlabel('Dataset', fontsize=12)
    ax.set_ylabel(error_metric_label, fontsize=12)
    ax.set_title(title, fontsize=15, pad=12)
    ax.grid(axis='y', linestyle='--', alpha=0.3)

    medians = combined_df.groupby('Dataset')['Error'].median()
    for idx, dataset_label in enumerate(combined_df['Dataset'].unique()):
        if dataset_label in medians:
            median_val = medians[dataset_label]
            ax.text(idx, median_val, f"Median: {median_val:.2f}",
                    ha='center', va='bottom', fontsize=10,
                    color=OPSIN_COLORS['gray_7_dark'], fontweight='bold')

    ax.set_ylim(bottom=0)
    plt.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    return fig


def plot_error_violin_panel(error_df, dataset_col='Dataset',
                            metrics=None,
                            palette=None,
                            figsize=(18, 6),
                            title="Error Distributions by Dataset",
                            output_path=None):
    """Plot multiple error metrics as side-by-side violin plots."""
    if metrics is None:
        metrics = [
            ('backbone_rmsd', 'Backbone RMSD (Å)'),
            ('pocket_rmsd', 'Binding Pocket RMSD (Å)'),
            ('retinal_rmsd', 'Retinal RMSD (Å)')
        ]

    filtered_df = error_df.copy()
    filtered_df = filtered_df.dropna(subset=[dataset_col])
    if filtered_df.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No error data available", ha='center', va='center', fontsize=12)
        ax.axis('off')
        if output_path:
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
        return fig

    available_metrics = [(col, label) for col, label in metrics if col in filtered_df.columns]
    if not available_metrics:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "Specified metrics not found", ha='center', va='center', fontsize=12)
        ax.axis('off')
        if output_path:
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
        return fig

    categories = filtered_df[dataset_col].dropna().unique().tolist()
    if palette is None:
        base_colors = [OPSIN_COLORS['warm_orange_medium'], OPSIN_COLORS['cold_blue_medium'], OPSIN_COLORS['utility_teal']]
        palette = {cat: base_colors[i % len(base_colors)] for i, cat in enumerate(categories)}

    fig, axes = plt.subplots(1, len(available_metrics), figsize=figsize, sharey=False)
    if len(available_metrics) == 1:
        axes = [axes]

    for ax, (metric_col, metric_label) in zip(axes, available_metrics):
        sns.violinplot(data=filtered_df, x=dataset_col, y=metric_col,
                       order=categories, palette=palette, cut=0,
                       inner='quartile', linewidth=1.1, saturation=0.85, ax=ax)
        sns.stripplot(data=filtered_df, x=dataset_col, y=metric_col,
                      order=categories, color=OPSIN_COLORS['gray_7_dark'],
                      size=3, alpha=0.45, jitter=True, ax=ax)

        ax.set_xlabel('Dataset', fontsize=12)
        ax.set_ylabel(metric_label, fontsize=12)
        ax.grid(axis='y', linestyle='--', alpha=0.3)
        medians = filtered_df.groupby(dataset_col)[metric_col].median()
        for idx, cat in enumerate(categories):
            if cat in medians and pd.notna(medians[cat]):
                ax.text(idx, medians[cat], f"Median: {medians[cat]:.2f}",
                        ha='center', va='bottom', fontsize=10,
                        color=OPSIN_COLORS['gray_7_dark'], fontweight='bold')
        ax.set_ylim(bottom=0)

    fig.suptitle(title, fontsize=16, y=1.02)
    fig.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')
    return fig


def plot_error_box_comparison(set_a_df: pd.DataFrame,
                              set_b_df: pd.DataFrame,
                              metrics=None,
                              dataset_labels=("Benchmark set", "Blind test set"),
                              output_path=None,
                              figure_size=(13, 6)):
    """Plot paired boxplots with muted colors for the key RMSD metrics."""
    metrics = metrics or ["backbone_rmsd", "pocket_rmsd", "retinal_rmsd"]
    metric_display = {
        "backbone_rmsd": "Cα RMSD",
        "pocket_rmsd": "Binding Pocket iRMSD",
        "retinal_rmsd": "Retinal L-RMSD",
    }

    set_a = set_a_df.copy()
    set_b = set_b_df.copy()
    if set_a.empty and set_b.empty:
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No error data available", ha='center', va='center', fontsize=12)
        ax.axis('off')
        if output_path:
            fig.savefig(output_path, dpi=300, bbox_inches='tight')
        return fig, pd.DataFrame()

    if 'protein' not in set_a.columns:
        set_a = set_a.copy()
        set_a['protein'] = [f"setA_{i}" for i in range(len(set_a))]
    if 'protein' not in set_b.columns:
        set_b = set_b.copy()
        set_b['protein'] = [f"setB_{i}" for i in range(len(set_b))]

    set_a['Dataset'] = dataset_labels[0]
    set_b['Dataset'] = dataset_labels[1]
    combined = pd.concat([set_a, set_b], ignore_index=True)
    combined = combined.dropna(subset=metrics)

    summary = combined.groupby('Dataset')[metrics].agg(['mean', 'median', 'std', 'min', 'max'])

    melted = combined.melt(
        id_vars=['protein', 'Dataset'],
        value_vars=metrics,
        var_name='RMSD_Type',
        value_name='RMSD_Value'
    )
    melted['RMSD_Type'] = melted['RMSD_Type'].map(metric_display).fillna(melted['RMSD_Type'])

    fig, ax = plt.subplots(figsize=figure_size, dpi=300)
    colors = [OPSIN_COLORS['muted_blue_gray'], OPSIN_COLORS['muted_coral']]
    positions_train = np.arange(1, len(metrics) * 3, 3)
    positions_val = positions_train + 1

    label_fontsize = 18
    annotation_fontsize = max(label_fontsize - 2, 12)
    annotation_color = OPSIN_COLORS.get('gray_7_dark', OPSIN_COLORS.get('black', 'black'))
    annotation_entries = []

    for idx, metric_name in enumerate([metric_display.get(m, m) for m in metrics]):
        subset = melted[melted['RMSD_Type'] == metric_name]
        train_vals = subset[subset['Dataset'] == dataset_labels[0]]['RMSD_Value'].values
        val_vals = subset[subset['Dataset'] == dataset_labels[1]]['RMSD_Value'].values

        bp_train = ax.boxplot([train_vals], positions=[positions_train[idx]], widths=0.6,
                              patch_artist=True, showmeans=True,
                              meanprops=dict(marker='D', markerfacecolor='white', markeredgecolor='black',
                                             markersize=4, markeredgewidth=0.8),
                              medianprops=dict(color='black', linewidth=1.2),
                              boxprops=dict(linewidth=0.8, edgecolor='black'),
                              whiskerprops=dict(linewidth=0.8, color='black'),
                              capprops=dict(linewidth=0.8, color='black'),
                              flierprops=dict(marker='o', markerfacecolor='white', markeredgecolor='black',
                                              markersize=5, markeredgewidth=0.8))

        bp_val = ax.boxplot([val_vals], positions=[positions_val[idx]], widths=0.6,
                             patch_artist=True, showmeans=True,
                             meanprops=dict(marker='D', markerfacecolor='white', markeredgecolor='black',
                                            markersize=4, markeredgewidth=0.8),
                             medianprops=dict(color='black', linewidth=1.2),
                             boxprops=dict(linewidth=0.8, edgecolor='black'),
                             whiskerprops=dict(linewidth=0.8, color='black'),
                             capprops=dict(linewidth=0.8, color='black'),
                             flierprops=dict(marker='o', markerfacecolor='white', markeredgecolor='black',
                                             markersize=5, markeredgewidth=0.8))

        if bp_train['boxes']:
            bp_train['boxes'][0].set_facecolor(colors[0])
            bp_train['boxes'][0].set_alpha(0.85)
        if bp_val['boxes']:
            bp_val['boxes'][0].set_facecolor(colors[1])
            bp_val['boxes'][0].set_alpha(0.85)

        for data_vals, pos in ((train_vals, positions_train[idx]), (val_vals, positions_val[idx])):
            if data_vals.size == 0:
                continue
            mean_val = float(np.mean(data_vals))
            max_val = float(np.max(data_vals))
            annotation_entries.append((pos, mean_val, max_val))

    tick_fontsize = label_fontsize
    ax.set_ylabel('Prediction error (Å)', fontsize=label_fontsize)
    tick_positions = (positions_train + positions_val) / 2
    tick_labels = [metric_display.get(m, m) for m in metrics]
    ax.set_xticks(tick_positions)
    ax.tick_params(axis='x', which='both', length=0, labelsize=tick_fontsize)
    ax.set_xticklabels(tick_labels, fontsize=label_fontsize)
    ax.tick_params(axis='y', labelsize=tick_fontsize)
    if len(positions_val):
        ax.set_xlim(0.5, positions_val[-1] + 1.5)

    ax.grid(False)
    ax.set_axisbelow(True)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(1.0)
    ax.spines['bottom'].set_linewidth(1.0)

    legend_labels = (
        f"{dataset_labels[0]} (N={set_a_df.shape[0]})",
        f"{dataset_labels[1]} (N={set_b_df.shape[0]})"
    )
    legend_elements = [
        Patch(facecolor=colors[0], edgecolor='black', linewidth=0.8, alpha=0.85, label=legend_labels[0]),
        Patch(facecolor=colors[1], edgecolor='black', linewidth=0.8, alpha=0.85, label=legend_labels[1])
    ]
    ax.legend(
        handles=legend_elements,
        loc='upper left',
        fontsize=label_fontsize,
        frameon=False,
        handlelength=1.5,
        borderpad=0.4,
        bbox_to_anchor=(0.02, 0.98),
        borderaxespad=0.0,
        bbox_transform=ax.transAxes
    )

    if annotation_entries:
        base_bottom, base_top = ax.get_ylim()
        span = base_top - base_bottom
        if span <= 0:
            span = 1.0
        offset = 0.04 * span
        top_needed = base_top
        for pos, mean_val, max_val in annotation_entries:
            y = max_val + offset
            ax.text(
                pos,
                y,
                f"μ={mean_val:.2f}",
                ha='center',
                va='bottom',
                fontsize=annotation_fontsize,
                color=annotation_color
            )
            top_needed = max(top_needed, y + 0.05 * span)
        ax.set_ylim(base_bottom, top_needed)

    fig.tight_layout()
    if output_path:
        fig.savefig(output_path, dpi=300, bbox_inches='tight')

    return fig, summary


def plot_conservation_around_x50(residue_table_data, conservation_scores_dict=None, figsize_tuple=(20,10)): # Renamed
    """Sequence logos around X.50 using AMINO_ACID_LOGO_COLORS and conservation scores."""
    if residue_table_data is None or residue_table_data.empty:
        fig, ax = plt.subplots(figsize=(8,5)); ax.text(0.5,0.5,"No data",ha='center',va='center'); return fig

    window_size = 4; fig, axes = plt.subplots(1,7,figsize=figsize_tuple,sharey=True)
    tm_cols_list = [c for c in residue_table_data.columns if str(c).count('.')==1 and str(c).split('.')[0].isdigit() and 1<=int(str(c).split('.')[0])<=7]
    filtered_res_table = residue_table_data[tm_cols_list] # Renamed
    if filtered_res_table.empty:
        for i in range(7): axes[i].text(0.5,0.5,f"No TM data H{i+1}",ha='center',va='center'); plt.tight_layout(); return fig

    cons_scores_local = conservation_scores_dict if conservation_scores_dict is not None else {} # Renamed
    if not cons_scores_local: # Calculate if not provided
        for pos_grn_str in filtered_res_table.columns: # Renamed
            aas = [str(cell)[0].upper() for cell in filtered_res_table[pos_grn_str].dropna() if isinstance(cell,str) and cell and cell[0]!='-' and cell[0].upper() in AMINO_ACID_LOGO_COLORS]
            if aas: counts=Counter(aas); total=sum(counts.values()); cons_scores_local[pos_grn_str] = (max(counts.values())/total if total >0 else 0)
            else: cons_scores_local[pos_grn_str] = 0.0

    for helix_idx_val in range(1,8): # Renamed
        ax = axes[helix_idx_val-1]; pivot_grn_id = f"{helix_idx_val}.50" # Renamed
        helix_grn_cols = sorted([c for c in filtered_res_table.columns if str(c).startswith(f"{helix_idx_val}.")], key=lambda p: float(str(p).split('.')[1])) # Renamed
        actual_pivot_grn_id = pivot_grn_id # Renamed
        try: pivot_idx_val = helix_grn_cols.index(pivot_grn_id) # Renamed
        except ValueError:
            closest = sorted([(abs(float(str(p).split('.')[1])-50.0),p) for p in helix_grn_cols])
            if closest: actual_pivot_grn_id = closest[0][1]; pivot_idx_val = helix_grn_cols.index(actual_pivot_grn_id)
            else: ax.text(0.5,0.5,f"No data H{helix_idx_val}",ha='center',va='center'); continue
        start_plot_idx_val = max(0, pivot_idx_val - window_size) # Renamed
        end_plot_idx_val = min(len(helix_grn_cols), pivot_idx_val + window_size + 1) # Renamed
        window_grn_ids = helix_grn_cols[start_plot_idx_val:end_plot_idx_val] # Renamed

        aa_alphabet_list = list(AMINO_ACID_LOGO_COLORS.keys()); aa_alphabet_list.remove('-') # Renamed
        logo_data_df_plot = pd.DataFrame(0.0, index=range(len(window_grn_ids)), columns=aa_alphabet_list) # Renamed

        for i_win_plot, grn_id_win in enumerate(window_grn_ids): # Renamed
            aas_at_pos_list = [str(cell)[0].upper() for cell in filtered_res_table[grn_id_win].dropna() if isinstance(cell,str) and cell and cell[0]!='-' and cell[0].upper() in aa_alphabet_list]
            if aas_at_pos_list:
                counts_map = Counter(aas_at_pos_list); total_val = sum(counts_map.values()) # Renamed
                for aa_code_val, num_val in counts_map.items(): logo_data_df_plot.loc[i_win_plot, aa_code_val] = num_val/total_val # Renamed

        plot_tick_labels = [int(float(str(p).split('.')[1]) - float(str(actual_pivot_grn_id).split('.')[1])) for p in window_grn_ids] # Renamed

        if not logo_data_df_plot.empty and logo_data_df_plot.sum().sum() > 0:
            try:
                logomaker.Logo(logo_data_df_plot, ax=ax, color_scheme=AMINO_ACID_LOGO_COLORS) # Use our scheme
                ax.set_title(f"Helix {helix_idx_val}", color=HELIX_NUMBER_COLORS.get(helix_idx_val, OPSIN_COLORS['black']), fontweight='bold')
                ax.set_xlabel("Offset from X.50", fontsize=10)
                ax.set_xticks(range(len(window_grn_ids))); ax.set_xticklabels(plot_tick_labels, fontsize=9)
                try:
                    pivot_in_win_idx_val = window_grn_ids.index(actual_pivot_grn_id) # Renamed
                    ax.axvline(x=pivot_in_win_idx_val, color=HELIX_NUMBER_COLORS.get(helix_idx_val, OPSIN_COLORS['gray_5_mid']), linestyle=':', alpha=0.6)
                except ValueError: pass
                y_lims_plot = ax.get_ylim(); text_y_score = y_lims_plot[0]-0.1*(y_lims_plot[1]-y_lims_plot[0]); text_y_aa = y_lims_plot[0]-0.25*(y_lims_plot[1]-y_lims_plot[0])
                for i_plot_val, grn_id_curr in enumerate(window_grn_ids): # Renamed
                    score = cons_scores_local.get(grn_id_curr,0)
                    if score > 0: ax.text(i_plot_val, text_y_score, f"{int(score*100)}%", ha='center',fontsize=8,rotation=90,va='top')
                    if score > 0.4 and not logo_data_df_plot.empty and i_plot_val < len(logo_data_df_plot) and logo_data_df_plot.loc[i_plot_val].sum() > 0:
                        max_aa_code = logo_data_df_plot.loc[i_plot_val].idxmax() # Renamed
                        ax.text(i_plot_val, text_y_aa, max_aa_code, ha='center', fontweight='bold', fontsize=10)
            except Exception as e: ax.text(0.5,0.5,"Logo Error",ha='center',va='center'); print(f"Logo Error H{helix_idx_val}: {e}")
        else: ax.text(0.5,0.5,f"No AA data H{helix_idx_val}",ha='center',va='center')

    for ax_curr_plot in axes: ax_curr_plot.set_ylim(-0.6, np.log2(len(aa_alphabet_list)) + 0.1) # Adjusted ylim
    axes[0].set_ylabel("Information (bits)", fontsize=12)
    plt.suptitle("AA Conservation around X.50 (TM Helices 1-7)", fontsize=14, y=0.98)
    plt.tight_layout(rect=[0,0.03,1,0.95]); return fig
