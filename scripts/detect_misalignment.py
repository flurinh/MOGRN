#!/usr/bin/env python3
"""
Detect misaligned residues in GRN alignment.

Uses the same alignment logic as the interactive visualization.
"""

import pickle
import sys
from pathlib import Path
from typing import Dict, Tuple

import numpy as np
import pandas as pd

# Project setup
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))  # Add project root for 'src' package

from src.visualize_alignment_grn import (
    extract_ca_coordinates_with_grn,
    apply_alignment_transformations,
)


def load_workflow_data(cache_dir: Path, chain_id: str = "A") -> Tuple[Dict, Dict]:
    """Load processed structures and alignment paths from cache."""
    cache_file = cache_dir / f"structure_comparison_{chain_id}.pkl"

    if cache_file.exists():
        print(f"[INFO] Loading from: {cache_file}")
        with open(cache_file, "rb") as f:
            data = pickle.load(f)
        return data.get('processed_structures', {}), data.get('alignment_paths', {})

    print("[ERROR] Cache file not found!")
    return {}, {}


def load_grn_table(grn_file: Path) -> pd.DataFrame:
    """Load GRN table."""
    if grn_file.exists():
        return pd.read_csv(grn_file, index_col=0)
    return pd.DataFrame()


def grn_sort_key(grn):
    """Sort GRN positions: 1.XX, 2.XX, ..., 7.XX"""
    parts = grn.split('.')
    if len(parts) != 2:
        return (99, 0)
    helix, pos = parts
    try:
        helix_num = int(helix)
    except ValueError:
        return (99, 0)
    try:
        pos_num = int(pos)
    except ValueError:
        pos_num = 0
    return (helix_num, pos_num)


def calculate_distances_from_average(
    aligned_structures: Dict,
    grn_df: pd.DataFrame
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Calculate distances from average position for each GRN.
    """
    # Build GRN -> coordinates mapping for each structure
    grn_coords = {}

    # Create case-insensitive lookup for GRN table index
    grn_index_lower = {idx.lower(): idx for idx in grn_df.index}

    for struct_id, data in aligned_structures.items():
        # Get GRN table ID (may be stored as grn_table_id, or use case-insensitive lookup)
        grn_id = data.get('grn_table_id', struct_id)
        if grn_id not in grn_df.index:
            # Try case-insensitive lookup
            if struct_id.lower() in grn_index_lower:
                grn_id = grn_index_lower[struct_id.lower()]
            else:
                continue

        coords = data.get('coords')
        grn_positions = data.get('grn_positions', data.get('grn', []))

        if coords is None or len(grn_positions) == 0:
            continue

        struct_grn_coords = {}
        for i, grn_pos in enumerate(grn_positions):
            if pd.notna(grn_pos) and i < len(coords):
                struct_grn_coords[grn_pos] = coords[i]

        grn_coords[struct_id] = struct_grn_coords

    # Collect all GRN positions
    all_grns = set()
    for struct_coords in grn_coords.values():
        all_grns.update(struct_coords.keys())

    # Filter to TM helices only (1-7)
    all_grns = [g for g in all_grns if '.' in g and g.split('.')[0].isdigit()
                and int(g.split('.')[0]) in range(1, 8)]
    all_grns = sorted(all_grns, key=grn_sort_key)
    all_structs = sorted(grn_coords.keys())

    print(f"[INFO] Processing {len(all_grns)} GRN positions across {len(all_structs)} structures")

    # Initialize distance matrix
    distance_matrix = pd.DataFrame(
        index=all_structs,
        columns=all_grns,
        dtype=float
    )
    distance_matrix[:] = np.nan

    coord_data = []
    stats_data = []

    for grn_pos in all_grns:
        coords_list = []
        struct_ids = []

        for struct_id in all_structs:
            if grn_pos in grn_coords.get(struct_id, {}):
                coords_list.append(grn_coords[struct_id][grn_pos])
                struct_ids.append(struct_id)

        if len(coords_list) < 2:
            continue

        coords_array = np.array(coords_list)
        avg_pos = np.mean(coords_array, axis=0)

        # Calculate distances from average
        distances = np.linalg.norm(coords_array - avg_pos, axis=1)

        for struct_id, dist, coords in zip(struct_ids, distances, coords_list):
            distance_matrix.loc[struct_id, grn_pos] = dist
            coord_data.append({
                'grn': grn_pos,
                'structure': struct_id,
                'distance': dist,
                'x': coords[0],
                'y': coords[1],
                'z': coords[2],
                'avg_x': avg_pos[0],
                'avg_y': avg_pos[1],
                'avg_z': avg_pos[2],
            })

        # Statistics
        mean_dist = np.mean(distances)
        std_dist = np.std(distances)
        max_dist = np.max(distances)
        max_struct = struct_ids[np.argmax(distances)]

        threshold = max(mean_dist + 2 * std_dist, 5.0)
        n_outliers = np.sum(distances > threshold)
        outlier_structs = [s for s, d in zip(struct_ids, distances) if d > threshold]

        stats_data.append({
            'grn': grn_pos,
            'n_structures': len(struct_ids),
            'mean_distance': mean_dist,
            'std_distance': std_dist,
            'max_distance': max_dist,
            'max_structure': max_struct,
            'threshold': threshold,
            'n_outliers': n_outliers,
            'outlier_structures': ';'.join(outlier_structs) if outlier_structs else '',
        })

    return distance_matrix, pd.DataFrame(stats_data), pd.DataFrame(coord_data)


def find_systematic_outliers(distance_matrix: pd.DataFrame, threshold: float = 5.0) -> pd.DataFrame:
    """Find structures that are systematic outliers."""
    outlier_counts = []

    for struct_id in distance_matrix.index:
        distances = distance_matrix.loc[struct_id].dropna()
        n_positions = len(distances)

        if n_positions == 0:
            continue

        n_outliers = (distances > threshold).sum()
        mean_dist = distances.mean()
        max_dist = distances.max()
        max_grn = distances.idxmax() if len(distances) > 0 else ''

        outlier_counts.append({
            'structure': struct_id,
            'n_positions': n_positions,
            'n_outliers': n_outliers,
            'outlier_fraction': n_outliers / n_positions if n_positions > 0 else 0,
            'mean_distance': mean_dist,
            'max_distance': max_dist,
            'max_grn': max_grn,
        })

    df = pd.DataFrame(outlier_counts)
    return df.sort_values('mean_distance', ascending=False)


def main():
    """Main function."""
    print("=" * 60)
    print("MISALIGNMENT DETECTION")
    print("(Using same alignment as interactive visualization)")
    print("=" * 60)

    # Paths
    cache_dir = PROJECT_ROOT / "opsin_output" / "cache"
    grn_file = PROJECT_ROOT / "opsin_output" / "curated_grn_postprocessed.csv"
    output_dir = PROJECT_ROOT / "opsin_output" / "misalignment_analysis"
    output_dir.mkdir(exist_ok=True)

    # Load data
    print("\n[STEP 1] Loading data...")
    processed_structures, alignment_paths = load_workflow_data(cache_dir)
    grn_df = load_grn_table(grn_file)
    print(f"  Structures in cache: {len(processed_structures)}")
    print(f"  Alignment paths: {len(alignment_paths)}")
    print(f"  Structures in GRN table: {len(grn_df)}")

    # Extract coordinates with GRN mapping (same as interactive plot)
    print("\n[STEP 2] Extracting CA coordinates with GRN mapping...")
    structures = extract_ca_coordinates_with_grn(processed_structures, grn_df, chain_id='A', use_helix_only=True)
    print(f"  Extracted {len(structures)} structures with coordinates")

    # Align structures to reference (same as interactive plot)
    print("\n[STEP 3] Aligning structures to reference...")
    reference_id = '6xl3'  # Use the global reference
    aligned_structures = apply_alignment_transformations(structures, alignment_paths, reference_id)
    print(f"  Aligned {len(aligned_structures)} structures")

    # Calculate distances
    print("\n[STEP 4] Calculating distances from average positions...")
    distance_matrix, stats_df, coord_df = calculate_distances_from_average(aligned_structures, grn_df)

    # Find outliers
    print("\n[STEP 5] Finding systematic outliers...")
    outlier_df = find_systematic_outliers(distance_matrix, threshold=5.0)

    # Save results
    print("\n[STEP 6] Saving results...")
    distance_matrix.to_csv(output_dir / "distance_matrix.csv")
    stats_df.to_csv(output_dir / "grn_statistics.csv", index=False)
    coord_df.to_csv(output_dir / "coordinate_details.csv", index=False)
    outlier_df.to_csv(output_dir / "structure_outliers.csv", index=False)
    print(f"  Saved to: {output_dir}")

    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)

    print("\n[GRN positions with most outliers (distance > 5Å)]")
    top_grn = stats_df[stats_df['n_outliers'] > 0].nlargest(15, 'n_outliers')
    if not top_grn.empty:
        print(top_grn[['grn', 'n_structures', 'mean_distance', 'max_distance', 'n_outliers', 'outlier_structures']].to_string())

    print("\n[Structures with highest mean distance from consensus]")
    print(outlier_df.head(20)[['structure', 'n_positions', 'mean_distance', 'max_distance', 'n_outliers', 'max_grn']].to_string())

    return distance_matrix, stats_df, outlier_df


if __name__ == "__main__":
    main()
