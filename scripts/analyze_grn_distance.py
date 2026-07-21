#!/usr/bin/env python3
"""
Calculate distances from GRN-annotated residues to retinal.

This script:
1. Loads structures from the processed structures cache (DataFrame format)
2. Uses GRN assignments from the canonical runtime grn_reference.csv
3. Calculates distances to retinal for each residue:
   - CA distance (alpha carbon to closest retinal atom)
   - Sidechain distance (closest sidechain atom to closest retinal atom)
4. Outputs distance tables with same structure as the GRN table

Output:
  - opsin_output/ca_distance_table_grn.csv
  - opsin_output/distance_table_grn.csv
"""

import pickle
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "src"))
sys.path.insert(0, str(PROJECT_ROOT / "protos" / "src"))


# =============================================================================
# Compatibility class for pickle loading
# =============================================================================

class DatasetCompat:
    """Compatibility class for loading pickled data."""
    def __init__(self, pdb_ids, data=None):
        self.pdb_ids = list(pdb_ids)
        self.data = data if data is not None else pd.DataFrame()

    def __reduce__(self):
        return (self.__class__, (self.pdb_ids, self.data))


# Backbone atom names (to exclude for sidechain distance)
BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'H', 'HA', 'HA2', 'HA3', 'OXT'}


def extract_residue_number(value: str) -> int | None:
    """Extract residue number from a GRN value like 'K266' or 'A123'."""
    if pd.isna(value) or value == '-':
        return None
    match = re.search(r'(-?\d+)', str(value))
    if match:
        return int(match.group(1))
    return None


def get_retinal_coords(df_ret: pd.DataFrame) -> np.ndarray | None:
    """
    Get all atom coordinates of retinal.

    Args:
        df_ret: DataFrame with retinal atoms

    Returns:
        Array of shape (N, 3) with retinal atom coordinates, or None if empty
    """
    if df_ret is None or df_ret.empty:
        return None

    coords = df_ret[['x', 'y', 'z']].values.astype(float)
    return coords if len(coords) > 0 else None


def get_residue_ca_coord(df: pd.DataFrame, res_num: int) -> np.ndarray | None:
    """Get CA coordinate for a specific residue."""
    mask = (df['auth_seq_id'] == res_num) & (df['atom_name'] == 'CA')
    atoms = df[mask]

    if atoms.empty:
        return None

    return atoms[['x', 'y', 'z']].values[0].astype(float)


def get_residue_sidechain_coords(df: pd.DataFrame, res_num: int) -> np.ndarray | None:
    """Get all sidechain atom coordinates for a specific residue."""
    mask = (df['auth_seq_id'] == res_num) & (~df['atom_name'].isin(BACKBONE_ATOMS))
    atoms = df[mask]

    if atoms.empty:
        return None

    return atoms[['x', 'y', 'z']].values.astype(float)


def min_distance(coords1: np.ndarray, coords2: np.ndarray) -> float:
    """Calculate minimum distance between two sets of coordinates."""
    if coords1.ndim == 1:
        coords1 = coords1.reshape(1, 3)

    # Compute pairwise distances
    diff = coords1[:, np.newaxis, :] - coords2[np.newaxis, :, :]
    distances = np.sqrt(np.sum(diff ** 2, axis=2))
    return float(np.min(distances))


def calculate_distances_for_structure(
    df: pd.DataFrame,
    df_ret: pd.DataFrame,
    grn_row: pd.Series,
    grn_columns: list
) -> tuple[dict, dict]:
    """
    Calculate CA and sidechain distances for all GRN positions in a structure.

    Args:
        df: DataFrame with all atoms
        df_ret: DataFrame with retinal atoms
        grn_row: Row from GRN table with residue assignments
        grn_columns: List of GRN column names

    Returns:
        Tuple of (ca_distances, sidechain_distances) dictionaries
    """
    ca_distances = {}
    sidechain_distances = {}

    # Get retinal coordinates
    retinal_coords = get_retinal_coords(df_ret)
    if retinal_coords is None:
        return ca_distances, sidechain_distances

    for grn_col in grn_columns:
        value = grn_row.get(grn_col)
        res_num = extract_residue_number(value)

        if res_num is None:
            continue

        # CA distance
        ca_coord = get_residue_ca_coord(df, res_num)
        if ca_coord is not None:
            ca_dist = min_distance(ca_coord, retinal_coords)
            ca_distances[grn_col] = round(ca_dist, 2)

        # Sidechain distance (closest atom)
        sc_coords = get_residue_sidechain_coords(df, res_num)
        if sc_coords is not None and len(sc_coords) > 0:
            sc_dist = min_distance(sc_coords, retinal_coords)
            sidechain_distances[grn_col] = round(sc_dist, 2)
        elif ca_coord is not None:
            # Fallback to CA for glycine or if no sidechain atoms found
            sidechain_distances[grn_col] = ca_distances.get(grn_col)

    return ca_distances, sidechain_distances


def main():
    print("=" * 60)
    print("GRN DISTANCE ANALYSIS")
    print("=" * 60)

    # Load GRN table
    grn_file = PROJECT_ROOT / "opsin_output" / "grn_reference.csv"
    print(f"\n[INFO] Loading GRN table: {grn_file}")

    grn_df = pd.read_csv(grn_file, index_col=0, dtype={0: str})
    grn_df.index = grn_df.index.astype(str)
    print(f"[INFO] Loaded {len(grn_df)} structures x {len(grn_df.columns)} columns")

    # Get all columns (GRN positions)
    grn_columns = list(grn_df.columns)

    # Load processed structures
    cache_file = PROJECT_ROOT / "opsin_output" / "cache" / "processed_structures_A.pkl"
    print(f"\n[INFO] Loading structures from: {cache_file}")

    if not cache_file.exists():
        print(f"[ERROR] Cache file not found: {cache_file}")
        sys.exit(1)

    with open(cache_file, 'rb') as f:
        cache_data = pickle.load(f)

    processed_structures = cache_data.get('processed_structures', {})
    print(f"[INFO] Loaded {len(processed_structures)} structures")

    # Initialize result DataFrames
    ca_distance_df = pd.DataFrame(index=grn_df.index, columns=grn_columns, dtype=float)
    sidechain_distance_df = pd.DataFrame(index=grn_df.index, columns=grn_columns, dtype=float)

    # Create lookup for processed structures (lowercase keys)
    ps_lookup = {k.lower(): k for k in processed_structures.keys()}

    # Process each structure
    print("\n[INFO] Calculating distances...")

    success_count = 0
    error_count = 0
    no_retinal_count = 0

    for idx, struct_id in enumerate(grn_df.index):
        if (idx + 1) % 20 == 0:
            print(f"[INFO] Processing {idx + 1}/{len(grn_df)}...")

        # Find matching structure in processed_structures
        struct_id_lower = struct_id.lower()

        # Try different ID formats
        ps_key = None

        # Try exact match (lowercase)
        if struct_id_lower in ps_lookup:
            ps_key = ps_lookup[struct_id_lower]
        else:
            # Try removing _model_0 suffix
            base_id = struct_id_lower.replace('_model_0', '')
            if base_id in ps_lookup:
                ps_key = ps_lookup[base_id]
            else:
                # Try extracting PDB ID (first 4 chars if alphanumeric)
                if len(struct_id_lower) >= 4:
                    pdb_prefix = struct_id_lower[:4]
                    if pdb_prefix in ps_lookup:
                        ps_key = ps_lookup[pdb_prefix]

        if ps_key is None:
            error_count += 1
            continue

        entry = processed_structures[ps_key]

        # Get DataFrames
        df = entry.get('df_norm')
        if df is None or (isinstance(df, pd.DataFrame) and df.empty):
            df = entry.get('df')
        df_ret = entry.get('df_ret')

        if df is None or (isinstance(df, pd.DataFrame) and df.empty):
            error_count += 1
            continue

        if df_ret is None or df_ret.empty:
            no_retinal_count += 1
            continue

        # Get GRN assignments for this structure
        grn_row = grn_df.loc[struct_id]

        # Calculate distances
        ca_distances, sc_distances = calculate_distances_for_structure(
            df, df_ret, grn_row, grn_columns
        )

        if not ca_distances:
            no_retinal_count += 1
            continue

        success_count += 1

        # Store results
        for grn_col, dist in ca_distances.items():
            ca_distance_df.loc[struct_id, grn_col] = dist

        for grn_col, dist in sc_distances.items():
            sidechain_distance_df.loc[struct_id, grn_col] = dist

    print(f"\n[INFO] Results:")
    print(f"  Successful: {success_count}")
    print(f"  No retinal found: {no_retinal_count}")
    print(f"  Structure not found: {error_count}")

    # Save results
    ca_output = PROJECT_ROOT / "opsin_output" / "ca_distance_table_grn.csv"
    sc_output = PROJECT_ROOT / "opsin_output" / "distance_table_grn.csv"

    print(f"\n[INFO] Saving CA distance table: {ca_output}")
    ca_distance_df.to_csv(ca_output)

    print(f"[INFO] Saving sidechain distance table: {sc_output}")
    sidechain_distance_df.to_csv(sc_output)

    # Statistics
    print("\n[INFO] Distance statistics (CA):")
    ca_valid = ca_distance_df.notna().sum().sum()
    print(f"  Valid distances: {ca_valid}")
    if ca_valid > 0:
        ca_mean = ca_distance_df.mean().mean()
        ca_min = ca_distance_df.min().min()
        ca_max = ca_distance_df.max().max()
        print(f"  Mean: {ca_mean:.2f} Å")
        print(f"  Range: {ca_min:.2f} - {ca_max:.2f} Å")

    print("\n[INFO] Distance statistics (Sidechain):")
    sc_valid = sidechain_distance_df.notna().sum().sum()
    print(f"  Valid distances: {sc_valid}")
    if sc_valid > 0:
        sc_mean = sidechain_distance_df.mean().mean()
        sc_min = sidechain_distance_df.min().min()
        sc_max = sidechain_distance_df.max().max()
        print(f"  Mean: {sc_mean:.2f} Å")
        print(f"  Range: {sc_min:.2f} - {sc_max:.2f} Å")

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
