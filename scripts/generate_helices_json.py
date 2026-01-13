#!/usr/bin/env python3
"""
Extend helix boundaries using phi/psi angles.

Approach:
1. Load helices_curated.json (correct helix assignments, but truncated boundaries)
2. For each structure, extend each helix boundary using phi/psi angles
3. Extend N-terminal and C-terminal until non-helical residue encountered
4. Output: {structure_id: {H1: [start, end], H2: [start, end], ...}}

This preserves the correct helix numbering from manual curation while
determining the true extent of each helix from phi/psi angles.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import numpy as np
import pandas as pd

# Project setup
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "protos" / "src"))

import protos
protos.set_data_path(str(PROJECT_ROOT / "data"))

from protos.processing.structure import StructureProcessor


# Alpha helix phi/psi ranges (in degrees)
# Core range for reliable helix detection
HELIX_PHI_RANGE = (-80, -40)  # typical: -60
HELIX_PSI_RANGE = (-60, -20)  # typical: -45

# Extended range for boundary detection (slightly permissive at helix ends)
# Made stricter to avoid including disordered residues
HELIX_PHI_RANGE_EXT = (-90, -35)   # was (-100, -30)
HELIX_PSI_RANGE_EXT = (-65, -15)   # was (-70, -10)


def calculate_dihedral(p1, p2, p3, p4):
    """Calculate dihedral angle between 4 points in degrees."""
    b1 = p2 - p1
    b2 = p3 - p2
    b3 = p4 - p3

    n1 = np.cross(b1, b2)
    n2 = np.cross(b2, b3)

    n1_norm = np.linalg.norm(n1)
    n2_norm = np.linalg.norm(n2)

    if n1_norm < 1e-6 or n2_norm < 1e-6:
        return np.nan

    n1 = n1 / n1_norm
    n2 = n2 / n2_norm

    m1 = np.cross(n1, b2 / np.linalg.norm(b2))

    x = np.dot(n1, n2)
    y = np.dot(m1, n2)

    return -np.degrees(np.arctan2(y, x))  # Negative for correct sign convention


def calculate_phi_psi_for_residues(df: pd.DataFrame) -> Dict[int, Tuple[float, float]]:
    """
    Calculate phi/psi angles for all residues from backbone coordinates.
    Returns {auth_seq_id: (phi, psi)}
    """
    atom_col = 'res_atom_name' if 'res_atom_name' in df.columns else 'atom_name'

    # Get backbone atoms
    backbone = df[df[atom_col].isin(['N', 'CA', 'C'])].copy()
    backbone = backbone.sort_values(['auth_seq_id', atom_col])

    # Group by residue
    residues = backbone.groupby('auth_seq_id')

    phi_psi = {}
    seq_ids = sorted(backbone['auth_seq_id'].unique())

    for i, seq_id in enumerate(seq_ids):
        try:
            res = residues.get_group(seq_id)
            n_atom = res[res[atom_col] == 'N'][['x', 'y', 'z']].values
            ca_atom = res[res[atom_col] == 'CA'][['x', 'y', 'z']].values
            c_atom = res[res[atom_col] == 'C'][['x', 'y', 'z']].values

            if len(n_atom) == 0 or len(ca_atom) == 0 or len(c_atom) == 0:
                continue

            n_atom = n_atom[0]
            ca_atom = ca_atom[0]
            c_atom = c_atom[0]

            phi = np.nan
            psi = np.nan

            # Calculate phi (need C from previous residue)
            if i > 0:
                prev_seq_id = seq_ids[i - 1]
                try:
                    prev_res = residues.get_group(prev_seq_id)
                    prev_c = prev_res[prev_res[atom_col] == 'C'][['x', 'y', 'z']].values
                    if len(prev_c) > 0:
                        phi = calculate_dihedral(prev_c[0], n_atom, ca_atom, c_atom)
                except:
                    pass

            # Calculate psi (need N from next residue)
            if i < len(seq_ids) - 1:
                next_seq_id = seq_ids[i + 1]
                try:
                    next_res = residues.get_group(next_seq_id)
                    next_n = next_res[next_res[atom_col] == 'N'][['x', 'y', 'z']].values
                    if len(next_n) > 0:
                        psi = calculate_dihedral(n_atom, ca_atom, c_atom, next_n[0])
                except:
                    pass

            phi_psi[seq_id] = (phi, psi)
        except:
            continue

    return phi_psi


def get_phi_psi_from_protos(df: pd.DataFrame) -> Dict[int, Tuple[float, float]]:
    """Get phi/psi from protos dataframe if available."""
    if 'phi' not in df.columns or 'psi' not in df.columns:
        return {}

    atom_col = 'res_atom_name' if 'res_atom_name' in df.columns else 'atom_name'
    ca_df = df[df[atom_col] == 'CA'].copy()

    phi_psi = {}
    for _, row in ca_df.iterrows():
        seq_id = int(row['auth_seq_id'])
        phi = row.get('phi')
        psi = row.get('psi')
        if pd.notna(phi) and pd.notna(psi):
            phi_psi[seq_id] = (phi, psi)

    return phi_psi


def get_phi_psi(df: pd.DataFrame) -> Dict[int, Tuple[float, float]]:
    """Get phi/psi values, calculating if not available from protos."""
    # Try protos first
    phi_psi = get_phi_psi_from_protos(df)

    # If not available, calculate
    if len(phi_psi) < 10:  # Threshold to detect missing data
        phi_psi = calculate_phi_psi_for_residues(df)

    return phi_psi


def is_helical(phi: float, psi: float, extended: bool = False) -> bool:
    """Check if phi/psi angles indicate alpha helix."""
    if pd.isna(phi) or pd.isna(psi):
        return False

    if extended:
        phi_range = HELIX_PHI_RANGE_EXT
        psi_range = HELIX_PSI_RANGE_EXT
    else:
        phi_range = HELIX_PHI_RANGE
        psi_range = HELIX_PSI_RANGE

    phi_ok = phi_range[0] <= phi <= phi_range[1]
    psi_ok = psi_range[0] <= psi <= psi_range[1]
    return phi_ok and psi_ok


def extend_helix(start: int, end: int,
                 phi_psi: Dict[int, Tuple[float, float]],
                 all_seq_ids: List[int],
                 next_helix_start: Optional[int] = None,
                 prev_helix_end: Optional[int] = None,
                 max_extension: int = 12) -> Tuple[int, int]:
    """
    Extend a helix boundary using phi/psi angles.

    Args:
        start, end: Current helix boundaries
        phi_psi: Dict of {seq_id: (phi, psi)}
        all_seq_ids: All sequence IDs in order
        next_helix_start: Start of next helix (don't extend past this)
        prev_helix_end: End of previous helix (don't extend before this)
        max_extension: Maximum number of residues to extend in each direction

    Returns:
        Extended (start, end) tuple
    """
    seq_set = set(all_seq_ids)

    new_start = start
    new_end = end

    # Extend N-terminal (up to max_extension residues)
    n_extended = 0
    while n_extended < max_extension:
        next_pos = new_start - 1
        if next_pos not in seq_set:
            break
        if prev_helix_end is not None and next_pos <= prev_helix_end:
            break
        if next_pos in phi_psi:
            phi, psi = phi_psi[next_pos]
            if is_helical(phi, psi, extended=True):
                new_start = next_pos
                n_extended += 1
            else:
                break
        else:
            break

    # Extend C-terminal (up to max_extension residues)
    c_extended = 0
    while c_extended < max_extension:
        next_pos = new_end + 1
        if next_pos not in seq_set:
            break
        if next_helix_start is not None and next_pos >= next_helix_start:
            break
        if next_pos in phi_psi:
            phi, psi = phi_psi[next_pos]
            if is_helical(phi, psi, extended=True):
                new_end = next_pos
                c_extended += 1
            else:
                break
        else:
            break

    return new_start, new_end


def resolve_helix_overlaps(helices: Dict[str, List[int]]) -> Dict[str, List[int]]:
    """
    Resolve overlaps between adjacent helices.

    When two helices overlap, split the overlapping region at the midpoint.
    """
    helix_nums = sorted(helices.keys(), key=lambda x: int(x))
    result = {h: list(helices[h]) for h in helix_nums}  # Copy

    for i in range(len(helix_nums) - 1):
        h_n = helix_nums[i]
        h_m = helix_nums[i + 1]

        n_end = result[h_n][1]
        m_start = result[h_m][0]

        if n_end >= m_start:
            # Overlap detected - resolve by splitting at midpoint
            midpoint = (n_end + m_start) // 2
            result[h_n][1] = midpoint
            result[h_m][0] = midpoint + 1

    return result


def process_structure(struct_id: str, df: pd.DataFrame,
                      curated_helices: Dict[str, List[int]]) -> Dict[str, List[int]]:
    """
    Process a single structure: extend helix boundaries using phi/psi.

    Args:
        struct_id: Structure identifier
        df: Structure dataframe
        curated_helices: Curated helix boundaries {helix_num: [start, end]}

    Returns:
        Extended helix boundaries {helix_num: [start, end]}
    """
    atom_col = 'res_atom_name' if 'res_atom_name' in df.columns else 'atom_name'

    # Get all sequence IDs
    ca_df = df[df[atom_col] == 'CA']
    all_seq_ids = sorted(ca_df['auth_seq_id'].unique())

    if not all_seq_ids:
        return curated_helices

    # Get phi/psi values
    phi_psi = get_phi_psi(df)

    if not phi_psi:
        return curated_helices

    # Sort helices by number
    helix_nums = sorted(curated_helices.keys(), key=lambda x: int(x))

    result = {}

    for i, helix_num in enumerate(helix_nums):
        start, end = curated_helices[helix_num]

        # Get boundaries of neighboring helices
        prev_helix_end = None
        next_helix_start = None

        if i > 0:
            prev_num = helix_nums[i - 1]
            prev_helix_end = curated_helices[prev_num][1]

        if i < len(helix_nums) - 1:
            next_num = helix_nums[i + 1]
            next_helix_start = curated_helices[next_num][0]

        # Extend helix
        new_start, new_end = extend_helix(
            start, end, phi_psi, all_seq_ids,
            next_helix_start, prev_helix_end
        )

        result[helix_num] = [new_start, new_end]

    # Resolve any overlaps
    result = resolve_helix_overlaps(result)

    return result


def main():
    """Main function."""
    print("=" * 60)
    print("EXTENDING HELIX BOUNDARIES USING PHI/PSI ANGLES")
    print("=" * 60)

    # Load curated helices
    curated_file = PROJECT_ROOT / "property" / "helices_curated.json"
    print(f"\n[INFO] Loading curated helices from: {curated_file}")

    with open(curated_file) as f:
        curated_helices = json.load(f)

    print(f"[INFO] Loaded {len(curated_helices)} structure definitions")

    # Initialize processor
    processor = StructureProcessor("helix_extend")

    # Get structure IDs from curated_grn.csv - these are the structures we need
    grn_file = PROJECT_ROOT / "opsin_output" / "curated_grn.csv"
    print(f"\n[INFO] Loading structure IDs from: {grn_file}")
    grn_df = pd.read_csv(grn_file, index_col=0)
    # Keep original IDs and create lowercase mapping for helices_curated.json lookup
    original_ids = list(grn_df.index)
    # Map lowercase -> original for loading structures with correct casing
    id_case_map = {sid.lower(): sid for sid in original_ids}
    all_structures = set(id_case_map.keys())  # lowercase for matching
    print(f"[INFO] Found {len(all_structures)} structures in curated_grn.csv")
    print(f"[INFO] Structures with curated helices: {len(curated_helices)}")

    # Process structures
    results = {}
    extended_count = 0
    unchanged_count = 0
    missing_count = 0

    for struct_id_lower in sorted(all_structures):
        if struct_id_lower not in curated_helices:
            print(f"  {struct_id_lower}... MISSING from curated")
            missing_count += 1
            continue

        # Use original casing for loading from protos, fall back to lowercase
        struct_id_original = id_case_map[struct_id_lower]
        print(f"  {struct_id_lower}...", end=" ", flush=True)

        try:
            # Try original casing first (for predictions), then lowercase (for experimental)
            df = processor.load_entity(struct_id_original)
            if df is None:
                df = processor.load_entity(struct_id_lower)
            if df is None:
                print("LOAD FAILED")
                results[struct_id_lower] = curated_helices[struct_id_lower]
                continue

            df = df.reset_index()

            # For predicted structures, set all chains to A
            if '_model_0' in struct_id_lower:
                df['auth_chain_id'] = 'A'

            if 'auth_chain_id' in df.columns:
                df = df[df['auth_chain_id'] == 'A']

            # Add res_atom_name alias
            if 'atom_name' in df.columns and 'res_atom_name' not in df.columns:
                df['res_atom_name'] = df['atom_name']

            # Process
            curated = curated_helices[struct_id_lower]
            extended = process_structure(struct_id_lower, df, curated)
            results[struct_id_lower] = extended

            # Compare
            changes = []
            for h in sorted(curated.keys(), key=int):
                if h in extended:
                    old_start, old_end = curated[h]
                    new_start, new_end = extended[h]
                    if old_start != new_start or old_end != new_end:
                        delta_n = old_start - new_start
                        delta_c = new_end - old_end
                        changes.append(f"H{h}:{delta_n:+d}/{delta_c:+d}")

            if changes:
                print(f"EXTENDED: {', '.join(changes)}")
                extended_count += 1
            else:
                print("OK (no change)")
                unchanged_count += 1

        except Exception as e:
            print(f"ERROR: {e}")
            results[struct_id_lower] = curated_helices[struct_id_lower]

    print(f"\n[INFO] Extended: {extended_count}, Unchanged: {unchanged_count}, Missing: {missing_count}")

    # Show statistics
    print("\n=== Extension Statistics ===")
    total_n_ext = 0
    total_c_ext = 0
    ext_counts = {h: {'n': 0, 'c': 0} for h in ['1', '2', '3', '4', '5', '6', '7']}

    for struct_id in results:
        if struct_id not in curated_helices:
            continue
        curated = curated_helices[struct_id]
        extended = results[struct_id]

        for h in curated:
            if h in extended:
                old_start, old_end = curated[h]
                new_start, new_end = extended[h]
                n_ext = old_start - new_start
                c_ext = new_end - old_end
                total_n_ext += n_ext
                total_c_ext += c_ext
                if n_ext > 0:
                    ext_counts[h]['n'] += 1
                if c_ext > 0:
                    ext_counts[h]['c'] += 1

    print(f"Total N-terminal extensions: {total_n_ext} residues")
    print(f"Total C-terminal extensions: {total_c_ext} residues")
    print("\nExtensions per helix:")
    for h in ['1', '2', '3', '4', '5', '6', '7']:
        print(f"  H{h}: {ext_counts[h]['n']} N-term, {ext_counts[h]['c']} C-term")

    # Save results
    output_file = PROJECT_ROOT / "property" / "helices_extended.json"
    with open(output_file, 'w') as f:
        json.dump(results, f, indent=2, sort_keys=True)
    print(f"\n[INFO] Saved to {output_file}")

    return results


if __name__ == "__main__":
    main()
