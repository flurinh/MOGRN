#!/usr/bin/env python3
"""
Postprocess curated GRN table to include all residues.

Logic:
1. Trust curated GRN assignments as ground truth
2. Extend helices to cover all phi/psi-detected helical residues (from helices_extended.json)
3. Add loops and tails for all non-helix residues
4. Validate: all residues present and sequential when sorted by GRN

Key principle: Within each helix, residues MUST be sequential.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
import pandas as pd
import re

# Project setup
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "protos" / "src"))

import protos
protos.set_data_path(str(PROJECT_ROOT / "data"))

from protos.processing.structure import StructureProcessor
from protos.processing.grn.grn_utils import (
    sort_grns_str,
    normalize_grn_format,
    validate_grn_string
)


def load_curated_grn() -> pd.DataFrame:
    """Load curated GRN table and normalize column names."""
    grn_file = PROJECT_ROOT / "opsin_output" / "curated_grn.csv"
    df = pd.read_csv(grn_file, index_col=0, dtype=str)

    # Normalize column names using protos (e.g., 1.4 -> 1.40)
    df.columns = [normalize_grn_format(str(c)) for c in df.columns]
    return df


def load_helices_extended() -> Dict[str, Dict[str, List[int]]]:
    """Load phi/psi-based helix definitions."""
    helix_file = PROJECT_ROOT / "property" / "helices_extended.json"
    with open(helix_file) as f:
        return json.load(f)


def parse_grn_value(val) -> Tuple[Optional[str], Optional[int]]:
    """Parse GRN value like 'A123' -> ('A', 123)."""
    if pd.isna(val) or val == '-':
        return (None, None)
    val_str = str(val)
    if len(val_str) < 2:
        return (None, None)
    try:
        aa = val_str[0]
        res_id = int(val_str[1:])
        return (aa, res_id)
    except:
        return (None, None)


def get_structure_residues(struct_id: str, processor: StructureProcessor) -> Dict[int, str]:
    """Get all protein residues for a structure: {res_id: aa}."""
    try:
        # Try original case first, then lowercase
        df = processor.load_entity(struct_id)
        if df is None:
            df = processor.load_entity(struct_id.lower())
        if df is None:
            return {}

        df = df.reset_index()
        if '_model_0' in struct_id:
            df['auth_chain_id'] = 'A'
        if 'auth_chain_id' in df.columns:
            df = df[df['auth_chain_id'] == 'A']

        atom_col = 'atom_name' if 'atom_name' in df.columns else 'res_atom_name'
        ca_df = df[df[atom_col] == 'CA']

        # Filter out non-protein
        if 'res_name3l' in ca_df.columns:
            ca_df = ca_df[~ca_df['res_name3l'].isin(['RET', 'LIG', 'LYR', 'HOH', 'WAT'])]

        res_map = {}
        for _, row in ca_df.iterrows():
            res_id = int(row['auth_seq_id'])
            aa = row.get('res_name1l', 'X')
            res_map[res_id] = aa

        return res_map
    except Exception as e:
        print(f"  [WARN] Could not load {struct_id}: {e}")
        return {}


def extract_helix_assignments(row: pd.Series, helix_num: int) -> List[Tuple[int, int]]:
    """
    Extract assigned (grn_pos, res_id) pairs for a helix from a row.
    Returns sorted list of (grn_position, residue_id) tuples.
    """
    assignments = []
    helix_prefix = f"{helix_num}."

    for col in row.index:
        col_str = str(col)
        # Match helix columns like "1.50" but not loop columns like "12.001"
        if re.match(rf'^{helix_num}\.(\d+)$', col_str):
            match = re.match(rf'^{helix_num}\.(\d+)$', col_str)
            pos = int(match.group(1))
            aa, res_id = parse_grn_value(row[col])
            if res_id is not None:
                assignments.append((pos, res_id))

    return sorted(assignments, key=lambda x: x[0])


def extend_helix(assignments: List[Tuple[int, int]],
                 helix_bounds: Tuple[int, int],
                 res_map: Dict[int, str]) -> Dict[int, Tuple[int, str]]:
    """
    Extend helix assignments to cover all helical residues.
    Maintains sequential residue order within the helix.

    Returns: Dict mapping grn_pos -> (res_id, aa)
    """
    if not assignments:
        return {}

    h_start, h_end = helix_bounds
    result = {}

    # Add existing assignments
    assigned_res_ids = set()
    for grn_pos, res_id in assignments:
        if res_id in res_map:
            result[grn_pos] = (res_id, res_map[res_id])
            assigned_res_ids.add(res_id)

    # Sort by residue ID
    sorted_by_res = sorted(assignments, key=lambda x: x[1])
    if not sorted_by_res:
        return result

    # Extend N-terminally
    first_grn, first_res = sorted_by_res[0]
    for res_id in range(h_start, first_res):
        if res_id in res_map and res_id not in assigned_res_ids:
            grn_pos = first_grn - (first_res - res_id)
            result[grn_pos] = (res_id, res_map[res_id])
            assigned_res_ids.add(res_id)

    # Extend C-terminally
    last_grn, last_res = sorted_by_res[-1]
    for res_id in range(last_res + 1, h_end + 1):
        if res_id in res_map and res_id not in assigned_res_ids:
            grn_pos = last_grn + (res_id - last_res)
            result[grn_pos] = (res_id, res_map[res_id])
            assigned_res_ids.add(res_id)

    # Fill internal gaps
    for i in range(len(sorted_by_res) - 1):
        grn1, res1 = sorted_by_res[i]
        grn2, res2 = sorted_by_res[i + 1]
        for res_id in range(res1 + 1, res2):
            if res_id in res_map and res_id not in assigned_res_ids:
                grn_pos = grn1 + (res_id - res1)
                result[grn_pos] = (res_id, res_map[res_id])
                assigned_res_ids.add(res_id)

    return result


def assign_flexible_regions(res_map: Dict[int, str],
                            helix_assignments: Dict[int, Set[int]],
                            helix_bounds: Dict[int, Tuple[int, int]]) -> Dict[str, Tuple[int, str]]:
    """
    Assign GRN positions to non-helix residues (loops and tails).
    Returns: Dict mapping grn_col -> (res_id, aa)
    """
    result = {}
    all_assigned = set()
    for res_ids in helix_assignments.values():
        all_assigned.update(res_ids)

    sorted_helices = sorted(helix_bounds.keys())
    if not sorted_helices:
        return result

    first_helix_start = helix_bounds[sorted_helices[0]][0]
    last_helix_end = helix_bounds[sorted_helices[-1]][1]

    for res_id, aa in res_map.items():
        if res_id in all_assigned:
            continue

        # N-terminal
        if res_id < first_helix_start:
            dist = first_helix_start - res_id
            grn_col = f"n.{dist}"
            result[grn_col] = (res_id, aa)
            continue

        # C-terminal
        if res_id > last_helix_end:
            dist = res_id - last_helix_end
            grn_col = f"c.{dist}"
            result[grn_col] = (res_id, aa)
            continue

        # Loop between helices
        for i in range(len(sorted_helices) - 1):
            h_n = sorted_helices[i]
            h_m = sorted_helices[i + 1]
            h_n_end = helix_bounds[h_n][1]
            h_m_start = helix_bounds[h_m][0]

            if h_n_end < res_id < h_m_start:
                dist_to_n = res_id - h_n_end
                dist_to_m = h_m_start - res_id

                if dist_to_n <= dist_to_m:
                    grn_col = f"{h_n}{h_m}.{dist_to_n:03d}"
                else:
                    grn_col = f"{h_m}{h_n}.{dist_to_m:03d}"

                result[grn_col] = (res_id, aa)
                break

    return result


def validate_structure(struct_id: str, row: pd.Series,
                       res_map: Dict[int, str], columns: List[str]) -> List[str]:
    """
    Validate that all residues are present and sequential.
    Returns list of error messages.
    """
    errors = []

    # Extract assigned residues in GRN order
    assigned = []
    for col in columns:
        if col in row.index:
            aa, res_id = parse_grn_value(row[col])
            if res_id is not None:
                assigned.append((col, res_id))

    # Check all residues assigned
    assigned_ids = {a[1] for a in assigned}
    missing = set(res_map.keys()) - assigned_ids
    if missing:
        errors.append(f"Missing residues: {sorted(missing)}")

    # Check sequential order (residue IDs should increase monotonically)
    res_sequence = [a[1] for a in assigned]
    for i in range(len(res_sequence) - 1):
        if res_sequence[i + 1] != res_sequence[i] + 1:
            errors.append(f"Non-sequential: res {res_sequence[i]} -> {res_sequence[i+1]}")
            if len(errors) > 5:  # Limit error count
                errors.append("... (more errors truncated)")
                break

    return errors


def postprocess_curated_grn():
    """Main postprocessing function."""
    print("=" * 70)
    print("POSTPROCESSING CURATED GRN TABLE")
    print("=" * 70)

    # Load data
    print("\n[INFO] Loading data...")
    curated_df = load_curated_grn()
    helices_ext = load_helices_extended()
    processor = StructureProcessor("postprocess")

    print(f"  Curated table: {len(curated_df)} structures x {len(curated_df.columns)} columns")
    print(f"  Helix definitions: {len(helices_ext)} structures")

    # Keep only helix columns (remove existing loops/tails)
    print("\n[INFO] Extracting helix-only columns...")
    helix_cols = []
    for col in curated_df.columns:
        valid, _ = validate_grn_string(str(col))
        if valid and re.match(r'^\d\.\d+$', str(col)):
            helix_cols.append(col)

    helix_only_df = curated_df[helix_cols].copy()
    print(f"  Helix-only table: {len(helix_only_df.columns)} columns")

    # Process each structure
    print("\n[INFO] Processing structures...")
    all_new_columns = set(helix_only_df.columns)  # Start with existing columns
    structure_data = {}
    structures_without_helix_def = []

    for struct_id in curated_df.index:
        helix_def = helices_ext.get(struct_id) or helices_ext.get(struct_id.lower())
        if not helix_def:
            # Keep original curated data for structures without helix definitions
            structures_without_helix_def.append(struct_id)
            row = curated_df.loc[struct_id]
            new_data = {}
            for col in curated_df.columns:
                val = row[col]
                if pd.notna(val) and val != '-':
                    col_normalized = normalize_grn_format(str(col))
                    new_data[col_normalized] = val
                    all_new_columns.add(col_normalized)
            structure_data[struct_id] = new_data
            continue

        res_map = get_structure_residues(struct_id, processor)
        if not res_map:
            continue

        row = helix_only_df.loc[struct_id]
        helix_bounds = {int(h): tuple(bounds) for h, bounds in helix_def.items()}
        helix_assignments = {}
        new_data = {}

        # Process helices
        for helix_num in range(1, 8):
            if helix_num not in helix_bounds:
                continue

            assignments = extract_helix_assignments(row, helix_num)
            extended = extend_helix(assignments, helix_bounds[helix_num], res_map)

            helix_assignments[helix_num] = set()
            for grn_pos, (res_id, aa) in extended.items():
                col = f"{helix_num}.{grn_pos}"
                col = normalize_grn_format(col)
                new_data[col] = f"{aa}{res_id}"
                helix_assignments[helix_num].add(res_id)
                all_new_columns.add(col)

        # Process flexible regions
        flex = assign_flexible_regions(res_map, helix_assignments, helix_bounds)
        for col, (res_id, aa) in flex.items():
            col = normalize_grn_format(col)
            new_data[col] = f"{aa}{res_id}"
            all_new_columns.add(col)

        structure_data[struct_id] = new_data

    structures_extended = len(structure_data) - len(structures_without_helix_def)
    print(f"  Processed {len(structure_data)} structures")
    print(f"    - Extended (with helix definitions): {structures_extended}")
    print(f"    - Preserved (no helix definitions): {len(structures_without_helix_def)}")
    print(f"  Total columns: {len(all_new_columns)}")

    # Build final dataframe with sorted columns
    print("\n[INFO] Building final table...")
    sorted_cols = sort_grns_str(list(all_new_columns))

    final_df = pd.DataFrame(index=curated_df.index, columns=sorted_cols)
    final_df = final_df.fillna('-')

    for struct_id, data in structure_data.items():
        for col, val in data.items():
            if col in final_df.columns:
                final_df.loc[struct_id, col] = val

    print(f"  Final table: {len(final_df)} x {len(final_df.columns)}")

    # Validate
    print("\n[INFO] Validating...")
    errors_by_struct = {}

    for struct_id in structure_data.keys():
        res_map = get_structure_residues(struct_id, processor)
        if not res_map:
            continue

        errors = validate_structure(struct_id, final_df.loc[struct_id], res_map, sorted_cols)
        if errors:
            errors_by_struct[struct_id] = errors

    if errors_by_struct:
        print(f"\n[WARN] {len(errors_by_struct)} structures with errors:")
        for struct_id, errors in list(errors_by_struct.items())[:10]:
            print(f"  {struct_id}:")
            for err in errors[:3]:
                print(f"    - {err}")
    else:
        print("  All structures valid!")

    # Save
    output_file = PROJECT_ROOT / "opsin_output" / "curated_grn_postprocessed.csv"
    final_df.to_csv(output_file)
    print(f"\n[INFO] Saved to: {output_file}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    helix_c = len([c for c in final_df.columns if re.match(r'^\d\.\d+$', str(c))])
    n_c = len([c for c in final_df.columns if str(c).startswith('n.')])
    c_c = len([c for c in final_df.columns if str(c).startswith('c.')])
    loop_c = len([c for c in final_df.columns if re.match(r'^\d\d\.', str(c))])
    print(f"  Helix: {helix_c}, N-term: {n_c}, C-term: {c_c}, Loop: {loop_c}")
    print(f"  Total: {len(final_df.columns)} columns")
    print(f"  Errors: {len(errors_by_struct)} structures")

    return final_df, errors_by_struct


if __name__ == "__main__":
    postprocess_curated_grn()
