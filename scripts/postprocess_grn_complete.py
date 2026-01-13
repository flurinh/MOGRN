#!/usr/bin/env python3
"""
Postprocess curated_grn.csv to include ALL residues.

Logic:
1. Trust curated GRN assignments as ground truth
2. Use helices_curated.json for helix boundaries
3. Extend helices to cover all residues within boundaries
4. Add loops between helices
5. Add N-terminal and C-terminal tails
6. Report: duplicates, missing residues
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
from protos.processing.grn.grn_utils import sort_grns_str, normalize_grn_format


def load_data():
    """Load curated GRN and helix definitions."""
    # Load curated GRN
    grn_file = PROJECT_ROOT / "opsin_output" / "curated_grn.csv"
    grn_df = pd.read_csv(grn_file, index_col=0, dtype=str)

    # Normalize column names
    grn_df.columns = [normalize_grn_format(str(c)) for c in grn_df.columns]

    # Load helix definitions
    helix_file = PROJECT_ROOT / "property" / "helices_curated.json"
    with open(helix_file) as f:
        helices = json.load(f)

    # Build case-insensitive lookup
    helix_lookup = {}
    for k, v in helices.items():
        helix_lookup[k] = v
        helix_lookup[k.lower()] = v

    return grn_df, helix_lookup


def get_structure_residues(struct_id: str, processor: StructureProcessor) -> Dict[int, str]:
    """Get all protein residues: {res_id: aa}"""
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

    if 'res_name3l' in ca_df.columns:
        # Keep LYR (retinal-bound lysine) but filter other non-protein
        ca_df = ca_df[~ca_df['res_name3l'].isin(['RET', 'LIG', 'HOH', 'WAT'])]

    res_map = {}
    for _, row in ca_df.iterrows():
        res_id = int(row['auth_seq_id'])
        aa = row.get('res_name1l', 'X')
        # LYR is lysine bound to retinal - treat as K
        if row.get('res_name3l') == 'LYR':
            aa = 'K'
        res_map[res_id] = aa
    return res_map


def parse_grn_value(val) -> Tuple[Optional[str], Optional[int]]:
    """Parse 'A123' -> ('A', 123)"""
    if pd.isna(val) or val == '-':
        return None, None
    val = str(val)
    if len(val) < 2:
        return None, None
    try:
        return val[0], int(val[1:])
    except:
        return None, None


def get_helix_assignments(row: pd.Series, helix_num: int) -> Dict[int, int]:
    """Extract {grn_position: res_id} for a helix from curated row."""
    assignments = {}
    pattern = re.compile(rf'^{helix_num}\.(\d+)$')

    for col in row.index:
        match = pattern.match(str(col))
        if match:
            pos = int(match.group(1))
            aa, res_id = parse_grn_value(row[col])
            if res_id is not None:
                assignments[pos] = res_id

    return assignments


def extend_helix_assignments(
    assignments: Dict[int, int],
    helix_bounds: Tuple[int, int],
    res_map: Dict[int, str]
) -> Dict[str, Tuple[int, str]]:
    """
    Extend helix to cover all residues within bounds.
    Returns {grn_position_str: (res_id, aa)}

    Uses insertion notation (e.g., "491" for position between 49 and 50)
    when a residue doesn't fit in the standard numbering.
    """
    if not assignments:
        return {}

    h_start, h_end = helix_bounds
    result = {}  # str -> (res_id, aa)

    # Add existing assignments (convert int keys to str)
    assigned_res = set()
    for grn_pos, res_id in assignments.items():
        if res_id in res_map:
            result[str(grn_pos)] = (res_id, res_map[res_id])
            assigned_res.add(res_id)

    if not result:
        return {}

    # Sort by residue ID to find anchor points
    sorted_by_res = sorted(assignments.items(), key=lambda x: x[1])
    first_grn, first_res = sorted_by_res[0]
    last_grn, last_res = sorted_by_res[-1]

    # Extend N-terminally (residues before first assigned, within helix bounds)
    for res_id in range(h_start, first_res):
        if res_id in res_map and res_id not in assigned_res:
            grn_pos = first_grn - (first_res - res_id)
            grn_key = str(grn_pos)
            if grn_key not in result:
                result[grn_key] = (res_id, res_map[res_id])
                assigned_res.add(res_id)

    # Extend C-terminally (residues after last assigned, within helix bounds)
    for res_id in range(last_res + 1, h_end + 1):
        if res_id in res_map and res_id not in assigned_res:
            grn_pos = last_grn + (res_id - last_res)
            grn_key = str(grn_pos)
            if grn_key not in result:
                result[grn_key] = (res_id, res_map[res_id])
                assigned_res.add(res_id)

    # Fill internal gaps - use insertion notation if position is taken
    for i in range(len(sorted_by_res) - 1):
        grn1, res1 = sorted_by_res[i]
        grn2, res2 = sorted_by_res[i + 1]

        missing_residues = [r for r in range(res1 + 1, res2) if r in res_map and r not in assigned_res]
        available_positions = grn2 - grn1 - 1  # positions between grn1 and grn2

        if len(missing_residues) <= available_positions:
            # Normal case: enough positions for missing residues
            for j, res_id in enumerate(missing_residues):
                grn_pos = grn1 + j + 1
                grn_key = str(grn_pos)
                if grn_key not in result:
                    result[grn_key] = (res_id, res_map[res_id])
                    assigned_res.add(res_id)
        else:
            # Insertion case: more residues than positions
            # First fill available positions
            for j in range(available_positions):
                if j < len(missing_residues):
                    res_id = missing_residues[j]
                    grn_pos = grn1 + j + 1
                    grn_key = str(grn_pos)
                    if grn_key not in result:
                        result[grn_key] = (res_id, res_map[res_id])
                        assigned_res.add(res_id)

            # Remaining residues need insertion positions (e.g., 491, 492 for insertions after 49)
            remaining = [r for r in missing_residues if r not in assigned_res]
            for j, res_id in enumerate(remaining):
                # Create insertion position: grn1 + "1", "2", etc.
                insertion_key = f"{grn1}{j+1}"
                result[insertion_key] = (res_id, res_map[res_id])
                assigned_res.add(res_id)

    return result


def assign_loops_and_tails(
    res_map: Dict[int, str],
    helix_assigned: Set[int],
    helix_bounds: Dict[int, Tuple[int, int]]
) -> Dict[str, Tuple[int, str]]:
    """
    Assign GRN positions to non-helix residues.
    Returns {grn_col: (res_id, aa)}
    """
    result = {}
    sorted_helices = sorted(helix_bounds.keys())

    if not sorted_helices:
        return result

    first_helix_start = helix_bounds[sorted_helices[0]][0]
    last_helix_end = helix_bounds[sorted_helices[-1]][1]

    for res_id, aa in res_map.items():
        if res_id in helix_assigned:
            continue

        # N-terminal tail
        if res_id < first_helix_start:
            dist = first_helix_start - res_id
            grn_col = f"n.{dist}"
            result[grn_col] = (res_id, aa)
            continue

        # C-terminal tail
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


def process_structure(
    struct_id: str,
    row: pd.Series,
    helix_bounds_raw: Dict[str, List[int]],
    res_map: Dict[int, str]
) -> Tuple[Dict[str, str], Dict[str, any]]:
    """
    Process a single structure.
    Returns (grn_data, issues)
    """
    # Convert helix bounds to int keys
    helix_bounds = {int(k): tuple(v) for k, v in helix_bounds_raw.items()}

    grn_data = {}
    helix_assigned = set()

    # Step 1: Process each helix
    for helix_num in sorted(helix_bounds.keys()):
        if helix_num not in helix_bounds:
            continue

        # Get existing assignments from curated
        assignments = get_helix_assignments(row, helix_num)

        # Extend to cover all residues in helix bounds
        extended = extend_helix_assignments(
            assignments, helix_bounds[helix_num], res_map
        )

        # Add to results
        for grn_pos, (res_id, aa) in extended.items():
            col = normalize_grn_format(f"{helix_num}.{grn_pos}")
            grn_data[col] = f"{aa}{res_id}"
            helix_assigned.add(res_id)

    # Step 2: Add loops and tails
    flex = assign_loops_and_tails(res_map, helix_assigned, helix_bounds)
    for col, (res_id, aa) in flex.items():
        col = normalize_grn_format(col)
        grn_data[col] = f"{aa}{res_id}"

    # Step 3: Check for issues
    all_assigned = set()
    duplicates = {}
    for col, val in grn_data.items():
        aa, res_id = parse_grn_value(val)
        if res_id is not None:
            if res_id in all_assigned:
                if res_id not in duplicates:
                    duplicates[res_id] = []
                duplicates[res_id].append(col)
            all_assigned.add(res_id)

    missing = set(res_map.keys()) - all_assigned

    issues = {
        'duplicates': duplicates,
        'missing': sorted(missing)
    }

    return grn_data, issues


def main():
    print("=" * 70)
    print("POSTPROCESSING CURATED GRN TABLE")
    print("=" * 70)

    # Load data
    print("\n[INFO] Loading data...")
    grn_df, helix_lookup = load_data()
    processor = StructureProcessor("postprocess")

    print(f"  Curated GRN: {len(grn_df)} structures")

    # Process structures
    print("\n[INFO] Processing structures...")
    all_data = {}
    all_issues = {}
    all_columns = set()

    for struct_id in grn_df.index:
        # Get helix bounds
        h_bounds = helix_lookup.get(struct_id) or helix_lookup.get(struct_id.lower())
        if not h_bounds:
            all_issues[struct_id] = {'error': 'No helix definition'}
            continue

        # Get structure residues
        res_map = get_structure_residues(struct_id, processor)
        if not res_map:
            all_issues[struct_id] = {'error': 'Cannot load structure'}
            continue

        # Process
        row = grn_df.loc[struct_id]
        grn_data, issues = process_structure(struct_id, row, h_bounds, res_map)

        all_data[struct_id] = grn_data
        all_columns.update(grn_data.keys())

        if issues['duplicates'] or issues['missing']:
            all_issues[struct_id] = issues

    print(f"  Processed: {len(all_data)} structures")
    print(f"  Columns: {len(all_columns)}")

    # Build final dataframe
    print("\n[INFO] Building final table...")
    sorted_cols = sort_grns_str(list(all_columns))

    final_df = pd.DataFrame(index=grn_df.index, columns=sorted_cols)
    final_df = final_df.fillna('-')

    for struct_id, data in all_data.items():
        for col, val in data.items():
            if col in final_df.columns:
                final_df.loc[struct_id, col] = val

    # Save
    output_file = PROJECT_ROOT / "opsin_output" / "curated_grn_postprocessed.csv"
    final_df.to_csv(output_file)
    print(f"  Saved to: {output_file}")

    # Report issues
    print("\n" + "=" * 70)
    print("VALIDATION REPORT")
    print("=" * 70)

    # Duplicates
    structs_with_dups = {s: i for s, i in all_issues.items()
                         if 'duplicates' in i and i['duplicates']}
    print(f"\n[DUPLICATES] {len(structs_with_dups)} structures have duplicate residue assignments:")
    for s, i in sorted(structs_with_dups.items()):
        print(f"  {s}:")
        for res_id, cols in i['duplicates'].items():
            print(f"    res {res_id} -> {cols}")

    # Missing
    structs_with_missing = {s: i for s, i in all_issues.items()
                            if 'missing' in i and i['missing']}
    print(f"\n[MISSING] {len(structs_with_missing)} structures have missing residues:")
    for s, i in sorted(structs_with_missing.items()):
        missing = i['missing']
        if len(missing) <= 10:
            print(f"  {s}: {missing}")
        else:
            print(f"  {s}: {len(missing)} residues - {missing[:5]}...{missing[-5:]}")

    # Errors
    structs_with_errors = {s: i for s, i in all_issues.items() if 'error' in i}
    if structs_with_errors:
        print(f"\n[ERRORS] {len(structs_with_errors)} structures failed:")
        for s, i in structs_with_errors.items():
            print(f"  {s}: {i['error']}")

    # Summary
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"  Total structures: {len(grn_df)}")
    print(f"  Successfully processed: {len(all_data)}")
    print(f"  With duplicates: {len(structs_with_dups)}")
    print(f"  With missing residues: {len(structs_with_missing)}")
    print(f"  Errors: {len(structs_with_errors)}")
    print(f"  Output columns: {len(sorted_cols)}")

    return final_df, all_issues


if __name__ == "__main__":
    main()
