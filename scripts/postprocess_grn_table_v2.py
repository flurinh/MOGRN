#!/usr/bin/env python3
"""
Postprocess GRN table - Version 2 (cleaner implementation).

Pipeline:
1. Load curated_grn.csv and helices_extended.json
2. For each structure:
   a. Get helix boundaries
   b. Get X.50 pivot positions from the GRN table
   c. Assign ALL helix residues based on offset from X.50
   d. Assign loop residues (between helices)
   e. Assign tail residues (N-tail before H1, C-tail after H7)
3. Each residue is assigned exactly once
"""

import json
import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import pandas as pd

# Project setup
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "protos" / "src"))

# Import protos
try:
    import protos
    if not hasattr(protos, '_data_path_set'):
        protos.set_data_path(str(PROJECT_ROOT / "data"))
        protos._data_path_set = True
    from protos.processing.structure import StructureProcessor
except Exception as e:
    try:
        from protos.processing.structure import StructureProcessor
    except:
        print(f"[ERROR] Could not import StructureProcessor: {e}")
        sys.exit(1)


def load_helices_extended() -> Dict[str, Dict[str, List[int]]]:
    """Load extended helix definitions."""
    helix_file = PROJECT_ROOT / "property" / "helices_extended.json"
    with open(helix_file) as f:
        return json.load(f)


def get_structure_residues(struct_id: str, processor: StructureProcessor,
                           id_case_map: Optional[Dict[str, str]] = None) -> Dict[int, str]:
    """
    Get all residues for a structure.

    Returns:
        Dict mapping auth_seq_id -> amino acid 1-letter code
    """
    try:
        # Try original casing first, then lowercase
        df = None
        if id_case_map and struct_id in id_case_map:
            df = processor.load_entity(id_case_map[struct_id])
        if df is None:
            df = processor.load_entity(struct_id)
        if df is None:
            return {}

        df = df.reset_index()

        # Handle chain filtering
        if '_model_0' in struct_id or (id_case_map and struct_id in id_case_map and '_model_0' in id_case_map[struct_id]):
            df['auth_chain_id'] = 'A'
        if 'auth_chain_id' in df.columns:
            df = df[df['auth_chain_id'] == 'A']

        # Get CA atoms for residue info
        atom_col = 'res_atom_name' if 'res_atom_name' in df.columns else 'atom_name'
        ca_df = df[df[atom_col] == 'CA'].copy()

        # Build mapping
        res_map = {}
        for _, row in ca_df.iterrows():
            seq_id = row['auth_seq_id']
            aa = row.get('res_name1l', '?')
            res_map[seq_id] = aa

        return res_map

    except Exception as e:
        return {}


def get_pivot_from_grn(struct_id: str, helix_num: int, grn_df: pd.DataFrame) -> Optional[int]:
    """
    Get the X.50 pivot residue number from the GRN table.

    First tries to get X.50 directly. If not available, infers from any
    other assigned position in that helix.

    Returns:
        The residue number at position X.50, or None if not found
    """
    if struct_id not in grn_df.index:
        return None

    row = grn_df.loc[struct_id]

    # First try direct X.50 lookup
    col_50 = f"{helix_num}.50"
    if col_50 in grn_df.columns:
        val = row[col_50]
        if pd.notna(val) and val != '-' and isinstance(val, str) and len(val) > 1:
            try:
                return int(val[1:])
            except ValueError:
                pass

    # If X.50 not available, infer from other positions
    for col in grn_df.columns:
        if not col.startswith(f"{helix_num}."):
            continue
        try:
            pos = int(col.split('.')[1])
        except:
            continue

        val = row[col]
        if pd.isna(val) or val == '-':
            continue
        if not isinstance(val, str) or len(val) < 2:
            continue

        try:
            res_num = int(val[1:])
            # Calculate what residue would be at X.50
            pivot = res_num + (50 - pos)
            return pivot
        except ValueError:
            continue

    return None


def assign_residues(struct_id: str,
                    helix_boundaries: Dict[str, List[int]],
                    pivots: Dict[int, int],
                    residue_map: Dict[int, str],
                    original_7_50: Optional[str] = None) -> Dict[str, str]:
    """
    Assign all residues to GRN positions.

    Each residue is assigned exactly once to one of:
    - n.XX (N-tail)
    - H.XX (helix H, position XX)
    - HK.XX (loop between helix H and K)
    - c.XX (C-tail)

    Special case: 7.50 uses the original value from curated_grn.csv because
    the Schiff base lysine is often modified (LYR, LIG, etc.) and would
    appear as 'X' if read from structure files.

    Returns:
        Dict mapping GRN position -> residue annotation (e.g., "K266")
    """
    assignments = {}
    assigned_residues = set()

    # Get all residue numbers sorted
    all_residues = sorted(residue_map.keys())
    if not all_residues:
        return assignments

    # First, assign helix residues
    for helix_str, (h_start, h_end) in sorted(helix_boundaries.items(), key=lambda x: int(x[0])):
        helix_num = int(helix_str)

        if helix_num not in pivots:
            continue

        pivot = pivots[helix_num]

        for res_num in range(h_start, h_end + 1):
            if res_num not in residue_map:
                continue
            if res_num in assigned_residues:
                continue

            # Calculate GRN position based on offset from pivot (X.50)
            grn_pos = 50 + (res_num - pivot)

            # Skip invalid positions
            if grn_pos < 0 or grn_pos > 99:
                continue

            col_name = f"{helix_num}.{grn_pos}"

            # Special case: use original 7.50 value (Schiff base lysine)
            if helix_num == 7 and grn_pos == 50 and original_7_50 is not None:
                assignments[col_name] = original_7_50
            else:
                aa = residue_map[res_num]
                assignments[col_name] = f"{aa}{res_num}"

            assigned_residues.add(res_num)

    # Get helix start/end for determining loops and tails
    helix_ranges = []
    for helix_str in sorted(helix_boundaries.keys(), key=int):
        h_start, h_end = helix_boundaries[helix_str]
        helix_ranges.append((int(helix_str), h_start, h_end))

    if not helix_ranges:
        return assignments

    # N-tail: residues before first helix
    first_helix_start = helix_ranges[0][1]
    n_tail_residues = [r for r in all_residues if r < first_helix_start and r not in assigned_residues]
    for i, res_num in enumerate(sorted(n_tail_residues, reverse=True)):
        # Number from helix inward: n.01 is closest to helix
        pos = i + 1
        col_name = f"n.{pos:02d}"
        aa = residue_map[res_num]
        assignments[col_name] = f"{aa}{res_num}"
        assigned_residues.add(res_num)

    # Loops: residues between consecutive helices
    for i in range(len(helix_ranges) - 1):
        h1_num, h1_start, h1_end = helix_ranges[i]
        h2_num, h2_start, h2_end = helix_ranges[i + 1]

        loop_residues = [r for r in all_residues
                        if h1_end < r < h2_start and r not in assigned_residues]

        for j, res_num in enumerate(sorted(loop_residues)):
            pos = j + 1
            col_name = f"{h1_num}{h2_num}.{pos:02d}"
            aa = residue_map[res_num]
            assignments[col_name] = f"{aa}{res_num}"
            assigned_residues.add(res_num)

    # C-tail: residues after last helix
    last_helix_end = helix_ranges[-1][2]
    c_tail_residues = [r for r in all_residues if r > last_helix_end and r not in assigned_residues]
    for i, res_num in enumerate(sorted(c_tail_residues)):
        pos = i + 1
        col_name = f"c.{pos:02d}"
        aa = residue_map[res_num]
        assignments[col_name] = f"{aa}{res_num}"
        assigned_residues.add(res_num)

    return assignments


def grn_sort_key(col: str) -> Tuple:
    """Sort key for GRN columns."""
    col_str = str(col)

    # N-terminal: n.XX (comes first, ordered from far to close)
    if col_str.startswith('n.'):
        try:
            dist = int(col_str.split('.')[1])
            return (0, 0, -dist)  # Negative for reverse order
        except:
            return (0, 0, 0)

    # TM helices (1-7): H.XX
    match = re.match(r'^(\d)\.(\d+)$', col_str)
    if match:
        helix = int(match.group(1))
        pos = int(match.group(2))
        return (helix, 0, pos)

    # Loops: HK.XX (between helix H and K)
    match = re.match(r'^(\d)(\d)\.(\d+)$', col_str)
    if match:
        h1 = int(match.group(1))
        h2 = int(match.group(2))
        pos = int(match.group(3))
        # Loops come after helix h1
        return (h1, 1, pos)

    # C-terminal: c.XX (comes last)
    if col_str.startswith('c.'):
        try:
            dist = int(col_str.split('.')[1])
            return (8, 0, dist)
        except:
            return (8, 0, 0)

    # Unknown format
    return (99, 0, 0)


def main():
    """Main function."""
    import argparse

    parser = argparse.ArgumentParser(description="Postprocess GRN table (v2)")
    parser.add_argument(
        "--input", "-i",
        default="opsin_output/curated_grn.csv",
        help="Input GRN table CSV"
    )
    parser.add_argument(
        "--output", "-o",
        default="opsin_output/curated_grn_postprocessed.csv",
        help="Output postprocessed GRN table CSV"
    )
    args = parser.parse_args()

    print("=" * 60)
    print("GRN TABLE POSTPROCESSING (v2)")
    print("=" * 60)

    # Resolve paths
    input_file = PROJECT_ROOT / args.input
    output_file = PROJECT_ROOT / args.output

    # Load extended helix definitions
    helices_ext = load_helices_extended()
    print(f"[INFO] Loaded helix definitions for {len(helices_ext)} structures")

    # Initialize processor
    processor = StructureProcessor("grn_postprocess_v2")

    # Load existing GRN table
    if not input_file.exists():
        print(f"[ERROR] GRN table not found: {input_file}")
        return

    print(f"[INFO] Loading GRN table: {input_file}")
    grn_df = pd.read_csv(input_file, index_col=0)
    print(f"[INFO] Original table: {len(grn_df)} structures x {len(grn_df.columns)} columns")

    # Create case map before normalizing
    id_case_map = {sid.lower(): sid for sid in grn_df.index}

    # Normalize structure IDs to lowercase
    grn_df.index = grn_df.index.str.lower()

    # Process each structure
    print("\n[INFO] Processing structures...")
    all_assignments = {}  # struct_id -> {grn_col -> value}
    all_columns = set()

    processed = 0
    skipped = 0

    for struct_id in grn_df.index:
        if struct_id not in helices_ext:
            skipped += 1
            continue

        # Get helix boundaries
        helix_boundaries = helices_ext[struct_id]

        # Get pivot positions from original GRN table
        pivots = {}
        for helix_num in range(1, 8):
            pivot = get_pivot_from_grn(struct_id, helix_num, grn_df)
            if pivot is not None:
                pivots[helix_num] = pivot

        if not pivots:
            print(f"  [WARN] {struct_id}: No pivot positions found, skipping")
            skipped += 1
            continue

        # Load structure residues
        residue_map = get_structure_residues(struct_id, processor, id_case_map)
        if not residue_map:
            print(f"  [WARN] {struct_id}: Could not load residues, skipping")
            skipped += 1
            continue

        # Get original 7.50 value (Schiff base lysine - often modified as LYR, LIG, etc.)
        original_7_50 = None
        if '7.50' in grn_df.columns:
            val = grn_df.loc[struct_id, '7.50']
            if pd.notna(val) and val != '-':
                original_7_50 = val

        # Assign all residues
        assignments = assign_residues(struct_id, helix_boundaries, pivots, residue_map, original_7_50)

        all_assignments[struct_id] = assignments
        all_columns.update(assignments.keys())
        processed += 1

    print(f"\n[INFO] Processed {processed} structures, skipped {skipped}")

    # Build output DataFrame
    print("\n[INFO] Building output table...")

    # Sort columns
    sorted_columns = sorted(all_columns, key=grn_sort_key)

    # Create DataFrame
    output_df = pd.DataFrame(index=list(all_assignments.keys()), columns=sorted_columns)
    output_df = output_df.fillna('-')

    for struct_id, assignments in all_assignments.items():
        for col, val in assignments.items():
            output_df.loc[struct_id, col] = val

    print(f"[INFO] Output table: {len(output_df)} structures x {len(output_df.columns)} columns")

    # Statistics
    print("\n[INFO] Column statistics:")
    n_tail_cols = [c for c in sorted_columns if c.startswith('n.')]
    c_tail_cols = [c for c in sorted_columns if c.startswith('c.')]
    loop_cols = [c for c in sorted_columns if re.match(r'^\d\d\.', c)]
    helix_cols = [c for c in sorted_columns if re.match(r'^\d\.', c)]

    print(f"  N-tail columns: {len(n_tail_cols)}")
    for h in range(1, 8):
        h_cols = [c for c in helix_cols if c.startswith(f"{h}.")]
        print(f"  Helix {h} columns: {len(h_cols)}")
        if h < 7:
            l_cols = [c for c in loop_cols if c.startswith(f"{h}{h+1}.")]
            print(f"  Loop {h}-{h+1} columns: {len(l_cols)}")
    print(f"  C-tail columns: {len(c_tail_cols)}")

    # Validate: check for duplicates
    print("\n[INFO] Validating assignments...")
    dup_count = 0
    for struct_id in output_df.index:
        residue_nums = []
        for col in output_df.columns:
            val = output_df.loc[struct_id, col]
            if pd.notna(val) and val != '-':
                match = re.search(r'(\d+)', str(val))
                if match:
                    residue_nums.append(int(match.group(1)))

        # Check for duplicates
        from collections import Counter
        counts = Counter(residue_nums)
        dups = [(num, cnt) for num, cnt in counts.items() if cnt > 1]
        if dups:
            dup_count += 1
            if dup_count <= 3:
                print(f"  [WARN] {struct_id}: duplicate residues: {dups[:5]}")

    if dup_count == 0:
        print("[INFO] No duplicate residues found!")
    else:
        print(f"[WARN] {dup_count} structures have duplicate residues")

    # Save
    print(f"\n[INFO] Saving to: {output_file}")
    output_df.to_csv(output_file)

    print("\n" + "=" * 60)
    print("POSTPROCESSING COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
