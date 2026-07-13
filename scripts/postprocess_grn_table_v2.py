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
from typing import Dict, List, Optional
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
    from protos.processing.grn.grn_utils import sort_grns_str
except Exception as e:
    try:
        from protos.processing.structure import StructureProcessor
        from protos.processing.grn.grn_utils import sort_grns_str
    except:
        print(f"[ERROR] Could not import protos modules: {e}")
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


def get_original_helix6_values(struct_id: str, grn_df: pd.DataFrame) -> tuple[Dict[int, str], bool]:
    """
    Extract original 6.42-6.50 values from the curated GRN table.

    These positions may contain valid gaps that should be preserved.
    We copy them first before any extension to avoid filling in gaps.

    Args:
        struct_id: Structure ID
        grn_df: Original GRN DataFrame

    Returns:
        Tuple of:
        - Dict mapping GRN position (42-50) -> value (e.g., "V181" or "-")
        - Boolean indicating if there's a gap in the range (requires position adjustment)
    """
    if struct_id not in grn_df.index:
        return {}, False

    row = grn_df.loc[struct_id]
    values = {}
    has_gap = False

    for grn_pos in range(42, 51):
        # Handle column naming: "6.42", "6.43", ..., "6.49", "6.5" (or "6.50")
        if grn_pos == 50:
            col = '6.50' if '6.50' in grn_df.columns else ('6.5' if '6.5' in grn_df.columns else None)
        else:
            col = f"6.{grn_pos}"

        if col and col in grn_df.columns:
            val = row[col]
            if pd.notna(val):
                values[grn_pos] = str(val)
                if val == '-':
                    has_gap = True

    return values, has_gap


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

    # First try direct X.50 lookup (handle both "7.50" and "7.5" formats)
    col_50 = f"{helix_num}.50"
    col_50_alt = f"{helix_num}.5"  # Some CSVs use "7.5" instead of "7.50"
    if col_50 in grn_df.columns or col_50_alt in grn_df.columns:
        col_50 = col_50 if col_50 in grn_df.columns else col_50_alt
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
                    original_7_50: Optional[str] = None,
                    original_helix6_values: Optional[Dict[int, str]] = None,
                    helix6_has_gap: bool = False) -> Dict[str, str]:
    """
    Assign all residues to GRN positions.

    Each residue is assigned exactly once to one of:
    - n.XX (N-tail)
    - H.XX (helix H, position XX)
    - HK.XX (loop between helix H and K)
    - c.XX (C-tail)

    Special cases using original values from curated_grn.csv:
    - 7.50: Schiff base lysine (often modified as LYR, LIG, etc.)
    - 6.42-6.50: May contain valid gaps that should be preserved

    The helix 6 positions 42-50 are copied FIRST from the original CSV,
    marking those residue numbers as assigned. For structures with gaps,
    the pivot calculation is adjusted by -1 to account for the shift.

    Returns:
        Dict mapping GRN position -> residue annotation (e.g., "K266")
    """
    assignments = {}
    assigned_residues = set()

    # Pre-populate helix 6 positions 42-50 from original CSV
    # This preserves valid gaps and must happen BEFORE extension
    if original_helix6_values:
        for grn_pos, val in original_helix6_values.items():
            col_name = f"6.{grn_pos:02d}"
            assignments[col_name] = val
            # Mark the residue number as assigned (if not a gap)
            if val != '-':
                match = re.search(r'(-?\d+)', str(val))
                if match:
                    res_num = int(match.group(1))
                    assigned_residues.add(res_num)

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

            # For helix 6 with gaps: the gap causes pivot-based calculation to be
            # off by 1 for ALL residues before/at the pre-populated range boundary.
            # Adjust by -1 for positions at or before position 42 (start of copied range)
            if helix_num == 6 and helix6_has_gap and grn_pos <= 42:
                grn_pos -= 1

            # Skip invalid positions
            if grn_pos < 0 or grn_pos > 99:
                continue

            col_name = f"{helix_num}.{grn_pos:02d}"

            # Skip if this column is already assigned
            if col_name in assignments:
                continue

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

    # N-tail: residues before first assigned helix
    # Find the first helix that has a pivot (was actually assigned)
    first_assigned_helix_idx = 0
    for idx, (h_num, h_start, h_end) in enumerate(helix_ranges):
        if h_num in pivots:
            first_assigned_helix_idx = idx
            break

    # All residues before the first assigned helix go to N-tail
    first_assigned_helix_start = helix_ranges[first_assigned_helix_idx][1]
    n_tail_residues = [r for r in all_residues if r < first_assigned_helix_start and r not in assigned_residues]
    for i, res_num in enumerate(sorted(n_tail_residues, reverse=True)):
        # Number from helix inward: n.01 is closest to helix
        pos = i + 1
        col_name = f"n.{pos:02d}"
        aa = residue_map[res_num]
        assignments[col_name] = f"{aa}{res_num}"
        assigned_residues.add(res_num)

    # Loops: residues between consecutive helices
    # Format: <closer helix><farther helix>.<distance to closer helix:03d>
    # Skip loops where either flanking helix has no pivot (those residues go to tails)
    for i in range(len(helix_ranges) - 1):
        h1_num, h1_start, h1_end = helix_ranges[i]
        h2_num, h2_start, h2_end = helix_ranges[i + 1]

        # Skip if either helix has no pivot - those residues should go to N-tail or C-tail
        if h1_num not in pivots or h2_num not in pivots:
            continue

        loop_residues = [r for r in all_residues
                        if h1_end < r < h2_start and r not in assigned_residues]

        for res_num in sorted(loop_residues):
            # Calculate distance to each flanking helix
            dist_to_h1 = res_num - h1_end
            dist_to_h2 = h2_start - res_num

            if dist_to_h1 <= dist_to_h2:
                # Closer to helix 1
                col_name = f"{h1_num}{h2_num}.{dist_to_h1:03d}"
            else:
                # Closer to helix 2
                col_name = f"{h2_num}{h1_num}.{dist_to_h2:03d}"

            aa = residue_map[res_num]
            assignments[col_name] = f"{aa}{res_num}"
            assigned_residues.add(res_num)

    # C-tail: residues after last assigned helix
    # Find the last helix that has a pivot (was actually assigned)
    last_assigned_helix_idx = len(helix_ranges) - 1
    for idx in range(len(helix_ranges) - 1, -1, -1):
        h_num, h_start, h_end = helix_ranges[idx]
        if h_num in pivots:
            last_assigned_helix_idx = idx
            break

    # All residues after the last assigned helix go to C-tail
    last_assigned_helix_end = helix_ranges[last_assigned_helix_idx][2]
    c_tail_residues = [r for r in all_residues if r > last_assigned_helix_end and r not in assigned_residues]
    for i, res_num in enumerate(sorted(c_tail_residues)):
        pos = i + 1
        col_name = f"c.{pos:02d}"
        aa = residue_map[res_num]
        assignments[col_name] = f"{aa}{res_num}"
        assigned_residues.add(res_num)

    return assignments


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
    # Read with index as string to prevent scientific notation parsing (e.g., "1e12" -> 1000000000000)
    grn_df = pd.read_csv(input_file, index_col=0, dtype={0: str})
    # Ensure index is string type
    grn_df.index = grn_df.index.astype(str)
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
        # Handle both "7.50" and "7.5" column names
        original_7_50 = None
        col_7_50 = '7.50' if '7.50' in grn_df.columns else ('7.5' if '7.5' in grn_df.columns else None)
        if col_7_50:
            val = grn_df.loc[struct_id, col_7_50]
            if pd.notna(val) and val != '-':
                original_7_50 = val

        # Get original helix 6 positions 42-50 (may contain valid gaps)
        # These must be copied BEFORE extension to preserve gaps
        original_helix6_values, helix6_has_gap = get_original_helix6_values(struct_id, grn_df)

        # Assign all residues
        assignments = assign_residues(struct_id, helix_boundaries, pivots, residue_map,
                                      original_7_50, original_helix6_values, helix6_has_gap)

        all_assignments[struct_id] = assignments
        all_columns.update(assignments.keys())
        processed += 1

    print(f"\n[INFO] Processed {processed} structures, skipped {skipped}")

    # Build output DataFrame
    print("\n[INFO] Building output table...")

    # Sort columns using protos GRN utilities
    sorted_columns = sort_grns_str(list(all_columns))

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

    # Validate sequences for gaps
    print("\n[INFO] Checking sequences for gaps...")
    structures_with_gaps = []
    unassigned_residues_count = 0

    for struct_id in output_df.index:
        # Build dict of {res_num: res_name} from the row (what we assigned)
        seq_dict = {}
        for col in output_df.columns:
            val = output_df.loc[struct_id, col]
            if pd.notna(val) and val != '-' and isinstance(val, str) and len(val) >= 2:
                # Parse value like "M1", "K266", "S-1" (negative residue IDs)
                aa = val[0]
                res_str = val[1:]
                try:
                    res_num = int(res_str)
                    seq_dict[res_num] = aa
                except ValueError:
                    continue

        if not seq_dict:
            continue

        # Check for gaps (missing residues between min and max)
        sorted_res_nums = sorted(seq_dict.keys())
        min_res = sorted_res_nums[0]
        max_res = sorted_res_nums[-1]

        # Find gaps (excluding the beginning - residues before min_res are OK to be missing)
        expected_residues = set(range(min_res, max_res + 1))

        # Special case: 8xx8 has no residue 220 (7.49) - this is a known gap, not an error
        if struct_id == '8xx8' and 220 in expected_residues:
            expected_residues.discard(220)
        actual_residues = set(seq_dict.keys())
        missing_residues = expected_residues - actual_residues

        if missing_residues:
            # Check if missing residues are in the original structure or truly missing
            original_residue_map = get_structure_residues(struct_id, processor, id_case_map)
            missing_from_structure = missing_residues - set(original_residue_map.keys())
            present_but_unassigned = missing_residues & set(original_residue_map.keys())

            if present_but_unassigned:
                unassigned_residues_count += len(present_but_unassigned)
                # Debug: show which residues are unassigned and why (with residue names)
                unassigned_with_names = [f"{original_residue_map.get(r, '?')}{r}" for r in sorted(present_but_unassigned)]
                print(f"  [DEBUG] {struct_id}: unassigned residues {unassigned_with_names}")
                # Check where these residues fall relative to helix boundaries
                if struct_id in helices_ext:
                    h_bounds = helices_ext[struct_id]
                    for res in sorted(present_but_unassigned)[:5]:
                        aa = original_residue_map.get(res, '?')
                        location = "unknown"
                        for h_str, (h_start, h_end) in h_bounds.items():
                            if h_start <= res <= h_end:
                                location = f"helix {h_str} ({h_start}-{h_end})"
                                break
                        # Check loops
                        sorted_helices = sorted(h_bounds.items(), key=lambda x: int(x[0]))
                        for i in range(len(sorted_helices) - 1):
                            h1_str, (h1_start, h1_end) = sorted_helices[i]
                            h2_str, (h2_start, h2_end) = sorted_helices[i + 1]
                            if h1_end < res < h2_start:
                                location = f"loop {h1_str}-{h2_str} (between {h1_end} and {h2_start})"
                                break
                        # Check tails
                        first_h = sorted_helices[0][1][0]
                        last_h = sorted_helices[-1][1][1]
                        if res < first_h:
                            location = f"N-tail (before {first_h})"
                        elif res > last_h:
                            location = f"C-tail (after {last_h})"
                        print(f"    res {aa}{res}: {location}")

            # Consolidate into ranges for cleaner output
            missing_sorted = sorted(missing_residues)
            ranges = []
            start = missing_sorted[0]
            end = start
            for res in missing_sorted[1:]:
                if res == end + 1:
                    end = res
                else:
                    ranges.append((start, end))
                    start = res
                    end = res
            ranges.append((start, end))

            # Determine gap type
            if len(missing_from_structure) == len(missing_residues):
                gap_type = "disordered"  # All missing from structure file
            elif len(present_but_unassigned) == len(missing_residues):
                gap_type = "UNASSIGNED"  # Bug: present in structure but not assigned
            else:
                gap_type = "mixed"

            structures_with_gaps.append((struct_id, len(missing_residues), ranges, gap_type,
                                         len(missing_from_structure), len(present_but_unassigned)))

    if structures_with_gaps:
        print(f"[INFO] {len(structures_with_gaps)} structures have gaps in their sequences:")
        for struct_id, gap_count, ranges, gap_type, disordered, unassigned in structures_with_gaps[:20]:
            range_strs = []
            for start, end in ranges[:5]:
                if start == end:
                    range_strs.append(str(start))
                else:
                    range_strs.append(f"{start}-{end}")
            suffix = "..." if len(ranges) > 5 else ""
            if gap_type == "disordered":
                print(f"  {struct_id}: {gap_count} disordered (not in structure) at [{', '.join(range_strs)}{suffix}]")
            elif gap_type == "UNASSIGNED":
                print(f"  [BUG] {struct_id}: {gap_count} UNASSIGNED (in structure but not assigned) at [{', '.join(range_strs)}{suffix}]")
            else:
                print(f"  {struct_id}: {gap_count} gaps ({disordered} disordered, {unassigned} unassigned) at [{', '.join(range_strs)}{suffix}]")
        if len(structures_with_gaps) > 20:
            print(f"  ... and {len(structures_with_gaps) - 20} more structures with gaps")

        if unassigned_residues_count > 0:
            print(f"\n[BUG] {unassigned_residues_count} residues are present in structures but not assigned!")
        else:
            print(f"\n[INFO] All gaps are disordered regions (missing from structure files) - this is expected")
    else:
        print("[INFO] All sequences are complete (no internal gaps)")

    print("\n" + "=" * 60)
    print("POSTPROCESSING COMPLETE")
    print("=" * 60)


if __name__ == "__main__":
    main()
