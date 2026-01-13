#!/usr/bin/env python3
"""
Generate helix boundaries based on GRN positions.

For each helix, uses positions X.41 to X.59 (±9 from X.50) to define
the helix boundaries in terms of residue numbers.

This is the inverse of helix extension - instead of extending based on
structure, we define helix ranges based on a fixed GRN window.

Output: property/helices_grn.json
"""

import json
import re
from pathlib import Path
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parent.parent

# GRN range for each helix: X.41 to X.59 (±9 from X.50)
HELIX_GRN_START = 41
HELIX_GRN_END = 59


def extract_residue_number(value: str) -> int | None:
    """Extract residue number from a GRN value like 'K266' or 'A123'."""
    if pd.isna(value) or value == '-':
        return None

    match = re.search(r'(-?\d+)', str(value))
    if match:
        return int(match.group(1))
    return None


def get_helix_range(row: pd.Series, helix_num: int, columns: list) -> tuple[int, int] | None:
    """
    Get the residue number range for a helix based on GRN positions.

    Args:
        row: DataFrame row for a structure
        helix_num: Helix number (1-7)
        columns: List of column names in the DataFrame

    Returns:
        Tuple of (start_residue, end_residue) or None if not found
    """
    start_col = f"{helix_num}.{HELIX_GRN_START}"
    end_col = f"{helix_num}.{HELIX_GRN_END}"

    start_res = None
    end_res = None

    # Find the first assigned residue in the range (for start)
    for pos in range(HELIX_GRN_START, HELIX_GRN_END + 1):
        col = f"{helix_num}.{pos}"
        if col in columns:
            res_num = extract_residue_number(row[col])
            if res_num is not None:
                start_res = res_num
                break

    # Find the last assigned residue in the range (for end)
    for pos in range(HELIX_GRN_END, HELIX_GRN_START - 1, -1):
        col = f"{helix_num}.{pos}"
        if col in columns:
            res_num = extract_residue_number(row[col])
            if res_num is not None:
                end_res = res_num
                break

    if start_res is not None and end_res is not None:
        return (min(start_res, end_res), max(start_res, end_res))

    return None


def main():
    print("=" * 60)
    print("GENERATING HELIX BOUNDARIES FROM GRN POSITIONS")
    print("=" * 60)
    print(f"Using GRN range: X.{HELIX_GRN_START} to X.{HELIX_GRN_END} (±9 from X.50)")

    # Load postprocessed GRN table
    grn_file = PROJECT_ROOT / "opsin_output" / "curated_grn_postprocessed.csv"
    print(f"\n[INFO] Loading GRN table: {grn_file}")

    df = pd.read_csv(grn_file, index_col=0)
    print(f"[INFO] Loaded {len(df)} structures x {len(df.columns)} columns")

    columns = list(df.columns)

    # Extract helix ranges for each structure
    helices_grn = {}

    print("\n[INFO] Extracting helix ranges...")

    for struct_id in df.index:
        row = df.loc[struct_id]
        struct_helices = {}

        for helix_num in range(1, 8):
            helix_range = get_helix_range(row, helix_num, columns)
            if helix_range is not None:
                struct_helices[str(helix_num)] = list(helix_range)

        if struct_helices:
            helices_grn[struct_id] = struct_helices

    print(f"[INFO] Extracted helix ranges for {len(helices_grn)} structures")

    # Statistics
    print("\n[INFO] Helix coverage statistics:")
    for helix_num in range(1, 8):
        count = sum(1 for s in helices_grn.values() if str(helix_num) in s)

        # Calculate average helix length
        lengths = []
        for s in helices_grn.values():
            if str(helix_num) in s:
                start, end = s[str(helix_num)]
                lengths.append(end - start + 1)

        avg_len = sum(lengths) / len(lengths) if lengths else 0
        print(f"  Helix {helix_num}: {count} structures, avg length: {avg_len:.1f} residues")

    # Save to JSON
    output_file = PROJECT_ROOT / "property" / "helices_grn.json"
    print(f"\n[INFO] Saving to: {output_file}")

    with open(output_file, 'w') as f:
        json.dump(helices_grn, f, indent=2)

    print("\n" + "=" * 60)
    print("COMPLETE")
    print("=" * 60)

    # Show a sample
    print("\n[INFO] Sample output (first 3 structures):")
    for struct_id in list(helices_grn.keys())[:3]:
        print(f"  {struct_id}:")
        for h, (start, end) in helices_grn[struct_id].items():
            print(f"    Helix {h}: {start}-{end} ({end-start+1} residues)")


if __name__ == "__main__":
    main()
