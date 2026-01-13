#!/usr/bin/env python3
"""
Postprocess GRN table to include all helical residues from all structures.

Uses helices_extended.json to identify all helical residues for each structure,
then extends the GRN table columns to ensure complete coverage.

Extension strategy:
- For each helix, find the X.50 position (middle/pivot)
- Extend N-terminally and C-terminally from X.50
- Add new GRN columns as needed to cover all helical residues
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional
import numpy as np
import pandas as pd
import re

# Project setup
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "protos" / "src"))

# Guard protos import to avoid re-initialization when imported as module
StructureProcessor = None
try:
    import protos
    # Only set data path if not already initialized
    if not hasattr(protos, '_data_path_set'):
        protos.set_data_path(str(PROJECT_ROOT / "data"))
        protos._data_path_set = True
    from protos.processing.structure import StructureProcessor
except Exception as e:
    # If protos is already initialized elsewhere, just import what we need
    try:
        from protos.processing.structure import StructureProcessor
    except:
        print(f"[WARNING] Could not import StructureProcessor: {e}")


def load_helices_extended() -> Dict[str, Dict[str, List[int]]]:
    """Load extended helix definitions."""
    helix_file = PROJECT_ROOT / "property" / "helices_extended.json"
    with open(helix_file) as f:
        return json.load(f)


def parse_grn_column(col: str) -> Tuple[Optional[int], Optional[int]]:
    """
    Parse a GRN column name to extract helix number and position.

    Args:
        col: Column name like "1.50", "2.41", etc.

    Returns:
        Tuple of (helix_num, position) or (None, None) if not a helix column
    """
    match = re.match(r'^(\d)\.(\d+)$', str(col))
    if match:
        return int(match.group(1)), int(match.group(2))
    return None, None


def get_structure_residue_map(struct_id: str, processor: StructureProcessor) -> Dict[int, str]:
    """
    Get mapping from auth_seq_id to amino acid for a structure.

    Args:
        struct_id: Structure identifier
        processor: StructureProcessor instance

    Returns:
        Dict mapping auth_seq_id -> amino acid 1-letter code
    """
    try:
        df = processor.load_entity(struct_id)
        if df is None:
            return {}

        df = df.reset_index()

        # Handle chain filtering
        if '_model_0' in struct_id:
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
        print(f"  [WARN] Could not load structure {struct_id}: {e}")
        return {}


def analyze_grn_table(msa_df: pd.DataFrame) -> Dict[str, Dict]:
    """
    Analyze current GRN table to find helix boundaries and pivot positions.

    Args:
        msa_df: MSA table with GRN columns

    Returns:
        Dict with helix info: {helix_num: {'min_pos': N, 'max_pos': M, 'columns': [...]}}
    """
    helix_info = {}

    for col in msa_df.columns:
        helix_num, position = parse_grn_column(col)
        if helix_num is not None:
            if helix_num not in helix_info:
                helix_info[helix_num] = {
                    'min_pos': position,
                    'max_pos': position,
                    'columns': []
                }
            helix_info[helix_num]['min_pos'] = min(helix_info[helix_num]['min_pos'], position)
            helix_info[helix_num]['max_pos'] = max(helix_info[helix_num]['max_pos'], position)
            helix_info[helix_num]['columns'].append((col, position))

    # Sort columns by position
    for helix_num in helix_info:
        helix_info[helix_num]['columns'].sort(key=lambda x: x[1])

    return helix_info


def get_structure_grn_mapping(struct_id: str, msa_df: pd.DataFrame) -> Dict[int, str]:
    """
    Get mapping from auth_seq_id to GRN position for a structure.

    Args:
        struct_id: Structure ID
        msa_df: MSA table

    Returns:
        Dict mapping auth_seq_id -> GRN position string
    """
    if struct_id not in msa_df.index:
        return {}

    row = msa_df.loc[struct_id]
    grn_map = {}

    for col in msa_df.columns:
        val = row[col]
        if pd.isna(val) or val == '-':
            continue

        # Parse value like "A123" to extract sequence ID
        if isinstance(val, str) and len(val) > 1:
            seq_id_str = val[1:]  # Remove amino acid letter
            try:
                seq_id = int(seq_id_str)
                grn_map[seq_id] = col
            except ValueError:
                continue

    return grn_map


def find_helix_pivot(struct_id: str, helix_num: int,
                     msa_df: pd.DataFrame, helix_info: Dict) -> Optional[Tuple[int, int]]:
    """
    Find the pivot point (X.50) for a structure's helix.

    Args:
        struct_id: Structure ID
        helix_num: Helix number (1-7)
        msa_df: MSA table
        helix_info: Helix analysis from analyze_grn_table

    Returns:
        Tuple of (auth_seq_id_at_50, offset_from_actual_50) or None
    """
    if struct_id not in msa_df.index:
        return None

    if helix_num not in helix_info:
        return None

    row = msa_df.loc[struct_id]

    # Find X.50 position
    col_50 = f"{helix_num}.50"
    if col_50 in msa_df.columns:
        val = row[col_50]
        if pd.notna(val) and val != '-' and isinstance(val, str) and len(val) > 1:
            try:
                seq_id = int(val[1:])
                return (seq_id, 0)
            except ValueError:
                pass

    # If X.50 is not present, find closest assigned position
    for col, pos in helix_info[helix_num]['columns']:
        val = row[col]
        if pd.notna(val) and val != '-' and isinstance(val, str) and len(val) > 1:
            try:
                seq_id = int(val[1:])
                offset = 50 - pos
                return (seq_id, offset)
            except ValueError:
                continue

    return None


def extend_grn_table(msa_df: pd.DataFrame,
                     helices_ext: Dict[str, Dict[str, List[int]]],
                     processor: StructureProcessor) -> pd.DataFrame:
    """
    Extend GRN table to include all helical residues.

    Args:
        msa_df: Original MSA table
        helices_ext: Extended helix definitions
        processor: StructureProcessor for loading structures

    Returns:
        Extended MSA table with additional GRN columns
    """
    print("[INFO] Analyzing current GRN table...")
    helix_info = analyze_grn_table(msa_df)

    for h_num, info in helix_info.items():
        print(f"  Helix {h_num}: positions {info['min_pos']}-{info['max_pos']} ({len(info['columns'])} columns)")

    # Track new columns needed per helix
    new_columns_needed = {h: {'n_term': [], 'c_term': []} for h in range(1, 8)}

    # Track values for new columns: {col_name: {struct_id: value}}
    new_column_values = {}

    # Track values to update in existing columns (fill gaps)
    existing_column_updates = {}  # {col_name: {struct_id: value}}

    print("\n[INFO] Processing structures...")
    processed = 0
    extended = 0
    gaps_filled = 0

    for struct_id in msa_df.index:
        if struct_id not in helices_ext:
            continue

        processed += 1
        struct_extended = False
        struct_gaps_filled = False

        # Get residue info for this structure
        res_map = get_structure_residue_map(struct_id, processor)
        if not res_map:
            continue

        struct_helices = helices_ext[struct_id]

        for helix_str, (h_start, h_end) in struct_helices.items():
            helix_num = int(helix_str)

            if helix_num not in helix_info:
                continue

            # Get pivot for this structure's helix
            pivot_result = find_helix_pivot(struct_id, helix_num, msa_df, helix_info)
            if pivot_result is None:
                continue

            pivot_seq_id, offset = pivot_result
            # Actual X.50 seq_id
            actual_50_seq = pivot_seq_id + offset

            # Check if helix extends beyond current GRN range
            current_min = helix_info[helix_num]['min_pos']
            current_max = helix_info[helix_num]['max_pos']

            # Map helix residues to GRN positions
            for seq_id in range(h_start, h_end + 1):
                if seq_id not in res_map:
                    continue

                # Calculate GRN position for this residue
                grn_pos = 50 + (seq_id - actual_50_seq)
                col_name = f"{helix_num}.{grn_pos}"
                aa = res_map[seq_id]
                value = f"{aa}{seq_id}"

                # Check if this position is outside current range
                if grn_pos < current_min:
                    # N-terminal extension needed
                    if grn_pos not in new_columns_needed[helix_num]['n_term']:
                        new_columns_needed[helix_num]['n_term'].append(grn_pos)

                    if col_name not in new_column_values:
                        new_column_values[col_name] = {}

                    new_column_values[col_name][struct_id] = value
                    struct_extended = True

                elif grn_pos > current_max:
                    # C-terminal extension needed
                    if grn_pos not in new_columns_needed[helix_num]['c_term']:
                        new_columns_needed[helix_num]['c_term'].append(grn_pos)

                    if col_name not in new_column_values:
                        new_column_values[col_name] = {}

                    new_column_values[col_name][struct_id] = value
                    struct_extended = True

                else:
                    # Position is within existing range - check if it's a gap
                    if col_name in msa_df.columns:
                        current_val = msa_df.loc[struct_id, col_name]
                        if pd.isna(current_val) or current_val == '-':
                            # This is a gap - fill it
                            if col_name not in existing_column_updates:
                                existing_column_updates[col_name] = {}
                            existing_column_updates[col_name][struct_id] = value
                            struct_gaps_filled = True

        if struct_extended:
            extended += 1
        if struct_gaps_filled:
            gaps_filled += 1

    print(f"\n[INFO] Processed {processed} structures")
    print(f"[INFO] {extended} structures need new columns (extensions)")
    print(f"[INFO] {gaps_filled} structures have gaps to fill in existing columns")

    # Report extensions needed
    total_new_cols = 0
    for h_num in range(1, 8):
        n_term = sorted(new_columns_needed[h_num]['n_term'], reverse=True)
        c_term = sorted(new_columns_needed[h_num]['c_term'])

        if n_term or c_term:
            print(f"  Helix {h_num}:")
            if n_term:
                print(f"    N-terminal: {h_num}.{n_term[-1]} to {h_num}.{n_term[0]} ({len(n_term)} new positions)")
                total_new_cols += len(n_term)
            if c_term:
                print(f"    C-terminal: {h_num}.{c_term[0]} to {h_num}.{c_term[-1]} ({len(c_term)} new positions)")
                total_new_cols += len(c_term)

    print(f"\n[INFO] Total new columns to add: {total_new_cols}")
    print(f"[INFO] Total gap fills in existing columns: {sum(len(v) for v in existing_column_updates.values())}")

    if total_new_cols == 0 and not existing_column_updates:
        print("[INFO] No extensions or gap fills needed, returning original table")
        return msa_df

    # Build extended dataframe
    print("[INFO] Building extended GRN table...")

    # Start with original columns
    extended_df = msa_df.copy()

    # First, fill gaps in existing columns
    gap_fill_count = 0
    for col_name, values in existing_column_updates.items():
        for struct_id, val in values.items():
            extended_df.loc[struct_id, col_name] = val
            gap_fill_count += 1

    print(f"[INFO] Filled {gap_fill_count} gaps in existing columns")

    # Add new columns using concat for better performance
    if new_column_values:
        new_cols_df = pd.DataFrame(index=extended_df.index, columns=list(new_column_values.keys()))
        new_cols_df = new_cols_df.fillna('-')

        for col_name, values in new_column_values.items():
            for struct_id, val in values.items():
                new_cols_df.loc[struct_id, col_name] = val

        extended_df = pd.concat([extended_df, new_cols_df], axis=1)

    # Sort columns in proper GRN order
    extended_df = sort_grn_columns(extended_df)

    print(f"[INFO] Extended table: {len(extended_df.columns)} columns (was {len(msa_df.columns)})")

    return extended_df


def sort_grn_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Sort columns in proper GRN order using protos GRN utilities."""
    try:
        from protos.processing.grn.grn_utils import sort_grns_str
        sorted_cols = sort_grns_str(list(df.columns))
        return df[sorted_cols]
    except ImportError:
        # Fallback to custom sorting if protos not available
        return _sort_grn_columns_fallback(df)


def validate_grn_columns(df: pd.DataFrame) -> Dict[str, Tuple[bool, str]]:
    """
    Validate all GRN column names in the DataFrame.

    Args:
        df: GRN table with column names to validate

    Returns:
        Dict mapping column name -> (is_valid, message)
    """
    try:
        from protos.processing.grn.grn_utils import validate_grn_string

        results = {}
        invalid_count = 0

        for col in df.columns:
            is_valid, message = validate_grn_string(str(col))
            results[col] = (is_valid, message)
            if not is_valid:
                invalid_count += 1

        if invalid_count > 0:
            print(f"[WARN] {invalid_count} invalid GRN column names found:")
            for col, (is_valid, msg) in results.items():
                if not is_valid:
                    print(f"  '{col}': {msg}")
        else:
            print(f"[INFO] All {len(df.columns)} GRN column names are valid")

        return results

    except ImportError:
        print("[WARN] protos GRN validation not available")
        return {}


def validate_grn_values(df: pd.DataFrame) -> Dict[str, List[Tuple[str, str]]]:
    """
    Validate GRN cell values in the DataFrame.

    Checks that values follow expected format: [A-Z][0-9]+ (e.g., K204, A123)
    or are gap markers ('-').

    Args:
        df: GRN table with values to validate

    Returns:
        Dict mapping struct_id -> list of (column, invalid_value) tuples
    """
    invalid_values = {}
    value_pattern = re.compile(r'^[A-Z]\d+$')

    for struct_id in df.index:
        try:
            row = df.loc[struct_id]
            # Handle case where loc returns a DataFrame (duplicate indices)
            if isinstance(row, pd.DataFrame):
                row = row.iloc[0]
            struct_invalid = []

            for col in df.columns:
                val = row[col]
                # Handle Series (shouldn't happen but be safe)
                if isinstance(val, pd.Series):
                    val = val.iloc[0] if len(val) > 0 else None

                if pd.isna(val):
                    continue
                if isinstance(val, str) and val == '-':
                    continue
                if isinstance(val, str) and value_pattern.match(val):
                    continue

                # Invalid value found
                struct_invalid.append((col, str(val)))

            if struct_invalid:
                invalid_values[struct_id] = struct_invalid
        except Exception as e:
            print(f"  [WARN] Error validating {struct_id}: {e}")
            continue

    if invalid_values:
        total_invalid = sum(len(v) for v in invalid_values.values())
        print(f"[WARN] {total_invalid} invalid GRN values found in {len(invalid_values)} structures:")
        for struct_id in list(invalid_values.keys())[:5]:
            for col, val in invalid_values[struct_id][:3]:
                print(f"  {struct_id}[{col}] = '{val}'")
    else:
        print(f"[INFO] All GRN values are valid")

    return invalid_values


def _sort_grn_columns_fallback(df: pd.DataFrame) -> pd.DataFrame:
    """Fallback sorting function if protos is not available."""

    def grn_sort_key(col):
        col_str = str(col)

        # N-terminal: n.XX
        if col_str.startswith('n.'):
            try:
                dist = int(col_str.split('.')[1])
                return (0, 0, -dist)  # Negative for reverse order (n.10 before n.1)
            except:
                return (0, 0, 0)

        # TM helices (1-7)
        match = re.match(r'^(\d)\.(\d+)$', col_str)
        if match:
            helix = int(match.group(1))
            pos = int(match.group(2))
            return (helix, 0, pos)

        # Loops (12.XX, 21.XX, 23.XX, 32.XX, etc.)
        match = re.match(r'^(\d)(\d)\.(\d+)$', col_str)
        if match:
            h1, h2 = int(match.group(1)), int(match.group(2))
            dist = int(match.group(3))
            # Place loops between their flanking helices
            if h1 < h2:
                # h1h2.XX: after helix h1 (e.g., 12.01 after helix 1)
                return (h1 + 0.3, h2, dist)
            else:
                # h2h1.XX: before helix h1 (e.g., 21.01 before helix 2)
                return (h1 - 0.3, h2, -dist)

        # C-terminal: c.XX
        if col_str.startswith('c.'):
            try:
                dist = int(col_str.split('.')[1])
                return (8, 0, dist)
            except:
                return (8, 0, 0)

        # Unknown
        return (9, 0, 0)

    sorted_cols = sorted(df.columns, key=grn_sort_key)
    return df[sorted_cols]


def analyze_coverage(df: pd.DataFrame, helices_ext: Dict[str, Dict[str, List[int]]]) -> Dict:
    """
    Analyze coverage of helical residues in the GRN table.

    Args:
        df: GRN table
        helices_ext: Extended helix definitions

    Returns:
        Coverage statistics
    """
    stats = {
        'total_helical_residues': 0,
        'assigned_residues': 0,
        'by_helix': {}
    }

    for helix_num in range(1, 8):
        stats['by_helix'][helix_num] = {
            'total': 0,
            'assigned': 0
        }

    for struct_id in df.index:
        if struct_id not in helices_ext:
            continue

        row = df.loc[struct_id]
        struct_helices = helices_ext[struct_id]

        for helix_str, (h_start, h_end) in struct_helices.items():
            helix_num = int(helix_str)
            helix_len = h_end - h_start + 1

            stats['total_helical_residues'] += helix_len
            stats['by_helix'][helix_num]['total'] += helix_len

            # Count assigned residues in this helix
            assigned = 0
            for col in df.columns:
                h_num, pos = parse_grn_column(col)
                if h_num == helix_num:
                    val = row[col]
                    if pd.notna(val) and val != '-':
                        assigned += 1

            stats['assigned_residues'] += assigned
            stats['by_helix'][helix_num]['assigned'] += assigned

    return stats


def main():
    """Main function."""
    print("=" * 60)
    print("GRN TABLE POSTPROCESSING - EXTEND HELIX COVERAGE")
    print("=" * 60)

    # Load extended helix definitions
    helices_ext = load_helices_extended()
    print(f"[INFO] Loaded helix definitions for {len(helices_ext)} structures")

    # Initialize processor
    processor = StructureProcessor("grn_postprocess")

    # Load existing GRN table
    grn_dir = PROJECT_ROOT / "opsin_output" / "global_reference_grn"
    msa_file = grn_dir / "msa_table_grn.csv"

    if not msa_file.exists():
        print(f"[ERROR] GRN table not found: {msa_file}")
        return

    print(f"[INFO] Loading GRN table: {msa_file}")
    msa_df = pd.read_csv(msa_file, index_col=0)
    print(f"[INFO] Original table: {len(msa_df)} structures x {len(msa_df.columns)} positions")

    # Analyze current coverage
    print("\n[INFO] Analyzing current coverage...")
    stats_before = analyze_coverage(msa_df, helices_ext)
    coverage_before = stats_before['assigned_residues'] / stats_before['total_helical_residues'] * 100 if stats_before['total_helical_residues'] > 0 else 0
    print(f"[INFO] Current coverage: {stats_before['assigned_residues']}/{stats_before['total_helical_residues']} ({coverage_before:.1f}%)")

    for h_num in range(1, 8):
        h_stats = stats_before['by_helix'][h_num]
        pct = h_stats['assigned'] / h_stats['total'] * 100 if h_stats['total'] > 0 else 0
        print(f"  Helix {h_num}: {h_stats['assigned']}/{h_stats['total']} ({pct:.1f}%)")

    # Extend GRN table
    print("\n" + "=" * 40)
    extended_df = extend_grn_table(msa_df, helices_ext, processor)
    print("=" * 40)

    # Analyze new coverage
    print("\n[INFO] Analyzing extended coverage...")
    stats_after = analyze_coverage(extended_df, helices_ext)
    coverage_after = stats_after['assigned_residues'] / stats_after['total_helical_residues'] * 100 if stats_after['total_helical_residues'] > 0 else 0
    print(f"[INFO] Extended coverage: {stats_after['assigned_residues']}/{stats_after['total_helical_residues']} ({coverage_after:.1f}%)")

    for h_num in range(1, 8):
        h_stats = stats_after['by_helix'][h_num]
        pct = h_stats['assigned'] / h_stats['total'] * 100 if h_stats['total'] > 0 else 0
        print(f"  Helix {h_num}: {h_stats['assigned']}/{h_stats['total']} ({pct:.1f}%)")

    # Calculate improvement
    improvement = coverage_after - coverage_before
    print(f"\n[INFO] Coverage improvement: +{improvement:.1f}%")

    # Anchor Schiff base lysine at 7.50
    print("\n" + "=" * 40)
    anchored_df = anchor_schiff_base_lysine(extended_df)
    print("=" * 40)

    # Squash gaps (fix register shifts far from .50 positions)
    print("\n" + "=" * 40)
    squashed_df = squash_gaps(anchored_df)
    print("=" * 40)

    # Detect missing residues
    missing = detect_missing_residues(squashed_df, helices_ext)

    # Expand to include flexible regions (loops and tails)
    print("\n" + "=" * 40)
    final_df = expand_to_flexible_regions(squashed_df, helices_ext, processor)
    print("=" * 40)

    # Validate GRN columns and values
    print("\n" + "=" * 40)
    print("[INFO] Validating GRN table...")
    validate_grn_columns(final_df)
    validate_grn_values(final_df)
    print("=" * 40)

    # Save extended and squashed table
    output_file = grn_dir / "msa_table_grn_extended.csv"
    final_df.to_csv(output_file)
    print(f"\n[INFO] Saved extended GRN table to: {output_file}")

    # Also process distance table
    dist_file = grn_dir / "distance_table_grn.csv"
    if dist_file.exists():
        print(f"\n[INFO] Processing distance table...")
        dist_df = pd.read_csv(dist_file, index_col=0)

        # Add new columns to distance table (with NaN values)
        new_cols = set(final_df.columns) - set(dist_df.columns)
        for col in new_cols:
            dist_df[col] = np.nan

        # Sort columns
        dist_df = sort_grn_columns(dist_df)

        dist_output = grn_dir / "distance_table_grn_extended.csv"
        dist_df.to_csv(dist_output)
        print(f"[INFO] Saved extended distance table to: {dist_output}")

    return squashed_df


def anchor_schiff_base_lysine(df: pd.DataFrame) -> pd.DataFrame:
    """
    Anchor the Schiff base lysine (K) at position 7.50.

    All microbial opsins require a lysine at position 7.50 for the retinal
    Schiff base linkage. This function identifies structures where the K
    is misaligned (e.g., at 7.49 instead of 7.50) and shifts all helix 7
    residues to correct the alignment.

    Args:
        df: GRN table with residue annotations

    Returns:
        DataFrame with helix 7 anchored on lysine at 7.50
    """
    print("\n[INFO] Anchoring Schiff base lysine at 7.50...")

    result_df = df.copy()
    structures_fixed = []
    structures_non_k = []

    for struct_id in result_df.index:
        row = result_df.loc[struct_id].copy()

        # Get all helix 7 columns with values for this structure
        h7_data = []
        for col in result_df.columns:
            h_num, pos = parse_grn_column(col)
            if h_num == 7:
                val = row[col]
                if pd.notna(val) and val != '-' and isinstance(val, str):
                    aa = val[0]  # First character is amino acid
                    h7_data.append((pos, col, val, aa))

        if not h7_data:
            continue

        # Check what's at 7.50
        current_7_50 = None
        lysine_positions = []
        for pos, col, val, aa in h7_data:
            if pos == 50:
                current_7_50 = (val, aa)
            if aa == 'K':
                lysine_positions.append((pos, col, val))

        # If 7.50 already has K, no change needed
        if current_7_50 and current_7_50[1] == 'K':
            continue

        # Find the lysine closest to position 50
        if not lysine_positions:
            # No lysine in helix 7 - this is unusual
            if current_7_50:
                structures_non_k.append((struct_id, current_7_50[0]))
            continue

        # Find K closest to 50
        closest_k = min(lysine_positions, key=lambda x: abs(x[0] - 50))
        k_pos, k_col, k_val = closest_k

        # Calculate shift needed
        shift = 50 - k_pos  # Positive = shift up, negative = shift down

        if shift == 0:
            continue  # Already at 50

        # Shift all helix 7 residues
        # First, clear all helix 7 values
        for pos, col, val, aa in h7_data:
            row[col] = '-'

        # Then, assign to new positions
        for pos, col, val, aa in h7_data:
            new_pos = pos + shift
            new_col = f"7.{new_pos}"
            if new_col not in result_df.columns:
                result_df[new_col] = '-'
            row[new_col] = val

        result_df.loc[struct_id] = row
        structures_fixed.append((struct_id, shift, k_val))

    print(f"[INFO] Fixed {len(structures_fixed)} structures with K anchor shift:")
    for struct_id, shift, k_val in structures_fixed[:10]:
        direction = "+" if shift > 0 else ""
        print(f"  {struct_id}: shifted {direction}{shift} to put {k_val} at 7.50")
    if len(structures_fixed) > 10:
        print(f"  ... and {len(structures_fixed) - 10} more")

    if structures_non_k:
        print(f"[WARN] {len(structures_non_k)} structures have non-K at 7.50:")
        for struct_id, val in structures_non_k[:5]:
            print(f"  {struct_id}: {val}")

    # Sort columns
    result_df = sort_grn_columns(result_df)

    return result_df


def squash_gaps(df: pd.DataFrame) -> pd.DataFrame:
    """
    Squash gaps in GRN table to fix register shifts far from .50 positions.

    Gap squashing rules:
    - Within ±5 of .50: gaps of size 1 are OK (allow minor alignment variations)
    - Outside ±5 of .50: ALL gaps should be squashed (even size 1)

    "Squashing" means shifting residues to be adjacent to their neighbors,
    eliminating artificial gaps from register shifts.

    Args:
        df: GRN table with residue annotations

    Returns:
        DataFrame with gaps squashed
    """
    print("\n[INFO] Squashing gaps in GRN table...")

    # Make a copy
    result_df = df.copy()

    # Track statistics
    total_gaps_squashed = 0
    structures_modified = set()

    # Process each structure (row)
    for struct_id in result_df.index:
        row = result_df.loc[struct_id].copy()

        # Process each helix (1-7)
        for helix_num in range(1, 8):
            # Get all columns for this helix
            helix_cols = []
            for col in result_df.columns:
                h_num, pos = parse_grn_column(col)
                if h_num == helix_num:
                    helix_cols.append((col, pos))

            if not helix_cols:
                continue

            # Sort by position
            helix_cols.sort(key=lambda x: x[1])

            # Find assigned positions (non-gap) for this structure
            assigned = []
            for col, pos in helix_cols:
                val = row[col]
                if pd.notna(val) and val != '-' and isinstance(val, str):
                    assigned.append((pos, col, val))

            if len(assigned) < 2:
                continue

            # Check for gaps that need squashing
            # Process from .50 outward in both directions
            # First, find the .50 position
            pos_50 = 50

            # Separate N-terminal (< 50) and C-terminal (> 50) positions
            n_term = [(p, c, v) for p, c, v in assigned if p < pos_50]
            c_term = [(p, c, v) for p, c, v in assigned if p > pos_50]
            at_50 = [(p, c, v) for p, c, v in assigned if p == pos_50]

            # Process N-terminal (from 50 downward)
            # Sort descending (closest to 50 first)
            n_term.sort(key=lambda x: -x[0])
            prev_pos = pos_50
            for i, (pos, col, val) in enumerate(n_term):
                gap_size = prev_pos - pos - 1

                if gap_size > 0:
                    dist_from_50 = abs(pos - pos_50)

                    should_squash = False
                    if dist_from_50 > 5:
                        # Outside ±5: always squash (even size 1 gaps)
                        should_squash = True
                    else:
                        # Within ±5: only squash if gap > 1
                        if gap_size > 1:
                            should_squash = True

                    if should_squash:
                        # Calculate new position (adjacent to previous)
                        new_pos = prev_pos - 1
                        new_col = f"{helix_num}.{new_pos}"

                        # Move the value to the new position
                        # Clear old position
                        row[col] = '-'
                        # Set new position (create column if needed)
                        if new_col not in result_df.columns:
                            result_df[new_col] = '-'
                        row[new_col] = val

                        total_gaps_squashed += 1
                        structures_modified.add(struct_id)
                        prev_pos = new_pos
                        continue

                prev_pos = pos

            # Process C-terminal (from 50 upward)
            # Sort ascending (closest to 50 first)
            c_term.sort(key=lambda x: x[0])
            prev_pos = pos_50
            for i, (pos, col, val) in enumerate(c_term):
                gap_size = pos - prev_pos - 1

                if gap_size > 0:
                    dist_from_50 = abs(pos - pos_50)

                    should_squash = False
                    if dist_from_50 > 5:
                        # Outside ±5: always squash (even size 1 gaps)
                        should_squash = True
                    else:
                        # Within ±5: only squash if gap > 1
                        if gap_size > 1:
                            should_squash = True

                    if should_squash:
                        # Calculate new position (adjacent to previous)
                        new_pos = prev_pos + 1
                        new_col = f"{helix_num}.{new_pos}"

                        # Move the value to the new position
                        # Clear old position
                        row[col] = '-'
                        # Set new position
                        if new_col not in result_df.columns:
                            result_df[new_col] = '-'
                        row[new_col] = val

                        total_gaps_squashed += 1
                        structures_modified.add(struct_id)
                        prev_pos = new_pos
                        continue

                prev_pos = pos

        # Update the row
        result_df.loc[struct_id] = row

    print(f"[INFO] Squashed {total_gaps_squashed} gaps in {len(structures_modified)} structures")

    # Sort columns
    result_df = sort_grn_columns(result_df)

    return result_df


def detect_missing_residues(df: pd.DataFrame, helices_ext: Dict[str, Dict[str, List[int]]]) -> Dict[str, List[Tuple[int, int, int]]]:
    """
    Detect residues that are in helix definitions but missing from the GRN table.

    Args:
        df: GRN table
        helices_ext: Extended helix definitions

    Returns:
        Dict mapping struct_id to list of (helix_num, start, end) missing ranges
    """
    print("\n[INFO] Detecting missing residues...")

    missing = {}

    for struct_id in df.index:
        if struct_id not in helices_ext:
            continue

        row = df.loc[struct_id]
        struct_missing = []

        # Get all assigned residue IDs for this structure
        assigned_res_ids = set()
        for col in df.columns:
            val = row[col]
            if pd.notna(val) and val != '-' and isinstance(val, str) and len(val) > 1:
                try:
                    res_id = int(val[1:])
                    assigned_res_ids.add(res_id)
                except ValueError:
                    pass

        # Check each helix
        for helix_str, (h_start, h_end) in helices_ext[struct_id].items():
            helix_num = int(helix_str)
            missing_in_helix = []

            for res_id in range(h_start, h_end + 1):
                if res_id not in assigned_res_ids:
                    missing_in_helix.append(res_id)

            if missing_in_helix:
                # Consolidate into ranges
                ranges = []
                start = missing_in_helix[0]
                end = start
                for res_id in missing_in_helix[1:]:
                    if res_id == end + 1:
                        end = res_id
                    else:
                        ranges.append((helix_num, start, end))
                        start = res_id
                        end = res_id
                ranges.append((helix_num, start, end))
                struct_missing.extend(ranges)

        if struct_missing:
            missing[struct_id] = struct_missing

    # Report
    total_missing = sum(len(v) for v in missing.values())
    print(f"[INFO] Found {total_missing} missing residue ranges in {len(missing)} structures")

    if missing:
        # Show first few examples
        print("[INFO] Examples of missing residues:")
        for struct_id in list(missing.keys())[:5]:
            for helix_num, start, end in missing[struct_id][:3]:
                if start == end:
                    print(f"  {struct_id}: helix {helix_num} residue {start}")
                else:
                    print(f"  {struct_id}: helix {helix_num} residues {start}-{end}")

    return missing


def expand_to_flexible_regions(df: pd.DataFrame, helices_ext: Dict[str, Dict[str, List[int]]],
                                processor: 'StructureProcessor') -> pd.DataFrame:
    """
    Expand GRN table to include non-helical residues (loops and tails).

    Notation:
    - N-terminal tail: n1.XX (counting from helix 1 start backwards)
    - Loop between helix N and N+1:
      - NM.XX if closer to helix N (e.g., 12.01, 12.02... from H1 end)
      - MN.XX if closer to helix M (e.g., 21.01, 21.02... from H2 start)
    - C-terminal tail: 7c.XX (counting from helix 7 end forwards)

    Args:
        df: GRN table with helix residues
        helices_ext: Extended helix definitions
        processor: StructureProcessor for loading structure data

    Returns:
        DataFrame with flexible region annotations added
    """
    print("\n[INFO] Expanding GRN table to include flexible regions (loops and tails)...")

    result_df = df.copy()

    # Track new columns to add
    new_columns_needed = set()

    # First pass: determine what columns we need
    struct_flexible_data = {}  # struct_id -> list of (grn_col, residue_value)

    for struct_id in df.index:
        if struct_id not in helices_ext:
            continue

        # Get helix boundaries for this structure (sorted by helix number)
        helix_bounds = {}
        for helix_str, (h_start, h_end) in helices_ext[struct_id].items():
            helix_num = int(helix_str)
            helix_bounds[helix_num] = (h_start, h_end)

        if not helix_bounds:
            continue

        # Get the sequence data for this structure
        try:
            struct_df = processor.load_entity(struct_id)
            if struct_df is None or (isinstance(struct_df, pd.DataFrame) and struct_df.empty):
                continue

            # Filter to chain A and CA atoms for residue identification
            chain_df = struct_df[struct_df['auth_chain_id'] == 'A'].copy()
            # Handle different column name conventions
            atom_col = 'atom_name' if 'atom_name' in chain_df.columns else 'res_atom_name'
            ca_df = chain_df[chain_df[atom_col] == 'CA'].drop_duplicates(subset=['auth_seq_id'])

            if ca_df.empty:
                continue

            # Get min/max residue IDs in the structure
            all_res_ids = sorted(ca_df['auth_seq_id'].unique())
            if not all_res_ids:
                continue

        except Exception as e:
            continue

        flexible_residues = []
        sorted_helices = sorted(helix_bounds.keys())

        # Get residues already assigned in the GRN table for this structure
        # This is more reliable than using helix boundaries from JSON
        already_assigned_res_ids = set()
        row = df.loc[struct_id]
        for col in df.columns:
            val = row[col]
            if pd.notna(val) and val != '-' and isinstance(val, str) and len(val) > 1:
                try:
                    # Extract residue ID from value like "M13" or "K204" or "S-1"
                    res_str = val[1:]
                    # Handle negative residue IDs (e.g., S-1)
                    if res_str.startswith('-'):
                        assigned_res_id = int(res_str)
                    else:
                        assigned_res_id = int(res_str)
                    already_assigned_res_ids.add(assigned_res_id)
                except ValueError:
                    pass

        # Process each residue in the structure
        for res_id in all_res_ids:
            # Skip if already assigned in GRN table
            if res_id in already_assigned_res_ids:
                continue  # Skip helix residues

            # Get residue info
            res_row = ca_df[ca_df['auth_seq_id'] == res_id]
            if res_row.empty:
                continue
            res_name = res_row['res_name3l'].iloc[0]
            if res_name in ['RET', 'LIG', 'LYR', 'HOH', 'WAT']:
                continue  # Skip non-protein

            # One letter code
            aa_map = {
                'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
                'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
                'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
                'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
            }
            aa = aa_map.get(res_name, 'X')

            # Determine which flexible region this residue belongs to
            grn_col = None

            # Check N-terminal (before first helix)
            first_helix = sorted_helices[0]
            first_helix_start = helix_bounds[first_helix][0]
            if res_id < first_helix_start:
                # N-terminal tail: n.X where X is distance from helix 1 start (no leading zeros)
                distance = first_helix_start - res_id
                grn_col = f"n.{distance}"

            # Check C-terminal (after last helix)
            last_helix = sorted_helices[-1]
            last_helix_end = helix_bounds[last_helix][1]
            if res_id > last_helix_end:
                # C-terminal tail: c.X where X is distance from helix 7 end (no leading zeros)
                distance = res_id - last_helix_end
                grn_col = f"c.{distance}"

            # Check loops between helices
            if grn_col is None:
                for i in range(len(sorted_helices) - 1):
                    h_n = sorted_helices[i]
                    h_m = sorted_helices[i + 1]
                    h_n_end = helix_bounds[h_n][1]
                    h_m_start = helix_bounds[h_m][0]

                    if h_n_end < res_id < h_m_start:
                        # In loop between helix N and M
                        dist_to_n = res_id - h_n_end
                        dist_to_m = h_m_start - res_id

                        if dist_to_n <= dist_to_m:
                            # Closer to helix N - use 3-digit format (e.g., 12.001)
                            grn_col = f"{h_n}{h_m}.{dist_to_n:03d}"
                        else:
                            # Closer to helix M - use 3-digit format (e.g., 21.001)
                            grn_col = f"{h_m}{h_n}.{dist_to_m:03d}"
                        break

            if grn_col:
                new_columns_needed.add(grn_col)
                flexible_residues.append((grn_col, f"{aa}{res_id}"))

        if flexible_residues:
            struct_flexible_data[struct_id] = flexible_residues

    # Add new columns efficiently using pd.concat to avoid fragmentation
    print(f"[INFO] Adding {len(new_columns_needed)} new columns for flexible regions")
    new_cols_to_add = [col for col in new_columns_needed if col not in result_df.columns]
    if new_cols_to_add:
        # Create new DataFrame with all new columns at once
        new_cols_df = pd.DataFrame(
            '-',
            index=result_df.index,
            columns=new_cols_to_add
        )
        result_df = pd.concat([result_df, new_cols_df], axis=1)

    # Second pass: fill in the values
    structures_expanded = 0
    total_residues_added = 0

    for struct_id, flexible_residues in struct_flexible_data.items():
        for grn_col, res_value in flexible_residues:
            result_df.loc[struct_id, grn_col] = res_value
            total_residues_added += 1
        structures_expanded += 1

    print(f"[INFO] Added {total_residues_added} flexible region residues in {structures_expanded} structures")

    # Sort columns
    result_df = sort_grn_columns(result_df)

    # Report statistics
    n_term_cols = [c for c in result_df.columns if c.startswith('n.')]
    c_term_cols = [c for c in result_df.columns if c.startswith('c.')]
    loop_cols = [c for c in result_df.columns if re.match(r'^\d\d\.', c)]

    print(f"[INFO] Flexible region columns:")
    print(f"  N-terminal (n.XX): {len(n_term_cols)} columns")
    print(f"  C-terminal (c.XX): {len(c_term_cols)} columns")
    print(f"  Loops (NM.XX): {len(loop_cols)} columns")

    return result_df


if __name__ == "__main__":
    main()
