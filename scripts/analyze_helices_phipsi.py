#!/usr/bin/env python3
"""
Analyze helix boundaries using phi/psi angles and sequence alignment.

This script:
1. Uses phi/psi angles with sliding window to detect helical regions
2. Uses sequence alignment to assign helix numbers (1-7)
3. Generates a new helices JSON with consistent boundaries

Alpha helix characteristics:
- phi ≈ -60° (range: -80° to -40°)
- psi ≈ -45° (range: -60° to -20°)
"""

import json
import sys
from pathlib import Path
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional

# Project setup
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "protos" / "src"))

import protos
protos.set_data_path(str(PROJECT_ROOT / "data"))

from protos.processing.structure import StructureProcessor


# Alpha helix phi/psi ranges (in degrees)
HELIX_PHI_RANGE = (-80, -40)  # typical: -60
HELIX_PSI_RANGE = (-60, -20)  # typical: -45

# Sliding window parameters
MIN_HELIX_LENGTH = 15  # Minimum consecutive helical residues
WINDOW_SIZE = 5  # Smoothing window for helix detection


def is_helical(phi: float, psi: float) -> bool:
    """Check if phi/psi angles indicate alpha helix."""
    if pd.isna(phi) or pd.isna(psi):
        return False

    phi_ok = HELIX_PHI_RANGE[0] <= phi <= HELIX_PHI_RANGE[1]
    psi_ok = HELIX_PSI_RANGE[0] <= psi <= HELIX_PSI_RANGE[1]

    return phi_ok and psi_ok


def detect_helical_regions(df: pd.DataFrame, window_size: int = WINDOW_SIZE) -> List[bool]:
    """
    Detect helical regions using phi/psi angles with smoothing.

    Returns list of booleans indicating helical state for each residue.
    """
    # Get unique residues in order
    if 'auth_seq_id' not in df.columns:
        return []

    # Get CA atoms with phi/psi
    ca_df = df[df['res_atom_name'] == 'CA'].copy() if 'res_atom_name' in df.columns else df[df['atom_name'] == 'CA'].copy()

    if ca_df.empty or 'phi' not in ca_df.columns or 'psi' not in ca_df.columns:
        return []

    ca_df = ca_df.sort_values('auth_seq_id')

    # Check each residue for helical angles
    helical = []
    for _, row in ca_df.iterrows():
        helical.append(is_helical(row.get('phi'), row.get('psi')))

    if not helical:
        return []

    # Apply smoothing window (majority vote)
    smoothed = []
    for i in range(len(helical)):
        start = max(0, i - window_size // 2)
        end = min(len(helical), i + window_size // 2 + 1)
        window = helical[start:end]
        smoothed.append(sum(window) > len(window) / 2)

    return smoothed, list(ca_df['auth_seq_id'])


def find_helix_segments(helical: List[bool], seq_ids: List[int],
                        min_length: int = MIN_HELIX_LENGTH) -> List[Tuple[int, int]]:
    """
    Find continuous helical segments above minimum length.

    Returns list of (start_seq_id, end_seq_id) tuples.
    """
    segments = []
    in_helix = False
    start_idx = 0

    for i, is_helix in enumerate(helical):
        if is_helix and not in_helix:
            # Start of helix
            in_helix = True
            start_idx = i
        elif not is_helix and in_helix:
            # End of helix
            in_helix = False
            length = i - start_idx
            if length >= min_length:
                segments.append((seq_ids[start_idx], seq_ids[i - 1]))

    # Handle helix extending to end
    if in_helix:
        length = len(helical) - start_idx
        if length >= min_length:
            segments.append((seq_ids[start_idx], seq_ids[-1]))

    return segments


# Known structures with non-7 helix counts (biologically valid)
NON_7TM_STRUCTURES = {
    'S13_Bin138_Proteo_SR': 6,  # Bacterial signal activator with 6 TM helices
    'S13_Bin138_Proteo_SR_model_0': 6,
}


def assign_helix_numbers(segments: List[Tuple[int, int]],
                         ref_helices: Dict[str, List[int]],
                         ref_alignment: Dict[int, int] = None,
                         struct_id: str = None) -> Dict[str, List[int]]:
    """
    Assign helix numbers (1-7) to detected segments based on reference.

    Uses sequence alignment to map detected helices to reference helix numbers.
    Most microbial opsins have 7 TM helices, but some variants have 6.
    """
    expected = NON_7TM_STRUCTURES.get(struct_id, 7)

    if len(segments) != expected:
        if len(segments) < expected:
            print(f"  [WARN] Found {len(segments)} helical segments, expected {expected}")
        elif len(segments) > expected:
            print(f"  [INFO] Found {len(segments)} segments, keeping first {expected}")

    # For now, assign by order (N-to-C terminal)
    # In a full implementation, we'd use sequence alignment
    result = {}
    for i, (start, end) in enumerate(segments[:expected]):
        helix_num = str(i + 1)
        result[helix_num] = [int(start), int(end)]

    return result


def analyze_structure(struct_id: str, processor: StructureProcessor,
                      ref_helices: Dict[str, List[int]] = None) -> Optional[Dict[str, List[int]]]:
    """
    Analyze a single structure to determine helix boundaries.
    """
    try:
        df = processor.load_entity(struct_id)
        if df is None or df.empty:
            return None

        df = df.reset_index()

        # Add res_atom_name alias if needed
        if 'atom_name' in df.columns and 'res_atom_name' not in df.columns:
            df['res_atom_name'] = df['atom_name']

        # Filter to chain A
        if 'auth_chain_id' in df.columns:
            df = df[df['auth_chain_id'] == 'A']

        if df.empty:
            return None

        # Detect helical regions
        result = detect_helical_regions(df)
        if not result or len(result) != 2:
            print(f"  [WARN] Could not detect helical regions for {struct_id}")
            return None

        helical, seq_ids = result

        # Find segments
        segments = find_helix_segments(helical, seq_ids)

        if not segments:
            print(f"  [WARN] No helical segments found for {struct_id}")
            return None

        # Assign helix numbers
        helices = assign_helix_numbers(segments, ref_helices, struct_id=struct_id)

        return helices

    except Exception as e:
        print(f"  [ERROR] Failed to analyze {struct_id}: {e}")
        return None


def compare_with_existing(new_helices: Dict[str, Dict],
                          existing_file: str = "property/helices_curated.json"):
    """Compare new helix assignments with existing ones."""

    existing_path = PROJECT_ROOT / existing_file
    if not existing_path.exists():
        print(f"[WARN] Existing file not found: {existing_file}")
        return

    with open(existing_path) as f:
        existing = json.load(f)

    print("\n=== Comparison with existing helices_curated.json ===\n")

    differences = []
    for struct_id in new_helices:
        if struct_id not in existing:
            continue

        new = new_helices[struct_id]
        old = existing[struct_id]

        for h in range(1, 8):
            h_key = str(h)
            if h_key in new and h_key in old:
                new_start, new_end = new[h_key]
                old_start, old_end = old[h_key]
                new_len = new_end - new_start + 1
                old_len = old_end - old_start + 1

                if abs(new_len - old_len) > 5:
                    differences.append({
                        'struct': struct_id,
                        'helix': h,
                        'old': old[h_key],
                        'new': new[h_key],
                        'old_len': old_len,
                        'new_len': new_len
                    })

    if differences:
        print(f"Found {len(differences)} significant differences (>5 residues):\n")
        for d in differences[:20]:
            print(f"  {d['struct']} H{d['helix']}: "
                  f"old={d['old']} ({d['old_len']}), "
                  f"new={d['new']} ({d['new_len']})")
    else:
        print("No significant differences found.")


def main():
    """Main analysis function."""
    print("=" * 60)
    print("HELIX BOUNDARY ANALYSIS USING PHI/PSI ANGLES")
    print("=" * 60)

    # Initialize processor
    processor = StructureProcessor("helix_analysis")

    # Get list of structures
    datasets = ["mo_exp_A", "mo_exp_B", "mo_pred_exp", "mo_pred_novel"]
    all_structures = []

    for dataset in datasets:
        if processor.dataset_manager.dataset_exists(dataset):
            structures = processor.get_dataset_entities(dataset)
            all_structures.extend(structures)
            print(f"[INFO] Found {len(structures)} structures in {dataset}")

    print(f"\n[INFO] Total structures to analyze: {len(all_structures)}")

    # Load reference helices (4fbz as reference)
    ref_helices = None
    existing_path = PROJECT_ROOT / "property" / "helices_curated.json"
    if existing_path.exists():
        with open(existing_path) as f:
            existing = json.load(f)
        if '4fbz' in existing:
            ref_helices = existing['4fbz']
            print(f"[INFO] Using 4fbz as reference: {ref_helices}")

    # Analyze structures
    new_helices = {}
    analyzed = 0
    failed = 0

    print("\n[INFO] Analyzing structures...")

    for struct_id in all_structures[:20]:  # Start with first 20 for testing
        print(f"  Analyzing {struct_id}...", end=" ")
        result = analyze_structure(struct_id, processor, ref_helices)

        if result:
            new_helices[struct_id] = result
            analyzed += 1
            print(f"OK - {len(result)} helices")
        else:
            failed += 1
            print("FAILED")

    print(f"\n[INFO] Analyzed: {analyzed}, Failed: {failed}")

    # Show sample results
    print("\n=== Sample Results ===\n")
    for struct_id in list(new_helices.keys())[:5]:
        print(f"{struct_id}:")
        for h in range(1, 8):
            h_key = str(h)
            if h_key in new_helices[struct_id]:
                start, end = new_helices[struct_id][h_key]
                length = end - start + 1
                print(f"  Helix {h}: {start}-{end} ({length} residues)")
        print()

    # Compare with existing
    compare_with_existing(new_helices)

    # Save results
    output_file = PROJECT_ROOT / "property" / "helices_phipsi.json"
    with open(output_file, 'w') as f:
        json.dump(new_helices, f, indent=2)
    print(f"\n[INFO] Saved results to {output_file}")

    return new_helices


if __name__ == "__main__":
    main()
