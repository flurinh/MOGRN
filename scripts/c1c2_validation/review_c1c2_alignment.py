#!/usr/bin/env python3
"""Inspect sequence and structure consistency for the C1C2 entry."""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple
import pandas as pd
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.Polypeptide import is_aa
from Bio.SeqUtils import seq1
from Bio.PDB.Superimposer import Superimposer
from Bio.PDB.Atom import Atom
from Bio.PDB.Residue import Residue
from Bio import pairwise2

DEFAULT_PROPERTY_FILE = Path("property/mo_exp.csv")
DEFAULT_CHAIN_ID = "A"
DEFAULT_SEARCH_DIRS = (
    Path("data/structure/mmcif"),
    Path("structures/mo_pred"),
    Path("structures/mo_exp"),
    Path("structures/hideaki_exp"),
    Path("structures/hideaki_pred"),
)


@dataclass
class StructureInfo:
    path: Path
    sequence: str
    residues: List[Residue]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Review the C1C2 experimental vs predicted structure alignment.")
    parser.add_argument("--property-file", type=Path, default=DEFAULT_PROPERTY_FILE,
                        help="Path to the property CSV (default: property/mo_exp.csv).")
    parser.add_argument("--target", default="C1C2",
                        help="Target identifier to look up (matches against name, opsin name, short_name, or pdb_id).")
    parser.add_argument("--chain", default=DEFAULT_CHAIN_ID,
                        help="Chain identifier to inspect (default: A).")
    parser.add_argument("--exp-path", type=Path,
                        help="Optional explicit path to the experimental structure (mmCIF).")
    parser.add_argument("--pred-path", type=Path,
                        help="Optional explicit path to the predicted structure (mmCIF).")
    parser.add_argument("--show-alignment", action="store_true",
                        help="Print the aligned sequences for manual inspection.")
    return parser.parse_args()


def load_property_row(property_file: Path, target: str) -> pd.Series:
    df = pd.read_csv(property_file)
    candidates = df[
        (df.get("name", "") == target)
        | (df.get("opsin name", "") == target)
        | (df.get("short_name", "") == target)
        | (df.get("display_name", "") == target)
        | (df.get("pdb_id", "") == target)
    ]
    if candidates.empty:
        raise ValueError(f"Could not locate property row for target '{target}' in {property_file}.")
    if len(candidates) > 1:
        print(f"[INFO] Found multiple property rows for '{target}'. Using the first match (index {candidates.index[0]}).")
    row = candidates.iloc[0]
    return row


def clean_short_name(name: str) -> str:
    cleaned = name.strip()
    return cleaned.replace(".", "").replace("+", "").replace("-", "_")


def discover_structure_path(name_candidates: Sequence[str], overrides: Optional[Path] = None) -> Path:
    if overrides:
        if overrides.exists():
            return overrides
        raise FileNotFoundError(f"Provided path does not exist: {overrides}")

    for candidate in name_candidates:
        filename = f"{candidate}.cif"
        for root in DEFAULT_SEARCH_DIRS:
            candidate_path = root / filename
            if candidate_path.exists():
                return candidate_path

    # Fallback: fuzzy search by substring
    hints = [cand for cand in name_candidates if cand]
    if hints:
        look_for = hints[0]
        for root in DEFAULT_SEARCH_DIRS:
            matches = list(root.glob(f"*{look_for}*.cif"))
            if matches:
                return matches[0]

    raise FileNotFoundError(f"Could not locate structure file for candidates: {name_candidates}")


def extract_sequence_info(path: Path, chain_id: str) -> StructureInfo:
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(path.stem, str(path))
    if chain_id not in structure[0]:
        raise ValueError(f"Chain '{chain_id}' not found in structure {path}")
    chain = structure[0][chain_id]
    residues: List[Residue] = []
    seq_chars: List[str] = []
    for residue in chain.get_residues():
        if not is_aa(residue, standard=True):
            continue
        try:
            aa = seq1(residue.get_resname())
        except Exception:
            aa = "X"
        residues.append(residue)
        seq_chars.append(aa)
    sequence = "".join(seq_chars)
    return StructureInfo(path=path, sequence=sequence, residues=residues)


def align_sequences(exp_seq: str, pred_seq: str) -> Tuple[str, str]:
    alignment = pairwise2.align.globalxx(exp_seq, pred_seq, one_alignment_only=True)[0]
    return alignment.seqA, alignment.seqB


def collect_atom_pairs(exp_info: StructureInfo, pred_info: StructureInfo, aligned_exp: str,
                       aligned_pred: str) -> Tuple[List[Atom], List[Atom], Dict[str, int]]:
    exp_index = -1
    pred_index = -1
    atoms_exp: List[Atom] = []
    atoms_pred: List[Atom] = []
    stats = {"skipped_missing_ca": 0, "paired_positions": 0}

    for aa_exp, aa_pred in zip(aligned_exp, aligned_pred):
        if aa_exp != "-":
            exp_index += 1
        if aa_pred != "-":
            pred_index += 1

        if aa_exp == "-" or aa_pred == "-":
            continue

        if exp_index >= len(exp_info.residues) or pred_index >= len(pred_info.residues):
            break

        exp_res = exp_info.residues[exp_index]
        pred_res = pred_info.residues[pred_index]

        if "CA" not in exp_res or "CA" not in pred_res:
            stats["skipped_missing_ca"] += 1
            continue

        atoms_exp.append(exp_res["CA"])
        atoms_pred.append(pred_res["CA"])
        stats["paired_positions"] += 1

    return atoms_exp, atoms_pred, stats


def compute_naive_pairs(exp_info: StructureInfo, pred_info: StructureInfo) -> Tuple[List[Atom], List[Atom], Dict[str, int]]:
    limit = min(len(exp_info.residues), len(pred_info.residues))
    atoms_exp: List[Atom] = []
    atoms_pred: List[Atom] = []
    skipped = 0
    for idx in range(limit):
        exp_res = exp_info.residues[idx]
        pred_res = pred_info.residues[idx]
        if "CA" not in exp_res or "CA" not in pred_res:
            skipped += 1
            continue
        atoms_exp.append(exp_res["CA"])
        atoms_pred.append(pred_res["CA"])
    return atoms_exp, atoms_pred, {"skipped_missing_ca": skipped, "paired_positions": len(atoms_exp)}


def summarize_alignment(aligned_exp: str, aligned_pred: str) -> Dict[str, int]:
    leading_exp = next((i for i, aa in enumerate(aligned_exp) if aa != "-"), 0)
    leading_pred = next((i for i, aa in enumerate(aligned_pred) if aa != "-"), 0)
    trailing_exp = next((i for i, aa in enumerate(reversed(aligned_exp)) if aa != "-"), 0)
    trailing_pred = next((i for i, aa in enumerate(reversed(aligned_pred)) if aa != "-"), 0)
    matches = sum(1 for a, b in zip(aligned_exp, aligned_pred) if a == b)
    mismatches = sum(1 for a, b in zip(aligned_exp, aligned_pred) if a != b and a != "-" and b != "-")
    gaps_exp = aligned_exp.count("-")
    gaps_pred = aligned_pred.count("-")
    return {
        "leading_gap_exp": leading_exp,
        "leading_gap_pred": leading_pred,
        "trailing_gap_exp": trailing_exp,
        "trailing_gap_pred": trailing_pred,
        "matches": matches,
        "mismatches": mismatches,
        "gaps_exp": gaps_exp,
        "gaps_pred": gaps_pred,
    }


def format_atom_count(stats: Dict[str, int]) -> str:
    return (
        f"paired CA atoms: {stats['paired_positions']}"
        f" (skipped missing CA: {stats['skipped_missing_ca']})"
    )


def main() -> None:
    args = parse_args()

    property_row = load_property_row(args.property_file, args.target)
    target_label = property_row.get("name", property_row.get("short_name", args.target))
    print(f"Analyzing property entry '{target_label}' (index {property_row.name})")
    print(f"  short_name: {property_row.get('short_name')}  |  display_name: {property_row.get('display_name')}")
    print(f"  pdb_id: {property_row.get('pdb_id')}  |  opsin name: {property_row.get('opsin name')}")

    cleaned_short = clean_short_name(str(property_row.get("short_name", "")))
    predicted_candidates = [
        f"{cleaned_short}_smile_model_0",
        f"{cleaned_short}_model_0",
        cleaned_short,
    ]

    if property_row.get("pdb_id") and isinstance(property_row.get("pdb_id"), str):
        experimental_candidates = [property_row.get("pdb_id").upper(), property_row.get("pdb_id").lower()]
    else:
        experimental_candidates = []

    exp_path = discover_structure_path(experimental_candidates, args.exp_path)
    pred_path = discover_structure_path(predicted_candidates, args.pred_path)

    print(f"Experimental structure: {exp_path}")
    print(f"Predicted structure:   {pred_path}")

    exp_info = extract_sequence_info(exp_path, args.chain)
    pred_info = extract_sequence_info(pred_path, args.chain)

    property_seq = str(property_row.get("sequence")) if pd.notna(property_row.get("sequence")) else ""
    structure_seq = str(property_row.get("structure_sequence")) if pd.notna(property_row.get("structure_sequence")) else ""

    print("\nSequence lengths (AA):")
    print(f"  property sequence:   {len(property_seq)} (empty -> {property_seq == ''})")
    print(f"  structure sequence:  {len(structure_seq)}")
    print(f"  experimental chain:  {len(exp_info.sequence)}")
    print(f"  predicted chain:     {len(pred_info.sequence)}")

    aligned_exp, aligned_pred = align_sequences(exp_info.sequence, pred_info.sequence)
    alignment_stats = summarize_alignment(aligned_exp, aligned_pred)

    print("\nAlignment summary (experimental vs predicted):")
    print(f"  leading gaps (exp/pred): {alignment_stats['leading_gap_exp']} / {alignment_stats['leading_gap_pred']}")
    print(f"  trailing gaps (exp/pred): {alignment_stats['trailing_gap_exp']} / {alignment_stats['trailing_gap_pred']}")
    print(f"  matches: {alignment_stats['matches']}  mismatches: {alignment_stats['mismatches']}")
    print(f"  gaps (exp/pred): {alignment_stats['gaps_exp']} / {alignment_stats['gaps_pred']}")

    naive_exp_atoms, naive_pred_atoms, naive_stats = compute_naive_pairs(exp_info, pred_info)
    if naive_exp_atoms and naive_pred_atoms:
        sup = Superimposer()
        sup.set_atoms(naive_exp_atoms, naive_pred_atoms)
        print("\nNaive RMSD (aligned by residue index, no trimming):")
        print(f"  RMSD: {sup.rms:.3f} Å using {format_atom_count(naive_stats)}")
    else:
        print("\nNaive RMSD: cannot compute (no overlapping CA atoms)")

    aligned_exp_atoms, aligned_pred_atoms, filtered_stats = collect_atom_pairs(exp_info, pred_info, aligned_exp, aligned_pred)
    if aligned_exp_atoms and aligned_pred_atoms:
        sup = Superimposer()
        sup.set_atoms(aligned_exp_atoms, aligned_pred_atoms)
        print("\nTrimmed RMSD (sequence alignment guided):")
        print(f"  RMSD: {sup.rms:.3f} Å using {format_atom_count(filtered_stats)}")
    else:
        print("\nTrimmed RMSD: cannot compute (no overlapping CA atoms after alignment)")

    if args.show_alignment:
        print("\nAligned sequences (experimental vs predicted):")
        width = 70
        for start in range(0, len(aligned_exp), width):
            end = start + width
            print(aligned_exp[start:end])
            print(aligned_pred[start:end])
            print()


if __name__ == "__main__":
    main()
