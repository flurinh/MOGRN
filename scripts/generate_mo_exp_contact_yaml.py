#!/usr/bin/env python3
"""Generate experimental opsin YAML configs with ligand-binding constraints.

This script loads the `mo_exp` structural dataset using the protos CIF processor,
identifies the Schiff base lysine for each structure, and writes YAML configs with
pocket constraints suitable for Boltz structure prediction runs.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq

# ---------------------------------------------------------------------------
# Project-relative imports (ensure protos package resolves to this repo copy)
# ---------------------------------------------------------------------------
REPO_ROOT = Path(__file__).resolve().parents[1]
PROTOS_SRC = REPO_ROOT / "protos" / "src"
if not PROTOS_SRC.exists():
    raise FileNotFoundError(f"Expected protos sources at {PROTOS_SRC}")

sys.path.insert(0, str(PROTOS_SRC))
sys.path.insert(0, str(REPO_ROOT))

from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.processing.structure.struct_utils import ALPHA_CARBON
from protos.io.paths.path_config import ProtosPaths

from src.lyr_processing import process_lyr_in_processor_data
from src.common_utils import find_retinal_within_cutoff

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEFAULT_DATASET_ID = "mo_exp"
DEFAULT_CHAIN_ID = "A"
RETINAL_CUTOFF = 6.0
MAX_SCHIFF_DISTANCE = 3.0  # Å

# ---------------------------------------------------------------------------
# Data containers
# ---------------------------------------------------------------------------

@dataclass
class SchiffBaseInfo:
    """Summary of the Schiff base lysine for a structure."""

    sequence: str
    lys_seq_index: int
    lys_auth_seq_id: Optional[int]
    lys_gen_seq_id: Optional[int]
    distance: float
    chain_id: str
    lys_residue: str
    lys_atom_name: str
    ret_atom_name: str


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def sanitize_name(raw: str) -> str:
    """Convert identifiers into filesystem-safe YAML stems."""
    if raw is None:
        raise ValueError("Cannot sanitise an empty name")
    token = str(raw).strip()
    if not token:
        raise ValueError("Name resolves to empty string")
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", token)
    safe = re.sub(r"_+", "_", safe).strip("_")
    if not safe:
        raise ValueError(f"Failed to sanitise name '{raw}'")
    return safe


def configure_protos_paths(data_root: Path) -> None:
    """Ensure protos uses the repository-local data directory."""
    data_root = data_root.resolve()
    os.environ["PROTOS_DATA_ROOT"] = str(data_root)
    os.environ["PROTOS_REF_DATA_ROOT"] = str(data_root)
    ProtosPaths(user_data_root=str(data_root), ref_data_root=str(data_root), create_dirs=True)


def load_structures(dataset_path: Path, data_root: Path) -> Tuple[CifBaseProcessor, List[str]]:
    """Load all structures listed in the dataset definition via CifBaseProcessor."""
    with dataset_path.open("r", encoding="utf-8") as handle:
        dataset = json.load(handle)

    pdb_ids: List[str] = dataset.get("pdb_ids", [])
    if not pdb_ids:
        raise ValueError(f"Dataset definition '{dataset_path}' does not list any PDB IDs")

    processor = CifBaseProcessor(
        name=f"{dataset.get('id', DEFAULT_DATASET_ID)}_processor",
        data_root=str(data_root.resolve()),
        processor_data_dir="structure",
    )
    processor.load_structures(pdb_ids, apply_dtypes=True, debug=False)
    process_lyr_in_processor_data(processor, retinal_res_name="RET")
    return processor, pdb_ids


def resolve_pdb_id(processor: CifBaseProcessor, pdb_id: str) -> str:
    """Return the canonical identifier used inside the processor for a PDB ID."""
    target = pdb_id.lower()
    for candidate in processor.pdb_ids:
        if candidate.lower() == target:
            return candidate
    raise KeyError(f"PDB ID '{pdb_id}' not found in loaded dataset")


def build_sequence_entries(df: pd.DataFrame) -> List[Dict[str, Optional[int]]]:
    """Collapse atomic records into ordered residue entries for one chain."""
    if df.empty:
        raise ValueError("Chain dataframe is empty")

    ca_df = df[df["res_atom_name"] == ALPHA_CARBON].copy()
    if ca_df.empty:
        raise ValueError("No alpha-carbon atoms found for chain")

    ca_df.sort_values(by=["gen_seq_id", "auth_seq_id", "atom_id"], inplace=True)

    entries: List[Dict[str, Optional[int]]] = []
    seen: set = set()

    for _, row in ca_df.iterrows():
        gen_seq = row["gen_seq_id"]
        auth_seq = row["auth_seq_id"]
        key = (
            int(gen_seq) if not pd.isna(gen_seq) else None,
            int(auth_seq) if not pd.isna(auth_seq) else None,
        )
        if key in seen:
            continue
        seen.add(key)

        res_name3l = str(row["res_name3l"]).strip().upper()
        res_name1l = str(row["res_name1l"]).strip().upper()
        if len(res_name1l) != 1 or res_name1l == "?":
            res_name1l = "X"

        entries.append(
            {
                "seq_index": len(entries) + 1,
                "gen_seq_id": key[0],
                "auth_seq_id": key[1],
                "res_name3l": res_name3l,
                "res_name1l": res_name1l,
            }
        )

    if not entries:
        raise ValueError("Failed to derive residue entries from CA atoms")
    return entries


def compute_schiff_base(
    processor: CifBaseProcessor,
    pdb_id: str,
    chain_preference: str = DEFAULT_CHAIN_ID,
) -> SchiffBaseInfo:
    """Identify the lysine forming the Schiff base and return sequence metadata."""
    canonical_id = resolve_pdb_id(processor, pdb_id)
    structure_df = processor.data[
        processor.data["pdb_id"].str.lower() == canonical_id.lower()
    ]
    if structure_df.empty:
        raise ValueError(f"No structural data for {pdb_id}")

    available_chains = structure_df["auth_chain_id"].unique().tolist()
    chain_id = chain_preference if chain_preference in available_chains else available_chains[0]

    chain_df = structure_df[
        (structure_df["auth_chain_id"] == chain_id)
        & (structure_df["group"] == "ATOM")
    ].copy()
    if chain_df.empty:
        raise ValueError(f"Chain {chain_id} for {pdb_id} has no ATOM records")

    entries = build_sequence_entries(chain_df)
    sequence = "".join(entry["res_name1l"] for entry in entries)

    lys_entries = [entry for entry in entries if entry["res_name3l"] == "LYS"]
    if not lys_entries:
        raise ValueError(f"No lysine residues found in chain {chain_id} of {pdb_id}")

    retinal_df = find_retinal_within_cutoff(
        structure_df,
        chain_df,
        cutoff=RETINAL_CUTOFF,
        retinal_name="RET",
    )
    if retinal_df.empty:
        raise ValueError(f"No retinal-like residues found in {pdb_id}")
    ret_coords = retinal_df[["x", "y", "z"]].to_numpy(dtype=float)

    best_distance = math.inf
    best_entry: Optional[Dict[str, Optional[int]]] = None
    best_nz_atom = None
    best_ret_atom = None

    for entry in lys_entries:
        auth_seq = entry["auth_seq_id"]
        residue_atoms = chain_df[
            (chain_df["res_name3l"] == "LYS")
            & (chain_df["auth_seq_id"] == auth_seq)
            & (chain_df["res_atom_name"] == "NZ")
        ]
        if residue_atoms.empty:
            continue
        nz_coords = residue_atoms[["x", "y", "z"]].to_numpy(dtype=float)
        distances = np.linalg.norm(nz_coords[:, None, :] - ret_coords[None, :, :], axis=-1)
        min_idx = np.unravel_index(np.argmin(distances), distances.shape)
        min_distance = float(distances[min_idx])
        if min_distance < best_distance:
            best_distance = min_distance
            best_entry = entry
            best_nz_atom = residue_atoms.iloc[min_idx[0]]
            best_ret_atom = retinal_df.iloc[min_idx[1]]

    if best_entry is None:
        raise ValueError(f"Could not locate NZ atoms for lysines in {pdb_id}")

    if not math.isfinite(best_distance) or best_distance > MAX_SCHIFF_DISTANCE:
        raise ValueError(
            f"Closest lys-retinal distance {best_distance:.2f}Å exceeds expected threshold "
            f"for {pdb_id} chain {chain_id}"
        )

    seq_index = int(best_entry["seq_index"])
    lys_letter = sequence[seq_index - 1] if 0 < seq_index <= len(sequence) else "?"

    lys_atom_name = (
        str(best_nz_atom["res_atom_name"]).strip().upper()
        if best_nz_atom is not None
        else "NZ"
    )
    ret_atom_name = (
        str(best_ret_atom["res_atom_name"]).strip().upper()
        if best_ret_atom is not None and pd.notna(best_ret_atom["res_atom_name"])
        else "C15"
    )

    return SchiffBaseInfo(
        sequence=sequence,
        lys_seq_index=seq_index,
        lys_auth_seq_id=best_entry["auth_seq_id"],
        lys_gen_seq_id=best_entry["gen_seq_id"],
        distance=best_distance,
        chain_id=chain_id,
        lys_residue=lys_letter,
        lys_atom_name=lys_atom_name,
        ret_atom_name=ret_atom_name,
    )


def normalise_pdb_id(raw: str) -> Optional[str]:
    """Standardise PDB identifiers from the property table."""
    if raw is None or (isinstance(raw, float) and math.isnan(raw)):
        return None
    token = str(raw).strip()
    if not token or token.lower() == "nan":
        return None
    token = token.replace(".0", "")
    token = token.upper()
    if len(token) != 4:
        # Some CSV entries may appear as scientific notation (e.g., 1E12)
        token = token[:4]
    return token


def build_yaml_doc(info: SchiffBaseInfo) -> CommentedMap:
    """Create the YAML document with ligand pocket and bond constraints."""
    doc = CommentedMap()
    doc["version"] = 1
    doc.yaml_add_eol_comment("Optional, defaults to 1", "version")

    sequences = CommentedSeq()

    protein_entry = CommentedMap()
    protein_entry["protein"] = CommentedMap()
    protein_entry["protein"]["id"] = "A"
    protein_entry["protein"]["sequence"] = info.sequence
    sequences.append(protein_entry)

    ligand_entry = CommentedMap()
    ligand_entry["ligand"] = CommentedMap()
    ligand_entry["ligand"]["id"] = "B"
    ligand_entry["ligand"]["ccd"] = "RSB"
    sequences.append(ligand_entry)

    doc["sequences"] = sequences

    constraints = CommentedSeq()

    bond_entry = CommentedMap()
    bond_entry["bond"] = CommentedMap()

    atom1 = CommentedSeq(["A", int(info.lys_seq_index), info.lys_atom_name or "NZ"])
    atom1.fa.set_flow_style()
    atom2 = CommentedSeq(["B", 1, info.ret_atom_name or "C15"])
    atom2.fa.set_flow_style()

    bond_entry["bond"]["atom1"] = atom1
    bond_entry["bond"]["atom2"] = atom2

    constraints.append(bond_entry)

    doc["constraints"] = constraints

    return doc


def write_yaml(output_path: Path, yaml_doc: CommentedMap) -> None:
    """Serialise YAML with stable formatting."""
    yaml = YAML()
    yaml.default_flow_style = False
    yaml.indent(mapping=2, sequence=4, offset=2)
    yaml.preserve_quotes = True
    yaml.width = 1000

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        yaml.dump(yaml_doc, handle)


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dataset-id",
        default=DEFAULT_DATASET_ID,
        help="Structure dataset identifier (default: mo_exp).",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=REPO_ROOT / "data",
        help="Root directory containing protos structural data (default: data/).",
    )
    parser.add_argument(
        "--property-csv",
        type=Path,
        default=REPO_ROOT / "property" / "mo_exp.csv",
        help="CSV providing short_name to PDB mappings (default: property/mo_exp.csv).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPO_ROOT / "yaml_configs" / "mo_exp_experimental_contact",
        help="Directory where augmented YAML configs will be written.",
    )
    parser.add_argument(
        "--chain-id",
        default=DEFAULT_CHAIN_ID,
        help="Preferred chain identifier when extracting sequences (default: A).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    dataset_path = args.data_root / "structure" / "structure_dataset" / "standard" / f"{args.dataset_id}.json"
    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset definition not found at {dataset_path}")

    if not args.property_csv.exists():
        raise FileNotFoundError(f"Property CSV not found at {args.property_csv}")

    configure_protos_paths(args.data_root)
    processor, pdb_ids = load_structures(dataset_path, args.data_root)

    property_df = pd.read_csv(args.property_csv)
    property_df = property_df[property_df["pdb_id"].notna()].copy()
    property_df["pdb_norm"] = property_df["pdb_id"].apply(normalise_pdb_id)
    property_df.dropna(subset=["pdb_norm", "short_name"], inplace=True)

    unique_pdbs = sorted({pid for pid in property_df["pdb_norm"] if pid})

    schiff_info: Dict[str, SchiffBaseInfo] = {}
    for pdb_id in unique_pdbs:
        info = compute_schiff_base(processor, pdb_id, chain_preference=args.chain_id)
        schiff_info[pdb_id] = info
        residue_check = "OK" if info.lys_residue == "K" else "WARNING"
        print(
            f"{pdb_id}: chain {info.chain_id}, Lys seq index {info.lys_seq_index}, "
            f"residue '{info.lys_residue}' [{residue_check}], distance {info.distance:.3f} Å"
        )
        if info.lys_residue != "K":
            print(
                f"  -> Expected lysine at position {info.lys_seq_index} but found '{info.lys_residue}'."
            )

    written = 0
    for _, row in property_df.iterrows():
        short_name = sanitize_name(row["short_name"])
        pdb_id = row["pdb_norm"]
        if pdb_id not in schiff_info:
            print(f"[WARNING] Skipping {short_name}: no Schiff base info for PDB {pdb_id}")
            continue

        info = schiff_info[pdb_id]
        if info.lys_residue != "K":
            raise ValueError(
                f"Sequence verification failed for {short_name} ({pdb_id}): "
                f"position {info.lys_seq_index} is '{info.lys_residue}', expected 'K'."
            )
        yaml_doc = build_yaml_doc(info)
        output_path = args.output_dir / f"{short_name}.yaml"
        write_yaml(output_path, yaml_doc)
        written += 1

    print(f"Wrote {written} YAML configs to {args.output_dir}")


if __name__ == "__main__":
    main()
