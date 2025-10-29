#!/usr/bin/env python3
"""Generate YAML configs for opsin sequences from FASTA files or structure datasets."""

from __future__ import annotations

import argparse
import json
import re
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict, Dict, Iterable, List, Optional, Sequence, Tuple

from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq

VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWY")
DEFAULT_LIGAND_SMILES = "CC=C(C)C=CC=C(C)C=CC1=C(CCCC1(C)C)C"
DEFAULT_DATA_ROOT = Path("data")
DEFAULT_DATASET_DIR = Path("structure/structure_dataset/standard")
DEFAULT_MMCIF_DIR = Path("structure/mmcif")
DEFAULT_SEARCH_DIRS = (
    Path("structures/mo_exp"),
    Path("structures/hideaki_exp"),
    Path("structures/mo_pred"),
    Path("structures/hideaki_pred"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create a YAML file per opsin sequence for use with the opsin workflow."
        )
    )
    parser.add_argument(
        "--fasta",
        type=Path,
        help="Path to a FASTA file containing opsin sequences (optional if --dataset is used).",
    )
    parser.add_argument(
        "--dataset",
        help="Structure dataset identifier (e.g., mo_exp) to extract sequences from experimental structures.",
    )
    parser.add_argument(
        "--data-root",
        type=Path,
        default=DEFAULT_DATA_ROOT,
        help="Root directory containing PROTOS data folders (default: data).",
    )
    parser.add_argument(
        "--chain",
        default="A",
        help="Chain identifier to extract from structures when using --dataset (default: A).",
    )
    parser.add_argument(
        "--property-file",
        type=Path,
        help="Optional CSV providing short_name mappings (default inferred from dataset if available).",
    )
    parser.add_argument(
        "--name-source",
        choices=["short_name", "pdb_id"],
        default="short_name",
        help="Field to use for YAML filenames when using --dataset (default: short_name).",
    )
    parser.add_argument(
        "--outputs-dir",
        type=Path,
        help="Directory where YAML configs will be written (default depends on mode).",
    )
    parser.add_argument(
        "--ligand-smiles",
        default=DEFAULT_LIGAND_SMILES,
        help="SMILES string to include for the ligand entry in each YAML file.",
    )
    return parser.parse_args()


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"FASTA file not found: {path}")

    records: List[Tuple[str, str]] = []
    header: str | None = None
    sequence_parts: List[str] = []

    with path.open("r", encoding="utf-8") as handle:
        for raw_line in handle:
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(sequence_parts)))
                header = line[1:].strip()
                sequence_parts = []
            else:
                sequence_parts.append(line.replace(" ", "").upper())

    if header is not None:
        records.append((header, "".join(sequence_parts)))

    if not records:
        raise ValueError(f"No sequences found in FASTA file: {path}")

    return records


def sanitise_name(raw_header: str) -> str:
    token = raw_header.split()[0]
    safe = re.sub(r"[^A-Za-z0-9_]+", "_", token).strip("_")
    if not safe:
        raise ValueError(f"Could not derive a valid file stem from FASTA header: {raw_header}")
    return safe


def validate_sequences(records: Iterable[Tuple[str, str]]) -> List[str]:
    problems: List[str] = []
    for header, sequence in records:
        invalid = sorted(set(sequence) - VALID_AMINO_ACIDS)
        if invalid:
            msg = (
                f"Sequence '{header}' contains non-standard amino acids: "
                f"{','.join(invalid)}"
            )
            problems.append(msg)
        if not sequence:
            problems.append(f"Sequence '{header}' is empty.")
    return problems


def load_property_mapping(property_file: Path) -> Dict[str, List[str]]:
    if not property_file or not property_file.exists():
        return {}

    import pandas as pd

    df = pd.read_csv(property_file)
    mapping: DefaultDict[str, List[str]] = defaultdict(list)
    for _, row in df.iterrows():
        pdb_id = str(row.get("pdb_id", "")).strip()
        short_name = str(row.get("short_name", "")).strip()
        if pdb_id and short_name:
            cleaned = sanitise_name(short_name)
            bucket = mapping[pdb_id.upper()]
            if cleaned and cleaned not in bucket:
                bucket.append(cleaned)
    return dict(mapping)


def discover_mmcif_path(pdb_id: str, data_root: Path) -> Path:
    candidates: Sequence[Path] = []
    pdb_upper = pdb_id.upper()
    pdb_lower = pdb_id.lower()

    mmcif_dir = data_root / DEFAULT_MMCIF_DIR
    candidates.append(mmcif_dir / f"{pdb_upper}.cif")
    candidates.append(mmcif_dir / f"{pdb_lower}.cif")

    for extra_root in DEFAULT_SEARCH_DIRS:
        candidates.append(extra_root / f"{pdb_upper}.cif")
        candidates.append(extra_root / f"{pdb_lower}.cif")

    for candidate in candidates:
        if candidate.exists():
            return candidate

    raise FileNotFoundError(f"Could not locate mmCIF file for {pdb_id} in known directories.")


def extract_chain_sequence(mmcif_path: Path, chain_id: str) -> str:
    from Bio.PDB.MMCIFParser import MMCIFParser
    from Bio.PDB.Polypeptide import is_aa
    from Bio.SeqUtils import seq1

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(mmcif_path.stem, str(mmcif_path))

    model = structure[0]
    if chain_id not in model:
        available = ", ".join(model.child_dict.keys())
        raise ValueError(f"Chain '{chain_id}' not found in {mmcif_path} (available: {available})")

    chain = model[chain_id]
    residues = [res for res in chain.get_residues() if is_aa(res, standard=True)]
    if not residues:
        raise ValueError(f"No standard amino acid residues found in {mmcif_path} chain {chain_id}")

    sequence_parts: List[str] = []
    for residue in residues:
        try:
            sequence_parts.append(seq1(residue.get_resname()))
        except Exception:
            sequence_parts.append("X")
    return "".join(sequence_parts)


def load_dataset_sequences(dataset_id: str, data_root: Path, chain_id: str,
                           property_mapping: Dict[str, str], name_source: str) -> List[Tuple[str, str]]:
    dataset_path = data_root / DEFAULT_DATASET_DIR / f"{dataset_id}.json"
    if not dataset_path.exists():
        raise FileNotFoundError(f"Dataset definition not found: {dataset_path}")

    with dataset_path.open("r", encoding="utf-8") as handle:
        dataset = json.load(handle)

    pdb_ids = dataset.get("pdb_ids")
    if not pdb_ids:
        raise ValueError(f"Dataset {dataset_id} does not list any pdb_ids")

    sequences: List[Tuple[str, str]] = []
    for pdb_id in pdb_ids:
        mmcif_path = discover_mmcif_path(pdb_id, data_root)
        sequence = extract_chain_sequence(mmcif_path, chain_id)

        if name_source == "short_name":
            short_names = property_mapping.get(pdb_id.upper(), [])
            if short_names:
                for candidate in short_names:
                    sequences.append((candidate, sequence))
                continue
            record_names = [sanitise_name(pdb_id)]
        else:
            record_names = [sanitise_name(pdb_id)]

        for record_name in record_names:
            sequences.append((record_name, sequence))
    return sequences


def build_yaml_doc(sequence: str, ligand_smiles: str) -> CommentedMap:
    doc = CommentedMap()
    doc["version"] = 1
    doc.yaml_add_eol_comment("Optional, defaults to 1", "version")

    sequences = CommentedSeq()

    protein_entry = CommentedMap()
    protein_entry["protein"] = CommentedMap()
    protein_entry["protein"]["id"] = "A"
    protein_entry["protein"]["sequence"] = sequence
    sequences.append(protein_entry)

    ligand_entry = CommentedMap()
    ligand_entry["ligand"] = CommentedMap()
    ligand_entry["ligand"]["id"] = "B"
    ligand_entry["ligand"]["smiles"] = ligand_smiles
    sequences.append(ligand_entry)

    doc["sequences"] = sequences
    return doc


def main() -> None:
    args = parse_args()

    if not args.fasta and not args.dataset:
        raise SystemExit("Please provide either --fasta or --dataset.")

    if args.dataset and args.property_file is None:
        default_property = Path("property/mo_exp.csv")
        if default_property.exists():
            args.property_file = default_property

    property_mapping = load_property_mapping(args.property_file) if args.dataset else {}

    if args.dataset:
        records = load_dataset_sequences(
            dataset_id=args.dataset,
            data_root=args.data_root,
            chain_id=args.chain,
            property_mapping=property_mapping,
            name_source=args.name_source,
        )
    else:
        records = read_fasta(args.fasta)
        validation_errors = validate_sequences(records)
        if validation_errors:
            joined = "\n".join(validation_errors)
            raise ValueError(f"Invalid sequences detected:\n{joined}")

    if args.output_dir is None:
        if args.dataset:
            args.output_dir = Path("yaml_configs") / f"{args.dataset}_experimental"
        else:
            args.output_dir = Path("yaml_configs/new_opsins")

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    yaml = YAML()
    yaml.default_flow_style = False
    yaml.indent(mapping=2, sequence=4, offset=2)
    yaml.preserve_quotes = True
    yaml.width = 1000

    from collections import defaultdict

    name_counts: Dict[str, int] = defaultdict(int)
    for header, sequence in records:
        base_name = sanitise_name(header)
        count = name_counts[base_name]
        name = base_name if count == 0 else f"{base_name}_{count+1}"
        name_counts[base_name] += 1

        yaml_doc = build_yaml_doc(sequence, args.ligand_smiles)
        output_path = output_dir / f"{name}.yaml"
        with output_path.open("w", encoding="utf-8") as handle:
            yaml.dump(yaml_doc, handle)
        print(f"Created YAML config: {output_path}")


if __name__ == "__main__":
    main()
