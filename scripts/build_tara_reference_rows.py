#!/usr/bin/env python3
"""Build and validate separate TARA R1/R2 reference rows.

R1 contains the manuscript-described TM5 pi-bulge.  Its hand-curated register
treats the helix as one ordinary position shorter at its N-terminal end: 5.40
is empty, W184--A188 occupy 5.41--5.45, the bulged I189 occupies 5.451, and
F190 remains at the downstream 5.46 anchor.  R2 retains the standard register
and has a gap only at the R1-specific column 5.451.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))
PROTOS_SRC = ROOT / "protos" / "src"
if PROTOS_SRC.exists():
    sys.path.insert(0, str(PROTOS_SRC))

from src.curated_grn_storage import synchronize_runtime_reference


TARA_A_TM5_REGISTER = {
    "5.40": "-",
    "5.41": "W184",
    "5.42": "M185",
    "5.43": "W186",
    "5.44": "Y187",
    "5.45": "A188",
    "5.451": "I189",
    "5.46": "F190",
}


def add_after(columns: list[str], existing: str, new: str) -> list[str]:
    result = list(columns)
    if new in result:
        return result
    result.insert(result.index(existing) + 1, new)
    return result


def alias_structure(processor, source_id: str, target_id: str) -> None:
    frame = processor.load_entity(source_id)
    if frame is None:
        raise ValueError(f"Missing split structure {source_id}")
    reset = frame.reset_index()
    reset["structure_id"] = target_id
    processor.save_entity(
        target_id,
        reset,
        metadata={
            "source": "tara_reference_alias",
            "parent_virtual_structure": source_id,
            "domain": target_id.rsplit("_", 1)[-1],
        },
    )


def apply_tara_a_tm5_register(row: pd.Series) -> pd.Series:
    """Apply the curated continuous TARA_A register around the 5.451 bulge."""

    corrected = row.copy()
    missing = set(TARA_A_TM5_REGISTER).difference(corrected.index)
    if missing:
        raise ValueError(f"TARA_A row is missing TM5 columns: {sorted(missing)}")
    for grn, residue in TARA_A_TM5_REGISTER.items():
        corrected[grn] = residue

    residues = [
        corrected[grn]
        for grn in TARA_A_TM5_REGISTER
        if corrected[grn] != "-"
    ]
    expected = ["W184", "M185", "W186", "Y187", "A188", "I189", "F190"]
    if residues != expected or len(residues) != len(set(residues)):
        raise ValueError(f"Invalid TARA_A TM5 register: {residues}")
    return corrected


def sync_runtime_reference(candidate: pd.DataFrame) -> str:
    """Validate generated TARA rows and restore exact authoritative copies."""

    source_path = (
        ROOT
        / "protos"
        / "src"
        / "protos"
        / "reference_data"
        / "grn"
        / "reference"
        / "type_I_opsins.csv"
    )
    authoritative = pd.read_csv(source_path, index_col=0, dtype=str).fillna("-")
    normalized = candidate.reindex(
        index=authoritative.index, columns=authoritative.columns, fill_value="-"
    ).fillna("-")
    if not normalized.equals(authoritative):
        differences = normalized.ne(authoritative)
        rows = differences.index[differences.any(axis=1)].tolist()
        raise ValueError(
            "Generated TARA candidate differs from authoritative ProtOS table "
            f"for rows: {rows}"
        )
    result = synchronize_runtime_reference(source_path, project_root=ROOT)
    return result["sha256"]


def main() -> int:
    output = ROOT / "opsin_output" / "dual_rhodopsins"
    source_path = (
        ROOT
        / "protos"
        / "src"
        / "protos"
        / "reference_data"
        / "grn"
        / "reference"
        / "type_I_opsins.csv"
    )
    annotations_path = output / "dual_domain_grn_annotations_local.csv"
    if not annotations_path.exists():
        raise FileNotFoundError("Run detect_annotate_dual_rhodopsins.py first")

    source = pd.read_csv(source_path, index_col=0, dtype=str)
    split = pd.read_csv(annotations_path, index_col=0, dtype=str)
    columns = add_after([str(column) for column in source.columns], "5.45", "5.451")
    candidate = source.drop(
        index=["7pl9", "TARA_A", "TARA_B"], errors="ignore"
    ).reindex(columns=columns, fill_value="-")

    tara_a = split.loc["7pl9_A"].reindex(columns, fill_value="-").copy()
    tara_b = split.loc["7pl9_B"].reindex(columns, fill_value="-").copy()
    # Confirm the uncurated Protos register before applying the manual shift.
    expected_r1_register = {
        "4.62": "T183",
        "5.40": "W184",
        "5.41": "M185",
        "5.42": "W186",
        "5.43": "Y187",
        "5.44": "A188",
        "5.45": "I189",
        "5.46": "F190",
    }
    assert all(tara_a[grn] == residue for grn, residue in expected_r1_register.items())
    tara_a = apply_tara_a_tm5_register(tara_a)
    tara_b["5.451"] = "-"
    candidate.loc["TARA_A"] = tara_a
    candidate.loc["TARA_B"] = tara_b
    candidate.index.name = source.index.name

    candidate_path = output / "type_I_with_tara_domains.csv"
    candidate.to_csv(candidate_path)
    candidate.to_csv(ROOT / "data" / "grn" / "reference" / "type_I_with_tara_domains.csv")

    import protos
    from protos.processing.grn import GRNProcessor
    from protos.processing.structure import StructureProcessor

    protos.set_data_path(str(ROOT / "data"))
    structures = StructureProcessor("tara_reference_rows")
    alias_structure(structures, "7pl9_A", "TARA_A")
    alias_structure(structures, "7pl9_B", "TARA_B")
    existing = (
        structures.get_dataset_entities("mo_dual_rhodopsin_domains")
        if structures.dataset_manager.dataset_exists("mo_dual_rhodopsin_domains")
        else []
    )
    entities = list(dict.fromkeys([*existing, "TARA_A", "TARA_B"]))
    structures.create_dataset(
        "mo_dual_rhodopsin_domains",
        entities,
        metadata={"source": "tandem_repeat_detection_and_tara_aliases"},
    )

    sequences = {}
    for structure_id in ["TARA_A", "TARA_B"]:
        frame = structures.load_entity(structure_id).reset_index()
        ca = frame[frame["atom_name"].astype(str).str.upper() == "CA"].sort_values("auth_seq_id")
        sequences[structure_id] = "".join(ca["res_name1l"].astype(str))
    grn = GRNProcessor("tara_reference_rows")
    roundtrip, summary = grn.annotate_sequences(
        sequences,
        reference_table="type_I_with_tara_domains",
        protein_family="mo",
        assign_unambiguous_insertions=True,
    )
    roundtrip.to_csv(output / "tara_reference_roundtrip.csv")
    with (output / "tara_reference_roundtrip_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)

    qc = []
    for structure_id in ["TARA_A", "TARA_B"]:
        expected = candidate.loc[structure_id]
        # Flexible-region expansion can emit duplicate loop labels.  Collapse
        # them for QC; the TM labels being validated are unique.
        observed_raw = roundtrip.loc[structure_id]
        observed = observed_raw.groupby(level=0, sort=False).first().reindex(columns, fill_value="-")
        tm_columns = [column for column in columns if column.split(".", 1)[0] in set("1234567")]
        differences = [column for column in tm_columns if expected[column] != observed[column]]
        qc.append(
            {
                "structure": structure_id,
                "selected_reference": summary["per_sequence"][structure_id]["reference"],
                "coverage": summary["per_sequence"][structure_id]["coverage"],
                "tm_roundtrip_differences": len(differences),
                "difference_columns": ";".join(differences),
                "grn_5_44": observed.get("5.44", "-"),
                "grn_5_45": observed.get("5.45", "-"),
                "grn_5_451": observed.get("5.451", "-"),
                "grn_5_46": observed.get("5.46", "-"),
            }
        )
    qc_table = pd.DataFrame(qc)
    qc_table.to_csv(output / "tara_reference_roundtrip_qc.csv", index=False)
    success = (qc_table["tm_roundtrip_differences"] == 0).all()
    digest = sync_runtime_reference(candidate) if success else "not-synchronized"
    print(qc_table.to_string(index=False))
    print(f"Candidate reference: {candidate_path}")
    print(f"Runtime reference SHA-256: {digest}")
    return 0 if success else 1


if __name__ == "__main__":
    raise SystemExit(main())
