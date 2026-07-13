#!/usr/bin/env python3
"""Build and validate separate TARA R1/R2 reference rows.

R1 contains the manuscript-described TM5 pi-bulge.  Broad structural consensus
places W184--W186 at 5.40--5.42, a gap at 5.43, Y187 at 5.44, A188 at 5.45,
the bulged I189 at 5.451, and F190 at the downstream 5.46 anchor.  R2 retains
the standard register and has a gap only at the R1-specific column 5.451.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
PROTOS_SRC = ROOT / "protos" / "src"
if PROTOS_SRC.exists():
    sys.path.insert(0, str(PROTOS_SRC))


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


def main() -> int:
    output = ROOT / "opsin_output" / "dual_rhodopsins"
    source_path = ROOT / "type_I.csv"
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
    # Restore the R1 pi-bulge register described in manuscript.md.
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
    tara_a["5.43"] = "-"
    tara_a["5.44"] = "Y187"
    tara_a["5.45"] = "A188"
    tara_a["5.451"] = "I189"
    tara_a["5.46"] = "F190"
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
    print(qc_table.to_string(index=False))
    print(f"Candidate reference: {candidate_path}")
    return 0 if (qc_table["tm_roundtrip_differences"] == 0).all() else 1


if __name__ == "__main__":
    raise SystemExit(main())
