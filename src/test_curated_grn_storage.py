"""Tests for the terminal curated-GRN structure overwrite."""

from __future__ import annotations

import hashlib
from pathlib import Path

import pandas as pd
from src.curated_grn_storage import (
    aliases_from_tandem_manifest,
    load_curated_structure_aliases,
    overwrite_frame_grns,
    overwrite_stored_structures_with_curated_grns,
    sequence_alignment_position_map,
    synchronize_runtime_reference,
)


def structure_frame(structure_id: str, residues: str) -> pd.DataFrame:
    rows = []
    atom_id = 1
    for sequence_number, residue in enumerate(residues, start=1):
        for atom_name in ("N", "CA"):
            rows.append(
                {
                    "structure_id": structure_id,
                    "atom_id": atom_id,
                    "auth_chain_id": "A",
                    "auth_seq_id": sequence_number,
                    "res_name1l": residue,
                    "atom_name": atom_name,
                    "grn": "provisional",
                }
            )
            atom_id += 1
    return pd.DataFrame(rows).set_index(["structure_id", "atom_id"])


class FakeProcessor:
    def __init__(self, frames: dict[str, pd.DataFrame]) -> None:
        self.frames = frames
        self.saved_metadata = {}

    def load_entity(self, structure_id: str) -> pd.DataFrame | None:
        return self.frames.get(structure_id)

    def save_entity(self, structure_id, frame, metadata=None) -> None:
        self.frames[structure_id] = frame.set_index(["structure_id", "atom_id"])
        self.saved_metadata[structure_id] = metadata


def test_overwrite_clears_provisional_grns_and_maps_curated_residues() -> None:
    row = pd.Series({"1.50": "A1", "1.51": "-", "1.52": "C2"})

    annotated, report = overwrite_frame_grns(
        structure_frame("TARA_A", "AC"),
        row,
        structure_id="TARA_A",
        reference_id="TARA_A",
    )

    assert set(annotated.loc[annotated["auth_seq_id"] == 1, "grn"]) == {"1.50"}
    assert set(annotated.loc[annotated["auth_seq_id"] == 2, "grn"]) == {"1.52"}
    assert report["stale_grn_atoms_removed"] == 4
    assert report["mapped_positions"] == 2
    assert report["missing_residue_count"] == 0


def test_authoritative_grn_overwrite_reports_but_keeps_sequence_variant() -> None:
    annotated, report = overwrite_frame_grns(
        structure_frame("TARA_A", "AC"),
        pd.Series({"1.50": "G1"}),
        structure_id="TARA_A",
        reference_id="TARA_A",
    )

    assert set(annotated.loc[annotated["auth_seq_id"] == 1, "grn"]) == {"1.50"}
    assert report["residue_mismatch_count"] == 1
    assert report["residue_mismatches"][0]["expected"] == "G"
    assert report["residue_mismatches"][0]["observed"] == ["A"]


def test_unknown_curated_residue_is_a_wildcard() -> None:
    annotated, report = overwrite_frame_grns(
        structure_frame("4kly", "M"),
        pd.Series({"1.50": "X1"}),
        structure_id="4kly",
        reference_id="4kly",
    )

    assert set(annotated["grn"]) == {"1.50"}
    assert report["residue_mismatch_count"] == 0


def test_sequence_alignment_maps_replacement_structure_numbering() -> None:
    frame = structure_frame("replacement", "QQACD")
    row = pd.Series({"1.50": "A1", "1.51": "C2", "1.52": "D3"})

    position_map, alignment = sequence_alignment_position_map(frame, row)
    annotated, report = overwrite_frame_grns(
        frame,
        row,
        structure_id="replacement",
        reference_id="source",
        position_map=position_map,
    )

    assert position_map == {1: 3, 2: 4, 3: 5}
    assert alignment["mismatches"] == []
    assert report["missing_residue_count"] == 0
    assert (
        annotated.loc[annotated["auth_seq_id"].eq(3), "grn"].unique().tolist()
        == ["1.50"]
    )


def test_current_structure_aliases_use_clean_protos_rows() -> None:
    aliases = load_curated_structure_aliases()

    assert aliases["9jws"]["reference_id"] == "hwmr_model_0"
    assert aliases["8rso"]["reference_id"] == "7avn"
    assert aliases["9j7w"]["reference_id"] == "knchr_j444_refine5"
    assert {
        aliases[structure_id]["position_mapping"]
        for structure_id in ("9jws", "8rso", "9j7w")
    } == {"sequence_alignment"}


def test_tandem_aliases_keep_a_and_b_as_separate_references() -> None:
    manifest = {
        "records": [
            {
                "structure": "Tara_RRB_model_0_A",
                "curated_reference_id": "TARA_A",
            },
            {
                "structure": "Tara_RRB_model_0_B",
                "curated_reference_id": "TARA_B",
            },
        ]
    }

    assert aliases_from_tandem_manifest(manifest) == {
        "Tara_RRB_model_0_A": "TARA_A",
        "Tara_RRB_model_0_B": "TARA_B",
    }


def test_persisted_tara_annotations_use_distinct_curated_rows(
    tmp_path, monkeypatch
) -> None:
    reference = tmp_path / "type_I_opsins.csv"
    pd.DataFrame(
        {
            "1.50": ["A1", "G1"],
            "1.51": ["C2", "T2"],
        },
        index=["TARA_A", "TARA_B"],
    ).to_csv(reference)
    processor = FakeProcessor(
        {
            "Tara_RRB_model_0_A": structure_frame("Tara_RRB_model_0_A", "AC"),
            "Tara_RRB_model_0_B": structure_frame("Tara_RRB_model_0_B", "GT"),
        }
    )
    monkeypatch.setattr(
        "src.curated_grn_storage.synchronize_runtime_reference",
        lambda source_path: {
            "source": str(source_path),
            "sha256": hashlib.sha256(source_path.read_bytes()).hexdigest(),
            "copies": [],
        },
    )

    summary = overwrite_stored_structures_with_curated_grns(
        processor,
        ["Tara_RRB_model_0_A", "Tara_RRB_model_0_B"],
        reference_path=reference,
        aliases={
            "Tara_RRB_model_0_A": "TARA_A",
            "Tara_RRB_model_0_B": "TARA_B",
        },
        output_dir=tmp_path / "output",
    )

    assert summary["annotated_structure_count"] == 2
    assert {
        row["structure_id"]: row["reference_id"] for row in summary["structures"]
    } == {
        "Tara_RRB_model_0_A": "TARA_A",
        "Tara_RRB_model_0_B": "TARA_B",
    }
    assert processor.saved_metadata["Tara_RRB_model_0_A"]["grn_reference_id"] == "TARA_A"
    assert processor.saved_metadata["Tara_RRB_model_0_B"]["grn_reference_id"] == "TARA_B"


def test_runtime_reference_copies_are_byte_identical(tmp_path) -> None:
    source = tmp_path / "source.csv"
    source.write_bytes(b",1.50\nTARA_A,A1\nTARA_B,G1\n")
    project = tmp_path / "project"

    result = synchronize_runtime_reference(source, project_root=project)

    assert result["sha256"] == hashlib.sha256(source.read_bytes()).hexdigest()
    for destination in result["copies"]:
        assert source.read_bytes() == Path(destination).read_bytes()
