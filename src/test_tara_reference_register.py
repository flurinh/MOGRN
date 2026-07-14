"""Regression tests for the hand-curated TARA TM5 register."""

from __future__ import annotations

import hashlib
import json

import pandas as pd

import scripts.build_tara_reference_rows as tara_builder
from scripts.audit_tara_tm5_consensus import HYPOTHESES
from scripts.build_tara_reference_rows import (
    TARA_A_TM5_REGISTER,
    apply_tara_a_tm5_register,
)


def test_tara_a_tm5_register_closes_internal_gap_and_keeps_insertion() -> None:
    automatic = pd.Series(
        {
            "5.40": "W184",
            "5.41": "M185",
            "5.42": "W186",
            "5.43": "Y187",
            "5.44": "A188",
            "5.45": "I189",
            "5.451": "-",
            "5.46": "F190",
        }
    )

    corrected = apply_tara_a_tm5_register(automatic)

    assert corrected[list(TARA_A_TM5_REGISTER)].to_dict() == TARA_A_TM5_REGISTER
    assert corrected["5.40"] == "-"
    assert corrected["5.43"] == "W186"
    assert corrected["5.451"] == "I189"
    assert [value for value in corrected if value != "-"] == [
        "W184",
        "M185",
        "W186",
        "Y187",
        "A188",
        "I189",
        "F190",
    ]


def test_curated_register_remains_an_explicit_structural_audit_hypothesis() -> None:
    assert HYPOTHESES["insertion_gap_5.40"] == {
        grn: int(residue[1:])
        for grn, residue in TARA_A_TM5_REGISTER.items()
        if residue != "-"
    }


def test_runtime_sync_restores_table_and_checksum_metadata(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(tara_builder, "ROOT", tmp_path)
    manifest = {
        "bundle_version": "old",
        "files": {"type_I_opsins.csv": "old"},
    }
    provenance = {
        "type_I": {"policy": "old", "source_table_sha256": "old"},
        "files": {"type_I_opsins.csv": "old"},
    }
    manifest_path = tmp_path / "data" / "grn" / "manifest.json"
    provenance_path = tmp_path / "data" / "grn" / "opsin_provenance.json"
    manifest_path.parent.mkdir(parents=True)
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    provenance_path.write_text(json.dumps(provenance), encoding="utf-8")
    candidate = pd.DataFrame({"5.451": ["I189", "-"]}, index=["TARA_A", "TARA_B"])

    digest = tara_builder.sync_runtime_reference(candidate)

    runtime = tmp_path / "data" / "grn" / "reference" / "type_I_opsins.csv"
    assert digest == hashlib.sha256(runtime.read_bytes()).hexdigest()
    assert pd.read_csv(runtime, index_col=0, dtype=str).fillna("-").equals(candidate)
    updated_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    updated_provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert updated_manifest["files"]["type_I_opsins.csv"] == digest
    assert updated_provenance["type_I"]["source_table_sha256"] == digest
    assert updated_provenance["files"]["type_I_opsins.csv"] == digest
