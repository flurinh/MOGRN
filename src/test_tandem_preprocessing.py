"""Focused tests for tandem-rhodopsin sequence preprocessing."""

from __future__ import annotations

import ast
import hashlib
import json
import random
import resource
import signal
import subprocess
import sys
from pathlib import Path
from urllib.parse import unquote

import pytest
import src.tandem_preprocessing as tandem_preprocessing

from src.tandem_preprocessing import (
    CuratedBoundaries,
    FastaRecord,
    MAX_AUTOMATIC_LENGTH,
    RepeatEvidence,
    cli_main,
    detect_repeat_evidence,
    load_curated_boundaries,
    preprocess_record,
    preprocess_records,
    publish_output_pair,
    read_fasta,
    sequence_complexity_rejection,
    stable_virtual_id,
    write_fasta,
)


KV_FIXTURE = Path(__file__).with_name("test_data") / "kv_rrb_regression.fasta"


def confident_evidence(**changes: object) -> RepeatEvidence:
    values = {
        "repeat_score": 150.0,
        "aligned_residues": 200,
        "repeat_identity": 0.35,
        "core_1_start": 50,
        "core_1_end": 250,
        "core_2_start": 350,
        "core_2_end": 550,
        "test_split": 300,
        "repeat_period": 300,
    }
    values.update(changes)
    return RepeatEvidence(**values)


def realistic_tandem_sequence() -> str:
    """Return deterministic, diverse homologous domains with realistic noise."""

    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    weights = [8, 2, 5, 6, 4, 7, 2, 6, 6, 9, 2, 5, 4, 5, 6, 5, 5, 7, 1, 4]
    rng = random.Random(17)
    domain_1 = "".join(rng.choices(amino_acids, weights=weights, k=250))
    domain_2 = "".join(
        amino_acids.replace(residue, "")[(position * 7) % 19]
        if position % 5 == 0
        else residue
        for position, residue in enumerate(domain_1)
    )
    flank = ("MSTNPKPQRITF" * 3)[:30]
    linker = "GQSPNATDKL" * 4
    return flank + domain_1 + linker + domain_2 + flank


def kv_fixture_records() -> dict[str, FastaRecord]:
    records: dict[str, FastaRecord] = {}
    for record in read_fasta(KV_FIXTURE):
        metadata = dict(
            field.split("=", 1)
            for field in record.record_id.split("|")[1:]
        )
        records[metadata["segment"]] = record
    return records


def limit_subprocess_file_size() -> None:
    """Make a manifest temp write fail after the smaller FASTA temp succeeds."""

    signal.signal(signal.SIGXFSZ, signal.SIG_IGN)
    resource.setrlimit(resource.RLIMIT_FSIZE, (4096, 4096))


def test_499_bypasses_without_calling_detector() -> None:
    calls = 0

    def detector(_sequence: str) -> None:
        nonlocal calls
        calls += 1

    result = preprocess_record("short", "A" * 499, detector=detector)

    assert calls == 0
    assert result.records == (FastaRecord("short", "A" * 499, ""),)
    assert result.manifest["status"] == "bypassed_too_short"
    assert result.manifest["eligible"] is False


def test_500_is_eligible_and_runs_detector() -> None:
    calls = 0

    def detector(sequence: str) -> None:
        nonlocal calls
        calls += 1
        assert len(sequence) == 500

    result = preprocess_record("boundary", "A" * 500, detector=detector)

    assert calls == 1
    assert result.manifest["eligible"] is True
    assert result.manifest["status"] == "pass_through_no_repeat"


def test_production_detector_accepts_realistic_synthetic_positive() -> None:
    sequence = realistic_tandem_sequence()

    evidence = detect_repeat_evidence(sequence)
    result = preprocess_record("realistic", sequence)

    assert evidence is not None
    assert evidence.ambiguous is False
    assert evidence.repeat_score >= 150
    assert evidence.aligned_residues >= 200
    assert evidence.repeat_identity >= 0.35
    assert result.manifest["status"] == "split_detected"
    assert len(result.records) == 2


def test_real_kv_fixture_accession_bounds_hashes_and_core_slices() -> None:
    expected = {
        "full": {
            "bounds": "1-1199",
            "length": 1199,
            "sha256": (
                "0e26463c085a28608660450e9436ce079b9da009f815ff299af57b66db2193c6"
            ),
        },
        "R1": {
            "bounds": "85-359",
            "length": 275,
            "sha256": (
                "fa62ea21c3f448ddb9588525b0c8b6bc4feba6bbb1bbccef0ecc8aa7901b3c04"
            ),
        },
        "R2": {
            "bounds": "472-749",
            "length": 278,
            "sha256": (
                "e407a9c1ed0f7d79234af6867ef82f2c351556e48d7a6d2e5a25b8c7afd3a013"
            ),
        },
    }
    assert hashlib.sha256(KV_FIXTURE.read_bytes()).hexdigest() == (
        "31beb248ee11319c552f6680738025804a1d0fa0d828aef39ad6d4a3ebca15e4"
    )
    records = kv_fixture_records()

    assert set(records) == set(expected)
    for segment, record in records.items():
        identifier, *fields = record.record_id.split("|")
        metadata = dict(field.split("=", 1) for field in fields)
        assert identifier == "Kv-RRB"
        assert metadata["protein"] == "USH44361.1"
        assert metadata["nucleotide"] == "MZ740268.1"
        assert metadata["coordinate_system"] == "protein_1based_inclusive"
        assert metadata["bounds"] == expected[segment]["bounds"]
        assert int(metadata["length"]) == expected[segment]["length"]
        assert len(record.sequence) == expected[segment]["length"]
        assert metadata["sha256"] == expected[segment]["sha256"]
        assert hashlib.sha256(record.sequence.encode()).hexdigest() == (
            expected[segment]["sha256"]
        )

    parent = records["full"].sequence
    assert parent[84:359] == records["R1"].sequence
    assert parent[471:749] == records["R2"].sequence


def test_real_kv_default_automatic_evidence_regression() -> None:
    parent = kv_fixture_records()["full"].sequence

    evidence = detect_repeat_evidence(parent)
    result = preprocess_record("Kv-RRB", parent)

    assert evidence is not None
    assert evidence.repeat_score == 157.5
    assert evidence.aligned_residues == 292
    assert evidence.repeat_identity == pytest.approx(0.5273972602739726)
    assert evidence.repeat_period == 387
    assert evidence.test_split == 450
    assert evidence.ambiguous is False
    assert (
        evidence.core_1_start,
        evidence.core_1_end,
        evidence.core_2_start,
        evidence.core_2_end,
    ) == (87, 378, 474, 768)
    assert result.manifest["status"] == "split_detected"
    assert [segment["parent_interval"] for segment in result.manifest["segments"]] == [
        [87, 426],
        [427, 768],
    ]


def test_real_kv_curated_tail_exclusion_and_exact_reconstruction() -> None:
    records = kv_fixture_records()
    parent = records["full"].sequence
    result = preprocess_record(
        "Kv-RRB",
        parent,
        description="USH44361.1 MZ740268.1 curated tandem",
        curated_boundaries=CuratedBoundaries.from_intervals((85, 359), (472, 749)),
    )

    assert result.manifest["status"] == "split_curated"
    assert result.manifest["ownership_cut"] == {
        "left_parent_position": 415,
        "right_parent_position": 416,
    }
    assert [segment["parent_interval"] for segment in result.manifest["segments"]] == [
        [85, 415],
        [416, 749],
    ]
    assert result.records[1].sequence == parent[415:749]
    assert result.records[1].sequence.endswith(records["R2"].sequence[-1])
    assert result.manifest["residual_regions"]["c_terminal"] == {
        "parent_interval": [750, 1199],
        "length": 450,
        "sequence": parent[749:],
    }

    reconstructed: list[str | None] = [None] * len(parent)
    by_id = {record.record_id: record for record in result.records}
    for residual in result.manifest["residual_regions"].values():
        if residual is not None:
            start, end = residual["parent_interval"]
            reconstructed[start - 1 : end] = residual["sequence"]
    for segment in result.manifest["segments"]:
        virtual = by_id[segment["virtual_id"]]
        for mapping in segment["local_to_parent"]:
            reconstructed[mapping["parent_position"] - 1] = virtual.sequence[
                mapping["local_position"] - 1
            ]
    assert all(residue is not None for residue in reconstructed)
    assert "".join(residue or "" for residue in reconstructed) == parent


@pytest.mark.parametrize("sequence", ["A" * 500, "AC" * 300, "AGP" * 200])
def test_production_detector_rejects_obvious_low_complexity(sequence: str) -> None:
    assert sequence_complexity_rejection(sequence) is not None
    assert detect_repeat_evidence(sequence) is None
    assert preprocess_record("low-complexity", sequence).manifest["status"] == (
        "pass_through_no_repeat"
    )


def test_production_detector_marks_native_multi_repeat_placement_ambiguous() -> None:
    amino_acids = "ACDEFGHIKLMNPQRSTVWY"
    base = "".join(random.Random(41).choices(amino_acids, k=200))
    sequence = base * 3

    evidence = detect_repeat_evidence(sequence)
    result = preprocess_record("three-copy", sequence)

    assert evidence is not None
    assert evidence.ambiguous is True
    assert "materially different near-tied" in (evidence.ambiguity_reason or "")
    assert result.manifest["status"] == "pass_through_ambiguous_evidence"
    assert len(result.records) == 1


def test_production_detector_does_not_split_unrelated_diverse_sequence() -> None:
    sequence = "".join(
        random.Random(73).choices("ACDEFGHIKLMNPQRSTVWY", k=700)
    )

    result = preprocess_record("diverse-negative", sequence)

    assert not result.manifest["status"].startswith("split_")


def test_bounded_automatic_detection_allows_curated_override() -> None:
    sequence = realistic_tandem_sequence() * 5
    assert len(sequence) > MAX_AUTOMATIC_LENGTH

    automatic = preprocess_record("long", sequence)
    curated = preprocess_record(
        "Kv-RRB-long",
        sequence,
        curated_boundaries=CuratedBoundaries.from_intervals((85, 359), (472, 749)),
    )

    assert automatic.manifest["status"] == "pass_through_automatic_length_limit"
    assert curated.manifest["status"] == "split_curated"
    assert len(curated.records) == 2


def test_nonhit_passes_through_unchanged() -> None:
    sequence = "aCdE" * 150
    result = preprocess_record("nonhit", sequence, detector=lambda _sequence: None)

    assert len(result.records) == 1
    assert result.records[0].record_id == "nonhit"
    assert result.records[0].sequence == sequence
    assert result.manifest["status"] == "pass_through_no_repeat"


def test_tab_header_pass_through_and_read_write_round_trip(tmp_path: Path) -> None:
    input_path = tmp_path / "tabbed.fasta"
    output_path = tmp_path / "roundtrip.fasta"
    description = "primary description\tsecondary field"
    input_path.write_text(">tabbed-id\t" + description + "\n" + "aCdE" * 100 + "\n")

    parsed = read_fasta(input_path)
    output, manifest = preprocess_records(parsed)
    write_fasta(output, output_path)
    reparsed = read_fasta(output_path)

    assert parsed == [FastaRecord("tabbed-id", "aCdE" * 100, description)]
    assert output == parsed
    assert reparsed == parsed
    assert manifest["records"][0]["parent_id"] == "tabbed-id"
    assert manifest["records"][0]["parent_description"] == description


def test_tab_header_curated_split_uses_clean_id_and_description(tmp_path: Path) -> None:
    input_path = tmp_path / "tabbed.fasta"
    description = "curated Kv metadata\taccession=USH44361.1"
    input_path.write_text(">Kv-tab\t" + description + "\n" + "A" * 800 + "\n")

    output, manifest = preprocess_records(
        read_fasta(input_path),
        curated_boundaries={
            "Kv-tab": CuratedBoundaries.from_intervals((85, 359), (472, 749))
        },
    )

    assert [record.record_id for record in output] == ["Kv-tab__R1", "Kv-tab__R2"]
    assert manifest["records"][0]["parent_id"] == "Kv-tab"
    assert manifest["records"][0]["parent_description"] == description
    assert all(
        segment["parent_description"] == description
        for segment in manifest["records"][0]["segments"]
    )


def test_tab_header_generated_id_collision_is_detected(tmp_path: Path) -> None:
    input_path = tmp_path / "tabbed-collision.fasta"
    input_path.write_text(
        ">sample\tparent description\n"
        + "A" * 800
        + "\n>sample__R1\tordinary child\n"
        + "C" * 100
        + "\n"
    )

    with pytest.raises(ValueError, match="output FASTA ID collision.*sample__R1"):
        preprocess_records(
            read_fasta(input_path),
            curated_boundaries={
                "sample": CuratedBoundaries.from_intervals(
                    (85, 359), (472, 749)
                )
            },
        )


def test_confident_hit_splits_at_exact_midpoint_with_reversible_mapping() -> None:
    sequence = "ACDEFGHIKLMNPQRSTVWY" * 30
    evidence = confident_evidence(
        core_1_start=51,
        core_1_end=240,
        core_2_start=361,
        core_2_end=550,
    )
    result = preprocess_record(
        "detected", sequence, detector=lambda _sequence: evidence
    )

    assert result.manifest["status"] == "split_detected"
    assert result.manifest["ownership_cut"] == {
        "left_parent_position": 300,
        "right_parent_position": 301,
    }
    assert [record.record_id for record in result.records] == [
        "detected__R1",
        "detected__R2",
    ]
    assert result.records[0].sequence == sequence[50:300]
    assert result.records[1].sequence == sequence[300:550]
    assert "|" not in result.records[0].sequence + result.records[1].sequence

    first_mapping = result.manifest["segments"][0]["local_to_parent"]
    second_mapping = result.manifest["segments"][1]["local_to_parent"]
    assert first_mapping[0] == {"local_position": 1, "parent_position": 51}
    assert first_mapping[-1] == {"local_position": 250, "parent_position": 300}
    assert second_mapping[0] == {"local_position": 1, "parent_position": 301}
    assert second_mapping[-1] == {"local_position": 250, "parent_position": 550}


def test_curated_kv_boundaries_yield_415_416_cut() -> None:
    sequence = "A" * 800
    override = CuratedBoundaries.from_intervals((85, 359), (472, 749))
    result = preprocess_record(
        "Kv-RRB",
        sequence,
        curated_boundaries=override,
        detector=lambda _sequence: (_ for _ in ()).throw(
            AssertionError("curated override should not need detection")
        ),
    )

    assert result.manifest["status"] == "split_curated"
    assert result.manifest["ownership_cut"] == {
        "left_parent_position": 415,
        "right_parent_position": 416,
    }
    assert [segment["parent_interval"] for segment in result.manifest["segments"]] == [
        [85, 415],
        [416, 749],
    ]
    assert [len(record.sequence) for record in result.records] == [331, 334]
    assert result.manifest["boundary_method"] == "curated_override"
    assert result.manifest["boundary_confidence"] == "curated"


def test_invalid_curated_boundaries_fail_closed() -> None:
    sequence = "A" * 800
    invalid = CuratedBoundaries.from_intervals((85, 500), (472, 749))

    result = preprocess_record("Kv-RRB", sequence, curated_boundaries=invalid)

    assert len(result.records) == 1
    assert result.records[0].record_id == "Kv-RRB"
    assert result.records[0].sequence == sequence
    assert result.manifest["status"] == "pass_through_invalid_curated_boundaries"
    assert "nonoverlapping" in result.manifest["error"]


def test_invalid_or_ambiguous_detected_evidence_passes_through() -> None:
    sequence = "A" * 800
    invalid = confident_evidence(core_1_end=500, core_2_start=472)
    invalid_result = preprocess_record(
        "invalid", sequence, detector=lambda _sequence: invalid
    )
    ambiguous = confident_evidence(
        ambiguous=True, ambiguity_reason="two incompatible repeat placements"
    )
    ambiguous_result = preprocess_record(
        "ambiguous", sequence, detector=lambda _sequence: ambiguous
    )

    assert invalid_result.manifest["status"] == (
        "pass_through_invalid_detected_boundaries"
    )
    assert ambiguous_result.manifest["status"] == (
        "pass_through_ambiguous_evidence"
    )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("repeat_score", float("nan")),
        ("repeat_score", float("inf")),
        ("repeat_score", True),
        ("aligned_residues", 200.5),
        ("aligned_residues", False),
        ("repeat_identity", 1.01),
        ("core_1_start", 0),
        ("core_2_end", 801),
        ("ambiguous", "false"),
        ("test_split", 900),
        ("repeat_period", -1),
    ],
)
def test_malformed_mapping_evidence_fails_closed(field: str, value: object) -> None:
    evidence: dict[str, object] = {
        "repeat_score": 150,
        "aligned_residues": 200,
        "repeat_identity": 0.35,
        "core_1_start": 85,
        "core_1_end": 359,
        "core_2_start": 472,
        "core_2_end": 749,
        "test_split": 415,
        "repeat_period": 387,
        "ambiguous": False,
    }
    evidence[field] = value

    result = preprocess_record(
        "malformed", "A" * 800, detector=lambda _sequence: evidence
    )

    assert result.manifest["status"] == "pass_through_detector_error"
    assert "error" in result.manifest
    assert len(result.records) == 1


def test_numeric_mapping_evidence_is_coerced_and_validated() -> None:
    evidence = {
        "repeat_score": "150",
        "aligned_residues": "200",
        "repeat_identity": "0.35",
        "core_1_start": "85",
        "core_1_end": "359",
        "core_2_start": "472",
        "core_2_end": "749",
        "test_split": "415",
        "repeat_period": "387",
        "ambiguous": False,
    }

    result = preprocess_record(
        "coerced", "A" * 800, detector=lambda _sequence: evidence
    )

    assert result.manifest["status"] == "split_detected"
    assert result.manifest["evidence"]["repeat_score"] == 150.0


def test_outer_residuals_are_preserved_and_tail_is_not_in_r2() -> None:
    sequence = "N" * 84 + "C" * 665 + "T" * 51
    override = CuratedBoundaries.from_intervals((85, 359), (472, 749))
    result = preprocess_record("Kv-tail", sequence, curated_boundaries=override)

    residuals = result.manifest["residual_regions"]
    assert residuals["n_terminal"] == {
        "parent_interval": [1, 84],
        "length": 84,
        "sequence": "N" * 84,
    }
    assert residuals["c_terminal"] == {
        "parent_interval": [750, 800],
        "length": 51,
        "sequence": "T" * 51,
    }
    assert result.records[1].sequence == sequence[415:749]
    assert "T" not in result.records[1].sequence


def test_virtual_ids_are_stable_and_retain_parent_and_domain() -> None:
    assert stable_virtual_id("Kv-RRB", "R1") == "Kv-RRB__R1"
    assert stable_virtual_id("Kv-RRB", "R2") == "Kv-RRB__R2"
    first = preprocess_record(
        "Kv-RRB",
        "A" * 800,
        curated_boundaries=CuratedBoundaries.from_intervals((85, 359), (472, 749)),
    )
    second = preprocess_record(
        "Kv-RRB",
        "A" * 800,
        curated_boundaries=CuratedBoundaries.from_intervals((85, 359), (472, 749)),
    )
    assert [record.record_id for record in first.records] == [
        record.record_id for record in second.records
    ]
    assert all(
        segment["parent_id"] == "Kv-RRB"
        for segment in first.manifest["segments"]
    )


def test_split_preserves_parent_description_in_manifest_and_output_metadata() -> None:
    description = "original accession metadata with spaces"
    result = preprocess_record(
        "described",
        realistic_tandem_sequence(),
        description=description,
    )

    assert result.manifest["parent_description"] == description
    for record, segment in zip(result.records, result.manifest["segments"]):
        assert segment["parent_description"] == description
        encoded = record.description.split("parent_description=", 1)[1]
        assert unquote(encoded) == description


def test_manifest_and_virtual_fasta_fully_reconstruct_parent_sequence() -> None:
    sequence = ("ACDEFGHIKLMNPQRSTVWY" * 60)[:1200]
    result = preprocess_record(
        "Kv-RRB",
        sequence,
        description="known Kv tandem with terminal tail",
        curated_boundaries=CuratedBoundaries.from_intervals((85, 359), (472, 749)),
    )
    reconstructed: list[str | None] = [None] * len(sequence)
    records_by_id = {record.record_id: record for record in result.records}

    for residual in result.manifest["residual_regions"].values():
        if residual is None:
            continue
        start, end = residual["parent_interval"]
        reconstructed[start - 1 : end] = residual["sequence"]
    for segment in result.manifest["segments"]:
        virtual = records_by_id[segment["virtual_id"]]
        for mapping in segment["local_to_parent"]:
            reconstructed[mapping["parent_position"] - 1] = virtual.sequence[
                mapping["local_position"] - 1
            ]

    assert "".join(residue or "" for residue in reconstructed) == sequence
    assert all(residue is not None for residue in reconstructed)


def test_batch_rejects_generated_and_pass_through_id_collision() -> None:
    records = [
        FastaRecord("sample", "A" * 800),
        FastaRecord("sample__R1", "C" * 100),
    ]
    overrides = {
        "sample": CuratedBoundaries.from_intervals((85, 359), (472, 749))
    }

    with pytest.raises(ValueError, match="output FASTA ID collision.*sample__R1"):
        preprocess_records(records, curated_boundaries=overrides)


def test_unmatched_override_is_rejected_before_detector_runs() -> None:
    calls = 0

    def detector(_sequence: str) -> None:
        nonlocal calls
        calls += 1

    with pytest.raises(ValueError, match="not found in FASTA: typo-id"):
        preprocess_records(
            [FastaRecord("present", "A" * 800)],
            detector=detector,
            curated_boundaries={
                "typo-id": CuratedBoundaries.from_intervals(
                    (85, 359), (472, 749)
                )
            },
        )
    assert calls == 0


@pytest.mark.parametrize(
    "payload",
    [
        {"Kv-RRB": 42},
        {"Kv-RRB": {"core_1": [85], "core_2": [472, 749]}},
        {"Kv-RRB": {"core_1": [85, "359"], "core_2": [472, 749]}},
        {"Kv-RRB": {"core_1": [85, 359], "core_2": [472, 749], "confidence": False}},
    ],
)
def test_malformed_curated_override_errors_are_contextual(
    tmp_path: Path,
    payload: dict[str, object],
) -> None:
    path = tmp_path / "bad-overrides.json"
    path.write_text(json.dumps(payload))

    with pytest.raises(ValueError, match="invalid curated override for 'Kv-RRB'"):
        load_curated_boundaries(path)


def test_malformed_curated_json_error_is_normalized(tmp_path: Path) -> None:
    path = tmp_path / "bad-overrides.json"
    path.write_text('{"Kv-RRB":')

    with pytest.raises(ValueError, match="invalid curated-boundary JSON at line"):
        load_curated_boundaries(path)


def test_cli_writes_expanded_fasta_and_json_manifest(tmp_path: Path) -> None:
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "expanded.fasta"
    manifest_path = tmp_path / "manifest.json"
    overrides_path = tmp_path / "overrides.json"
    input_fasta.write_text(">Kv-RRB known record\n" + "A" * 800 + "\n")
    overrides_path.write_text(
        json.dumps({"Kv-RRB": {"core_1": [85, 359], "core_2": [472, 749]}})
    )

    exit_code = cli_main(
        [
            str(input_fasta),
            str(output_fasta),
            str(manifest_path),
            "--curated-boundaries",
            str(overrides_path),
        ]
    )

    assert exit_code == 0
    output = output_fasta.read_text()
    assert output.count(">") == 2
    assert ">Kv-RRB__R1 " in output
    assert ">Kv-RRB__R2 " in output
    manifest = json.loads(manifest_path.read_text())
    assert manifest["schema"] == "mogrn.tandem-preprocessing/v1"
    assert manifest["records"][0]["status"] == "split_curated"


@pytest.mark.parametrize(
    "option",
    [
        ["--scan-step", "0"],
        ["--scan-step", "-1"],
        ["--min-repeat-score", "nan"],
        ["--min-repeat-score", "0"],
        ["--min-aligned-residues", "0"],
        ["--min-repeat-identity", "0"],
        ["--min-repeat-identity", "1.1"],
    ],
)
def test_cli_rejects_invalid_global_configuration(
    tmp_path: Path,
    option: list[str],
) -> None:
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "output.fasta"
    manifest = tmp_path / "manifest.json"
    input_fasta.write_text(">record\n" + realistic_tandem_sequence() + "\n")

    with pytest.raises(SystemExit) as error:
        cli_main(
            [str(input_fasta), str(output_fasta), str(manifest), *option]
        )

    assert error.value.code == 2
    assert not output_fasta.exists()
    assert not manifest.exists()


@pytest.mark.parametrize("alias", ["outputs", "input-output", "input-manifest"])
def test_cli_rejects_aliased_paths_without_overwriting(
    tmp_path: Path,
    alias: str,
) -> None:
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "output.fasta"
    manifest = tmp_path / "manifest.json"
    original = ">record\n" + "A" * 499 + "\n"
    input_fasta.write_text(original)
    if alias == "outputs":
        manifest = output_fasta
    elif alias == "input-output":
        output_fasta = input_fasta
    else:
        manifest = input_fasta

    with pytest.raises(SystemExit) as error:
        cli_main([str(input_fasta), str(output_fasta), str(manifest)])

    assert error.value.code == 2
    assert input_fasta.read_text() == original


def test_cli_requires_one_output_destination_directory(tmp_path: Path) -> None:
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "fasta" / "output.fasta"
    manifest = tmp_path / "metadata" / "manifest.json"
    input_fasta.write_text(">record\n" + "A" * 499 + "\n")

    with pytest.raises(SystemExit) as error:
        cli_main([str(input_fasta), str(output_fasta), str(manifest)])

    assert error.value.code == 2
    assert not output_fasta.exists()
    assert not manifest.exists()


def test_cli_collision_fails_before_writing_outputs(tmp_path: Path) -> None:
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "output.fasta"
    manifest = tmp_path / "manifest.json"
    overrides = tmp_path / "overrides.json"
    input_fasta.write_text(
        ">sample\n" + "A" * 800 + "\n>sample__R1\n" + "C" * 100 + "\n"
    )
    overrides.write_text(
        json.dumps({"sample": {"core_1": [85, 359], "core_2": [472, 749]}})
    )

    with pytest.raises(SystemExit) as error:
        cli_main(
            [
                str(input_fasta),
                str(output_fasta),
                str(manifest),
                "--curated-boundaries",
                str(overrides),
            ]
        )

    assert error.value.code == 2
    assert not output_fasta.exists()
    assert not manifest.exists()


def test_executable_cli_handles_production_positive_and_low_complexity(
    tmp_path: Path,
) -> None:
    root = Path(__file__).resolve().parent.parent
    script = root / "scripts" / "preprocess_tandem_rhodopsins.py"
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "output.fasta"
    manifest_path = tmp_path / "manifest.json"
    input_fasta.write_text(
        ">positive original metadata\n"
        + realistic_tandem_sequence()
        + "\n>low_complexity\n"
        + "A" * 500
        + "\n"
    )

    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            str(input_fasta),
            str(output_fasta),
            str(manifest_path),
        ],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 0, completed.stderr
    output_records = output_fasta.read_text()
    assert output_records.count(">") == 3
    assert ">positive__R1 " in output_records
    assert ">positive__R2 " in output_records
    assert ">low_complexity\n" in output_records
    entries = {
        entry["parent_id"]: entry
        for entry in json.loads(manifest_path.read_text())["records"]
    }
    assert entries["positive"]["status"] == "split_detected"
    assert entries["low_complexity"]["status"] == "pass_through_no_repeat"


def test_executable_cli_returns_two_for_invalid_configuration(tmp_path: Path) -> None:
    root = Path(__file__).resolve().parent.parent
    script = root / "scripts" / "preprocess_tandem_rhodopsins.py"
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "output.fasta"
    manifest_path = tmp_path / "manifest.json"
    input_fasta.write_text(">record\n" + realistic_tandem_sequence() + "\n")

    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            str(input_fasta),
            str(output_fasta),
            str(manifest_path),
            "--scan-step",
            "0",
        ],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 2
    assert "scan step must be a positive integer" in completed.stderr
    assert not output_fasta.exists()
    assert not manifest_path.exists()


@pytest.mark.parametrize("preexisting", [False, True])
def test_executable_serialization_failure_leaves_pair_untouched(
    tmp_path: Path,
    preexisting: bool,
) -> None:
    root = Path(__file__).resolve().parent.parent
    script = root / "scripts" / "preprocess_tandem_rhodopsins.py"
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "output.fasta"
    manifest_path = tmp_path / "manifest.json"
    input_fasta.write_text(">record\n" + "A" * 499 + "\n")
    if preexisting:
        output_fasta.write_text("original FASTA artifact\n")
        manifest_path.write_text("original manifest artifact\n")

    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            str(input_fasta),
            str(output_fasta),
            str(manifest_path),
        ],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
        preexec_fn=limit_subprocess_file_size,
    )

    assert completed.returncode == 2, completed.stderr
    if preexisting:
        assert output_fasta.read_text() == "original FASTA artifact\n"
        assert manifest_path.read_text() == "original manifest artifact\n"
    else:
        assert not output_fasta.exists()
        assert not manifest_path.exists()
    assert not list(tmp_path.glob(".output.fasta.*"))
    assert not list(tmp_path.glob(".manifest.json.*"))


@pytest.mark.parametrize("preexisting_fasta", [False, True])
def test_executable_second_publish_failure_rolls_back_pair(
    tmp_path: Path,
    preexisting_fasta: bool,
) -> None:
    root = Path(__file__).resolve().parent.parent
    script = root / "scripts" / "preprocess_tandem_rhodopsins.py"
    input_fasta = tmp_path / "input.fasta"
    output_fasta = tmp_path / "output.fasta"
    manifest_path = tmp_path / "manifest.json"
    input_fasta.write_text(">record\n" + "A" * 499 + "\n")
    if preexisting_fasta:
        output_fasta.write_text("original FASTA artifact\n")
    manifest_path.mkdir()
    sentinel = manifest_path / "sentinel"
    sentinel.write_text("do not modify\n")

    completed = subprocess.run(
        [
            sys.executable,
            str(script),
            str(input_fasta),
            str(output_fasta),
            str(manifest_path),
        ],
        cwd=root,
        text=True,
        capture_output=True,
        check=False,
    )

    assert completed.returncode == 2, completed.stderr
    if preexisting_fasta:
        assert output_fasta.read_text() == "original FASTA artifact\n"
    else:
        assert not output_fasta.exists()
    assert manifest_path.is_dir()
    assert sentinel.read_text() == "do not modify\n"
    assert not list(tmp_path.glob(".output.fasta.*"))
    assert not list(tmp_path.glob(".manifest.json.*"))


def test_second_publish_failure_restores_both_preexisting_files(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_fasta = tmp_path / "output.fasta"
    manifest_path = tmp_path / "manifest.json"
    output_fasta.write_text("original FASTA artifact\n")
    manifest_path.write_text("original manifest artifact\n")
    original_replace = tandem_preprocessing.os.replace

    def fail_manifest_temp_publish(source: object, destination: object) -> None:
        source_path = Path(source)
        destination_path = Path(destination)
        if source_path.suffix == ".tmp" and destination_path == manifest_path:
            raise OSError("forced second publish failure")
        original_replace(source, destination)

    monkeypatch.setattr(
        tandem_preprocessing.os,
        "replace",
        fail_manifest_temp_publish,
    )

    with pytest.raises(OSError, match="forced second publish failure"):
        publish_output_pair(
            [FastaRecord("record", "A" * 499)],
            {"schema": "test"},
            output_fasta,
            manifest_path,
        )

    assert output_fasta.read_text() == "original FASTA artifact\n"
    assert manifest_path.read_text() == "original manifest artifact\n"
    assert not list(tmp_path.glob(".output.fasta.*"))
    assert not list(tmp_path.glob(".manifest.json.*"))


def test_failed_restore_retains_exact_original_backup(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    output_fasta = tmp_path / "output.fasta"
    manifest_path = tmp_path / "manifest.json"
    original_fasta = b"original FASTA artifact\n"
    original_manifest = b"original manifest artifact\n"
    output_fasta.write_bytes(original_fasta)
    manifest_path.write_bytes(original_manifest)
    original_replace = tandem_preprocessing.os.replace

    def fail_manifest_publish_and_restore(
        source: object,
        destination: object,
    ) -> None:
        source_path = Path(source)
        destination_path = Path(destination)
        if source_path.suffix == ".tmp" and destination_path == manifest_path:
            raise OSError("forced second publish failure")
        if source_path.suffix == ".backup" and destination_path == manifest_path:
            raise OSError("forced manifest restore failure")
        original_replace(source, destination)

    monkeypatch.setattr(
        tandem_preprocessing.os,
        "replace",
        fail_manifest_publish_and_restore,
    )

    with pytest.raises(RuntimeError, match="rollback was incomplete"):
        publish_output_pair(
            [FastaRecord("record", "A" * 499)],
            {"schema": "test"},
            output_fasta,
            manifest_path,
        )

    assert output_fasta.read_bytes() == original_fasta
    assert not manifest_path.exists()
    retained_backups = list(tmp_path.glob(".manifest.json.*.backup"))
    assert len(retained_backups) == 1
    assert retained_backups[0].read_bytes() == original_manifest
    assert not list(tmp_path.glob(".output.fasta.*"))
    assert not list(tmp_path.glob(".manifest.json.*.tmp"))


def test_new_preprocessing_code_has_no_protos_dependency() -> None:
    root = Path(__file__).resolve().parent.parent
    paths = [
        root / "src" / "tandem_preprocessing.py",
        root / "scripts" / "preprocess_tandem_rhodopsins.py",
    ]
    for path in paths:
        tree = ast.parse(path.read_text())
        imported_names = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        }
        imported_names.update(
            node.module or ""
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
        )
        assert not any(
            name == "protos" or name.startswith("protos.")
            for name in imported_names
        )
