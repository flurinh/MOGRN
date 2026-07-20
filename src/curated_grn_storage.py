"""Persist the authoritative ProtOS type-I GRNs on registered structures.

The structural alignment stages may calculate provisional GRNs, but the pulled
ProtOS ``type_I_opsins.csv`` table is the final curated source of truth.  This
module deliberately clears existing structure GRNs before applying that table
and saving the registered PKL entities.
"""

from __future__ import annotations

import hashlib
import json
import re
import shutil
from pathlib import Path
from typing import Any, Iterable, Mapping

import pandas as pd
from Bio.Align import PairwiseAligner


PROJECT_ROOT = Path(__file__).resolve().parent.parent
CURATED_REFERENCE = (
    PROJECT_ROOT
    / "protos"
    / "src"
    / "protos"
    / "reference_data"
    / "grn"
    / "reference"
    / "type_I_opsins.csv"
)
CURATED_STRUCTURE_ALIASES = (
    PROJECT_ROOT / "src" / "resources" / "curated_structure_aliases.json"
)
CURATED_CELL = re.compile(r"^([A-Z])(-?\d+)$")
MODIFIED_RESIDUES = {
    "LYR": "K",
    "MSE": "M",
}


def sha256_file(path: Path) -> str:
    """Return the SHA-256 digest of *path*."""

    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_curated_reference(path: Path = CURATED_REFERENCE) -> pd.DataFrame:
    """Load and minimally validate the authoritative curated table."""

    if not path.is_file():
        raise FileNotFoundError(f"Curated ProtOS GRN table not found: {path}")
    table = pd.read_csv(path, index_col=0, dtype=str).fillna("-")
    table.index = table.index.astype(str)
    table.columns = table.columns.astype(str)
    if table.index.has_duplicates:
        raise ValueError("Curated ProtOS GRN table has duplicate structure IDs")
    missing_tara = {"TARA_A", "TARA_B"}.difference(table.index)
    if missing_tara:
        raise ValueError(f"Curated ProtOS GRN table is missing {sorted(missing_tara)}")
    return table


def synchronize_runtime_reference(
    source_path: Path = CURATED_REFERENCE,
    *,
    project_root: Path = PROJECT_ROOT,
) -> dict[str, Any]:
    """Copy the authoritative table byte-for-byte to derived runtime locations."""

    source_path = source_path.resolve()
    source_digest = sha256_file(source_path)
    destinations = [
        project_root / "data" / "grn" / "reference" / "type_I_opsins.csv",
        project_root / "opsin_output" / "grn_reference.csv",
    ]
    paper_copy = (
        project_root / "opsin_output" / "paper_figures" / "type_I_opsins.csv"
    )
    if paper_copy.parent.exists():
        destinations.append(paper_copy)

    copied = []
    for destination in destinations:
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.resolve() != source_path:
            shutil.copyfile(source_path, destination)
        destination_digest = sha256_file(destination)
        if destination_digest != source_digest:
            raise RuntimeError(f"Curated-reference copy failed verification: {destination}")
        copied.append(str(destination))

    return {"source": str(source_path), "sha256": source_digest, "copies": copied}


def aliases_from_tandem_manifest(manifest: Mapping[str, Any]) -> dict[str, str]:
    """Return configured virtual-structure to curated-reference aliases."""

    aliases: dict[str, str] = {}
    for record in manifest.get("records", []):
        structure_id = str(record.get("structure", "")).strip()
        reference_id = str(record.get("curated_reference_id", "")).strip()
        if structure_id and reference_id:
            aliases[structure_id] = reference_id
    return aliases


def load_curated_structure_aliases(
    path: Path = CURATED_STRUCTURE_ALIASES,
) -> dict[str, dict[str, Any]]:
    """Load replacements whose residue numbering must be sequence-reconciled."""

    if not path.is_file():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    aliases = payload.get("aliases", {})
    if not isinstance(aliases, dict):
        raise ValueError("Curated structure aliases must be a JSON object")
    for structure_id, spec in aliases.items():
        if not isinstance(spec, dict) or not str(spec.get("reference_id", "")).strip():
            raise ValueError(f"Invalid curated structure alias for {structure_id!r}")
        if spec.get("position_mapping") not in {None, "direct", "sequence_alignment"}:
            raise ValueError(
                f"Unsupported position_mapping for curated alias {structure_id!r}"
            )
    return aliases


def _reference_mapping(row: pd.Series) -> dict[int, tuple[str, str]]:
    """Map author residue numbers to ``(GRN, expected amino acid)``."""

    mapping: dict[int, tuple[str, str]] = {}
    for grn, cell in row.items():
        match = CURATED_CELL.fullmatch(str(cell).strip())
        if not match:
            continue
        residue = match.group(1)
        sequence_number = int(match.group(2))
        if sequence_number in mapping:
            previous = mapping[sequence_number][0]
            raise ValueError(
                f"Curated residue {sequence_number} occurs at both {previous} and {grn}"
            )
        mapping[sequence_number] = (str(grn), residue)
    return mapping


def sequence_alignment_position_map(
    frame: pd.DataFrame,
    row: pd.Series,
    *,
    chain_id: str = "A",
) -> tuple[dict[int, int], dict[str, Any]]:
    """Map curated source positions to observed author positions by sequence."""

    reset = frame.reset_index()
    chain = reset.loc[
        reset["auth_chain_id"].astype(str).eq(str(chain_id))
        & reset["res_name1l"].fillna("").astype(str).str.strip().ne("")
    ].copy()
    if "atom_name" in chain.columns:
        ca = chain.loc[chain["atom_name"].astype(str).str.upper().eq("CA")]
        if not ca.empty:
            chain = ca
    chain["auth_seq_numeric"] = pd.to_numeric(chain["auth_seq_id"], errors="coerce")
    chain = (
        chain.dropna(subset=["auth_seq_numeric"])
        .sort_values("auth_seq_numeric")
        .drop_duplicates("auth_seq_numeric")
    )
    observed_positions = chain["auth_seq_numeric"].astype(int).tolist()
    observed_sequence = "".join(
        chain["res_name1l"].astype(str).str.strip().str.upper()
    )

    reference = _reference_mapping(row)
    reference_positions = sorted(reference)
    if not reference_positions:
        raise ValueError("Curated reference row has no populated positions")
    expected_range = list(
        range(reference_positions[0], reference_positions[-1] + 1)
    )
    if reference_positions != expected_range:
        raise ValueError(
            "Sequence-aligned curated aliases require a contiguous reference row"
        )
    reference_sequence = "".join(reference[position][1] for position in expected_range)

    aligner = PairwiseAligner(
        mode="global",
        match_score=2,
        mismatch_score=-2,
        open_gap_score=-5,
        extend_gap_score=-0.5,
    )
    alignment = aligner.align(reference_sequence, observed_sequence)[0]
    position_map: dict[int, int] = {}
    identities = 0
    mismatches = []
    for target_start, target_end, query_start, query_end in zip(
        alignment.coordinates[0][:-1],
        alignment.coordinates[0][1:],
        alignment.coordinates[1][:-1],
        alignment.coordinates[1][1:],
    ):
        if target_end - target_start != query_end - query_start:
            continue
        for target_index, query_index in zip(
            range(target_start, target_end),
            range(query_start, query_end),
        ):
            source_position = expected_range[target_index]
            observed_position = observed_positions[query_index]
            position_map[source_position] = observed_position
            expected = reference_sequence[target_index]
            observed = observed_sequence[query_index]
            if expected == observed or expected == "X":
                identities += 1
            else:
                mismatches.append(
                    {
                        "reference_seq_id": source_position,
                        "auth_seq_id": observed_position,
                        "expected": expected,
                        "observed": observed,
                    }
                )

    report = {
        "method": "global_sequence_alignment",
        "chain_id": str(chain_id),
        "reference_positions": len(reference_positions),
        "observed_residues": len(observed_positions),
        "mapped_positions": len(position_map),
        "identity_positions": identities,
        "mismatches": mismatches,
        "missing_reference_positions": sorted(set(reference_positions) - set(position_map)),
    }
    return position_map, report


def overwrite_frame_grns(
    frame: pd.DataFrame,
    row: pd.Series,
    *,
    structure_id: str,
    reference_id: str,
    chain_id: str = "A",
    position_map: Mapping[int, int] | None = None,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Clear provisional GRNs and apply one curated reference row."""

    reset = frame.reset_index()
    required = {"auth_chain_id", "auth_seq_id", "res_name1l"}
    missing_columns = required.difference(reset.columns)
    if missing_columns:
        raise ValueError(
            f"{structure_id} lacks columns needed for curated GRNs: "
            f"{sorted(missing_columns)}"
        )

    stale_grn_atoms = 0
    if "grn" in reset.columns:
        stale_grn_atoms = int(
            reset["grn"].fillna("").astype(str).str.strip().ne("").sum()
        )
    reset["grn"] = ""

    chain_mask = reset["auth_chain_id"].astype(str) == str(chain_id)
    numeric_ids = pd.to_numeric(reset["auth_seq_id"], errors="coerce")
    atom_column = "atom_name" if "atom_name" in reset.columns else None
    residue_view = reset.loc[chain_mask]
    if atom_column is not None:
        ca_view = residue_view[
            residue_view[atom_column].astype(str).str.upper() == "CA"
        ]
        if not ca_view.empty:
            residue_view = ca_view

    mapping = _reference_mapping(row)
    missing_residues: list[int] = []
    residue_mismatches: list[dict[str, Any]] = []
    mapped_positions = 0
    for sequence_number, (grn, expected_residue) in mapping.items():
        structure_sequence_number = (
            position_map.get(sequence_number)
            if position_map is not None
            else sequence_number
        )
        if structure_sequence_number is None:
            missing_residues.append(sequence_number)
            continue
        residue_mask = chain_mask & numeric_ids.eq(structure_sequence_number)
        if not residue_mask.any():
            missing_residues.append(sequence_number)
            continue
        observed = sorted(
            {
                str(value).strip().upper()
                for value in reset.loc[residue_mask, "res_name1l"].dropna()
                if str(value).strip()
            }
        )
        if observed == ["X"] and "res_name3l" in reset.columns:
            modified = {
                MODIFIED_RESIDUES.get(str(value).strip().upper())
                for value in reset.loc[residue_mask, "res_name3l"].dropna()
            }
            modified.discard(None)
            if modified:
                observed = sorted(modified)
        # ``X`` in a curated cell records a known source position whose
        # residue identity was unresolved; it is not a literal amino acid.
        if expected_residue != "X" and observed and observed != [expected_residue]:
            residue_mismatches.append(
                {
                    "reference_seq_id": sequence_number,
                    "auth_seq_id": structure_sequence_number,
                    "grn": grn,
                    "expected": expected_residue,
                    "observed": observed,
                }
            )
        reset.loc[residue_mask, "grn"] = grn
        mapped_positions += 1

    annotated_atoms = int(reset["grn"].astype(str).str.strip().ne("").sum())
    report = {
        "structure_id": structure_id,
        "reference_id": reference_id,
        "chain_id": str(chain_id),
        "curated_positions": len(mapping),
        "mapped_positions": mapped_positions,
        "missing_residues": missing_residues,
        "missing_residue_count": len(missing_residues),
        "residue_mismatches": residue_mismatches,
        "residue_mismatch_count": len(residue_mismatches),
        "annotated_atoms": annotated_atoms,
        "stale_grn_atoms_removed": stale_grn_atoms,
    }
    return reset, report


def overwrite_stored_structures_with_curated_grns(
    processor: Any,
    structure_ids: Iterable[str],
    *,
    reference_path: Path = CURATED_REFERENCE,
    aliases: Mapping[str, str | Mapping[str, Any]] | None = None,
    chain_id: str = "A",
    output_dir: Path | None = None,
) -> dict[str, Any]:
    """Overwrite and persist GRNs for every structure covered by the table."""

    reference_path = Path(reference_path)
    table = load_curated_reference(reference_path)
    reference_lookup = {name.casefold(): name for name in table.index}
    aliases = dict(aliases or {})
    alias_lookup = {name.casefold(): target for name, target in aliases.items()}
    reports: list[dict[str, Any]] = []
    skipped: list[str] = []

    for structure_id in dict.fromkeys(str(value) for value in structure_ids):
        alias = alias_lookup.get(structure_id.casefold())
        alias_spec: Mapping[str, Any] = {}
        if isinstance(alias, Mapping):
            alias_spec = alias
            reference_id = str(alias.get("reference_id", "")).strip()
        else:
            reference_id = alias
        if not reference_id:
            reference_id = reference_lookup.get(structure_id.casefold())
        if reference_id is None:
            skipped.append(structure_id)
            continue
        if reference_id not in table.index:
            raise ValueError(
                f"Curated alias for {structure_id!r} points to missing row "
                f"{reference_id!r}"
            )
        frame = processor.load_entity(structure_id)
        if frame is None:
            raise ValueError(f"Registered structure could not be loaded: {structure_id}")
        position_map = None
        alignment_report = None
        if alias_spec.get("position_mapping") == "sequence_alignment":
            position_map, alignment_report = sequence_alignment_position_map(
                frame,
                table.loc[reference_id],
                chain_id=chain_id,
            )
        annotated, report = overwrite_frame_grns(
            frame,
            table.loc[reference_id],
            structure_id=structure_id,
            reference_id=reference_id,
            chain_id=chain_id,
            position_map=position_map,
        )
        if alias_spec:
            report["alias"] = dict(alias_spec)
        if alignment_report is not None:
            report["position_alignment"] = alignment_report
        processor.save_entity(
            structure_id,
            annotated,
            metadata={
                "grn_source": str(reference_path.resolve()),
                "grn_source_sha256": sha256_file(reference_path),
                "grn_reference_id": reference_id,
                "grn_policy": "authoritative_curated_overwrite",
            },
        )
        reports.append(report)

    reference_sync = synchronize_runtime_reference(reference_path)
    summary = {
        "policy": "authoritative_curated_overwrite",
        "reference": reference_sync,
        "annotated_structure_count": len(reports),
        "skipped_structure_count": len(skipped),
        "skipped_structures": skipped,
        "structures": reports,
    }
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        report_path = output_dir / "curated_structure_grn_overwrite.json"
        report_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        pd.DataFrame(reports).to_csv(
            output_dir / "curated_structure_grn_overwrite.csv", index=False
        )
    return summary
