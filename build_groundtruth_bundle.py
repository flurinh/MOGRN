#!/usr/bin/env python3
"""Build the canonical, self-contained GRN opsin ground-truth Zenodo bundle."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import io
import json
import re
import shutil
import stat
import subprocess
import sys
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd


ROOT = Path(__file__).resolve().parent
PROTOS_SRC = ROOT / "protos" / "src"
if str(PROTOS_SRC) not in sys.path:
    sys.path.insert(0, str(PROTOS_SRC))

from src.curated_grn_storage import (  # noqa: E402
    CURATED_CELL,
    CURATED_REFERENCE,
    MODIFIED_RESIDUES,
    load_curated_structure_aliases,
    overwrite_frame_grns,
    sequence_alignment_position_map,
    sha256_file,
)


PROPERTY_COLUMNS = (
    "short_name",
    "name",
    "species",
    "domain",
    "function",
    "function_detail",
    "pdb_id",
    "method",
    "resolution",
    "reference",
    "reference_year",
    "length",
    "uniprot_id",
    "sequence",
)
FUNCTION_CARD_COLUMNS = (
    "family",
    "protein function",
    "protein example",
    "grn",
    "function",
)
STRUCTURE_COLUMNS = (
    "atom_id",
    "atom_name",
    "element",
    "group",
    "auth_chain_id",
    "auth_seq_id",
    "res_name3l",
    "res_name1l",
    "grn",
    "x",
    "y",
    "z",
    "occupancy",
    "b_factor",
    "pdb_auth_seq_id",
)
MANIFEST_COLUMNS = (
    "structure_id",
    "structure_type",
    "n_atoms",
    "n_grn_residues",
    "file",
)
ARCHIVE_TIMESTAMP = (1980, 1, 1, 0, 0, 0)
YEAR_PATTERN = re.compile(r"\b(19|20)\d{2}\b")
PDB_PATTERN = re.compile(r"\b([0-9][A-Za-z0-9]{3})\b")
HIDEAKI_PATTERN = re.compile(r"^(.+?)_J\d+_refine\d+$", re.IGNORECASE)


def _clean(value: Any) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def _sequence(value: Any) -> str:
    return re.sub(r"\s+", "", _clean(value)).upper()


def _year(value: Any) -> str:
    match = YEAR_PATTERN.search(_clean(value))
    return match.group(0) if match else ""


def _pdb(value: Any) -> str:
    match = PDB_PATTERN.search(_clean(value))
    return match.group(1).lower() if match else ""


def _identifier(value: Any) -> str:
    return re.sub(r"_+", "_", _clean(value).casefold()).rstrip("_")


def source_record(path: Path, *, rows: int | None = None, columns=None) -> dict:
    record = {
        "file": str(path.resolve()),
        "sha256": sha256_file(path),
    }
    if rows is not None:
        record["rows"] = int(rows)
    if columns is not None:
        record["columns"] = [str(column) for column in columns]
    return record


def normalize_property(path: Path) -> pd.DataFrame:
    """Normalize the ST5 workbook to the canonical public schema."""

    source = pd.read_excel(path, dtype=str)
    required = {
        "display_name",
        "source (species)",
        "Rhodopsin Type (Microbial)",
        "molecular_function",
        "molecular_function_advanced",
        "PDB ID",
        "method",
        "resolution",
        "reference_discovery",
        "sequence",
        "short_name",
    }
    missing = required.difference(source.columns)
    if missing:
        raise ValueError(f"Property workbook lacks columns: {sorted(missing)}")

    sequence = source["sequence"].map(_sequence)
    normalized = pd.DataFrame(
        {
            "short_name": source["short_name"].map(_clean),
            "name": source["display_name"].map(_clean),
            "species": source["source (species)"].map(_clean),
            "domain": source["Rhodopsin Type (Microbial)"].map(_clean),
            "function": source["molecular_function"].map(_clean),
            "function_detail": source["molecular_function_advanced"].map(_clean),
            "pdb_id": source["PDB ID"].map(_pdb),
            "method": source["method"].map(_clean),
            "resolution": source["resolution"].map(_clean),
            "reference": source["reference_discovery"].map(_clean),
            "reference_year": source["reference_discovery"].map(_year),
            "length": sequence.map(lambda value: len(value) if value else ""),
            "uniprot_id": "",
            "sequence": sequence,
        },
        columns=PROPERTY_COLUMNS,
    )
    if normalized["short_name"].eq("").any():
        raise ValueError("Property workbook has empty short_name values")
    if normalized["short_name"].map(_identifier).duplicated().any():
        duplicates = normalized.loc[
            normalized["short_name"].map(_identifier).duplicated(False), "short_name"
        ].tolist()
        raise ValueError(f"Property workbook has duplicate short names: {duplicates}")
    return normalized


def load_function_cards(path: Path) -> pd.DataFrame:
    cards = pd.read_csv(path, dtype=str, keep_default_na=False)
    if tuple(cards.columns) != FUNCTION_CARD_COLUMNS:
        raise ValueError(f"Unexpected function-card columns: {list(cards.columns)}")
    cards = cards.loc[
        cards.apply(lambda row: any(_clean(value) for value in row), axis=1)
    ].copy()
    for column in FUNCTION_CARD_COLUMNS:
        cards[column] = cards[column].map(_clean)
    if cards[list(FUNCTION_CARD_COLUMNS)].eq("").any(axis=None):
        raise ValueError("Function cards contain empty required fields")
    if cards.duplicated(list(FUNCTION_CARD_COLUMNS)).any():
        raise ValueError("Function cards contain duplicate rows")
    return cards


def reconcile_function_card_sources(
    csv_cards: pd.DataFrame,
    workbook_path: Path,
) -> dict[str, Any]:
    """Compare the generated CSV with the non-empty rows in its V2 workbook."""

    workbook = pd.read_excel(workbook_path, dtype=str, keep_default_na=False)
    if tuple(workbook.columns) != FUNCTION_CARD_COLUMNS:
        return {
            "matching": False,
            "reason": "column_mismatch",
            "csv_columns": list(csv_cards.columns),
            "workbook_columns": list(workbook.columns),
        }
    workbook = workbook.loc[
        workbook.apply(lambda row: any(_clean(value) for value in row), axis=1)
    ].copy()
    for column in FUNCTION_CARD_COLUMNS:
        workbook[column] = workbook[column].map(_clean)
    workbook = workbook.reset_index(drop=True)
    csv_cards = csv_cards.reset_index(drop=True)
    differences = []
    for row_index in range(max(len(csv_cards), len(workbook))):
        for column in FUNCTION_CARD_COLUMNS:
            csv_value = (
                csv_cards.at[row_index, column] if row_index < len(csv_cards) else None
            )
            workbook_value = (
                workbook.at[row_index, column] if row_index < len(workbook) else None
            )
            if csv_value != workbook_value:
                differences.append(
                    {
                        "row": row_index + 2,
                        "column": column,
                        "csv": csv_value,
                        "workbook": workbook_value,
                    }
                )
    return {
        "matching": not differences,
        "csv_rows": len(csv_cards),
        "workbook_nonempty_rows": len(workbook),
        "cell_differences": differences,
    }


def git_revision(path: Path) -> str:
    """Return the checked-out revision containing *path*, if available."""

    try:
        return subprocess.check_output(
            ["git", "-C", str(path), "rev-parse", "HEAD"],
            text=True,
            stderr=subprocess.DEVNULL,
        ).strip()
    except (OSError, subprocess.CalledProcessError):
        return ""


def normalize_grn_label(value: str) -> str:
    value = _clean(value)
    match = re.fullmatch(r"([1-7])\.(\d+)", value)
    if not match:
        return value
    suffix = match.group(2)
    if len(suffix) < 2:
        suffix = suffix.zfill(2)
    return f"{match.group(1)}.{suffix}"


def property_join_map(
    property_table: pd.DataFrame,
    reference_ids: Iterable[str],
    structure_mapping_path: Path,
) -> tuple[dict[str, int], dict[str, Any]]:
    """Resolve every reference ID to exactly one normalized property row."""

    aliases: dict[str, set[int]] = {}

    def add(alias: str, row_index: int) -> None:
        alias = _identifier(alias)
        if alias:
            aliases.setdefault(alias, set()).add(int(row_index))

    for row_index, row in property_table.iterrows():
        short_name = row["short_name"]
        add(short_name, row_index)
        add(f"{short_name}_model_0", row_index)
        add(row["pdb_id"], row_index)

    structure_mapping = json.loads(structure_mapping_path.read_text(encoding="utf-8"))
    for experimental_id, predicted_id in structure_mapping.items():
        predicted_matches = aliases.get(_identifier(predicted_id), set())
        for row_index in predicted_matches:
            add(experimental_id, row_index)

    # The two curated domains join to the single parent property record.
    tara_parent = aliases.get(_identifier("7pl9"), set())
    for row_index in tara_parent:
        add("TARA_A", row_index)
        add("TARA_B", row_index)

    resolved: dict[str, int] = {}
    diagnostics = []
    for structure_id in reference_ids:
        candidate_keys = {_identifier(structure_id)}
        hideaki = HIDEAKI_PATTERN.match(structure_id)
        if hideaki:
            candidate_keys.add(_identifier(hideaki.group(1)))
        candidates: set[int] = set()
        for key in candidate_keys:
            candidates.update(aliases.get(key, set()))
        diagnostics.append(
            {
                "structure_id": structure_id,
                "candidate_keys": sorted(candidate_keys),
                "property_rows": sorted(candidates),
            }
        )
        if len(candidates) != 1:
            continue
        resolved[structure_id] = next(iter(candidates))

    failures = [
        row for row in diagnostics if len(row["property_rows"]) != 1
    ]
    return resolved, {"rows": diagnostics, "failures": failures}


def _observed_residue(frame: pd.DataFrame, mask: pd.Series) -> list[str]:
    observed = sorted(
        {
            _clean(value).upper()
            for value in frame.loc[mask, "res_name1l"]
            if _clean(value)
        }
    )
    if observed == ["X"] and "res_name3l" in frame.columns:
        modified = {
            MODIFIED_RESIDUES.get(_clean(value).upper())
            for value in frame.loc[mask, "res_name3l"]
        }
        modified.discard(None)
        if modified:
            observed = sorted(modified)
    return observed


def validate_structure_reference(
    frame: pd.DataFrame,
    reference_row: pd.Series,
    structure_id: str,
) -> dict[str, Any]:
    """Reconcile residue number, identity, and persisted GRN label cell by cell."""

    numeric_ids = pd.to_numeric(frame["auth_seq_id"], errors="coerce")
    chain_mask = frame["auth_chain_id"].astype(str).eq("A")
    missing = []
    residue_mismatches = []
    grn_mismatches = []
    for grn, cell in reference_row.items():
        match = CURATED_CELL.fullmatch(_clean(cell))
        if not match:
            continue
        expected_residue = match.group(1)
        sequence_number = int(match.group(2))
        mask = chain_mask & numeric_ids.eq(sequence_number)
        if not mask.any():
            missing.append({"grn": str(grn), "auth_seq_id": sequence_number})
            continue
        observed = _observed_residue(frame, mask)
        # A curated ``X<position>`` cell identifies a mapped position with an
        # unresolved residue identity, so any observed residue is compatible.
        if expected_residue != "X" and observed and observed != [expected_residue]:
            residue_mismatches.append(
                {
                    "grn": str(grn),
                    "auth_seq_id": sequence_number,
                    "expected": expected_residue,
                    "observed": observed,
                }
            )
        observed_grns = sorted(
            {
                _clean(value)
                for value in frame.loc[mask, "grn"]
                if _clean(value)
            }
        )
        if observed_grns != [str(grn)]:
            grn_mismatches.append(
                {
                    "auth_seq_id": sequence_number,
                    "expected": str(grn),
                    "observed": observed_grns,
                }
            )
    return {
        "structure_id": structure_id,
        "missing_reference_residues": missing,
        "residue_identity_mismatches": residue_mismatches,
        "grn_mismatches": grn_mismatches,
    }


def normalize_replacement_positions(
    frame: pd.DataFrame,
    position_map: dict[int, int],
    *,
    chain_id: str = "A",
) -> pd.DataFrame:
    """Expose curated source positions while retaining original PDB numbering."""

    normalized = frame.copy()
    normalized["pdb_auth_seq_id"] = normalized["auth_seq_id"]
    chain_mask = normalized["auth_chain_id"].astype(str).eq(str(chain_id))
    numeric_ids = pd.to_numeric(normalized["auth_seq_id"], errors="coerce")
    reverse_map = {observed: source for source, observed in position_map.items()}
    if len(reverse_map) != len(position_map):
        raise ValueError("Replacement alignment maps multiple source residues together")

    # Unaligned construct tags and unresolved residues have no corresponding
    # curated source position. Their original numbering remains available in
    # pdb_auth_seq_id, but cannot masquerade as a reference position.
    normalized.loc[chain_mask, "auth_seq_id"] = ""
    for observed_position, source_position in reverse_map.items():
        mask = chain_mask & numeric_ids.eq(observed_position)
        normalized.loc[mask, "auth_seq_id"] = source_position
    return normalized


def add_original_position_column(frame: pd.DataFrame) -> pd.DataFrame:
    """Add an explicit original-position column to non-replacement structures."""

    result = frame.copy()
    result["pdb_auth_seq_id"] = result["auth_seq_id"]
    return result


def gzip_csv(frame: pd.DataFrame, path: Path) -> None:
    csv_bytes = frame.to_csv(index=False, lineterminator="\n").encode("utf-8")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as compressed:
            compressed.write(csv_bytes)


def write_json(path: Path, payload: Any) -> None:
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def zip_bundle(bundle_dir: Path, output_zip: Path) -> str:
    temporary = output_zip.with_suffix(output_zip.suffix + ".tmp")
    output_zip.parent.mkdir(parents=True, exist_ok=True)
    try:
        with zipfile.ZipFile(
            temporary,
            mode="w",
            compression=zipfile.ZIP_DEFLATED,
            allowZip64=True,
        ) as archive:
            for path in sorted(
                (item for item in bundle_dir.rglob("*") if item.is_file()),
                key=lambda item: item.relative_to(bundle_dir.parent).as_posix(),
            ):
                name = path.relative_to(bundle_dir.parent).as_posix()
                info = zipfile.ZipInfo(name, date_time=ARCHIVE_TIMESTAMP)
                info.compress_type = zipfile.ZIP_DEFLATED
                info.external_attr = (stat.S_IFREG | 0o644) << 16
                with path.open("rb") as source, archive.open(info, "w") as target:
                    shutil.copyfileobj(source, target, length=1024 * 1024)
        temporary.replace(output_zip)
    finally:
        temporary.unlink(missing_ok=True)
    digest = sha256_file(output_zip)
    output_zip.with_suffix(output_zip.suffix + ".sha256").write_text(
        f"{digest}  {output_zip.name}\n", encoding="utf-8"
    )
    return digest


def build_bundle(args: argparse.Namespace) -> dict[str, Any]:
    args.output_zip.unlink(missing_ok=True)
    args.output_zip.with_suffix(args.output_zip.suffix + ".sha256").unlink(
        missing_ok=True
    )
    if args.output_dir.exists():
        shutil.rmtree(args.output_dir)
    if args.reports_dir.exists():
        shutil.rmtree(args.reports_dir)
    args.output_dir.mkdir(parents=True)
    args.reports_dir.mkdir(parents=True)
    structures_dir = args.output_dir / "structures"
    structures_dir.mkdir()

    property_table = normalize_property(args.property)
    reference = pd.read_csv(
        args.grn_reference, index_col=0, dtype=str, keep_default_na=False
    ).fillna("-")
    reference.index = reference.index.astype(str)
    reference.columns = reference.columns.astype(str)
    if reference.index.has_duplicates:
        raise ValueError("GRN reference contains duplicate structure_id values")
    if len(reference) != 130 or {"TARA_A", "TARA_B"}.difference(reference.index):
        raise ValueError(
            "Expected 130 curated reference rows including TARA_A and TARA_B"
        )
    structure_aliases = load_curated_structure_aliases()
    release_id_by_reference = {}
    for release_id, spec in structure_aliases.items():
        source_reference_id = str(spec["reference_id"])
        if source_reference_id in release_id_by_reference:
            raise ValueError(
                f"Multiple release structures replace {source_reference_id}"
            )
        release_id_by_reference[source_reference_id] = release_id
    unknown_aliases = set(release_id_by_reference).difference(reference.index)
    if unknown_aliases:
        raise ValueError(
            f"Release aliases target missing reference rows: {sorted(unknown_aliases)}"
        )
    canonical_reference = reference.rename(index=release_id_by_reference)
    if canonical_reference.index.has_duplicates:
        raise ValueError("Canonicalized GRN reference has duplicate structure IDs")
    cards = load_function_cards(args.function_cards)
    function_card_reconciliation = reconcile_function_card_sources(
        cards, args.function_cards_workbook
    )
    grn_config = json.loads(args.grn_config.read_text(encoding="utf-8"))
    if not {"standard", "strict"}.issubset(grn_config.get("mo", {})):
        raise ValueError("GRN config lacks mo.standard or mo.strict")

    property_table.to_csv(args.output_dir / "property.csv", index=False)
    canonical_reference.to_csv(
        args.output_dir / "grn_reference.csv",
        index=True,
        index_label="structure_id",
        lineterminator="\n",
    )
    shutil.copyfile(args.grn_config, args.output_dir / "grn_config.json")
    cards.to_csv(args.output_dir / "function_cards.csv", index=False, lineterminator="\n")

    import protos
    from protos.processing.structure import StructureProcessor

    protos.set_data_path(str(args.data_root))
    processor = StructureProcessor("groundtruth_bundle")
    datasets = ["mo_exp_A", "mo_exp_B", "mo_pred_exp", "mo_pred_novel"]
    registered = []
    for dataset in datasets:
        registered.extend(processor.get_dataset_entities(dataset))
    registered_lookup = {identifier.casefold(): identifier for identifier in registered}

    structure_rows = []
    structure_components = {}
    reconciliations = []
    for source_reference_id, row in reference.iterrows():
        reference_id = release_id_by_reference.get(
            source_reference_id, source_reference_id
        )
        stored_id = registered_lookup.get(reference_id.casefold())
        if stored_id is None:
            raise ValueError(f"No registered structure for reference row {reference_id}")
        frame = processor.load_entity(stored_id)
        if frame is None:
            raise ValueError(f"Registered structure failed to load: {stored_id}")
        alias_spec = structure_aliases.get(reference_id)
        position_map = None
        position_alignment = None
        if alias_spec is not None:
            position_map, position_alignment = sequence_alignment_position_map(
                frame,
                row,
                chain_id="A",
            )
        annotated, overwrite_report = overwrite_frame_grns(
            frame,
            row,
            structure_id=stored_id,
            reference_id=source_reference_id,
            chain_id="A",
            position_map=position_map,
        )
        if position_map is not None:
            annotated = normalize_replacement_positions(annotated, position_map)
            overwrite_report["alias"] = dict(alias_spec)
            overwrite_report["position_alignment"] = position_alignment
        else:
            annotated = add_original_position_column(annotated)
        validation = validate_structure_reference(annotated, row, reference_id)
        validation["source_reference_id"] = source_reference_id
        validation["overwrite"] = overwrite_report
        reconciliations.append(validation)

        missing_columns = set(STRUCTURE_COLUMNS).difference(annotated.columns)
        if missing_columns:
            raise ValueError(
                f"{reference_id} lacks canonical structure columns: "
                f"{sorted(missing_columns)}"
            )
        canonical = annotated.loc[:, STRUCTURE_COLUMNS].copy()
        relative_file = f"structures/{reference_id}.csv.gz"
        gzip_csv(canonical, args.output_dir / relative_file)
        structure_components[relative_file] = {
            "sha256": sha256_file(args.output_dir / relative_file),
            "bytes": (args.output_dir / relative_file).stat().st_size,
            "rows": len(canonical),
            "columns": list(STRUCTURE_COLUMNS),
            "provenance": {
                "registered_structure_id": stored_id,
                "grn_reference_structure_id": source_reference_id,
                "canonical_structure_id": reference_id,
                "grn_policy": "rewritten_from_grn_reference.csv",
            },
        }
        residue_count = (
            canonical.loc[
                canonical["grn"].map(_clean).ne(""),
                ["auth_chain_id", "auth_seq_id"],
            ]
            .drop_duplicates()
            .shape[0]
        )
        structure_rows.append(
            {
                "structure_id": reference_id,
                "structure_type": (
                    "predicted" if "_model_0" in reference_id.casefold() else "experimental"
                ),
                "n_atoms": len(canonical),
                "n_grn_residues": residue_count,
                "file": relative_file,
            }
        )

    manifest = pd.DataFrame(structure_rows, columns=MANIFEST_COLUMNS)
    manifest.to_csv(structures_dir / "manifest.csv", index=False, lineterminator="\n")

    joins, join_report = property_join_map(
        property_table, canonical_reference.index, args.structure_mapping
    )
    card_grns = {
        normalize_grn_label(value) for value in cards["grn"].astype(str)
    }
    reference_grns = set(canonical_reference.columns)
    structural_exception_mask = cards["family"].eq("Structural exceptions")
    invalid_card_grns = sorted(
        {
            cards.at[index, "grn"]
            for index in cards.index
            if normalize_grn_label(cards.at[index, "grn"]) not in reference_grns
            and not structural_exception_mask.at[index]
        }
    )

    reconciliation_failures = []
    for row in reconciliations:
        predicted = "_model_0" in row["structure_id"].casefold()
        if (
            row["residue_identity_mismatches"]
            or row["grn_mismatches"]
            or (predicted and row["missing_reference_residues"])
        ):
            reconciliation_failures.append(row)
    validation = {
        "property_rows": len(property_table),
        "grn_reference_rows": len(canonical_reference),
        "structure_manifest_rows": len(manifest),
        "function_cards": len(cards),
        "function_card_families": int(cards["family"].nunique()),
        "reference_structure_ids_unique": bool(reference.index.is_unique),
        "manifest_structure_ids_unique": bool(manifest["structure_id"].is_unique),
        "manifest_files_unique": bool(manifest["file"].is_unique),
        "property_joined_reference_rows": len(joins),
        "property_join_failures": join_report["failures"],
        "invalid_function_card_grns": invalid_card_grns,
        "function_card_sources": function_card_reconciliation,
        "structure_reference_failures": reconciliation_failures,
        "experimental_missing_reference_residues": [
            {
                "structure_id": row["structure_id"],
                "source_reference_id": row["source_reference_id"],
                "missing_reference_residues": row["missing_reference_residues"],
            }
            for row in reconciliations
            if "_model_0" not in row["structure_id"].casefold()
            and row["missing_reference_residues"]
        ],
    }
    validation["passed"] = bool(
        len(canonical_reference) == 130
        and len(manifest) == 130
        and len(cards) == 84
        and cards["family"].nunique() == 13
        and not join_report["failures"]
        and not invalid_card_grns
        and function_card_reconciliation["matching"]
        and not reconciliation_failures
    )
    write_json(args.reports_dir / "validation.json", validation)
    write_json(
        args.reports_dir / "reconciliation.json",
        {
            "function_cards": function_card_reconciliation,
            "property_joins": join_report,
            "structure_reference": reconciliations,
        },
    )

    tabular_components = {
        "property.csv": (len(property_table), list(PROPERTY_COLUMNS)),
        "grn_reference.csv": (
            len(canonical_reference),
            ["structure_id", *list(canonical_reference.columns)],
        ),
        "function_cards.csv": (len(cards), list(FUNCTION_CARD_COLUMNS)),
        "structures/manifest.csv": (len(manifest), list(MANIFEST_COLUMNS)),
    }
    components = {}
    for relative_file, (rows, columns) in tabular_components.items():
        path = args.output_dir / relative_file
        components[relative_file] = {
            "sha256": sha256_file(path),
            "bytes": path.stat().st_size,
            "rows": rows,
            "columns": columns,
        }
    config_path = args.output_dir / "grn_config.json"
    components["grn_config.json"] = {
        "sha256": sha256_file(config_path),
        "bytes": config_path.stat().st_size,
        "provenance": "byte_copy_of_pulled_protos_config",
    }
    components.update(structure_components)
    bundle = {
        "schema": "grn_opsins_groundtruth/v1",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "authoritative_grn_policy": (
            "All structure GRN labels are rewritten from the pulled ProtOS "
            "type_I_opsins.csv table; embedded labels are ignored."
        ),
        "counts": {
            "property_records": len(property_table),
            "grn_reference_rows": len(canonical_reference),
            "structure_files": len(manifest),
            "function_cards": len(cards),
            "function_card_families": int(cards["family"].nunique()),
        },
        "sources": {
            "property": source_record(
                args.property,
                rows=len(property_table),
                columns=property_table.columns,
            ),
            "grn_reference": source_record(
                args.grn_reference,
                rows=len(reference),
                columns=reference.columns,
            ),
            "grn_config": source_record(args.grn_config),
            "function_cards": source_record(
                args.function_cards,
                rows=len(cards),
                columns=cards.columns,
            ),
            "function_cards_workbook": source_record(args.function_cards_workbook),
            "structure_registry": {
                "path": str(args.data_root.resolve()),
                "datasets": datasets,
                "mogrn_revision": git_revision(ROOT),
                "protos_revision": git_revision(ROOT / "protos"),
            },
        },
        "canonical_structure_aliases": {
            release_id: {
                **dict(spec),
                "source_cells_preserved": bool(
                    canonical_reference.loc[release_id].equals(
                        reference.loc[str(spec["reference_id"])]
                    )
                ),
            }
            for release_id, spec in structure_aliases.items()
        },
        "validation": validation,
        "components": components,
    }
    write_json(args.output_dir / "bundle.json", bundle)

    if not validation["passed"]:
        raise RuntimeError(
            f"Bundle validation failed; inspect {args.reports_dir / 'validation.json'}"
        )

    archive_digest = zip_bundle(args.output_dir, args.output_zip)
    return {
        "bundle_dir": str(args.output_dir),
        "archive": str(args.output_zip),
        "archive_sha256": archive_digest,
        "counts": bundle["counts"],
        "validation": validation,
        "reports": str(args.reports_dir),
    }


def parse_args() -> argparse.Namespace:
    sibling = ROOT.parent / "grn_opsins"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--property", type=Path, default=ROOT / "mo_exp_ST5_HEK1.xlsx")
    parser.add_argument("--grn-reference", type=Path, default=CURATED_REFERENCE)
    parser.add_argument(
        "--grn-config",
        type=Path,
        default=PROTOS_SRC / "protos/reference_data/grn/configs/config.json",
    )
    parser.add_argument(
        "--function-cards",
        type=Path,
        default=sibling / "backend/data/grn_function_cards_v2.csv",
    )
    parser.add_argument(
        "--function-cards-workbook",
        type=Path,
        default=sibling / "grn_function_cards_v2_ST1_UTF8_HEK2.xlsx",
    )
    parser.add_argument("--data-root", type=Path, default=ROOT / "data")
    parser.add_argument(
        "--structure-mapping",
        type=Path,
        default=ROOT / "opsin_output/structure_mapping.json",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=ROOT / "zenodo_upload/grn_opsins_groundtruth",
    )
    parser.add_argument(
        "--output-zip",
        type=Path,
        default=ROOT / "zenodo_upload/grn_opsins_groundtruth.zip",
    )
    parser.add_argument(
        "--reports-dir",
        type=Path,
        default=ROOT / "zenodo_upload/grn_opsins_groundtruth_reports",
    )
    return parser.parse_args()


def main() -> int:
    result = build_bundle(parse_args())
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
