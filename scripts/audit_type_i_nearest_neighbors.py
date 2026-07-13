#!/usr/bin/env python3
"""Audit curated Type-I GRNs against each row's nearest alignment neighbour.

The CSV is treated as the curated alignment to audit.  Neighbours are selected
only from amino-acid identity in populated TM columns.  Each target structure is
then aligned to its neighbour as a complete seven-helix CA bundle with Protos.
The resulting coordinate frame is used to measure same-GRN errors, gap runs, and
whether a nearby GRN in the neighbour fits better than the assigned GRN.

No experimental/predicted categories are used by this analysis.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Iterable

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(__file__).resolve().parent.parent
PROTOS_SRC = PROJECT_ROOT / "protos" / "src"
if PROTOS_SRC.exists():
    sys.path.insert(0, str(PROTOS_SRC))

DEFAULT_DATASETS = [
    "mo_exp_A",
    "mo_exp_B",
    "mo_pred_exp",
    "mo_pred_novel",
    "mo_dual_rhodopsin_domains",
]
TM_PATTERN = re.compile(r"^([1-7])\.(\d+)$")
CELL_PATTERN = re.compile(r"^([A-Za-z])(-?\d+)$")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def populated(value: Any) -> bool:
    return pd.notna(value) and str(value).strip() not in {"", "-"}


def is_insertion_grn(label: str) -> bool:
    """Return whether a TM label is an explicit insertion coordinate (e.g. 5.451)."""
    match = TM_PATTERN.fullmatch(str(label))
    return bool(match and int(match.group(2)) > 99)


def residue_letter(value: Any) -> str | None:
    match = CELL_PATTERN.fullmatch(str(value).strip()) if populated(value) else None
    return match.group(1).upper() if match else None


def residue_number(value: Any) -> int | None:
    match = CELL_PATTERN.fullmatch(str(value).strip()) if populated(value) else None
    return int(match.group(2)) if match else None


def find_nearest_neighbours(
    table: pd.DataFrame,
    tm_columns: list[str],
    min_overlap: int,
    top_k: int,
) -> pd.DataFrame:
    """Rank alignment-table neighbours by common-column sequence identity."""
    letters = np.empty((len(table), len(tm_columns)), dtype=object)
    letters[:] = None
    for i, (_, row) in enumerate(table[tm_columns].iterrows()):
        letters[i] = [residue_letter(value) for value in row]

    present = letters != None  # noqa: E711 - intentional object-array comparison
    records: list[dict[str, Any]] = []
    ids = [str(value) for value in table.index]
    for i, target in enumerate(ids):
        ranked: list[dict[str, Any]] = []
        for j, neighbour in enumerate(ids):
            if i == j:
                continue
            common = present[i] & present[j]
            overlap = int(common.sum())
            if overlap < min_overlap:
                continue
            matches = int((letters[i, common] == letters[j, common]).sum())
            target_n = int(present[i].sum())
            neighbour_n = int(present[j].sum())
            coverage = overlap / max(min(target_n, neighbour_n), 1)
            union = int((present[i] | present[j]).sum())
            ranked.append(
                {
                    "target": target,
                    "neighbour": neighbour,
                    "identity": matches / overlap,
                    "matches": matches,
                    "common_residues": overlap,
                    "target_tm_residues": target_n,
                    "neighbour_tm_residues": neighbour_n,
                    "common_coverage": coverage,
                    "gap_disagreements": int((present[i] ^ present[j]).sum()),
                    "gap_agreement_jaccard": overlap / max(union, 1),
                }
            )
        ranked.sort(
            key=lambda x: (
                x["identity"],
                x["common_coverage"],
                x["common_residues"],
                -x["gap_disagreements"],
                x["neighbour"],
            ),
            reverse=True,
        )
        for rank, item in enumerate(ranked[:top_k], start=1):
            records.append({"rank": rank, **item})
    return pd.DataFrame(records)


def load_structure_id_map(processor: Any, datasets: Iterable[str]) -> dict[str, str]:
    """Return a case-insensitive ID map across the requested Protos datasets."""
    result: dict[str, str] = {}
    collisions: dict[str, set[str]] = defaultdict(set)
    for dataset in datasets:
        if not processor.dataset_manager.dataset_exists(dataset):
            raise ValueError(f"Protos dataset {dataset!r} does not exist")
        for entity in processor.get_dataset_entities(dataset):
            key = str(entity).casefold()
            collisions[key].add(str(entity))
            result.setdefault(key, str(entity))
    ambiguous = {key: values for key, values in collisions.items() if len(values) > 1}
    if ambiguous:
        raise ValueError(f"Case-insensitive structure ID collisions: {ambiguous}")
    return result


def load_ca_lookup(processor: Any, structure_id: str, preferred_chain: str) -> tuple[str, dict[int, dict[str, Any]]]:
    """Load one structure and map author residue number to its CA coordinate."""
    frame = processor.load_entity(structure_id)
    if frame is None:
        raise ValueError("load_entity returned None")
    frame = frame.reset_index()
    atom_col = "atom_name" if "atom_name" in frame.columns else "res_atom_name"
    ca = frame[frame[atom_col].astype(str).str.upper() == "CA"].copy()
    if ca.empty:
        raise ValueError("no CA atoms")
    counts = ca["auth_chain_id"].astype(str).value_counts()
    chain = preferred_chain if preferred_chain in counts.index else str(counts.index[0])
    ca = ca[ca["auth_chain_id"].astype(str) == chain]
    lookup: dict[int, dict[str, Any]] = {}
    for row in ca.itertuples(index=False):
        number = int(getattr(row, "auth_seq_id"))
        if number in lookup:
            continue
        lookup[number] = {
            "aa_structure": str(getattr(row, "res_name1l", "?")).upper(),
            "x": float(getattr(row, "x")),
            "y": float(getattr(row, "y")),
            "z": float(getattr(row, "z")),
        }
    return chain, lookup


def map_table_rows_to_coordinates(
    table: pd.DataFrame,
    tm_columns: list[str],
    processor: Any,
    id_map: dict[str, str],
    preferred_chain: str,
) -> tuple[dict[str, pd.DataFrame], pd.DataFrame, pd.DataFrame]:
    """Map curated CSV cells directly to CA coordinates by author residue ID."""
    mapped: dict[str, pd.DataFrame] = {}
    qc: list[dict[str, Any]] = []
    mismatches: list[dict[str, Any]] = []
    for number, (table_id, row) in enumerate(table[tm_columns].iterrows(), start=1):
        table_id = str(table_id)
        entity_id = id_map.get(table_id.casefold())
        if entity_id is None:
            qc.append({"structure": table_id, "status": "structure_not_registered"})
            continue
        try:
            chain, ca_lookup = load_ca_lookup(processor, entity_id, preferred_chain)
        except Exception as exc:
            qc.append({"structure": table_id, "entity_id": entity_id, "status": "load_failed", "error": str(exc)})
            continue

        records: list[dict[str, Any]] = []
        populated_cells = 0
        identity_mismatches = 0
        missing_coordinates = 0
        for column_index, grn in enumerate(tm_columns):
            value = row[grn]
            if not populated(value):
                continue
            populated_cells += 1
            author_id = residue_number(value)
            aa_table = residue_letter(value)
            ca = ca_lookup.get(author_id) if author_id is not None else None
            if ca is None:
                missing_coordinates += 1
                continue
            mismatch = aa_table != ca["aa_structure"] and ca["aa_structure"] not in {"?", "X"}
            identity_mismatches += int(mismatch)
            if mismatch:
                mismatches.append(
                    {
                        "structure": table_id,
                        "entity_id": entity_id,
                        "chain": chain,
                        "grn": grn,
                        "auth_seq_id": author_id,
                        "table_residue": aa_table,
                        "structure_residue": ca["aa_structure"],
                        "table_cell": str(value),
                        "structure_cell": f"{ca['aa_structure']}{author_id}",
                    }
                )
            records.append(
                {
                    "structure": table_id,
                    "entity_id": entity_id,
                    "chain": chain,
                    "grn": grn,
                    "helix": int(grn.split(".", 1)[0]),
                    "column_index": column_index,
                    "auth_seq_id": author_id,
                    "aa_table": aa_table,
                    "aa_structure": ca["aa_structure"],
                    "residue_identity_mismatch": mismatch,
                    "x": ca["x"],
                    "y": ca["y"],
                    "z": ca["z"],
                }
            )
        mapped[table_id] = pd.DataFrame(records)
        qc.append(
            {
                "structure": table_id,
                "entity_id": entity_id,
                "chain": chain,
                "status": "ok",
                "populated_tm_cells": populated_cells,
                "mapped_tm_coordinates": len(records),
                "missing_tm_coordinates": missing_coordinates,
                "residue_identity_mismatches": identity_mismatches,
            }
        )
        if number % 25 == 0 or number == len(table):
            print(f"[coordinates] {number}/{len(table)}")
    return mapped, pd.DataFrame(qc), pd.DataFrame(mismatches)


def ordered_bundle(table: pd.DataFrame) -> pd.DataFrame:
    """Order mapped coordinates as seven physically ordered helix fragments."""
    if table.empty:
        return table
    return table.sort_values(["helix", "auth_seq_id", "column_index"]).drop_duplicates(
        ["helix", "auth_seq_id"], keep="first"
    )


def structural_pair_audit(
    target: str,
    neighbour: str,
    coordinate_tables: dict[str, pd.DataFrame],
    tm_columns: list[str],
    local_radius: int,
    shift_error_floor: float,
    shift_improvement_floor: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Align a TM pair and calculate same-GRN and nearby-GRN residuals."""
    from protos.analysis.structure.alignment import get_structure_alignment

    target_table = coordinate_tables[target]
    neighbour_table = coordinate_tables[neighbour]
    target_bundle = ordered_bundle(target_table)
    neighbour_bundle = ordered_bundle(neighbour_table)
    fixed = neighbour_bundle[["x", "y", "z"]].to_numpy(float)
    mobile = target_bundle[["x", "y", "z"]].to_numpy(float)
    rotation, translation, path, rmsd = get_structure_alignment(fixed, mobile)
    transformed_all = target_table[["x", "y", "z"]].to_numpy(float) @ np.asarray(rotation) + np.asarray(translation)
    transformed = target_table.copy()
    transformed[["aligned_x", "aligned_y", "aligned_z"]] = transformed_all

    target_by_grn = transformed.set_index("grn")
    neighbour_by_grn = neighbour_table.set_index("grn")
    helix_columns: dict[int, list[str]] = {
        helix: [column for column in tm_columns if column.startswith(f"{helix}.")]
        for helix in range(1, 8)
    }
    records: list[dict[str, Any]] = []
    for grn in sorted(set(target_by_grn.index) & set(neighbour_by_grn.index), key=tm_columns.index):
        target_row = target_by_grn.loc[grn]
        neighbour_row = neighbour_by_grn.loc[grn]
        target_coord = target_row[["aligned_x", "aligned_y", "aligned_z"]].to_numpy(float)
        same_coord = neighbour_row[["x", "y", "z"]].to_numpy(float)
        same_error = float(np.linalg.norm(target_coord - same_coord))
        helix = int(str(grn).split(".", 1)[0])
        columns = helix_columns[helix]
        position = columns.index(grn)
        candidates: list[tuple[float, str, int]] = []
        for neighbour_position in range(max(0, position - local_radius), min(len(columns), position + local_radius + 1)):
            candidate_grn = columns[neighbour_position]
            if candidate_grn not in neighbour_by_grn.index:
                continue
            candidate_coord = neighbour_by_grn.loc[candidate_grn, ["x", "y", "z"]].to_numpy(float)
            candidates.append((float(np.linalg.norm(target_coord - candidate_coord)), candidate_grn, neighbour_position - position))
        best_error, best_grn, shift = min(candidates) if candidates else (math.nan, "", 0)
        improvement = same_error - best_error
        possible_shift = bool(
            best_grn != grn
            and same_error >= shift_error_floor
            and improvement >= shift_improvement_floor
        )
        records.append(
            {
                "target": target,
                "neighbour": neighbour,
                "grn": grn,
                "helix": helix,
                "target_cell": f"{target_row.aa_table}{int(target_row.auth_seq_id)}",
                "neighbour_cell": f"{neighbour_row.aa_table}{int(neighbour_row.auth_seq_id)}",
                "same_grn_error": same_error,
                "best_local_neighbour_grn": best_grn,
                "best_local_error": best_error,
                "shift_columns": shift,
                "shift_improvement": improvement,
                "possible_shift": possible_shift,
            }
        )
    path_a, path_b = path
    qc = {
        "target": target,
        "neighbour": neighbour,
        "status": "ok",
        "bundle_rmsd": float(rmsd),
        "aligned_bundle_residues": len(path_a),
        "target_bundle_residues": len(target_bundle),
        "neighbour_bundle_residues": len(neighbour_bundle),
        "target_alignment_coverage": len(path_b) / max(len(target_bundle), 1),
        "neighbour_alignment_coverage": len(path_a) / max(len(neighbour_bundle), 1),
    }
    return pd.DataFrame(records), qc


def gap_disagreements(table: pd.DataFrame, nearest: pd.DataFrame, tm_columns: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Record individual gap disagreements and collapse them into runs."""
    cells: list[dict[str, Any]] = []
    runs: list[dict[str, Any]] = []
    for pair in nearest.itertuples(index=False):
        target_row = table.loc[pair.target]
        neighbour_row = table.loc[pair.neighbour]
        pair_cells: list[dict[str, Any]] = []
        for index, grn in enumerate(tm_columns):
            target_present = populated(target_row[grn])
            neighbour_present = populated(neighbour_row[grn])
            if target_present == neighbour_present:
                continue
            pair_cells.append(
                {
                    "target": pair.target,
                    "neighbour": pair.neighbour,
                    "grn": grn,
                    "helix": int(grn.split(".", 1)[0]),
                    "column_index": index,
                    "disagreement": "target_only" if target_present else "neighbour_only",
                    "target_cell": target_row[grn],
                    "neighbour_cell": neighbour_row[grn],
                }
            )
        cells.extend(pair_cells)
        for direction in ("target_only", "neighbour_only"):
            subset = [cell for cell in pair_cells if cell["disagreement"] == direction]
            for helix in range(1, 8):
                helix_cells = sorted((cell for cell in subset if cell["helix"] == helix), key=lambda x: x["column_index"])
                current: list[dict[str, Any]] = []
                for cell in helix_cells:
                    if current and cell["column_index"] != current[-1]["column_index"] + 1:
                        runs.append(make_gap_run(current, direction, table, tm_columns))
                        current = []
                    current.append(cell)
                if current:
                    runs.append(make_gap_run(current, direction, table, tm_columns))
    return pd.DataFrame(cells), pd.DataFrame(runs)


def curated_internal_gap_inventory(
    table: pd.DataFrame, tm_columns: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Inventory dashes bounded by assigned residues in the same TM helix."""
    cells: list[dict[str, Any]] = []
    runs: list[dict[str, Any]] = []
    for structure, row in table.iterrows():
        for helix in range(1, 8):
            # A dash in an insertion-only column is the normal state, not a
            # deletion.  Standard gap inventory therefore ignores columns
            # such as 5.451 while the populated insertion is reported by its
            # own annotation column.
            columns = [
                column for column in tm_columns
                if column.startswith(f"{helix}.") and not is_insertion_grn(column)
            ]
            occupied = [index for index, grn in enumerate(columns) if populated(row[grn])]
            if not occupied:
                continue
            first, last = min(occupied), max(occupied)
            gap_positions = [index for index in range(first, last + 1) if not populated(row[columns[index]])]
            current: list[int] = []
            for index in gap_positions:
                if current and index != current[-1] + 1:
                    runs.append(make_curated_gap_run(str(structure), helix, current, columns, row))
                    current = []
                current.append(index)
                cells.append(
                    {
                        "structure": str(structure),
                        "helix": helix,
                        "grn": columns[index],
                        "left_grn": columns[index - 1] if index > first else "",
                        "left_cell": row[columns[index - 1]] if index > first else "",
                        "right_grn": columns[index + 1] if index < last else "",
                        "right_cell": row[columns[index + 1]] if index < last else "",
                    }
                )
            if current:
                runs.append(make_curated_gap_run(str(structure), helix, current, columns, row))
    cell_table = pd.DataFrame(cells)
    run_table = pd.DataFrame(runs)
    if cell_table.empty:
        return cell_table, run_table, pd.DataFrame()
    summaries: list[dict[str, Any]] = []
    for grn, group in cell_table.groupby("grn", sort=False):
        structures = sorted(group["structure"].astype(str))
        summaries.append(
            {
                "grn": grn,
                "helix": int(grn.split(".", 1)[0]),
                "n_structures_with_internal_gap": len(structures),
                "structures": ";".join(structures),
            }
        )
    summary = pd.DataFrame(summaries)
    summary["tm_column_order"] = summary["grn"].map(tm_columns.index)
    summary = summary.sort_values("tm_column_order").drop(columns="tm_column_order")
    return cell_table, run_table, summary


def make_curated_gap_run(
    structure: str,
    helix: int,
    positions: list[int],
    columns: list[str],
    row: pd.Series,
) -> dict[str, Any]:
    start, end = positions[0], positions[-1]
    return {
        "structure": structure,
        "helix": helix,
        "start_grn": columns[start],
        "end_grn": columns[end],
        "length": len(positions),
        "left_anchor_grn": columns[start - 1],
        "left_anchor_cell": row[columns[start - 1]],
        "right_anchor_grn": columns[end + 1],
        "right_anchor_cell": row[columns[end + 1]],
    }


def make_gap_run(run: list[dict[str, Any]], direction: str, table: pd.DataFrame, tm_columns: list[str]) -> dict[str, Any]:
    first = run[0]
    target = first["target"]
    neighbour = first["neighbour"]
    helix = first["helix"]
    helix_cols = [column for column in tm_columns if column.startswith(f"{helix}.")]
    union_positions = [
        index for index, column in enumerate(helix_cols)
        if populated(table.loc[target, column]) or populated(table.loc[neighbour, column])
    ]
    run_positions = [helix_cols.index(cell["grn"]) for cell in run]
    location = "internal"
    if union_positions and min(run_positions) == min(union_positions):
        location = "n_terminal"
    if union_positions and max(run_positions) == max(union_positions):
        location = "c_terminal" if location == "internal" else "both_terminals"
    return {
        "target": target,
        "neighbour": neighbour,
        "helix": helix,
        "disagreement": direction,
        "start_grn": run[0]["grn"],
        "end_grn": run[-1]["grn"],
        "length": len(run),
        "location": location,
    }


def contiguous_error_regions(errors: pd.DataFrame, tm_columns: list[str], error_floor: float) -> pd.DataFrame:
    """Collapse adjacent high-error or shift-flagged residues into regions."""
    regions: list[dict[str, Any]] = []
    for (target, neighbour, helix), group in errors.groupby(["target", "neighbour", "helix"]):
        columns = [column for column in tm_columns if column.startswith(f"{helix}.")]
        suspect = group[(group["same_grn_error"] >= error_floor) | group["possible_shift"]].copy()
        suspect["position"] = suspect["grn"].map(columns.index)
        suspect = suspect.sort_values("position")
        current: list[Any] = []
        for row in suspect.itertuples(index=False):
            if current and row.position != current[-1].position + 1:
                regions.append(make_error_region(target, neighbour, helix, current))
                current = []
            current.append(row)
        if current:
            regions.append(make_error_region(target, neighbour, helix, current))
    return pd.DataFrame(regions)


def coherent_shift_regions(
    errors: pd.DataFrame,
    tm_columns: list[str],
    minimum_residues: int,
    maximum_gap: int,
) -> pd.DataFrame:
    """Find repeated nearby-GRN offsets across a local helix segment."""
    regions: list[dict[str, Any]] = []
    flagged = errors[errors["possible_shift"]].copy()
    for (target, neighbour, helix, shift), group in flagged.groupby(
        ["target", "neighbour", "helix", "shift_columns"]
    ):
        columns = [column for column in tm_columns if column.startswith(f"{helix}.")]
        group["position"] = group["grn"].map(columns.index)
        group = group.sort_values("position")
        current: list[Any] = []
        for row in group.itertuples(index=False):
            if current and row.position - current[-1].position > maximum_gap:
                if len(current) >= minimum_residues:
                    regions.append(make_coherent_shift_region(target, neighbour, helix, int(shift), current))
                current = []
            current.append(row)
        if len(current) >= minimum_residues:
            regions.append(make_coherent_shift_region(target, neighbour, helix, int(shift), current))
    return pd.DataFrame(regions)


def make_coherent_shift_region(target: str, neighbour: str, helix: int, shift: int, rows: list[Any]) -> dict[str, Any]:
    return {
        "target": target,
        "neighbour": neighbour,
        "helix": helix,
        "shift_columns": shift,
        "start_grn": rows[0].grn,
        "end_grn": rows[-1].grn,
        "span_columns": int(rows[-1].position - rows[0].position + 1),
        "flagged_residues": len(rows),
        "mean_same_grn_error": float(np.mean([row.same_grn_error for row in rows])),
        "mean_shift_improvement": float(np.mean([row.shift_improvement for row in rows])),
        "max_shift_improvement": float(np.max([row.shift_improvement for row in rows])),
    }


def make_error_region(target: str, neighbour: str, helix: int, rows: list[Any]) -> dict[str, Any]:
    return {
        "target": target,
        "neighbour": neighbour,
        "helix": helix,
        "start_grn": rows[0].grn,
        "end_grn": rows[-1].grn,
        "length": len(rows),
        "mean_same_grn_error": float(np.mean([row.same_grn_error for row in rows])),
        "max_same_grn_error": float(np.max([row.same_grn_error for row in rows])),
        "n_possible_shifts": int(sum(row.possible_shift for row in rows)),
    }


def summarize(
    nearest: pd.DataFrame,
    errors: pd.DataFrame,
    gaps: pd.DataFrame,
    gap_runs: pd.DataFrame,
    fits: pd.DataFrame,
    regions: pd.DataFrame,
    coherent_shifts: pd.DataFrame,
) -> pd.DataFrame:
    records: list[dict[str, Any]] = []
    nearest_index = nearest.set_index("target")
    fit_index = fits.set_index("target")
    for target, group in errors.groupby("target"):
        n = nearest_index.loc[target]
        f = fit_index.loc[target]
        target_gaps = gaps[gaps["target"] == target] if not gaps.empty else gaps
        target_runs = gap_runs[gap_runs["target"] == target] if not gap_runs.empty else gap_runs
        target_regions = regions[regions["target"] == target] if not regions.empty else regions
        target_coherent = coherent_shifts[coherent_shifts["target"] == target] if not coherent_shifts.empty else coherent_shifts
        records.append(
            {
                "target": target,
                "nearest_neighbour": n["neighbour"],
                "sequence_identity": n["identity"],
                "common_tm_residues": n["common_residues"],
                "gap_disagreements": len(target_gaps),
                "internal_gap_runs": int((target_runs["location"] == "internal").sum()) if not target_runs.empty else 0,
                "max_gap_run": int(target_runs["length"].max()) if not target_runs.empty else 0,
                "bundle_rmsd": f["bundle_rmsd"],
                "aligned_bundle_residues": f["aligned_bundle_residues"],
                "n_common_grns": len(group),
                "median_same_grn_error": group["same_grn_error"].median(),
                "p95_same_grn_error": group["same_grn_error"].quantile(0.95),
                "max_same_grn_error": group["same_grn_error"].max(),
                "n_errors_gt_3": int((group["same_grn_error"] > 3).sum()),
                "n_errors_gt_5": int((group["same_grn_error"] > 5).sum()),
                "n_possible_shifts": int(group["possible_shift"].sum()),
                "max_shift_improvement": group["shift_improvement"].max(),
                "n_suspect_regions": len(target_regions),
                "max_suspect_region": int(target_regions["length"].max()) if not target_regions.empty else 0,
                "n_coherent_shift_regions": len(target_coherent),
                "coherent_shift_residues": int(target_coherent["flagged_residues"].sum()) if not target_coherent.empty else 0,
                "max_coherent_shift_span": int(target_coherent["span_columns"].max()) if not target_coherent.empty else 0,
            }
        )
    result = pd.DataFrame(records)
    if result.empty:
        return result
    # Shift evidence is most specific to numbering mistakes; generic structural
    # divergence and gap differences remain visible but rank below it.
    result["audit_score"] = (
        8 * result["coherent_shift_residues"]
        + 4 * result["n_coherent_shift_regions"]
        + result["n_possible_shifts"]
        + 2 * result["internal_gap_runs"]
        + result["n_errors_gt_5"]
        + 0.25 * result["n_errors_gt_3"]
        + 0.5 * result["max_suspect_region"]
    )
    return result.sort_values(
        ["audit_score", "coherent_shift_residues", "n_possible_shifts", "p95_same_grn_error"], ascending=False
    )


def write_report(output: Path, metadata: dict[str, Any], summary: pd.DataFrame, shifts: pd.DataFrame, gap_runs: pd.DataFrame) -> None:
    top = summary.head(20)
    report_columns = [
        "target", "nearest_neighbour", "sequence_identity", "bundle_rmsd",
        "n_coherent_shift_regions", "coherent_shift_residues", "internal_gap_runs",
        "p95_same_grn_error", "audit_score",
    ]
    header = "| " + " | ".join(report_columns) + " |"
    separator = "| " + " | ".join("---" for _ in report_columns) + " |"
    rows = []
    for values in top[report_columns].itertuples(index=False, name=None):
        formatted = [f"{value:.3f}" if isinstance(value, float) else str(value) for value in values]
        rows.append("| " + " | ".join(formatted) + " |")
    markdown_table = "\n".join([header, separator, *rows])
    lines = [
        "# Type-I nearest-neighbour GRN audit",
        "",
        "This audit compares every curated row with its closest TM-sequence row; it does not classify structures as experimental or predicted.",
        "",
        f"- Rows audited: {metadata['rows_audited']} / {metadata['reference_rows']}",
        f"- Successful TM-bundle pair alignments: {metadata['successful_pair_alignments']}",
        f"- Possible local GRN shifts: {metadata['possible_shifts']}",
        f"- Coherent multi-residue shift regions: {metadata['coherent_shift_regions']}",
        f"- Table/structure residue-letter mismatches: {metadata['table_structure_residue_mismatches']}",
        f"- Gap-disagreement runs: {metadata['gap_runs']}",
        f"- Curated internal gap cells/runs: {metadata['curated_internal_gap_cells']} / {metadata['curated_internal_gap_runs']}",
        "",
        "## Highest-priority rows",
        "",
        markdown_table,
        "",
        "## Interpretation",
        "",
        "`possible_shift` means the target CA is substantially closer to a nearby GRN in its neighbour than to the identically labelled GRN after whole-bundle alignment. Treat this as a visual-inspection candidate, not an automatic correction. Gap runs expose insertion/deletion and helix-boundary disagreements. A high bundle RMSD with few shift flags more likely reflects genuine structural divergence than a table error.",
        "",
        "## Files",
        "",
        "- `protein_audit_summary.csv`: ranked rows",
        "- `possible_grn_shifts.csv`: residue-level shift candidates",
        "- `coherent_shift_regions.csv`: repeated register shifts (highest-specificity candidates)",
        "- `suspect_regions.csv`: contiguous high-error regions",
        "- `gap_runs.csv` and `gap_disagreements.csv`: insertion/deletion patterns",
        "- `curated_internal_gap_runs.csv`: gaps bounded by residues, with squash anchors",
        "- `curated_internal_gap_positions.csv`: GRN-wise list of structures containing each gap",
        "- `per_grn_pair_errors.csv`: all same-GRN residuals",
        "- `nearest_neighbors.csv`: top alignment-table neighbours",
        "- `tm_bundle_pair_alignment_qc.csv`: structural-fit quality",
        "- `table_structure_residue_mismatches.csv`: direct residue-letter discrepancies",
    ]
    (output / "AUDIT_REPORT.md").write_text("\n".join(lines) + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-csv", type=Path, default=PROJECT_ROOT / "type_I.csv")
    parser.add_argument("--output-dir", type=Path, default=PROJECT_ROOT / "opsin_output" / "type_i_nearest_neighbor_audit")
    parser.add_argument("--datasets", nargs="+", default=DEFAULT_DATASETS)
    parser.add_argument("--chain-id", default="A")
    parser.add_argument("--min-overlap", type=int, default=100)
    parser.add_argument("--top-k", type=int, default=3)
    parser.add_argument("--local-radius", type=int, default=5)
    parser.add_argument("--shift-error-floor", type=float, default=3.0)
    parser.add_argument("--shift-improvement-floor", type=float, default=1.5)
    parser.add_argument("--region-error-floor", type=float, default=3.0)
    parser.add_argument("--coherent-shift-min-residues", type=int, default=3)
    parser.add_argument("--coherent-shift-max-gap", type=int, default=2)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table = pd.read_csv(args.reference_csv, index_col=0, dtype=str)
    table.index = table.index.map(str)
    if not table.index.is_unique:
        raise ValueError("Reference table index is not unique")
    tm_columns = [str(column) for column in table.columns if TM_PATTERN.fullmatch(str(column))]
    if not tm_columns:
        raise ValueError("No helix 1-7 GRN columns found")

    all_neighbours = find_nearest_neighbours(table, tm_columns, args.min_overlap, args.top_k)
    all_neighbours.to_csv(args.output_dir / "nearest_neighbors.csv", index=False)
    nearest = all_neighbours[all_neighbours["rank"] == 1].copy()
    if len(nearest) != len(table):
        missing = sorted(set(table.index) - set(nearest["target"]))
        raise ValueError(f"No qualifying neighbour for: {missing}")

    import protos
    from protos.processing.structure import StructureProcessor

    protos.set_data_path(str(PROJECT_ROOT / "data"))
    processor = StructureProcessor("type_i_nearest_neighbor_audit")
    id_map = load_structure_id_map(processor, args.datasets)
    coordinate_tables, mapping_qc, residue_mismatches = map_table_rows_to_coordinates(
        table, tm_columns, processor, id_map, args.chain_id
    )
    mapping_qc.to_csv(args.output_dir / "structure_mapping_qc.csv", index=False)
    residue_mismatches.to_csv(args.output_dir / "table_structure_residue_mismatches.csv", index=False)

    errors: list[pd.DataFrame] = []
    fits: list[dict[str, Any]] = []
    for number, pair in enumerate(nearest.itertuples(index=False), start=1):
        if pair.target not in coordinate_tables or pair.neighbour not in coordinate_tables:
            fits.append({"target": pair.target, "neighbour": pair.neighbour, "status": "coordinates_missing"})
            continue
        try:
            pair_errors, fit = structural_pair_audit(
                pair.target,
                pair.neighbour,
                coordinate_tables,
                tm_columns,
                args.local_radius,
                args.shift_error_floor,
                args.shift_improvement_floor,
            )
            errors.append(pair_errors)
            fits.append(fit)
        except Exception as exc:
            fits.append({"target": pair.target, "neighbour": pair.neighbour, "status": "alignment_failed", "error": str(exc)})
        if number % 20 == 0 or number == len(nearest):
            print(f"[pair alignment] {number}/{len(nearest)}")

    error_table = pd.concat(errors, ignore_index=True) if errors else pd.DataFrame()
    fit_table = pd.DataFrame(fits)
    error_table.to_csv(args.output_dir / "per_grn_pair_errors.csv", index=False)
    fit_table.to_csv(args.output_dir / "tm_bundle_pair_alignment_qc.csv", index=False)
    shifts = error_table[error_table["possible_shift"]].sort_values(
        ["shift_improvement", "same_grn_error"], ascending=False
    )
    shifts.to_csv(args.output_dir / "possible_grn_shifts.csv", index=False)

    gap_cells, gap_runs = gap_disagreements(table, nearest, tm_columns)
    gap_cells.to_csv(args.output_dir / "gap_disagreements.csv", index=False)
    gap_runs.sort_values(["length", "target"], ascending=[False, True]).to_csv(args.output_dir / "gap_runs.csv", index=False)
    internal_gap_cells, internal_gap_runs, internal_gap_positions = curated_internal_gap_inventory(
        table, tm_columns
    )
    internal_gap_cells.to_csv(args.output_dir / "curated_internal_gap_cells.csv", index=False)
    internal_gap_runs.to_csv(args.output_dir / "curated_internal_gap_runs.csv", index=False)
    internal_gap_positions.to_csv(args.output_dir / "curated_internal_gap_positions.csv", index=False)
    regions = contiguous_error_regions(error_table, tm_columns, args.region_error_floor)
    regions.sort_values(["n_possible_shifts", "length", "max_same_grn_error"], ascending=False).to_csv(
        args.output_dir / "suspect_regions.csv", index=False
    )
    coherent_shifts = coherent_shift_regions(
        error_table,
        tm_columns,
        args.coherent_shift_min_residues,
        args.coherent_shift_max_gap,
    )
    if not coherent_shifts.empty:
        coherent_shifts = coherent_shifts.sort_values(
            ["flagged_residues", "mean_shift_improvement"], ascending=False
        )
    coherent_shifts.to_csv(args.output_dir / "coherent_shift_regions.csv", index=False)
    successful_fits = fit_table[fit_table["status"] == "ok"]
    successful_targets = set(successful_fits["target"])
    summary = summarize(
        nearest[nearest["target"].isin(successful_targets)],
        error_table,
        gap_cells,
        gap_runs,
        successful_fits,
        regions,
        coherent_shifts,
    )
    summary.to_csv(args.output_dir / "protein_audit_summary.csv", index=False)

    metadata = {
        "reference_csv": str(args.reference_csv.resolve()),
        "reference_sha256": sha256_file(args.reference_csv),
        "reference_rows": len(table),
        "tm_columns": len(tm_columns),
        "rows_with_mapped_structures": len(coordinate_tables),
        "table_structure_residue_mismatches": len(residue_mismatches),
        "rows_audited": len(summary),
        "successful_pair_alignments": len(successful_fits),
        "possible_shifts": len(shifts),
        "coherent_shift_regions": len(coherent_shifts),
        "gap_disagreements": len(gap_cells),
        "gap_runs": len(gap_runs),
        "curated_internal_gap_cells": len(internal_gap_cells),
        "curated_internal_gap_runs": len(internal_gap_runs),
        "structures_with_curated_internal_gaps": int(internal_gap_cells["structure"].nunique()) if not internal_gap_cells.empty else 0,
        "parameters": {
            "min_overlap": args.min_overlap,
            "top_k": args.top_k,
            "local_radius": args.local_radius,
            "shift_error_floor": args.shift_error_floor,
            "shift_improvement_floor": args.shift_improvement_floor,
            "region_error_floor": args.region_error_floor,
            "coherent_shift_min_residues": args.coherent_shift_min_residues,
            "coherent_shift_max_gap": args.coherent_shift_max_gap,
        },
        "method": "nearest TM-sequence row; Protos whole-TM-bundle CE alignment; same/nearby-GRN CA residuals",
        "category_comparison": "none",
    }
    with (args.output_dir / "run_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2)
    write_report(args.output_dir, metadata, summary, shifts, gap_runs)
    print(json.dumps(metadata, indent=2))
    return 0 if len(summary) == len(table) else 1


if __name__ == "__main__":
    raise SystemExit(main())
