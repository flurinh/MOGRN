#!/usr/bin/env python3
"""Test alternative TARA_A TM5 registers against the full structural panel."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "protos" / "src"))

from scripts.audit_type_i_nearest_neighbors import (  # noqa: E402
    DEFAULT_DATASETS,
    find_nearest_neighbours,
    load_structure_id_map,
    map_table_rows_to_coordinates,
    ordered_bundle,
)


def build_hypotheses() -> dict[str, dict[str, int]]:
    ordinary = [f"5.{position}" for position in range(40, 46)]
    hypotheses = {
        "compressed_no_insertion": {
            **dict(zip(ordinary, range(184, 190))),
            "5.46": 190,
        },
        "dense_insertion": {
            **dict(zip(ordinary, range(183, 189))),
            "5.451": 189,
            "5.46": 190,
        },
    }
    # With I189 designated as the insertion and F190 fixed at the downstream
    # anchor, W184--A188 occupy five of the six standard columns. Test every
    # possible gap rather than assuming its location from sequence alone.
    for gap in ordinary:
        populated_columns = [grn for grn in ordinary if grn != gap]
        hypotheses[f"insertion_gap_{gap}"] = {
            **dict(zip(populated_columns, range(184, 189))),
            "5.451": 189,
            "5.46": 190,
        }
    return hypotheses


HYPOTHESES = build_hypotheses()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--reference-csv",
        type=Path,
        default=(
            ROOT
            / "protos/src/protos/reference_data/grn/reference/type_I_opsins.csv"
        ),
    )
    parser.add_argument(
        "--output-dir", type=Path,
        default=ROOT / "opsin_output" / "dual_rhodopsins" / "tm5_consensus",
    )
    parser.add_argument("--max-rmsd", type=float, default=3.0)
    parser.add_argument("--min-coverage", type=float, default=0.70)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    import protos
    from protos.analysis.structure.alignment import get_structure_alignment
    from protos.processing.structure import StructureProcessor

    protos.set_data_path(str(ROOT / "data"))
    table = pd.read_csv(args.reference_csv, index_col=0, dtype=str).fillna("-")
    tm_columns = [
        str(column) for column in table.columns
        if str(column).split(".", 1)[0] in set("1234567")
        and str(column).split(".", 1)[1].isdigit()
    ]
    processor = StructureProcessor("tara_tm5_consensus")
    id_map = load_structure_id_map(processor, DEFAULT_DATASETS)
    mapped, mapping_qc, _ = map_table_rows_to_coordinates(
        table, tm_columns, processor, id_map, preferred_chain="A"
    )
    mapping_qc.to_csv(args.output_dir / "mapping_qc.csv", index=False)

    target = mapped["TARA_A"]
    target_bundle = ordered_bundle(target)
    target_coords = target_bundle[["x", "y", "z"]].to_numpy(float)
    target_by_auth = target.drop_duplicates("auth_seq_id").set_index("auth_seq_id")
    fits: list[dict] = []
    errors: list[dict] = []
    correspondences: list[dict] = []

    for neighbour, neighbour_table in mapped.items():
        if neighbour == "TARA_A" or neighbour_table.empty:
            continue
        fixed_table = ordered_bundle(neighbour_table)
        fixed = fixed_table[["x", "y", "z"]].to_numpy(float)
        try:
            rotation, translation, path, rmsd = get_structure_alignment(fixed, target_coords)
        except Exception as exc:
            fits.append({"neighbour": neighbour, "status": "failed", "error": str(exc)})
            continue
        path_fixed, path_target = np.asarray(path[0]), np.asarray(path[1])
        coverage = len(path_target) / max(len(target_bundle), 1)
        fits.append({
            "neighbour": neighbour, "status": "ok", "bundle_rmsd": float(rmsd),
            "aligned_residues": len(path_target), "target_coverage": coverage,
        })
        if rmsd > args.max_rmsd or coverage < args.min_coverage:
            continue
        transformed = (
            target_by_auth[["x", "y", "z"]].to_numpy(float) @ np.asarray(rotation)
            + np.asarray(translation)
        )
        transformed_by_auth = dict(zip(target_by_auth.index.astype(int), transformed))
        neighbour_by_grn = neighbour_table.drop_duplicates("grn").set_index("grn")

        for hypothesis, register in HYPOTHESES.items():
            for grn, auth_seq_id in register.items():
                if grn not in neighbour_by_grn.index or auth_seq_id not in transformed_by_auth:
                    continue
                neighbour_coord = neighbour_by_grn.loc[grn, ["x", "y", "z"]].to_numpy(float)
                error = float(np.linalg.norm(transformed_by_auth[auth_seq_id] - neighbour_coord))
                errors.append({
                    "neighbour": neighbour, "hypothesis": hypothesis, "grn": grn,
                    "tara_auth_seq_id": auth_seq_id, "ca_error": error,
                    "bundle_rmsd": float(rmsd), "target_coverage": coverage,
                })

        fixed_rows = fixed_table.iloc[path_fixed].reset_index(drop=True)
        target_rows = target_bundle.iloc[path_target].reset_index(drop=True)
        for target_row, fixed_row in zip(target_rows.itertuples(), fixed_rows.itertuples()):
            auth_id = int(target_row.auth_seq_id)
            if 180 <= auth_id <= 195 and int(fixed_row.helix) == 5:
                correspondences.append({
                    "neighbour": neighbour, "tara_auth_seq_id": auth_id,
                    "neighbour_grn": fixed_row.grn,
                    "neighbour_auth_seq_id": int(fixed_row.auth_seq_id),
                    "bundle_rmsd": float(rmsd), "target_coverage": coverage,
                })

    fit_df = pd.DataFrame(fits).sort_values(["status", "bundle_rmsd"], na_position="last")
    error_df = pd.DataFrame(errors)
    correspondence_df = pd.DataFrame(correspondences)
    fit_df.to_csv(args.output_dir / "all_pair_fits.csv", index=False)
    error_df.to_csv(args.output_dir / "per_hypothesis_ca_errors.csv", index=False)
    correspondence_df.to_csv(args.output_dir / "alignment_path_correspondences.csv", index=False)

    neighbour_scores = (
        error_df.groupby(["hypothesis", "neighbour"])
        .agg(mean_ca_error=("ca_error", "mean"), median_ca_error=("ca_error", "median"), positions=("ca_error", "size"))
        .reset_index()
    )
    summary = (
        neighbour_scores.groupby("hypothesis")
        .agg(
            structures=("neighbour", "nunique"),
            median_of_structure_means=("mean_ca_error", "median"),
            mean_of_structure_means=("mean_ca_error", "mean"),
            median_of_structure_medians=("median_ca_error", "median"),
        )
        .sort_values("median_of_structure_means")
    )
    neighbour_scores.to_csv(args.output_dir / "per_structure_hypothesis_scores.csv", index=False)
    summary.to_csv(args.output_dir / "hypothesis_summary.csv")

    if not correspondence_df.empty:
        path_consensus = (
            correspondence_df.groupby(["tara_auth_seq_id", "neighbour_grn"])
            .size().rename("count").reset_index()
            .sort_values(["tara_auth_seq_id", "count"], ascending=[True, False])
        )
        path_consensus.to_csv(args.output_dir / "alignment_path_consensus.csv", index=False)
    print(summary.to_string())
    print(f"Accepted structural comparisons: {error_df['neighbour'].nunique()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
