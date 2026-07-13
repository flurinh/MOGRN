#!/usr/bin/env python3
"""Detect tandem rhodopsin repeats, split them, and annotate each scaffold.

Long single chains are screened for a strong non-overlapping internal repeat.
Detected repeats are split at their inferred repeat period into two virtual
structures.  The already curated scaffold in ``type_I.csv`` is converted into
a single-domain reference, after which Protos annotates both virtual domains.
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner


ROOT = Path(__file__).resolve().parent.parent
PROTOS_SRC = ROOT / "protos" / "src"
if PROTOS_SRC.exists():
    sys.path.insert(0, str(PROTOS_SRC))

DATASETS = ["mo_exp_A", "mo_exp_B", "mo_pred_exp", "mo_pred_novel"]
TM_PATTERN = re.compile(r"^([1-7])\.(\d+)$")


def residue_table(frame: pd.DataFrame, chain: str) -> pd.DataFrame:
    frame = frame.reset_index()
    atom_col = "atom_name" if "atom_name" in frame.columns else "res_atom_name"
    ca = frame[(frame[atom_col].astype(str).str.upper() == "CA") & (frame["auth_chain_id"].astype(str) == chain)].copy()
    ca = ca.sort_values(["label_seq_id", "auth_seq_id"]).drop_duplicates(
        ["auth_seq_id", "insertion"], keep="first"
    ).reset_index(drop=True)
    ca["ordinal"] = np.arange(1, len(ca) + 1)
    return ca


def chain_sequence(ca: pd.DataFrame) -> str:
    return "".join(ca["res_name1l"].fillna("X").astype(str).str.upper())


def repeat_aligner() -> PairwiseAligner:
    aligner = PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -6.0
    aligner.extend_gap_score = -0.5
    return aligner


def best_tandem_repeat(sequence: str, step: int = 25) -> dict[str, Any] | None:
    """Return the strongest repeat spanning opposite sides of a split point."""
    if len(sequence) < 500:
        return None
    aligner = repeat_aligner()
    best: dict[str, Any] | None = None
    for cut in range(250, len(sequence) - 249, step):
        alignment = aligner.align(sequence[:cut], sequence[cut:])[0]
        blocks_a = alignment.aligned[0]
        blocks_b = alignment.aligned[1]
        aligned = int(sum(min(a1 - a0, b1 - b0) for (a0, a1), (b0, b1) in zip(blocks_a, blocks_b)))
        if aligned < 100:
            continue
        matches = 0
        offsets: list[int] = []
        for (a0, a1), (b0, b1) in zip(blocks_a, blocks_b):
            n = min(a1 - a0, b1 - b0)
            for delta in range(n):
                apos = int(a0 + delta)
                bpos = int(cut + b0 + delta)
                matches += sequence[apos] == sequence[bpos]
                offsets.append(bpos - apos)
        record = {
            "repeat_score": float(alignment.score),
            "aligned_residues": aligned,
            "repeat_identity": matches / max(aligned, 1),
            "test_split": cut,
            "repeat_period": int(round(float(np.median(offsets)))),
            "domain_a_aligned_start": int(blocks_a[0][0]) + 1,
            "domain_a_aligned_end": int(blocks_a[-1][1]),
            "domain_b_aligned_start": cut + int(blocks_b[0][0]) + 1,
            "domain_b_aligned_end": cut + int(blocks_b[-1][1]),
        }
        if best is None or record["repeat_score"] > best["repeat_score"]:
            best = record
    return best


def detect_candidates(processor: Any, datasets: list[str], args: argparse.Namespace) -> tuple[pd.DataFrame, dict[str, dict[str, Any]]]:
    ids: list[str] = []
    for dataset in datasets:
        ids.extend(processor.get_dataset_entities(dataset))
    ids = list(dict.fromkeys(ids))
    records: list[dict[str, Any]] = []
    inputs: dict[str, dict[str, Any]] = {}
    for number, structure_id in enumerate(ids, start=1):
        frame = processor.load_entity(structure_id)
        if frame is None:
            continue
        reset = frame.reset_index()
        atom_col = "atom_name" if "atom_name" in reset.columns else "res_atom_name"
        ca_all = reset[reset[atom_col].astype(str).str.upper() == "CA"]
        counts = ca_all["auth_chain_id"].astype(str).value_counts()
        if counts.empty:
            continue
        chain = str(counts.index[0])
        ca = residue_table(frame, chain)
        sequence = chain_sequence(ca)
        repeat = best_tandem_repeat(sequence, args.scan_step)
        detected = bool(
            repeat
            and repeat["repeat_score"] >= args.min_repeat_score
            and repeat["aligned_residues"] >= args.min_aligned_residues
            and repeat["repeat_identity"] >= args.min_repeat_identity
        )
        record = {
            "structure": structure_id,
            "chain": chain,
            "chain_residues": len(sequence),
            "detected_dual_rhodopsin": detected,
            **(repeat or {}),
        }
        records.append(record)
        if detected:
            inputs[structure_id] = {"frame": frame, "chain": chain, "ca": ca, "sequence": sequence, **repeat}
        if number % 25 == 0 or number == len(ids):
            print(f"[scan] {number}/{len(ids)}")
    return pd.DataFrame(records), inputs


def build_single_domain_reference(
    curated: pd.DataFrame,
    anchor_id: str,
    anchor_ca: pd.DataFrame,
    domain_b_start: int,
    domain_b_end: int,
) -> pd.DataFrame:
    auth_to_ordinal = dict(zip(anchor_ca["auth_seq_id"].astype(int), anchor_ca["ordinal"].astype(int)))
    row: dict[str, str] = {}
    for grn, value in curated.loc[anchor_id].items():
        if not isinstance(value, str) or value == "-":
            row[str(grn)] = "-"
            continue
        try:
            ordinal = auth_to_ordinal[int(value[1:])]
        except (KeyError, ValueError):
            row[str(grn)] = "-"
            continue
        row[str(grn)] = (
            f"{value[0]}{ordinal - domain_b_start + 1}"
            if domain_b_start <= ordinal <= domain_b_end
            else "-"
        )
    return pd.DataFrame([row], index=[f"{anchor_id}_B"])


def split_virtual_structure(
    processor: Any,
    parent_id: str,
    info: dict[str, Any],
    label: str,
    start: int,
    end: int,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    frame = info["frame"].reset_index()
    chain = info["chain"]
    ca = info["ca"]
    selected = ca[(ca["ordinal"] >= start) & (ca["ordinal"] <= end)].copy()
    selected["local_position"] = selected["ordinal"] - start + 1
    key_to_local = {
        (int(row.auth_seq_id), str(row.insertion or "")): int(row.local_position)
        for row in selected.itertuples()
    }
    insertion = frame["insertion"].fillna("").astype(str)
    keys = list(zip(frame["auth_seq_id"].fillna(-10**9).astype(int), insertion))
    local = pd.Series([key_to_local.get(key) for key in keys], index=frame.index)
    subset = frame[(frame["auth_chain_id"].astype(str) == chain) & local.notna()].copy()
    subset["parent_auth_seq_id"] = subset["auth_seq_id"]
    subset["parent_label_seq_id"] = subset["label_seq_id"]
    subset["auth_seq_id"] = local.loc[subset.index].astype(int)
    subset["label_seq_id"] = local.loc[subset.index].astype(int)
    virtual_id = f"{parent_id}_{label}"
    subset["structure_id"] = virtual_id
    processor.save_entity(
        virtual_id,
        subset,
        metadata={
            "source": "dual_rhodopsin_split",
            "parent_structure": parent_id,
            "parent_chain": chain,
            "parent_ordinal_start": start,
            "parent_ordinal_end": end,
        },
    )
    mapping = selected[["local_position", "ordinal", "auth_seq_id", "label_seq_id", "res_name1l"]].copy()
    mapping.insert(0, "virtual_structure", virtual_id)
    mapping.insert(1, "parent_structure", parent_id)
    mapping.insert(2, "domain", label)
    return subset, mapping


def annotation_to_parent_ids(annotations: pd.DataFrame, mappings: pd.DataFrame) -> pd.DataFrame:
    result = annotations.copy()
    lookups = {
        virtual: group.set_index("local_position")
        for virtual, group in mappings.groupby("virtual_structure")
    }
    for virtual in result.index:
        lookup = lookups[virtual]
        for grn, value in result.loc[virtual].items():
            if not isinstance(value, str) or value == "-":
                continue
            local = int(value[1:])
            if local in lookup.index:
                result.at[virtual, grn] = f"{value[0]}{int(lookup.loc[local, 'auth_seq_id'])}"
    return result


def tm_bundle_qc(processor: Any, annotations: pd.DataFrame, reference_id: str) -> pd.DataFrame:
    from protos.analysis.structure.alignment import get_structure_alignment

    bundles: dict[str, np.ndarray] = {}
    for structure_id, row in annotations.iterrows():
        frame = processor.load_entity(structure_id).reset_index()
        ca = frame[frame["atom_name"].astype(str).str.upper() == "CA"].set_index("auth_seq_id")
        records: list[tuple[int, int, np.ndarray]] = []
        for column_index, (grn, value) in enumerate(row.items()):
            if not TM_PATTERN.fullmatch(str(grn)) or not isinstance(value, str) or value == "-":
                continue
            position = int(value[1:])
            if position not in ca.index:
                continue
            item = ca.loc[position]
            records.append((int(str(grn).split(".")[0]), position, item[["x", "y", "z"]].to_numpy(float)))
        records.sort(key=lambda item: (item[0], item[1]))
        bundles[structure_id] = np.vstack([item[2] for item in records])

    fixed = bundles[reference_id]
    output: list[dict[str, Any]] = []
    for structure_id, mobile in bundles.items():
        if structure_id == reference_id:
            rmsd, aligned = 0.0, len(fixed)
        else:
            _rotation, _translation, path, rmsd = get_structure_alignment(fixed, mobile)
            aligned = len(path[0])
        row = annotations.loc[structure_id]
        occupied_tm = [grn for grn, value in row.items() if TM_PATTERN.fullmatch(str(grn)) and value != "-"]
        output.append(
            {
                "structure": structure_id,
                "tm_bundle_rmsd": float(rmsd),
                "aligned_tm_residues": aligned,
                "annotated_tm_residues": len(occupied_tm),
                "annotated_helices": len({str(grn).split(".")[0] for grn in occupied_tm}),
                "grn_7_50": row.get("7.50", "-"),
            }
        )
    return pd.DataFrame(output)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference-csv", type=Path, default=ROOT / "type_I.csv")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "opsin_output" / "dual_rhodopsins")
    parser.add_argument("--anchor-id", default="7pl9")
    parser.add_argument("--scan-step", type=int, default=25)
    parser.add_argument("--min-repeat-score", type=float, default=150.0)
    parser.add_argument("--min-aligned-residues", type=int, default=200)
    parser.add_argument("--min-repeat-identity", type=float, default=0.35)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    import protos
    from protos.processing.grn import GRNProcessor
    from protos.processing.structure import StructureProcessor

    protos.set_data_path(str(ROOT / "data"))
    processor = StructureProcessor("dual_rhodopsin_domains")
    detector, candidates = detect_candidates(processor, DATASETS, args)
    detector.to_csv(args.output_dir / "dual_rhodopsin_detection.csv", index=False)
    if not candidates:
        print("No tandem rhodopsin candidates detected")
        return 1
    if args.anchor_id not in candidates:
        raise ValueError(f"Curated anchor {args.anchor_id!r} was not detected")

    anchor = candidates[args.anchor_id]
    period = int(anchor["repeat_period"])
    domain_a = (1, period)
    domain_b = (period + 1, 2 * period)
    curated = pd.read_csv(args.reference_csv, index_col=0, dtype=str)
    reference = build_single_domain_reference(curated, args.anchor_id, anchor["ca"], *domain_b)
    reference_path = ROOT / "data" / "grn" / "reference" / "type_I_dual_domain.csv"
    reference.to_csv(reference_path)
    reference.to_csv(args.output_dir / "type_I_dual_domain_reference.csv")

    sequences: dict[str, str] = {}
    mapping_parts: list[pd.DataFrame] = []
    virtual_ids: list[str] = []
    boundary_records: list[dict[str, Any]] = []
    for parent_id, info in candidates.items():
        # Candidates with the same repeat period can use the curated domain reference.
        if abs(int(info["repeat_period"]) - period) > 10:
            continue
        for label, (start, end) in {"A": domain_a, "B": domain_b}.items():
            virtual_id = f"{parent_id}_{label}"
            _frame, mapping = split_virtual_structure(processor, parent_id, info, label, start, end)
            mapping_parts.append(mapping)
            virtual_ids.append(virtual_id)
            sequences[virtual_id] = info["sequence"][start - 1 : end]
            boundary_records.append(
                {
                    "parent_structure": parent_id,
                    "virtual_structure": virtual_id,
                    "domain": label,
                    "parent_chain": info["chain"],
                    "parent_ordinal_start": start,
                    "parent_ordinal_end": end,
                    "domain_residues": end - start + 1,
                    "repeat_period": info["repeat_period"],
                }
            )
    processor.create_dataset(
        "mo_dual_rhodopsin_domains",
        virtual_ids,
        metadata={"source": "tandem_repeat_detection", "reference": "type_I_dual_domain"},
    )
    mappings = pd.concat(mapping_parts, ignore_index=True)
    mappings.to_csv(args.output_dir / "virtual_to_parent_residue_mapping.csv", index=False)
    pd.DataFrame(boundary_records).to_csv(args.output_dir / "domain_boundaries.csv", index=False)

    grn = GRNProcessor("dual_rhodopsin_domains")
    annotations, summary = grn.annotate_sequences(
        sequences,
        reference_table="type_I_dual_domain",
        protein_family="mo",
        assign_unambiguous_insertions=True,
    )
    annotations.to_csv(args.output_dir / "dual_domain_grn_annotations_local.csv")
    parent_annotations = annotation_to_parent_ids(annotations, mappings)
    parent_annotations.to_csv(args.output_dir / "dual_domain_grn_annotations_parent_ids.csv")
    with (args.output_dir / "annotation_summary.json").open("w") as handle:
        json.dump(summary, handle, indent=2)

    qc = tm_bundle_qc(processor, annotations, f"{args.anchor_id}_B")
    qc.to_csv(args.output_dir / "tm_bundle_alignment_qc.csv", index=False)
    metadata = {
        "registered_structures_scanned": len(detector),
        "long_chains_scored": int(detector["repeat_score"].notna().sum()),
        "dual_rhodopsin_parents": sorted(candidates),
        "virtual_structures": virtual_ids,
        "reference_anchor": f"{args.anchor_id}_B",
        "repeat_period": period,
        "domain_a_parent_ordinals": domain_a,
        "domain_b_parent_ordinals": domain_b,
        "all_virtual_structures_have_7_helices": bool((qc["annotated_helices"] == 7).all()),
        "all_virtual_structures_have_7_50": bool((qc["grn_7_50"] != "-").all()),
        "dataset": "mo_dual_rhodopsin_domains",
    }
    with (args.output_dir / "run_metadata.json").open("w") as handle:
        json.dump(metadata, handle, indent=2)
    print(json.dumps(metadata, indent=2))
    print(qc.to_string(index=False))
    return 0 if metadata["all_virtual_structures_have_7_helices"] and metadata["all_virtual_structures_have_7_50"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
