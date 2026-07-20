#!/usr/bin/env python3
"""Apply explicit GRN register corrections from manual visual inspection.

The corrections in this module are curated observations, not automatic
alignment decisions.  Each transformation checks its expected before/after
state so reruns are idempotent and unexpected table drift fails loudly.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
PSCHR2 = "pschr2_model_0"
HULACCR1 = "hulaccr1_model_0"


def _residue_number(value: str) -> int:
    match = re.fullmatch(r"[A-Z](\d+)", value)
    if match is None:
        raise ValueError(f"Expected a numbered residue, found {value!r}")
    return int(match.group(1))


def _expect_consecutive(values: list[str], start: int, context: str) -> None:
    observed = [_residue_number(value) for value in values]
    expected = list(range(start, start + len(values)))
    if observed != expected:
        raise ValueError(f"Unexpected {context} register: {observed}; expected {expected}")


def _curate_pschr2_helix_1(row: pd.Series) -> pd.Series:
    """Shift PsChR2 H1 one GRN later with no internal register gap."""

    corrected = row.copy()
    source_columns = [f"1.{position}" for position in range(30, 64)]
    target_columns = [f"1.{position}" for position in range(32, 64)]
    target_values = [corrected[column] for column in target_columns]

    if all(value != "-" for value in target_values):
        _expect_consecutive(target_values, 2, "curated PsChR2 helix-1")
    else:
        source_values = [
            corrected[column]
            for column in source_columns
            if corrected[column] != "-"
        ]
        by_number = {_residue_number(value): value for value in source_values}
        if sorted(by_number) != list(range(2, 34)):
            raise ValueError(
                "Unexpected PsChR2 helix-1 residues: "
                f"{sorted(by_number)}; expected 2--33"
            )

        for column in source_columns:
            corrected[column] = "-"
        for target, residue_number in zip(target_columns, range(2, 34)):
            corrected[target] = by_number[residue_number]

    expected_anchors = {
        "1.30": "-",
        "1.31": "-",
        "1.32": "T2",
        "1.44": "W14",
        "1.45": "Q15",
        "1.46": "L16",
        "1.47": "A17",
        "1.63": "W33",
    }
    if any(corrected[grn] != residue for grn, residue in expected_anchors.items()):
        raise ValueError("PsChR2 did not receive the requested continuous H1 shift")
    return corrected


def _curate_hulaccr1_helix_5(row: pd.Series) -> pd.Series:
    """Shift the complete HulaCCR1 helix-5 register one GRN earlier."""

    corrected = row.copy()
    helix_columns = [f"5.{position}" for position in range(42, 73)]
    loop_columns = ["56.001", "56.002", "56.003"]
    curated_values = [corrected[column] for column in helix_columns] + [
        corrected[column] for column in loop_columns[:2]
    ]

    if (
        corrected["5.451"] == "-"
        and corrected["56.003"] == "-"
        and all(value != "-" for value in curated_values)
    ):
        _expect_consecutive(curated_values, 219, "curated HulaCCR1 helix-5")
    else:
        source_columns = helix_columns + ["5.451"] + loop_columns
        source_values = [
            corrected[column]
            for column in source_columns
            if corrected[column] != "-"
        ]
        by_number = {_residue_number(value): value for value in source_values}
        if sorted(by_number) != list(range(219, 252)):
            raise ValueError(
                "Unexpected HulaCCR1 helix-5 residues: "
                f"{sorted(by_number)}; expected 219--251"
            )

        for target, residue_number in zip(helix_columns, range(219, 250)):
            corrected[target] = by_number[residue_number]
        corrected["5.451"] = "-"
        corrected["56.001"] = by_number[250]
        corrected["56.002"] = by_number[251]
        corrected["56.003"] = "-"

    expected_anchors = {
        "5.42": "K219",
        "5.43": "S220",
        "5.44": "T221",
        "5.45": "T222",
        "5.451": "-",
        "5.46": "D223",
        "5.71": "L248",
        "5.72": "K249",
    }
    if any(corrected[grn] != residue for grn, residue in expected_anchors.items()):
        raise ValueError("HulaCCR1 did not receive the requested continuous H5 shift")
    return corrected


def apply_visual_curations(table: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with all currently approved visual curations applied."""

    missing_rows = {PSCHR2, HULACCR1}.difference(table.index)
    if missing_rows:
        raise ValueError(f"Missing curated rows: {sorted(missing_rows)}")

    corrected = table.copy()
    corrected.loc[PSCHR2] = _curate_pschr2_helix_1(corrected.loc[PSCHR2])
    corrected.loc[HULACCR1] = _curate_hulaccr1_helix_5(corrected.loc[HULACCR1])
    return corrected


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=(
            ROOT
            / "protos/src/protos/reference_data/grn/reference/type_I_opsins.csv"
        ),
    )
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    output = args.output or args.input
    table = pd.read_csv(args.input, index_col=0, dtype=str).fillna("-")
    corrected = apply_visual_curations(table)
    changed = int(table.ne(corrected).to_numpy().sum())

    temporary = output.with_suffix(output.suffix + ".tmp")
    corrected.to_csv(temporary)
    temporary.replace(output)
    print(f"Applied {changed} cell changes to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
