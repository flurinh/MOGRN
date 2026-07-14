"""Tests for GRN corrections from manual visual inspection."""

from __future__ import annotations

import pandas as pd

from scripts.apply_visual_grn_curations import (
    HULACCR1,
    PSCHR2,
    apply_visual_curations,
)


def _fixture() -> pd.DataFrame:
    columns = (
        [f"1.{position}" for position in range(30, 64)]
        + ["5.451"]
        + [f"5.{position}" for position in range(42, 73)]
        + ["56.001", "56.002", "56.003"]
    )
    table = pd.DataFrame("-", index=[PSCHR2, HULACCR1, "untouched"], columns=columns)

    ps_residues = [
        "T2", "N3", "G4", "A5", "Q6", "T7", "A8", "S9",
        "N10", "V11", "L12", "Q13", "W14", "Q15", "L16",
        "A17", "A18", "G19", "F20", "S21", "I22", "L23", "L24",
        "L25", "M26", "F27", "Y28", "A29", "Y30", "Q31", "T32",
        "W33",
    ]
    table.loc[PSCHR2, [f"1.{position}" for position in range(31, 63)]] = ps_residues

    hula_residues = [
        "K219", "S220", "T221", "T222", "D223", "V224", "Y225",
        "L226", "L227", "F228",
        "M229", "L230", "S231", "A232", "I233", "A234", "M235",
        "V236", "V237", "L238", "Y239", "S240", "L241", "M242",
        "L243", "Y244", "G245", "V246", "A247", "L248", "K249",
        "W250", "D251",
    ]
    source_columns = (
        [f"5.{position}" for position in range(43, 73)]
        + ["56.001", "56.002", "56.003"]
    )
    table.loc[HULACCR1, source_columns] = hula_residues
    return table


def test_visual_curations_apply_exact_bounded_shifts() -> None:
    source = _fixture()
    corrected = apply_visual_curations(source)

    assert corrected.loc[PSCHR2, "1.30"] == "-"
    assert corrected.loc[PSCHR2, "1.31"] == "-"
    assert corrected.loc[PSCHR2, "1.32"] == "T2"
    assert corrected.loc[PSCHR2, "1.44"] == "W14"
    assert corrected.loc[PSCHR2, "1.45"] == "Q15"
    assert corrected.loc[PSCHR2, "1.46"] == "L16"
    assert corrected.loc[PSCHR2, "1.47"] == "A17"
    assert corrected.loc[PSCHR2, "1.63"] == "W33"
    assert corrected.loc[PSCHR2, [f"1.{position}" for position in range(32, 64)]].ne("-").all()

    assert corrected.loc[HULACCR1, "5.42"] == "K219"
    assert corrected.loc[HULACCR1, "5.43"] == "S220"
    assert corrected.loc[HULACCR1, "5.44"] == "T221"
    assert corrected.loc[HULACCR1, "5.45"] == "T222"
    assert corrected.loc[HULACCR1, "5.451"] == "-"
    assert corrected.loc[HULACCR1, "5.46"] == "D223"
    assert corrected.loc[HULACCR1, "5.71"] == "L248"
    assert corrected.loc[HULACCR1, "5.72"] == "K249"
    assert corrected.loc[HULACCR1, ["56.001", "56.002", "56.003"]].tolist() == [
        "W250", "D251", "-",
    ]
    assert corrected.loc["untouched"].equals(source.loc["untouched"])


def test_visual_curations_are_idempotent() -> None:
    corrected = apply_visual_curations(_fixture())
    assert apply_visual_curations(corrected).equals(corrected)


def test_visual_curations_remove_intermediate_hula_5_451_inference() -> None:
    intermediate = apply_visual_curations(_fixture())
    intermediate.loc[HULACCR1, "5.42"] = "-"
    intermediate.loc[HULACCR1, "5.43"] = "K219"
    intermediate.loc[HULACCR1, "5.44"] = "S220"
    intermediate.loc[HULACCR1, "5.45"] = "T221"
    intermediate.loc[HULACCR1, "5.451"] = "T222"

    corrected = apply_visual_curations(intermediate)

    assert corrected.loc[HULACCR1, "5.42"] == "K219"
    assert corrected.loc[HULACCR1, "5.45"] == "T222"
    assert corrected.loc[HULACCR1, "5.451"] == "-"


def test_visual_curations_replace_intermediate_pschr2_gap() -> None:
    intermediate = _fixture()
    values = intermediate.loc[
        PSCHR2, [f"1.{position}" for position in range(31, 63)]
    ].tolist()
    intermediate.loc[PSCHR2, [f"1.{position}" for position in range(30, 64)]] = "-"
    intermediate.loc[PSCHR2, [f"1.{position}" for position in range(30, 45)]] = values[:15]
    intermediate.loc[PSCHR2, [f"1.{position}" for position in range(46, 63)]] = values[15:]

    corrected = apply_visual_curations(intermediate)

    assert corrected.loc[PSCHR2, "1.32"] == "T2"
    assert corrected.loc[PSCHR2, "1.44"] == "W14"
    assert corrected.loc[PSCHR2, "1.45"] == "Q15"
    assert corrected.loc[PSCHR2, "1.46"] == "L16"
    assert corrected.loc[PSCHR2, "1.47"] == "A17"
    assert corrected.loc[PSCHR2, "1.63"] == "W33"
