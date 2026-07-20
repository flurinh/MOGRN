"""Load the authoritative ST5 metadata without initializing ProtOS."""

from __future__ import annotations

import re
from pathlib import Path

import pandas as pd


def load_st5_property_data(
    property_file: Path,
    legacy_property_file: Path,
) -> pd.DataFrame:
    """Load ST5 properties and add operational historical split metadata."""

    property_file = Path(property_file)
    legacy_property_file = Path(legacy_property_file)
    print(f"[INFO] Loading property data from: {property_file}")

    if not property_file.exists():
        raise FileNotFoundError(f"Property file not found: {property_file}")

    if property_file.suffix.casefold() in {".xlsx", ".xls"}:
        frame = pd.read_excel(property_file, dtype=str)
    else:
        frame = pd.read_csv(property_file, dtype={"PDB ID": str})
    print(f"[INFO] Loaded {len(frame)} entries from property file")

    frame.columns = frame.columns.str.strip()
    short_names = frame["short_name"].fillna("").astype(str).str.strip()
    if short_names.duplicated().any():
        raise ValueError("Property table has duplicate short_name values")

    legacy_splits = {}
    if legacy_property_file.is_file() and property_file != legacy_property_file:
        legacy = pd.read_csv(legacy_property_file, dtype=str).fillna("")
        legacy_splits = {
            str(row["short_name"]).strip(): str(row["dataset_split"]).strip()
            for _, row in legacy.iterrows()
            if str(row.get("short_name", "")).strip()
        }

    def fix_pdb_id(value):
        if pd.isna(value):
            return value
        text = str(value).strip()
        if "E+" in text.upper() or "E-" in text.upper():
            match = re.match(r"(\d+)\.?\d*[Ee][+]?(\d+)", text)
            if match:
                return f"{match.group(1)}E{match.group(2)}"
        match = re.search(r"\b([0-9][A-Za-z0-9]{3})\b", text)
        return match.group(1).upper() if match else ""

    if "PDB ID" in frame.columns:
        frame["PDB ID"] = frame["PDB ID"].apply(fix_pdb_id)

    frame["experimentally_determined"] = (
        frame["PDB ID"].fillna("").astype(str).str.strip().ne("").astype(int)
    )
    frame["dataset_split"] = (
        frame["short_name"].map(legacy_splits).fillna("unknown").astype(str).str.strip()
    )
    new_experimental = (
        frame["experimentally_determined"].eq(1)
        & ~frame["dataset_split"].isin({"A", "B"})
    )
    frame.loc[new_experimental, "dataset_split"] = "B"

    frame["exp_structure_id"] = frame["PDB ID"].apply(
        lambda value: (
            str(value).strip().lower()
            if pd.notna(value) and str(value).strip()
            else ""
        )
    )
    frame["pred_structure_id"] = frame["short_name"].apply(
        lambda value: (
            f"{str(value).strip()}_model_0"
            if pd.notna(value) and str(value).strip()
            else ""
        )
    )

    exp_count = int(frame["experimentally_determined"].eq(1).sum())
    pred_count = int(frame["experimentally_determined"].eq(0).sum())
    set_a_count = int(
        (
            frame["experimentally_determined"].eq(1)
            & frame["dataset_split"].eq("A")
        ).sum()
    )
    set_b_count = int(
        (
            frame["experimentally_determined"].eq(1)
            & frame["dataset_split"].eq("B")
        ).sum()
    )
    print(f"[INFO] Total entries: {len(frame)}")
    print(f"[INFO] Experimental entries: {exp_count}")
    print(f"[INFO]   - Set A (pre-Sept 2021): {set_a_count}")
    print(f"[INFO]   - Set B (post-Sept 2021): {set_b_count}")
    print(f"[INFO] Novel/undetermined entries: {pred_count}")

    return frame
