#!/usr/bin/env python3
"""Augment opsin property data with RCSB release metadata.

The script reads the input CSV (default: property/mo_exp.csv), queries the RCSB
REST API for each available PDB identifier, and adds the columns:
  - reference_date: ISO release date from RCSB (initial_release_date fallback)
  - reference_year: integer year (from release date or textual fallback)
  - dataset_split: Set A/B tag using the fixed cutoff date 2021-09-30

If the release date cannot be recovered from RCSB, the script optionally falls
back to parsing the existing "reference" column for a four-digit year.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple
from datetime import datetime, date

import pandas as pd

RCSB_ENTRY_URL = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
CUTOFF_DATE_STR = "2021-09-30"
CUTOFF_DATE = datetime.strptime(CUTOFF_DATE_STR, '%Y-%m-%d').date()
YEAR_PATTERN = re.compile(r"(19|20)\d{2}")


def fetch_rcsb_release(pdb_id: str, timeout: float) -> Tuple[Optional[str], Optional[date]]:
    """Return (release_date_iso, release_date_obj) for a PDB identifier."""
    url = RCSB_ENTRY_URL.format(pdb_id=pdb_id.lower())
    request = urllib.request.Request(url, headers={"User-Agent": "MOGRN-data-updater"})

    try:
        with urllib.request.urlopen(request, timeout=timeout) as handle:
            payload = json.load(handle)
    except urllib.error.HTTPError as exc:  # 404, etc.
        sys.stderr.write(f"[WARN] RCSB lookup failed for {pdb_id}: {exc.code}\n")
        return None, None
    except urllib.error.URLError as exc:
        sys.stderr.write(f"[WARN] Network error for {pdb_id}: {exc}\n")
        return None, None

    accession = payload.get("rcsb_accession_info", {})
    release_date = accession.get("initial_release_date") or accession.get("deposit_date")

    if not release_date:
        return None, None

    date_obj = None
    try:
        date_obj = datetime.strptime(release_date[:10], '%Y-%m-%d').date()
    except (TypeError, ValueError):
        pass
    return release_date, date_obj


def extract_pdb_ids(raw_value: object) -> Iterable[str]:
    """Yield cleaned PDB identifiers from a CSV field."""
    if raw_value is None or (isinstance(raw_value, float) and pd.isna(raw_value)):
        return []

    candidates = str(raw_value).replace(";", ",").split(",")
    for candidate in candidates:
        candidate = candidate.strip()
        if len(candidate) == 4 and candidate.isalnum():
            yield candidate.upper()


def parse_year_from_reference(reference: object) -> Optional[int]:
    if reference is None or (isinstance(reference, float) and pd.isna(reference)):
        return None
    match = YEAR_PATTERN.search(str(reference))
    if not match:
        return None
    try:
        return int(match.group(0))
    except ValueError:
        return None


def load_hideaki_ids(path: Optional[Path]) -> Dict[str, None]:
    if path is None or not path.exists():
        return {}
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}

    if isinstance(payload, dict) and "pdb_ids" in payload:
        return {pdb_id.upper(): None for pdb_id in payload["pdb_ids"]}
    return {}


def assign_dataset_split(release_date: Optional[date], fallback_year: Optional[int], force_b: bool) -> str:
    if force_b:
        return "B"
    if release_date is not None:
        return "A" if release_date <= CUTOFF_DATE else "B"
    if fallback_year is not None:
        # Assume end-of-year when only year is available
        fallback_date = date(fallback_year, 12, 31)
        return "A" if fallback_date <= CUTOFF_DATE else "B"
    return "unknown"


def main() -> None:
    parser = argparse.ArgumentParser(description="Add release metadata to opsin property CSV")
    parser.add_argument("--input", default="property/mo_exp.csv", type=Path)
    parser.add_argument("--outputs", type=Path, default=None,
                        help="Path for augmented CSV (defaults to overwriting input)")
    parser.add_argument("--hideaki-dataset", type=Path, default=Path("data/structure/structure_dataset/standard/hideaki_exp.json"),
                        help="JSON file listing hideaki experimental PDB IDs")
    parser.add_argument("--sleep", type=float, default=0.2,
                        help="Delay between RCSB requests to avoid throttling")
    parser.add_argument("--timeout", type=float, default=10.0,
                        help="HTTP timeout per RCSB request (seconds)")
    parser.add_argument("--no-reference-fallback", action="store_true",
                        help="Disable parsing publication year from the reference column")

    args = parser.parse_args()

    if not args.input.exists():
        parser.error(f"Input file not found: {args.input}")

    df = pd.read_csv(args.input)
    hideaki_ids = load_hideaki_ids(args.hideaki_dataset.resolve())

    cache: Dict[str, Tuple[Optional[str], Optional[date]]] = {}
    reference_dates = []
    reference_years = []
    dataset_split = []

    for idx, row in df.iterrows():
        pdb_candidates = list(extract_pdb_ids(row.get("PDB ID")))
        if not pdb_candidates:
            pdb_candidates = list(extract_pdb_ids(row.get("pdb_id")))

        release_date = None
        release_date_obj: Optional[date] = None

        for pdb_id in pdb_candidates:
            if pdb_id not in cache:
                cache[pdb_id] = fetch_rcsb_release(pdb_id, timeout=args.timeout)
                if args.sleep:
                    time.sleep(args.sleep)
            release_date, release_date_obj = cache[pdb_id]
            if release_date and release_date_obj:
                break

        release_year_value: Optional[int] = release_date_obj.year if release_date_obj else None
        if release_year_value is None and release_date:
            try:
                release_year_value = int(str(release_date)[:4])
            except ValueError:
                release_year_value = None

        fallback_year = release_year_value
        if fallback_year is None and not args.no_reference_fallback:
            fallback_year = parse_year_from_reference(row.get("reference"))

        reference_dates.append(release_date)
        reference_years.append(fallback_year)

        pdb_for_split = pdb_candidates[0] if pdb_candidates else None
        force_b = bool(pdb_for_split and pdb_for_split.upper() in hideaki_ids)
        dataset_split.append(assign_dataset_split(release_date_obj, fallback_year, force_b))

    df["reference_date"] = reference_dates
    df["reference_year"] = pd.Series(reference_years, dtype='Int64')
    df["dataset_split"] = pd.Series(dataset_split, dtype='string')

    output_path = args.output or args.input
    df.to_csv(output_path, index=False)
    print(f"Wrote augmented dataset to {output_path}")


if __name__ == "__main__":
    main()
