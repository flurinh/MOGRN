"""External preprocessing for tandem-domain structure entities.

The standard MOGRN pipeline expects one rhodopsin domain per structure ID.
This module splits configured long-chain parents before dataset analysis and
registers each domain as an ordinary structure entity.  Protos is used only
through its public load/save API; no annotation or analysis logic is changed.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
import pandas as pd


CONFIG_SCHEMA = "mogrn.tandem-structure-domains/v1"
MANIFEST_SCHEMA = "mogrn.tandem-structure-preprocessing/v1"
DEFAULT_MINIMUM_PARENT_RESIDUES = 500
RETINAL_NAMES = frozenset({"RET", "LIG"})


@dataclass(frozen=True)
class TandemStructureResult:
    """Registered virtual entities and their parent replacements."""

    replacements: dict[str, tuple[str, ...]]
    entity_parents: dict[str, str]
    records: tuple[dict[str, Any], ...]
    config_sha256: str

    def as_manifest(self) -> dict[str, Any]:
        return {
            "schema": MANIFEST_SCHEMA,
            "config_sha256": self.config_sha256,
            "replacements": {
                parent: list(children)
                for parent, children in self.replacements.items()
            },
            "entity_parents": self.entity_parents,
            "records": list(self.records),
        }


def preprocess_registered_tandem_structures(
    processor: Any,
    available_structure_ids: Sequence[str],
    config_path: str | Path,
    *,
    retinal_cutoff: float = 6.0,
) -> TandemStructureResult:
    """Split configured registered parents and save their virtual domains."""

    config_path = Path(config_path)
    raw_config = config_path.read_bytes()
    config = json.loads(raw_config)
    _validate_config(config)
    minimum_length = int(
        config.get("minimum_parent_residues", DEFAULT_MINIMUM_PARENT_RESIDUES)
    )
    available_lookup = {
        str(structure_id).lower(): str(structure_id)
        for structure_id in available_structure_ids
    }

    replacements: dict[str, tuple[str, ...]] = {}
    entity_parents: dict[str, str] = {}
    records: list[dict[str, Any]] = []
    for requested_parent, parent_config in config["parents"].items():
        parent_id = available_lookup.get(str(requested_parent).lower())
        if parent_id is None:
            continue
        frame = processor.load_entity(parent_id)
        if frame is None:
            raise ValueError(f"registered tandem parent {parent_id!r} could not be loaded")
        reset = frame.reset_index()
        chain = str(parent_config["chain"])
        residues = _ordered_ca_residues(reset, chain)
        if len(residues) < minimum_length:
            raise ValueError(
                f"configured tandem parent {parent_id!r} has {len(residues)} "
                f"residues, below the {minimum_length}-residue screen"
            )

        children: list[str] = []
        for domain in parent_config["domains"]:
            child_id = str(domain["id"])
            start = int(domain["ordinal_start"])
            end = int(domain["ordinal_end"])
            if not 1 <= start <= end <= len(residues):
                raise ValueError(
                    f"domain interval {start}-{end} for {child_id!r} is "
                    f"outside parent {parent_id!r} (1-{len(residues)})"
                )
            virtual_frame, split_report = split_structure_domain(
                reset,
                residues,
                parent_id=parent_id,
                child_id=child_id,
                chain=chain,
                ordinal_start=start,
                ordinal_end=end,
                retinal_cutoff=retinal_cutoff,
            )
            processor.save_entity(
                child_id,
                virtual_frame,
                metadata={
                    "source": "mogrn_tandem_structure_preprocessing",
                    "parent_structure": parent_id,
                    "parent_chain": chain,
                    "parent_ordinal_start": start,
                    "parent_ordinal_end": end,
                    "domain": str(domain["label"]),
                    "schema": MANIFEST_SCHEMA,
                },
            )
            children.append(child_id)
            entity_parents[child_id] = parent_id
            records.append(
                {
                    "parent_structure": parent_id,
                    "structure": child_id,
                    "domain": str(domain["label"]),
                    "curated_reference_id": domain.get("curated_reference_id"),
                    "parent_chain": chain,
                    "parent_ordinal_start": start,
                    "parent_ordinal_end": end,
                    "helix_definitions": domain.get("helices", {}),
                    **split_report,
                }
            )
        replacements[parent_id] = tuple(children)

    return TandemStructureResult(
        replacements=replacements,
        entity_parents=entity_parents,
        records=tuple(records),
        config_sha256=hashlib.sha256(raw_config).hexdigest(),
    )


def split_structure_domain(
    frame: pd.DataFrame,
    residues: pd.DataFrame,
    *,
    parent_id: str,
    child_id: str,
    chain: str,
    ordinal_start: int,
    ordinal_end: int,
    retinal_cutoff: float,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Return one locally renumbered domain and its nearest bound retinal."""

    selected = residues[
        residues["ordinal"].between(ordinal_start, ordinal_end)
    ].copy()
    selected["local_position"] = selected["ordinal"] - ordinal_start + 1
    residue_keys = {
        (_integer(row.auth_seq_id), _insertion(row.insertion)): int(row.local_position)
        for row in selected.itertuples()
    }
    frame_keys = list(
        zip(
            frame["auth_seq_id"].map(_integer_or_none),
            frame["insertion"].map(_insertion),
        )
    )
    local_positions = pd.Series(
        [residue_keys.get(key) for key in frame_keys], index=frame.index
    )
    protein_mask = (
        (frame["auth_chain_id"].astype(str) == chain)
        & (frame["group"].astype(str) == "ATOM")
        & local_positions.notna()
    )
    protein = frame.loc[protein_mask].copy()
    protein["parent_auth_seq_id"] = protein["auth_seq_id"]
    protein["parent_label_seq_id"] = protein["label_seq_id"]
    protein["auth_seq_id"] = local_positions.loc[protein.index].astype(int)
    protein["label_seq_id"] = local_positions.loc[protein.index].astype(int)

    ligand, ligand_distance = _nearest_bound_retinal(
        frame,
        protein,
        cutoff=retinal_cutoff,
    )
    parts = [protein]
    if not ligand.empty:
        parts.append(ligand.copy())
    virtual = pd.concat(parts, ignore_index=True, sort=False)
    virtual["structure_id"] = child_id
    if "pdb_id" in virtual.columns:
        virtual["pdb_id"] = child_id

    atom_col = "atom_name" if "atom_name" in virtual.columns else "res_atom_name"
    ca_count = int(
        (
            (virtual["group"].astype(str) == "ATOM")
            & (virtual[atom_col].astype(str).str.upper() == "CA")
        ).sum()
    )
    expected = ordinal_end - ordinal_start + 1
    if ca_count != expected:
        raise ValueError(
            f"split {child_id!r} contains {ca_count} CA residues; expected {expected}"
        )
    return virtual, {
        "parent_residues": len(residues),
        "domain_residues": ca_count,
        "retinal_atoms": len(ligand),
        "retinal_min_distance": ligand_distance,
    }


def apply_structure_replacements(
    structure_ids: Sequence[str],
    replacements: Mapping[str, Sequence[str]],
) -> list[str]:
    """Replace configured parent IDs while preserving order and uniqueness."""

    lookup = {str(parent).lower(): list(children) for parent, children in replacements.items()}
    result: list[str] = []
    seen: set[str] = set()
    for structure_id in structure_ids:
        expanded = lookup.get(str(structure_id).lower(), [str(structure_id)])
        for child in expanded:
            if child not in seen:
                result.append(child)
                seen.add(child)
    return result


def expand_structure_mapping(
    structure_mapping: Mapping[str, str],
    replacements: Mapping[str, Sequence[str]],
) -> dict[str, str]:
    """Expand parent experimental/predicted pairs into paired domains."""

    lookup = {str(parent).lower(): list(children) for parent, children in replacements.items()}
    expanded: dict[str, str] = {}
    for experimental, predicted in structure_mapping.items():
        experimental_children = lookup.get(str(experimental).lower(), [str(experimental)])
        predicted_children = lookup.get(str(predicted).lower(), [str(predicted)])
        if len(experimental_children) != len(predicted_children):
            raise ValueError(
                f"cannot pair tandem replacements for {experimental!r} and "
                f"{predicted!r}: {len(experimental_children)} vs "
                f"{len(predicted_children)} entities"
            )
        expanded.update(zip(experimental_children, predicted_children))
    return expanded


def dataset_signature(
    dataset_contents: Mapping[str, Sequence[str]],
    *,
    chain_id: str,
    retinal_name: str,
    preprocessing_signature: str | None = None,
) -> str:
    """Return a stable cache signature for the actual analysis inputs."""

    payload = {
        "schema": MANIFEST_SCHEMA,
        "chain_id": chain_id,
        "retinal_name": retinal_name,
        "preprocessing_signature": preprocessing_signature,
        "datasets": {
            name: list(entities) for name, entities in sorted(dataset_contents.items())
        },
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(encoded).hexdigest()


def update_tandem_helix_definitions(
    helix_file: str | Path,
    records: Sequence[Mapping[str, Any]],
) -> None:
    """Merge configured local helix intervals into the generated helix file."""

    helix_file = Path(helix_file)
    definitions = {}
    if helix_file.exists():
        with helix_file.open() as handle:
            definitions = json.load(handle)
    for record in records:
        helices = record.get("helix_definitions", {})
        if helices:
            definitions[str(record["structure"])] = helices
    helix_file.parent.mkdir(parents=True, exist_ok=True)
    with helix_file.open("w") as handle:
        json.dump(definitions, handle, indent=2)


def _ordered_ca_residues(frame: pd.DataFrame, chain: str) -> pd.DataFrame:
    _require_columns(
        frame,
        {
            "auth_chain_id",
            "auth_seq_id",
            "label_seq_id",
            "insertion",
            "group",
        },
    )
    atom_col = "atom_name" if "atom_name" in frame.columns else "res_atom_name"
    ca = frame[
        (frame["auth_chain_id"].astype(str) == chain)
        & (frame["group"].astype(str) == "ATOM")
        & (frame[atom_col].astype(str).str.upper() == "CA")
    ].copy()
    ca["_label_sort"] = pd.to_numeric(ca["label_seq_id"], errors="coerce")
    ca["_auth_sort"] = pd.to_numeric(ca["auth_seq_id"], errors="coerce")
    ca["_insertion_sort"] = ca["insertion"].map(_insertion)
    ca.sort_values(
        ["_label_sort", "_auth_sort", "_insertion_sort"], inplace=True
    )
    ca.drop_duplicates(["auth_seq_id", "_insertion_sort"], inplace=True)
    ca.reset_index(drop=True, inplace=True)
    ca["ordinal"] = np.arange(1, len(ca) + 1)
    return ca


def _nearest_bound_retinal(
    frame: pd.DataFrame,
    protein: pd.DataFrame,
    *,
    cutoff: float,
) -> tuple[pd.DataFrame, float | None]:
    residue_col = "res_name3l" if "res_name3l" in frame.columns else "res_name"
    candidates = frame[
        (frame["group"].astype(str) != "ATOM")
        & frame[residue_col].astype(str).str.upper().isin(RETINAL_NAMES)
    ].copy()
    if candidates.empty or protein.empty:
        return candidates.iloc[0:0], None
    group_columns = [
        column
        for column in (
            "auth_chain_id",
            "auth_seq_id",
            "label_seq_id",
            "insertion",
            residue_col,
        )
        if column in candidates.columns
    ]
    protein_coords = protein[["x", "y", "z"]].to_numpy(float)
    best_group = None
    best_distance = float("inf")
    for _key, group in candidates.groupby(group_columns, dropna=False):
        ligand_coords = group[["x", "y", "z"]].to_numpy(float)
        distances = np.linalg.norm(
            ligand_coords[:, None, :] - protein_coords[None, :, :], axis=2
        )
        distance = float(distances.min())
        if distance < best_distance:
            best_group = group
            best_distance = distance
    if best_group is None or best_distance > cutoff:
        return candidates.iloc[0:0], None
    return best_group.copy(), best_distance


def _validate_config(config: Mapping[str, Any]) -> None:
    if config.get("schema") != CONFIG_SCHEMA:
        raise ValueError(
            f"tandem structure config schema must be {CONFIG_SCHEMA!r}"
        )
    if not isinstance(config.get("parents"), Mapping):
        raise ValueError("tandem structure config must contain a parents mapping")
    for parent, entry in config["parents"].items():
        if not isinstance(entry, Mapping) or not entry.get("chain"):
            raise ValueError(f"tandem parent {parent!r} requires a chain")
        domains = entry.get("domains")
        if not isinstance(domains, list) or len(domains) != 2:
            raise ValueError(f"tandem parent {parent!r} requires exactly two domains")
        intervals = []
        for domain in domains:
            required = {"id", "label", "ordinal_start", "ordinal_end", "helices"}
            if not isinstance(domain, Mapping) or not required.issubset(domain):
                raise ValueError(f"invalid domain entry for tandem parent {parent!r}")
            start = int(domain["ordinal_start"])
            end = int(domain["ordinal_end"])
            if start < 1 or end < start:
                raise ValueError(f"invalid domain interval {start}-{end} for {parent!r}")
            helices = domain["helices"]
            if not isinstance(helices, Mapping) or set(helices) != {
                str(number) for number in range(1, 8)
            }:
                raise ValueError(
                    f"domain {domain['id']!r} requires helix definitions 1-7"
                )
            intervals.append((start, end))
        if intervals[0][1] >= intervals[1][0]:
            raise ValueError(f"tandem domains overlap for parent {parent!r}")


def _integer(value: Any) -> int:
    result = _integer_or_none(value)
    if result is None:
        raise ValueError(f"expected an integer residue identifier, found {value!r}")
    return result


def _integer_or_none(value: Any) -> int | None:
    if pd.isna(value):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _insertion(value: Any) -> str:
    return "" if pd.isna(value) else str(value).strip()


def _require_columns(frame: pd.DataFrame, required: set[str]) -> None:
    missing = sorted(required.difference(frame.columns))
    if missing:
        raise ValueError(f"structure frame is missing required columns: {missing}")
