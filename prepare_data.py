#!/usr/bin/env python3
"""
Data preparation script for MOGRN (Microbial Opsin Generic Residue Numbering).

This script:
1. Initializes protos with the project data directory
2. Loads the current ST5 opsin property workbook
3. Registers available structure files (CIF) with protos
4. Creates 4 datasets:
   - mo_exp_A: Experimental structures pre-Sept 2021 (within Boltz training window)
   - mo_exp_B: Experimental structures post-Sept 2021 (true test set)
   - mo_pred_exp: Boltz predictions of experimental structures
   - mo_pred_novel: Boltz predictions of experimentally undetermined structures
5. Creates mapping between experimental PDB IDs and predicted structure IDs

Usage:
    python prepare_data.py [--rebuild]
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd

# =============================================================================
# Project paths
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parent
DATA_ROOT = PROJECT_ROOT / "data"
PROPERTY_DIR = PROJECT_ROOT / "property"
OUTPUT_DIR = PROJECT_ROOT / "opsin_output"
SRC_DIR = PROJECT_ROOT / "src"

# Structure directories (CIF files are spread across subdirectories)
# NOTE: clustered_mo is excluded - not part of this analysis
STRUCTURE_DIRS = [
    PROJECT_ROOT / "structures" / "mo_pred",       # Predicted structures (~121)
    PROJECT_ROOT / "structures" / "hideaki_exp",   # Hideaki experimental (8)
    PROJECT_ROOT / "structures" / "hideaki_pred",  # Hideaki predicted (8)
    DATA_ROOT / "structure" / "mmcif",             # Downloaded from RCSB (experimental)
]

# Add protos to path
PROTOS_SRC = PROJECT_ROOT / "protos" / "src"
if PROTOS_SRC.exists():
    sys.path.insert(0, str(PROTOS_SRC))

# =============================================================================
# Initialize protos BEFORE importing processors
# =============================================================================

import protos
protos.set_data_path(str(DATA_ROOT))

from protos.processing.structure import StructureProcessor
from protos.io.ingest.structure_loader import StructureLoader
from src.tandem_structure_preprocessing import (
    apply_structure_replacements,
    expand_structure_mapping,
    preprocess_registered_tandem_structures,
    update_tandem_helix_definitions,
)
from src.property_data import load_st5_property_data

# =============================================================================
# Configuration
# =============================================================================

PROPERTY_FILE = PROJECT_ROOT / "mo_exp_ST5_HEK1.xlsx"
LEGACY_PROPERTY_FILE = PROPERTY_DIR / "mo_exp_ST1.csv"
TANDEM_STRUCTURE_CONFIG = SRC_DIR / "resources" / "tandem_structure_domains.json"
TANDEM_STRUCTURE_MANIFEST = OUTPUT_DIR / "tandem_structure_preprocessing.json"

# Dataset definitions based on the split logic:
#
# EXPERIMENTAL STRUCTURES (from mo_exp.csv rows with experimentally_determined=1):
#   - mo_exp_A: Set A = pre-Sept 2021 (within Boltz training window)
#   - mo_exp_B: Set B = post-Sept 2021 (true test set) + Hideaki experimental
#
# PREDICTED STRUCTURES (Boltz-1 predictions):
#   - mo_pred_exp: Predictions of experimental structures (63 from mo_exp.csv + 8 Hideaki = 71)
#   - mo_pred_novel: Predictions of novel opsins (134 - 71 = 63, minus 5 missing files = 58)
#
# Note: Hideaki structures are in mo_exp.csv with short_name like 'A1ACR1' but marked
# as experimentally_determined=0 (incorrectly). They map to hideaki_pred files with
# naming like 'A1ACR1_J318_refine3_model_0'.

# Hideaki short_names in mo_exp.csv - these map to hideaki_exp/hideaki_pred (different naming)
HIDEAKI_SHORT_NAMES = {'A1ACR1', 'ChroME2s', 'CnChR2', 'CoChR', 'KnChR', 'R2ACR', 'TsChR', 'bReaChES'}


def hideaki_short_name(structure_id: str) -> str:
    """Return the ST5 short name encoded in a collaborator structure ID."""

    return str(structure_id).split("_J", 1)[0]


def is_hideaki_entry(row) -> bool:
    """Check if a row is a Hideaki entry based on short_name."""
    short_name = str(row.get("short_name", "")).strip()
    return short_name in HIDEAKI_SHORT_NAMES


def is_experimental_entry(row) -> bool:
    """Check if entry is experimental (has experimentally_determined=1 OR is Hideaki)."""
    exp_determined = row.get("experimentally_determined", None)
    return (exp_determined == 1) or is_hideaki_entry(row)


DATASET_CONFIGS = {
    "mo_exp_A": {
        "description": "Experimental microbial opsin structures - Set A (pre-Sept 2021, within Boltz training)",
        "filter": lambda row: (
            row.get("experimentally_determined", 0) == 1 and
            row.get("dataset_split", "") == "A"
        ),
        "use_pdb_id": True,  # Use lowercase PDB ID as structure_id
        "include_hideaki": False,
    },
    "mo_exp_B": {
        "description": "Experimental microbial opsin structures - Set B (post-Sept 2021, true test set) + Hideaki exp",
        "filter": lambda row: (
            row.get("experimentally_determined", 0) == 1 and
            row.get("dataset_split", "") == "B"
        ),
        "use_pdb_id": True,
        "include_hideaki": True,  # Hideaki exp structures are part of Set B
        "hideaki_source": "hideaki_exp",
    },
    "mo_pred_exp": {
        "description": "Boltz-1 predictions of experimental microbial opsin structures (63 + 8 Hideaki = 71)",
        # Filter: experimentally_determined=1 OR is Hideaki entry
        "filter": is_experimental_entry,
        "use_pdb_id": False,  # Use short_name + _model_0
        "include_hideaki": True,  # Include hideaki predictions (different naming)
        "hideaki_source": "hideaki_pred",
        "exclude_hideaki_from_filter": True,  # Don't add short_name_model_0 for Hideaki
    },
    "mo_pred_novel": {
        "description": "Boltz-1 predictions of novel microbial opsins (no experimental counterpart, ~58)",
        # Filter: NOT experimental (excludes experimentally_determined=1 AND Hideaki entries)
        "filter": lambda row: not is_experimental_entry(row),
        "use_pdb_id": False,  # Use short_name + _model_0
        "include_hideaki": False,  # Hideaki predictions are in mo_pred_exp
    },
}


# =============================================================================
# Property Data Loading
# =============================================================================

def load_property_data(property_file: Path = PROPERTY_FILE) -> pd.DataFrame:
    """Load ST5 properties and add the operational dataset split metadata."""
    return load_st5_property_data(property_file, LEGACY_PROPERTY_FILE)


def create_structure_mapping(
    property_df: pd.DataFrame,
    available_structures: Dict[str, Path]
) -> Dict[str, str]:
    """Create mapping between experimental PDB IDs and predicted structure IDs.

    Maps: pdb_id (lowercase) -> short_name_model_0

    Args:
        property_df: Property dataframe with exp_structure_id and pred_structure_id columns
        available_structures: Dict of available structure IDs

    Returns:
        Dict mapping pdb_id -> pred_structure_id (for all experimental entries with predictions)
    """
    available_set = set(available_structures.keys())
    mapping = {}

    # Create mappings for experimental structures from property file
    # Map: pdb_id -> short_name_model_0
    exp_df = property_df[property_df["experimentally_determined"] == 1]

    for _, row in exp_df.iterrows():
        exp_id = row["exp_structure_id"]  # lowercase pdb_id
        pred_id = row["pred_structure_id"]  # short_name_model_0

        # Only require that the prediction file exists
        if exp_id and pred_id and pred_id in available_set:
            mapping[exp_id] = pred_id

    print(f"[INFO] Created {len(mapping)} PDB ID -> prediction mappings from mo_exp.csv")

    # Also add Hideaki structure mappings
    # These use different naming: A1ACR1_J318_refine3 -> A1ACR1_J318_refine3_model_0
    hideaki_exp_dir = PROJECT_ROOT / "structures" / "hideaki_exp"
    hideaki_pred_dir = PROJECT_ROOT / "structures" / "hideaki_pred"

    hideaki_count = 0
    property_by_short = property_df.set_index("short_name", drop=False)
    if hideaki_exp_dir.exists() and hideaki_pred_dir.exists():
        for exp_file in hideaki_exp_dir.glob("*.cif"):
            exp_id = exp_file.stem
            pred_id = exp_id + "_model_0"
            short_name = hideaki_short_name(exp_id)
            property_row = (
                property_by_short.loc[short_name]
                if short_name in property_by_short.index
                else None
            )
            current_exp_id = (
                str(property_row["exp_structure_id"]).strip()
                if property_row is not None
                else ""
            )
            if current_exp_id and current_exp_id in available_set:
                if pred_id in available_set:
                    mapping[current_exp_id] = pred_id
                continue
            if pred_id in available_set:
                mapping[exp_id] = pred_id
                hideaki_count += 1

    print(f"[INFO] Added {hideaki_count} Hideaki structure mappings")
    print(f"[INFO] Total structure mapping pairs: {len(mapping)}")
    return mapping


# =============================================================================
# Structure Download and Discovery
# =============================================================================

def download_experimental_structures(
    loader: StructureLoader,
    property_df: pd.DataFrame,
    available_structures: Dict[str, Path],
) -> Dict[str, bool]:
    """Download missing experimental structures from RCSB.

    Args:
        loader: StructureLoader instance
        property_df: Property dataframe with exp_structure_id column
        available_structures: Currently available structures

    Returns:
        Dict mapping pdb_id to success status
    """
    available_set = set(available_structures.keys())

    # Get all experimental PDB IDs that need to be downloaded
    exp_df = property_df[property_df["experimentally_determined"] == 1]
    pdb_ids_needed = []

    for _, row in exp_df.iterrows():
        pdb_id = row["exp_structure_id"]
        if pdb_id and pdb_id not in available_set:
            pdb_ids_needed.append(pdb_id)

    if not pdb_ids_needed:
        print("[INFO] All experimental structures already available")
        return {}

    print(f"\n[INFO] Downloading {len(pdb_ids_needed)} experimental structures from RCSB...")

    results = {}
    for pdb_id in pdb_ids_needed:
        try:
            # Download from RCSB and register
            loader.download_and_register(
                pdb_id,
                source="rcsb",
                metadata={"source": "rcsb", "type": "experimental"},
            )
            results[pdb_id] = True
            print(f"  [OK] Downloaded: {pdb_id}")
        except Exception as e:
            results[pdb_id] = False
            print(f"  [WARN] Failed to download {pdb_id}: {e}")

    success_count = sum(1 for v in results.values() if v)
    print(f"[INFO] Downloaded {success_count}/{len(pdb_ids_needed)} structures from RCSB")

    return results


def scan_available_structures(structure_dirs: List[Path] = None) -> Dict[str, Path]:
    """Scan all structure directories for available CIF files.

    Returns:
        Dict mapping structure_id to file path
    """
    if structure_dirs is None:
        structure_dirs = STRUCTURE_DIRS

    print(f"[INFO] Scanning for structures in {len(structure_dirs)} directories...")

    structures = {}
    for struct_dir in structure_dirs:
        if not struct_dir.exists():
            print(f"[WARN] Directory does not exist: {struct_dir}")
            continue

        cif_files = list(struct_dir.glob("*.cif"))
        for cif_file in cif_files:
            structures[cif_file.stem] = cif_file

        print(f"[INFO]   {struct_dir.name}: {len(cif_files)} CIF files")

    print(f"[INFO] Total: {len(structures)} CIF files found")
    return structures


def register_structures(
    processor: StructureProcessor,
    loader: StructureLoader,
    structure_paths: Dict[str, Path],
    force: bool = False
) -> Dict[str, bool]:
    """Register all structures with protos."""
    print(f"\n[INFO] Registering {len(structure_paths)} structures...")

    results = {}
    registered = 0
    skipped = 0

    for struct_id, cif_path in structure_paths.items():
        # Check if already registered
        if not force and processor.entity_registry.find_entity(struct_id, "structure"):
            results[struct_id] = True
            skipped += 1
            continue

        try:
            loader.download_and_register(
                str(cif_path),
                name=struct_id,
                source="local",
                metadata={
                    "source": "mogrn",
                    "file": str(cif_path),
                    "source_dir": cif_path.parent.name,
                },
            )
            results[struct_id] = True
            registered += 1
        except Exception as e:
            print(f"[WARN] Failed to register {struct_id}: {e}")
            results[struct_id] = False

    print(f"[INFO] Registered {registered} new, skipped {skipped} existing")

    return results


# =============================================================================
# Dataset Creation
# =============================================================================

def create_datasets(
    processor: StructureProcessor,
    property_df: pd.DataFrame,
    available_structures: Dict[str, Path],
    force: bool = False,
    structure_replacements: Dict[str, Tuple[str, ...]] = None,
) -> Dict[str, int]:
    """Create the 4 datasets based on property filters."""
    print("\n" + "=" * 60)
    print("CREATING DATASETS")
    print("=" * 60)

    available_set = set(available_structures.keys())
    structure_replacements = structure_replacements or {}
    results = {}
    created_datasets = {}  # Store content of created datasets for complement logic

    for dataset_name, config in DATASET_CONFIGS.items():
        print(f"\n[INFO] Processing dataset: {dataset_name}")

        # Build content list
        content = []
        filter_func = config.get("filter")

        if filter_func is not None:
            use_pdb_id = config.get("use_pdb_id", False)
            exclude_hideaki = config.get("exclude_hideaki_from_filter", False)

            for _, row in property_df.iterrows():
                row_dict = row.to_dict()
                if not filter_func(row_dict):
                    continue

                # Skip Hideaki entries when building content from filter
                # (they'll be added via include_hideaki with correct naming)
                if exclude_hideaki and is_hideaki_entry(row_dict):
                    continue

                # Select appropriate structure ID based on dataset type
                if use_pdb_id:
                    struct_id = row["exp_structure_id"]
                else:
                    struct_id = row["pred_structure_id"]

                if struct_id and struct_id in available_set:
                    content.append(struct_id)

        # Add Hideaki structures if configured
        if config.get("include_hideaki", False):
            hideaki_source = config.get("hideaki_source", "")
            hideaki_dir = PROJECT_ROOT / "structures" / hideaki_source
            if hideaki_dir.exists():
                hideaki_ids = [f.stem for f in hideaki_dir.glob("*.cif")]
                property_by_short = property_df.set_index("short_name", drop=False)
                for struct_id in hideaki_ids:
                    short_name = hideaki_short_name(
                        struct_id.removesuffix("_model_0")
                    )
                    if (
                        hideaki_source == "hideaki_exp"
                        and short_name in property_by_short.index
                    ):
                        current_exp_id = str(
                            property_by_short.loc[short_name, "exp_structure_id"]
                        ).strip()
                        if current_exp_id and current_exp_id in available_set:
                            continue
                    if struct_id in available_set and struct_id not in content:
                        content.append(struct_id)
                print(f"[INFO]   Added {len(hideaki_ids)} Hideaki structures from {hideaki_source}")

        content = apply_structure_replacements(content, structure_replacements)
        content = [structure_id for structure_id in content if structure_id in available_set]

        if not content:
            print(f"[WARN] No structures found for dataset: {dataset_name}")
            results[dataset_name] = 0
            continue

        # Reconcile existing datasets when preprocessing changes entity IDs.
        if processor.dataset_manager.dataset_exists(dataset_name):
            existing = processor.dataset_manager.get_dataset_entities(dataset_name)
            if not force and list(existing) == content:
                print(f"[INFO] Dataset is current with {len(existing)} entities")
                results[dataset_name] = len(existing)
                created_datasets[dataset_name] = existing
                continue
            processor.dataset_manager.delete_dataset(dataset_name)
            print(f"[INFO] Replacing stale dataset: {dataset_name}")

        # Create dataset
        try:
            processor.create_dataset(
                dataset_name,
                content,
                {"description": config["description"], "source": "mogrn"},
            )
            print(f"[OK] Created dataset: {dataset_name} ({len(content)} structures)")
            results[dataset_name] = len(content)
            created_datasets[dataset_name] = content
        except Exception as e:
            print(f"[ERROR] Failed to create dataset {dataset_name}: {e}")
            results[dataset_name] = 0

    return results


# =============================================================================
# Main
# =============================================================================

def main(
    rebuild: bool = False,
    download: bool = False,
    preprocess_tandem: bool = True,
):
    """Main entry point for data preparation.

    Args:
        rebuild: Force rebuild of existing datasets
        download: Download missing experimental structures from RCSB
    """
    print("\n" + "=" * 60)
    print("MOGRN DATA PREPARATION")
    print("=" * 60)
    print(f"Project root: {PROJECT_ROOT}")
    print(f"Data root: {DATA_ROOT}")
    print(f"Protos data path: {protos.get_data_path()}")
    print("=" * 60 + "\n")

    # Create processor and loader
    print("[INFO] Initializing StructureProcessor...")
    processor = StructureProcessor("mogrn")
    loader = StructureLoader(processor=processor)

    # Load property data
    property_df = load_property_data()

    # Scan available structures (returns dict: structure_id -> file_path)
    available_structures = scan_available_structures()

    # Download missing experimental structures from RCSB if requested
    if download:
        download_results = download_experimental_structures(loader, property_df, available_structures)
        # Re-scan to include newly downloaded structures
        if any(download_results.values()):
            # Also scan the protos mmcif directory for downloaded files
            protos_mmcif = DATA_ROOT / "structure" / "mmcif"
            if protos_mmcif.exists():
                for cif_file in protos_mmcif.glob("*.cif"):
                    if cif_file.stem not in available_structures:
                        available_structures[cif_file.stem] = cif_file
            print(f"[INFO] Updated available structures: {len(available_structures)}")

    # Register structures
    register_structures(processor, loader, available_structures, force=rebuild)

    tandem_result = None
    structure_replacements = {}
    entity_parents = {}
    if preprocess_tandem:
        print("\n[INFO] Preprocessing configured tandem-domain structures...")
        tandem_result = preprocess_registered_tandem_structures(
            processor,
            list(available_structures),
            TANDEM_STRUCTURE_CONFIG,
        )
        structure_replacements = tandem_result.replacements
        entity_parents = tandem_result.entity_parents
        for child_id in entity_parents:
            available_structures[child_id] = TANDEM_STRUCTURE_CONFIG
        OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
        with TANDEM_STRUCTURE_MANIFEST.open("w") as handle:
            json.dump(tandem_result.as_manifest(), handle, indent=2)
        update_tandem_helix_definitions(
            PROPERTY_DIR / "helices_grn.json",
            tandem_result.records,
        )
        print(
            f"[INFO] Registered {len(entity_parents)} virtual tandem domains; "
            f"manifest: {TANDEM_STRUCTURE_MANIFEST}"
        )

    # Create datasets
    dataset_counts = create_datasets(
        processor,
        property_df,
        available_structures,
        force=rebuild,
        structure_replacements=structure_replacements,
    )

    # Create and save structure mapping (exp_id -> pred_id)
    structure_mapping = create_structure_mapping(property_df, available_structures)
    structure_mapping = expand_structure_mapping(
        structure_mapping, structure_replacements
    )

    # Save mapping to JSON
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    mapping_file = OUTPUT_DIR / "structure_mapping.json"
    with open(mapping_file, "w") as f:
        json.dump(structure_mapping, f, indent=2)
    print(f"\n[INFO] Saved structure mapping to: {mapping_file}")

    # Validate datasets and print detailed summary
    print("\n" + "=" * 60)
    print("DATASET VALIDATION")
    print("=" * 60)

    # Build property lookup by various IDs
    property_by_pdb = {}  # pdb_id -> row
    property_by_pred = {}  # pred_structure_id -> row
    property_by_short = {}  # short_name -> row

    for _, row in property_df.iterrows():
        exp_id = row.get("exp_structure_id", "")
        pred_id = row.get("pred_structure_id", "")
        short_name = str(row.get("short_name", "")).strip()

        if exp_id:
            property_by_pdb[exp_id] = row
        if pred_id:
            property_by_pred[pred_id] = row
        if short_name:
            property_by_short[short_name] = row

    # Hideaki mapping: hideaki file stem -> short_name in property file
    hideaki_exp_dir = PROJECT_ROOT / "structures" / "hideaki_exp"
    hideaki_to_short = {}
    if hideaki_exp_dir.exists():
        for f in hideaki_exp_dir.glob("*.cif"):
            # e.g., A1ACR1_J318_refine3 -> A1ACR1
            stem = f.stem
            prefix = stem.split("_")[0]
            if prefix in HIDEAKI_SHORT_NAMES:
                hideaki_to_short[stem] = prefix
                hideaki_to_short[stem + "_model_0"] = prefix

    dataset_details = {}

    for dataset_name in ["mo_exp_A", "mo_exp_B", "mo_pred_exp", "mo_pred_novel"]:
        if not processor.dataset_manager.dataset_exists(dataset_name):
            print(f"\n[WARN] Dataset {dataset_name} does not exist")
            continue

        entities = processor.dataset_manager.get_dataset_entities(dataset_name)
        n_entities = len(entities)

        # Check prediction mapping (for experimental datasets)
        mapped_to_pred = 0
        mapped_to_props = 0

        for entity_id in entities:
            # Check prediction mapping
            if entity_id in structure_mapping:
                mapped_to_pred += 1
            elif "_model_0" in entity_id:
                # Predictions don't need prediction mapping
                mapped_to_pred += 1

            # Check property mapping
            has_props = False
            property_entity_id = entity_parents.get(entity_id, entity_id)
            if property_entity_id in property_by_pdb:
                has_props = True
            elif property_entity_id in property_by_pred:
                has_props = True
            elif entity_id in hideaki_to_short:
                short = hideaki_to_short[entity_id]
                if short in property_by_short:
                    has_props = True
            # For hideaki predictions, strip _model_0 and check
            elif entity_id.endswith("_model_0"):
                base = entity_id[:-8]  # Remove _model_0
                if base in hideaki_to_short:
                    short = hideaki_to_short[base]
                    if short in property_by_short:
                        has_props = True

            if has_props:
                mapped_to_props += 1

        # Determine dataset label
        if dataset_name == "mo_exp_A":
            label = "Dataset 1a (exp Set A)"
        elif dataset_name == "mo_exp_B":
            label = "Dataset 1b (exp Set B + Hideaki)"
        elif dataset_name == "mo_pred_exp":
            label = "Dataset 2a (pred of exp)"
        else:
            label = "Dataset 2b (pred novel)"

        print(f"\n{label}: {n_entities} entities - {mapped_to_pred}/{n_entities} mapped to predictions - {mapped_to_props}/{n_entities} mapped to properties")

        # Report unmapped entities
        if mapped_to_props < n_entities:
            unmapped = []
            for entity_id in entities:
                has_props = False
                property_entity_id = entity_parents.get(entity_id, entity_id)
                if property_entity_id in property_by_pdb:
                    has_props = True
                elif property_entity_id in property_by_pred:
                    has_props = True
                elif entity_id in hideaki_to_short:
                    short = hideaki_to_short[entity_id]
                    if short in property_by_short:
                        has_props = True
                elif entity_id.endswith("_model_0"):
                    base = entity_id[:-8]
                    if base in hideaki_to_short:
                        short = hideaki_to_short[base]
                        if short in property_by_short:
                            has_props = True
                if not has_props:
                    unmapped.append(entity_id)
            if unmapped:
                print(f"  [WARN] Unmapped to properties: {unmapped[:10]}{'...' if len(unmapped) > 10 else ''}")

        dataset_details[dataset_name] = {
            "entities": n_entities,
            "mapped_to_pred": mapped_to_pred,
            "mapped_to_props": mapped_to_props,
        }

    # Final summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Property file: {PROPERTY_FILE.name}")
    print(f"Total structures available: {len(available_structures)}")
    print(f"Property entries: {len(property_df)}")
    print(f"Structure mapping pairs: {len(structure_mapping)}")
    print("=" * 60 + "\n")

    return {
        "processor": processor,
        "property_df": property_df,
        "available_structures": available_structures,
        "dataset_counts": dataset_counts,
        "structure_mapping": structure_mapping,
        "tandem_structure_preprocessing": (
            tandem_result.as_manifest() if tandem_result is not None else None
        ),
        "dataset_details": dataset_details,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Prepare MOGRN data and register datasets")
    parser.add_argument("--rebuild", action="store_true", help="Force rebuild of existing datasets")
    parser.add_argument("--download", action="store_true", help="Download missing experimental structures from RCSB")
    parser.add_argument(
        "--skip-tandem-preprocessing",
        action="store_true",
        help="Do not replace configured long-chain tandem structures",
    )
    args = parser.parse_args()

    main(
        rebuild=args.rebuild,
        download=args.download,
        preprocess_tandem=not args.skip_tandem_preprocessing,
    )
