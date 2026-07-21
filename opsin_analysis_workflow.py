#!/usr/bin/env python3
"""
Opsin analysis workflow for MOGRN (Microbial Opsin Generic Residue Numbering).

This workflow:
1. Uses prepare_data.py to register structures and create datasets
2. Loads structures from registered datasets using protos API
3. Runs analysis pipeline:
   - Calculate structure errors (experimental vs predicted)
   - Helix annotation and alignment
   - Structure comparison (RMSD matrix)
   - GRN assignment

Usage:
    python opsin_analysis_workflow.py [options]
"""

import argparse
import json
import os
import pickle
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np
import pandas as pd
from tqdm import tqdm

# =============================================================================
# Project paths
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parent
DATA_ROOT = PROJECT_ROOT / "data"
OUTPUT_DIR = PROJECT_ROOT / "opsin_output"
PROPERTY_DIR = PROJECT_ROOT / "property"

# Add protos to path
PROTOS_SRC = PROJECT_ROOT / "protos" / "src"
if PROTOS_SRC.exists():
    sys.path.insert(0, str(PROTOS_SRC))

# =============================================================================
# Initialize protos
# =============================================================================

import protos
protos.set_data_path(str(DATA_ROOT))

from protos.processing.structure import StructureProcessor

# =============================================================================
# Import analysis modules
# =============================================================================

from src.error_analysis import calculate_structure_errors
from src.structure_comparison import (
    compare_structures,
    create_unified_structure_mapping,
)
from src.helix_analysis import load_and_apply_helix_annotations
from src.assign_grns import align_and_assign_grn
from src.lyr_processing import process_lyr_in_processor_data, standardize_retinal_naming
from src.property_data import load_st5_property_data
from src.tandem_structure_preprocessing import dataset_signature


# =============================================================================
# Compatibility class for error_analysis (must be at module level for pickling)
# =============================================================================

class DatasetCompat:
    """Compatibility class to provide pdb_ids and data for legacy error_analysis code."""
    def __init__(self, pdb_ids, data=None):
        self.pdb_ids = list(pdb_ids)
        self.data = data if data is not None else pd.DataFrame()

    def format_data_types(self):
        """Format data types for compatibility with legacy code."""
        if self.data.empty:
            return
        # Ensure coordinate columns are numeric
        for col in ['x', 'y', 'z']:
            if col in self.data.columns:
                self.data[col] = pd.to_numeric(self.data[col], errors='coerce')


# =============================================================================
# Structure Loading (using protos API)
# =============================================================================

def load_structure_mapping() -> Dict[str, str]:
    """Load the structure mapping from prepare_data.py output."""
    mapping_file = OUTPUT_DIR / "structure_mapping.json"
    if mapping_file.exists():
        with open(mapping_file) as f:
            return json.load(f)
    return {}


def load_tandem_structure_manifest() -> Dict[str, Any]:
    """Load the external preprocessing manifest written by prepare_data.py."""
    manifest_file = OUTPUT_DIR / "tandem_structure_preprocessing.json"
    if not manifest_file.exists():
        return {}
    with manifest_file.open() as handle:
        return json.load(handle)


def load_property_data(property_file: Path = None) -> pd.DataFrame:
    """Load the authoritative ST5 workbook with workflow split metadata."""

    property_file = property_file or PROJECT_ROOT / "mo_exp_ST5_HEK1.xlsx"
    return load_st5_property_data(
        property_file,
        PROPERTY_DIR / "mo_exp_ST1.csv",
    )


def load_structures_from_datasets(
    datasets: list,
    chain_id: str = "A",
    retinal_name: str = "RET",
    retinal_cutoff: float = 6.0,
    use_cache: bool = True,
    cache_dir: Path = None,
) -> Dict[str, Any]:
    """
    Load structures from protos datasets using the new StructureProcessor API.

    Args:
        datasets: List of dataset names to load
        chain_id: Chain ID to filter by
        retinal_name: Name of retinal residue
        retinal_cutoff: Distance cutoff for retinal selection
        use_cache: Whether to use cached data
        cache_dir: Directory for cache files

    Returns:
        Dictionary containing processed structures and processor
    """
    print("\n" + "=" * 60)
    print("LOADING STRUCTURES FROM DATASETS")
    print("=" * 60)

    if cache_dir is None:
        cache_dir = OUTPUT_DIR / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    # Cache file for processed structures
    cache_file = cache_dir / f"processed_structures_{chain_id}.pkl"

    # Resolve current dataset contents before consulting the cache.  A changed
    # preprocessing split must invalidate every downstream stage automatically.
    processor = StructureProcessor("mogrn_workflow")
    dataset_contents = {}
    for dataset_name in datasets:
        if not processor.dataset_manager.dataset_exists(dataset_name):
            continue
        dataset_contents[dataset_name] = processor.get_dataset_entities(dataset_name)
    input_signature = dataset_signature(
        dataset_contents,
        chain_id=chain_id,
        retinal_name=retinal_name,
        preprocessing_signature=load_tandem_structure_manifest().get(
            "config_sha256"
        ),
    )

    # Try loading from cache
    if use_cache and cache_file.exists():
        print(f"[INFO] Loading from cache: {cache_file}")
        try:
            with open(cache_file, "rb") as f:
                cached_data = pickle.load(f)
            if (
                "processed_structures" in cached_data
                and cached_data.get("_dataset_signature") == input_signature
            ):
                print(f"[INFO] Loaded {len(cached_data['processed_structures'])} structures from cache")
                return cached_data
            print("[INFO] Ignoring stale processed-structure cache")
        except Exception as e:
            print(f"[WARN] Failed to load cache: {e}")

    all_structures = {}

    for dataset_name in datasets:
        print(f"\n[INFO] Loading dataset: {dataset_name}")

        if not processor.dataset_manager.dataset_exists(dataset_name):
            print(f"[WARN] Dataset not found: {dataset_name}")
            continue

        # Get structure IDs from the dataset snapshot used for the cache key.
        structure_ids = dataset_contents[dataset_name]
        print(f"[INFO] Found {len(structure_ids)} structures in {dataset_name}")

        # Load structures using new API
        loaded_count = 0
        lyr_count = 0
        lig_count = 0
        for struct_id in structure_ids:
            try:
                df = processor.load_entity(struct_id)
                if df is not None:
                    # Reset index to get structure_id as column
                    df_reset = df.reset_index()
                    # Add res_atom_name alias for backward compatibility
                    if "atom_name" in df_reset.columns and "res_atom_name" not in df_reset.columns:
                        df_reset["res_atom_name"] = df_reset["atom_name"]
                    # Add pdb_id alias for backward compatibility
                    if "structure_id" in df_reset.columns and "pdb_id" not in df_reset.columns:
                        df_reset["pdb_id"] = df_reset["structure_id"]

                    # CRITICAL: Standardize retinal naming EARLY (before any caching)
                    # This converts LYR → LYS + RET and LIG → RET
                    had_lyr = (df_reset['res_name3l'] == 'LYR').any() if 'res_name3l' in df_reset.columns else False
                    had_lig = (df_reset['res_name3l'] == 'LIG').any() if 'res_name3l' in df_reset.columns else False
                    df_reset = standardize_retinal_naming(df_reset, retinal_res_name=retinal_name, verbose=False)
                    if had_lyr:
                        lyr_count += 1
                    if had_lig:
                        lig_count += 1

                    all_structures[struct_id] = {
                        "dataset": dataset_name,
                        "df": df_reset,
                    }
                    loaded_count += 1
            except Exception as e:
                print(f"[WARN] Failed to load {struct_id}: {e}")

        print(f"[INFO] Loaded {loaded_count}/{len(structure_ids)} structures for {dataset_name}")
        if lyr_count > 0 or lig_count > 0:
            print(f"[INFO] Retinal standardized: {lyr_count} LYR→LYS+RET, {lig_count} LIG→RET")

    if not all_structures:
        print("[ERROR] No structures loaded!")
        return {"processed_structures": {}, "processor": processor, "datasets": dataset_contents}

    # Process structures: filter by chain and retinal
    print(f"\n[INFO] Processing {len(all_structures)} structures...")

    processed_structures = {}

    for struct_id, struct_data in tqdm(all_structures.items(), desc="Processing"):
        df = struct_data["df"]

        if df.empty:
            continue

        # For predicted structures (_model_0), set all chains to A
        # (Boltz puts protein on chain A but LIG on chain B)
        if "_model_0" in struct_id:
            df = df.copy()
            df["auth_chain_id"] = chain_id

        # NOTE: LYR→LYS+RET and LIG→RET conversion is now done earlier
        # in the loading phase (standardize_retinal_naming)

        # Filter by chain
        df_chain = df[df["auth_chain_id"] == chain_id].copy()

        if df_chain.empty:
            continue

        # Find retinal (now consistently named RET after standardization)
        df_ret = df_chain[df_chain["res_name3l"] == retinal_name].copy()

        # Filter to keep only ATOM records and retinal HETATM
        # Note: LYR is now converted to LYS (ATOM) + RET (HETATM)
        is_atom = df_chain["group"] == "ATOM"
        is_ret = (df_chain["group"] == "HETATM") & (df_chain["res_name3l"] == retinal_name)
        df_filtered = df_chain[is_atom | is_ret].copy()

        # Create normalized dataframe
        df_norm = df_filtered.copy()
        for coord in ["x", "y", "z"]:
            if coord in df_norm.columns:
                df_norm[coord] = pd.to_numeric(df_norm[coord], errors="coerce")

        # Extract CA atoms - check for both res_atom_name and atom_name columns
        atom_name_col = "res_atom_name" if "res_atom_name" in df_norm.columns else "atom_name"
        df_ca_norm = df_norm[df_norm[atom_name_col] == "CA"].copy()

        # Determine structure type
        is_predicted = "_model_0" in struct_id or "_pred" in struct_id

        processed_structures[struct_id] = {
            "chain_id": chain_id,
            "dataset": struct_data["dataset"],
            "df": df_filtered,
            "df_norm": df_norm,
            "df_ca_norm": df_ca_norm,
            "df_ret": df_ret,
            "structure_type": "predicted" if is_predicted else "experimental",
        }

    print(f"[INFO] Processed {len(processed_structures)} structures")

    # Load structure mapping
    structure_mapping = load_structure_mapping()
    print(f"[INFO] Loaded {len(structure_mapping)} exp-pred structure pairs")

    # Create compatibility objects for error_analysis (expects cp_mo_exp, etc.)
    # Map dataset names to legacy processor names
    # mo_exp_A/B -> cp_mo_exp, mo_pred_exp/novel -> cp_mo_pred
    # hide_exp -> cp_hide_exp, hide_pred -> cp_hide_pred
    mo_exp_ids = set()
    mo_pred_ids = set()
    hide_exp_ids = set()
    hide_pred_ids = set()

    for ds_name, ids in dataset_contents.items():
        if ds_name in ["mo_exp_A", "mo_exp_B"]:
            mo_exp_ids.update(ids)
        elif ds_name in ["mo_pred_exp", "mo_pred_novel"]:
            mo_pred_ids.update(ids)
        elif ds_name == "hide_exp":
            hide_exp_ids.update(ids)
        elif ds_name == "hide_pred":
            hide_pred_ids.update(ids)

    # Build combined DataFrames for each group
    def build_combined_df(pdb_ids):
        """Combine structure DataFrames into one with pdb_id column."""
        dfs = []
        for pdb_id in pdb_ids:
            if pdb_id in processed_structures and 'df' in processed_structures[pdb_id]:
                df = processed_structures[pdb_id]['df'].copy()
                df['pdb_id'] = pdb_id
                dfs.append(df)
        return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()

    mo_exp_df = build_combined_df(mo_exp_ids)
    mo_pred_df = build_combined_df(mo_pred_ids)
    hide_exp_df = build_combined_df(hide_exp_ids)
    hide_pred_df = build_combined_df(hide_pred_ids)

    cp_mo_exp = DatasetCompat(mo_exp_ids, mo_exp_df)
    cp_mo_pred = DatasetCompat(mo_pred_ids, mo_pred_df)
    cp_hide_exp = DatasetCompat(hide_exp_ids, hide_exp_df)
    cp_hide_pred = DatasetCompat(hide_pred_ids, hide_pred_df)

    print(f"[INFO] Dataset compatibility: mo_exp={len(mo_exp_ids)}, mo_pred={len(mo_pred_ids)}, hide_exp={len(hide_exp_ids)}, hide_pred={len(hide_pred_ids)}")

    result = {
        "processed_structures": processed_structures,
        "structure_mapping": structure_mapping,
        "processor": processor,
        "datasets": dataset_contents,
        # Legacy compatibility for error_analysis
        "cp_mo_exp": cp_mo_exp,
        "cp_mo_pred": cp_mo_pred,
        "cp_hide_exp": cp_hide_exp,
        "cp_hide_pred": cp_hide_pred,
        "_dataset_signature": input_signature,
    }

    # Save to cache
    if use_cache:
        try:
            print(f"[INFO] Saving to cache: {cache_file}")
            with open(cache_file, "wb") as f:
                pickle.dump(result, f)
        except Exception as e:
            print(f"[WARN] Failed to save cache: {e}")

    return result


# =============================================================================
# Main Workflow
# =============================================================================

def run_opsin_analysis_workflow(
    output_dir: Path = None,
    visualize: bool = False,
    use_cache: bool = True,
    chain_id: str = "A",
    retinal_name: str = "RET",
    retinal_cutoff: float = 6.0,
    global_ref_override: str = None,
    helices_file: str = "property/helices_grn.json",
    skip_prepare: bool = False,
    datasets: list = None,
) -> Dict[str, Any]:
    """
    Run the full opsin analysis workflow.

    Args:
        output_dir: Directory for output files
        visualize: Whether to generate visualizations
        use_cache: Whether to use cached data
        chain_id: Chain ID to use for analysis
        retinal_name: Name of retinal residue
        retinal_cutoff: Distance cutoff for retinal selection
        global_ref_override: Override for global reference structure
        helices_file: Path to helix boundaries JSON
        skip_prepare: Skip running prepare_data.py
        datasets: List of dataset names to load (default: all 4)

    Returns:
        Dictionary with analysis results
    """
    print("\n" + "=" * 80)
    print("RUNNING OPSIN ANALYSIS WORKFLOW")
    print("=" * 80 + "\n")

    # Set up directories
    if output_dir is None:
        output_dir = OUTPUT_DIR
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    cache_dir = output_dir / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    print(f"Output directory: {output_dir}")
    print(f"Cache directory: {cache_dir}")
    print(f"Protos data path: {protos.get_data_path()}")

    # Step 0: Run prepare_data.py if needed
    if not skip_prepare:
        print("\n" + "-" * 60)
        print("Step 0: Preparing data and registering structures")
        print("-" * 60)

        from prepare_data import main as prepare_main
        prepare_result = prepare_main(rebuild=False, download=False)
        print(f"[OK] Data preparation complete")

    # Define default datasets
    if datasets is None:
        datasets = ["mo_exp_A", "mo_exp_B", "mo_pred_exp", "mo_pred_novel"]

    # Step 1: Load structures from datasets
    print("\n" + "-" * 60)
    print("Step 1: Loading structures from datasets")
    print("-" * 60)

    data = load_structures_from_datasets(
        datasets=datasets,
        chain_id=chain_id,
        retinal_name=retinal_name,
        retinal_cutoff=retinal_cutoff,
        use_cache=use_cache,
        cache_dir=cache_dir,
    )

    if not data.get("processed_structures"):
        print("[ERROR] No structures loaded. Exiting.")
        return data

    # Step 2: Load property data
    print("\n" + "-" * 60)
    print("Step 2: Loading property data")
    print("-" * 60)

    property_df = load_property_data()
    print(f"[INFO] Loaded {len(property_df)} property entries")
    tandem_manifest = load_tandem_structure_manifest()

    # Create property lookup for structures
    property_data = {
        "properties": {},
        "structure_mapping": data.get("structure_mapping", {}),
    }

    if not property_df.empty:
        for _, row in property_df.iterrows():
            row_dict = row.to_dict()

            # Extract key properties
            props = {
                "domain": row_dict.get("Rhodopsin Type (Microbial)", "Unknown"),
                "molecular_function": row_dict.get("molecular_function", "Unknown"),
                "experimentally_determined": row_dict.get("experimentally_determined", 0),
            }

            # Map to both PDB ID and short_name variants
            pdb_id = row_dict.get("PDB ID", row_dict.get("pdb_id", ""))
            short_name = row_dict.get("short_name", "")

            if pd.notna(pdb_id) and str(pdb_id).strip():
                pdb_id_clean = str(pdb_id).strip().lower()
                property_data["properties"][pdb_id_clean] = props

            if pd.notna(short_name) and str(short_name).strip():
                pred_id = str(short_name).strip() + "_model_0"
                property_data["properties"][pred_id] = props

    # Virtual tandem domains inherit the metadata of the parent structure.
    property_lookup = {
        str(key).lower(): key for key in property_data["properties"]
    }
    for child_id, parent_id in tandem_manifest.get("entity_parents", {}).items():
        parent_key = property_lookup.get(str(parent_id).lower())
        if parent_key is not None:
            property_data["properties"][child_id] = property_data["properties"][parent_key]

    # Update processed structures with properties
    for struct_id, struct_data in data["processed_structures"].items():
        if struct_id in property_data["properties"]:
            struct_data["properties"] = property_data["properties"][struct_id]

    # Step 3: Calculate structure errors
    print("\n" + "-" * 60)
    print("Step 3: Calculating structure errors (experimental vs predicted)")
    print("-" * 60)

    errors_cache = cache_dir / f"structure_errors_{chain_id}.pkl"
    input_signature = data.get("_dataset_signature")

    if use_cache and errors_cache.exists():
        print(f"[INFO] Loading errors from cache: {errors_cache}")
        try:
            with open(errors_cache, "rb") as f:
                errors_data = pickle.load(f)
            if errors_data.get("_dataset_signature") == input_signature:
                data.update(errors_data)
            else:
                print("[INFO] Ignoring stale structure-error cache")
                errors_data = None
        except Exception as e:
            print(f"[WARN] Failed to load cache: {e}")
            errors_data = None
    else:
        errors_data = None

    if errors_data is None:
        try:
            errors_data = calculate_structure_errors(data, output_dir=str(output_dir), visualize=visualize)
            errors_data["_dataset_signature"] = input_signature
            data.update(errors_data)

            if use_cache:
                with open(errors_cache, "wb") as f:
                    pickle.dump(errors_data, f)
                print(f"[INFO] Saved errors to cache")
        except Exception as e:
            print(f"[WARN] Error calculation failed: {e}")

    # Step 4: Load helix annotations from helices_grn.json
    # (Generated by scripts/generate_helices_grn.py from the postprocessed GRN table)
    print("\n" + "-" * 60)
    print("Step 4: Loading helix annotations from helices_grn.json")
    print("-" * 60)

    try:
        helix_data = load_and_apply_helix_annotations(data)
        data.update(helix_data)
    except Exception as e:
        print(f"[WARN] Helix annotation failed: {e}")

    # Step 5: Structure comparison (RMSD matrix)
    print("\n" + "-" * 60)
    print("Step 5: Comparing structures (RMSD matrix)")
    print("-" * 60)

    comparison_cache = cache_dir / f"structure_comparison_{chain_id}.pkl"

    if use_cache and comparison_cache.exists():
        print(f"[INFO] Loading comparison from cache: {comparison_cache}")
        try:
            with open(comparison_cache, "rb") as f:
                comparison_data = pickle.load(f)
            if comparison_data.get("_dataset_signature") == input_signature:
                data.update(comparison_data)
            else:
                print("[INFO] Ignoring stale structure-comparison cache")
                comparison_data = None
        except Exception as e:
            print(f"[WARN] Failed to load cache: {e}")
            comparison_data = None
    else:
        comparison_data = None

    if comparison_data is None:
        try:
            comparison_data = compare_structures(
                data, str(output_dir), visualize=visualize
            )
            comparison_data["_dataset_signature"] = input_signature
            data.update(comparison_data)

            if use_cache:
                with open(comparison_cache, "wb") as f:
                    pickle.dump(comparison_data, f)
                print(f"[INFO] Saved comparison to cache")
        except Exception as e:
            print(f"[WARN] Structure comparison failed: {e}")

    # Step 6: GRN assignment
    print("\n" + "-" * 60)
    print("Step 6: Assigning Generic Residue Numbers (GRNs)")
    print("-" * 60)

    grn_cache = cache_dir / f"grn_assignment_{chain_id}.pkl"

    if use_cache and grn_cache.exists():
        print(f"[INFO] Loading GRN data from cache: {grn_cache}")
        try:
            with open(grn_cache, "rb") as f:
                grn_data = pickle.load(f)
            if grn_data.get("_dataset_signature") == input_signature:
                data.update(grn_data)
            else:
                print("[INFO] Ignoring stale GRN cache")
                grn_data = None
        except Exception as e:
            print(f"[WARN] Failed to load cache: {e}")
            grn_data = None
    else:
        grn_data = None

    if grn_data is None:
        try:
            grn_data = align_and_assign_grn(
                data,
                str(output_dir),
                visualize=visualize,
                global_ref_override=global_ref_override,
                helices_file=helices_file,
            )
            grn_data["_dataset_signature"] = input_signature
            data.update(grn_data)

            if use_cache:
                with open(grn_cache, "wb") as f:
                    pickle.dump(grn_data, f)
                print(f"[INFO] Saved GRN data to cache")
        except Exception as e:
            print(f"[WARN] GRN assignment failed: {e}")

    # Save analysis summary
    print("\n" + "-" * 60)
    print("Saving analysis summary")
    print("-" * 60)

    try:
        summary = {
            "datasets": datasets,
            "structures_count": len(data.get("processed_structures", {})),
            "exp_pred_pairs": len(data.get("structure_mapping", {})),
            "chain_id": chain_id,
            "timestamp": pd.Timestamp.now().strftime("%Y-%m-%d %H:%M:%S"),
            "cache_enabled": use_cache,
            "tandem_domain_entities": sorted(
                tandem_manifest.get("entity_parents", {})
            ),
            "grn_stage": "raw_uncurated_baseline",
            "curated_application_command": "python apply_curated_grns.py",
        }

        with open(output_dir / "analysis_summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        print(f"[OK] Summary saved to {output_dir / 'analysis_summary.json'}")
    except Exception as e:
        print(f"[WARN] Could not save summary: {e}")

    print("\n" + "=" * 80)
    print(f"ANALYSIS COMPLETE. Results saved to {output_dir}")
    print("=" * 80)

    return data


# =============================================================================
# CLI Entry Point
# =============================================================================

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run opsin analysis workflow")

    parser.add_argument(
        "--output-dir", type=str, default=None,
        help="Directory for output files"
    )
    parser.add_argument(
        "--visualize", action="store_true",
        help="Generate visualizations"
    )
    parser.add_argument(
        "--no-cache", action="store_true",
        help="Disable caching"
    )
    parser.add_argument(
        "--chain-id", type=str, default="A",
        help="Chain ID to analyze (default: A)"
    )
    parser.add_argument(
        "--global-ref", type=str, default="7bmh",
        help="Global reference structure (default: 7bmh)"
    )
    parser.add_argument(
        "--helices-file", type=str, default="property/helices_grn.json",
        help="Path to helix boundaries JSON"
    )
    parser.add_argument(
        "--skip-prepare", action="store_true",
        help="Skip running prepare_data.py"
    )
    parser.add_argument(
        "--datasets", type=str, nargs="+",
        default=None,
        help="Specific datasets to load"
    )

    args = parser.parse_args()

    results = run_opsin_analysis_workflow(
        output_dir=Path(args.output_dir) if args.output_dir else None,
        visualize=args.visualize,
        use_cache=not args.no_cache,
        chain_id=args.chain_id,
        global_ref_override=args.global_ref,
        helices_file=args.helices_file,
        skip_prepare=args.skip_prepare,
        datasets=args.datasets,
    )
