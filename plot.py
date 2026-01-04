#!/usr/bin/env python3
"""
Visualization script for MOGRN (Microbial Opsin Generic Residue Numbering).

This script generates publication-quality figures from workflow outputs:
- Opsin overview plots
- RMSD heatmaps and clustering
- Distance plots (all-atom and CA)
- Helix logo plots
- Property analysis
- Interactive GRN alignment visualization

Usage:
    python plot.py [--input-dir DIR] [--output-dir DIR] [--only TYPE]
"""

import argparse
import json
import os
import pickle
import sys
import traceback
from collections import Counter
from pathlib import Path
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

# =============================================================================
# Project paths
# =============================================================================

PROJECT_ROOT = Path(__file__).resolve().parent
OUTPUT_DIR = PROJECT_ROOT / "opsin_output"
PROPERTY_DIR = PROJECT_ROOT / "property"
FIGURES_DIR = OUTPUT_DIR / "paper_figures"

# Add src to path
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

# Add protos to path
PROTOS_SRC = PROJECT_ROOT / "protos" / "src"
if PROTOS_SRC.exists() and str(PROTOS_SRC) not in sys.path:
    sys.path.insert(0, str(PROTOS_SRC))

# =============================================================================
# Imports from src modules
# =============================================================================

from src.opsin_color_scheme import OPSIN_COLORS
from src.visualization_functions import (
    compute_rmsd_metrics,
    create_opsin_overview_plot,
    visualize_rmsd_matrix_improved,
    plot_distances_with_std,
    plot_helix_logo_plots,
    plot_error_box_comparison,
)

try:
    from src.visualize_alignment_grn import (
        create_opsin_visualization_from_workflow,
        create_opsin_visualization_from_workflow_b,
    )
    HAS_INTERACTIVE = True
except ImportError:
    HAS_INTERACTIVE = False
    print("[WARN] Interactive visualization module not available")

try:
    from src.property_mapping import create_unified_property_mapper
    HAS_PROPERTY_MAPPER = True
except ImportError:
    HAS_PROPERTY_MAPPER = False
    print("[WARN] Property mapping module not available")

try:
    from protos.processing.grn.grn_utils import sort_grns_str, get_tm_residues
except ImportError:
    print("[WARN] GRN utilities not available")
    def sort_grns_str(grns):
        return sorted(grns)
    def get_tm_residues(grns):
        return [g for g in grns if "." in g]


# =============================================================================
# Data Loading Functions
# =============================================================================

def load_cached_data(cache_path: Path, description: str = "data") -> Optional[dict]:
    """Load pickled data from cache file."""
    if cache_path.exists():
        print(f"[INFO] Loading {description} from: {cache_path}")
        try:
            with open(cache_path, "rb") as f:
                return pickle.load(f)
        except Exception as e:
            print(f"[ERROR] Failed to load {description}: {e}")
    else:
        print(f"[WARN] Cache file not found: {cache_path}")
    return None


def load_structure_mapping(output_dir: Path) -> Dict[str, str]:
    """Load structure mapping from prepare_data.py output."""
    mapping_file = output_dir / "structure_mapping.json"
    if mapping_file.exists():
        with open(mapping_file) as f:
            return json.load(f)
    return {}


def load_property_data(property_file: Path = None) -> pd.DataFrame:
    """Load property data from CSV."""
    if property_file is None:
        property_file = PROPERTY_DIR / "mo_exp.csv"

    if not property_file.exists():
        print(f"[WARN] Property file not found: {property_file}")
        return pd.DataFrame()

    df = pd.read_csv(property_file)
    df.columns = df.columns.str.strip()

    # Clean function annotations
    if "molecular_function" in df.columns:
        df["molecular_function"] = df["molecular_function"].apply(
            lambda x: str(x).replace("?", "").strip() if pd.notna(x) else ""
        )

    return df


def load_workflow_data(input_dir: Path, chain_id: str = "A") -> dict:
    """Load all cached workflow data."""
    data = {}
    cache_dir = input_dir / "cache"

    cache_files = {
        "processed_structures": f"processed_structures_{chain_id}.pkl",
        "structure_errors": f"structure_errors_{chain_id}.pkl",
        "helix_annotations": f"helix_annotations_{chain_id}.pkl",
        "structure_comparison": f"structure_comparison_{chain_id}.pkl",
        "grn_assignment": f"grn_assignment_{chain_id}.pkl",
    }

    for key, filename in cache_files.items():
        cache_path = cache_dir / filename
        cached = load_cached_data(cache_path, key)
        if cached:
            if isinstance(cached, dict):
                data.update(cached)
            else:
                data[key] = cached

    # Load structure mapping
    data["structure_mapping"] = load_structure_mapping(input_dir)

    # Load GRN tables if available
    grn_tables_path = input_dir / "tree_based_grn" / "grn_tables_data.pkl"
    if not grn_tables_path.exists():
        grn_tables_path = input_dir / "grn_tables_data.pkl"

    grn_data = load_cached_data(grn_tables_path, "GRN tables")
    if grn_data:
        data.update(grn_data)

    return data


# =============================================================================
# Property Analysis Visualizations
# =============================================================================

def create_property_distribution_figure(
    experimental_props: list,
    predicted_only_props: list,
    combined_props: list,
    natural_domains: list,
    output_dir: Path,
):
    """Create property distribution heatmaps."""
    from src.opsin_color_scheme import RMSD_WHITE_TO_DARKGRAY_CMAP

    all_functions = sorted(set(p["molecular_function"] for p in combined_props))
    all_domains = sorted(natural_domains)

    fig = plt.figure(figsize=(20, 14))
    gs = fig.add_gridspec(3, 3, height_ratios=[2, 2, 1], hspace=0.3, wspace=0.25)

    # Combined heatmap
    ax1 = fig.add_subplot(gs[0, :])
    df_combined = pd.DataFrame(combined_props)
    pivot_combined = pd.crosstab(df_combined["molecular_function"], df_combined["domain"])
    pivot_combined = pivot_combined.reindex(index=all_functions, columns=all_domains, fill_value=0)

    annot_combined = pivot_combined.astype(str)
    annot_combined[pivot_combined == 0] = ""

    sns.heatmap(
        pivot_combined,
        annot=annot_combined,
        fmt="",
        cmap=RMSD_WHITE_TO_DARKGRAY_CMAP,
        linewidths=0.5,
        linecolor="gray",
        cbar_kws={"label": "Count"},
        vmin=0,
        ax=ax1,
    )
    ax1.set_xlabel("Domain", fontsize=14)
    ax1.set_ylabel("Molecular Function", fontsize=14)
    ax1.set_title("Combined Dataset (Experimental + Predicted-only)", fontsize=16)

    # Experimental heatmap
    ax2a = fig.add_subplot(gs[1, 0])
    if experimental_props:
        df_exp = pd.DataFrame(experimental_props)
        pivot_exp = pd.crosstab(df_exp["molecular_function"], df_exp["domain"])
        pivot_exp = pivot_exp.reindex(index=all_functions, columns=all_domains, fill_value=0)
        annot_exp = pivot_exp.astype(str)
        annot_exp[pivot_exp == 0] = ""
        sns.heatmap(
            pivot_exp,
            annot=annot_exp,
            fmt="",
            cmap=RMSD_WHITE_TO_DARKGRAY_CMAP,
            linewidths=0.5,
            ax=ax2a,
        )
        ax2a.set_title(f"Experimental Only (n={len(experimental_props)})")

    # Predicted-only heatmap
    ax2b = fig.add_subplot(gs[1, 1])
    if predicted_only_props:
        df_pred = pd.DataFrame(predicted_only_props)
        pivot_pred = pd.crosstab(df_pred["molecular_function"], df_pred["domain"])
        pivot_pred = pivot_pred.reindex(index=all_functions, columns=all_domains, fill_value=0)
        annot_pred = pivot_pred.astype(str)
        annot_pred[pivot_pred == 0] = ""
        sns.heatmap(
            pivot_pred,
            annot=annot_pred,
            fmt="",
            cmap="YlOrRd",
            linewidths=0.5,
            ax=ax2b,
        )
        ax2b.set_title(f"Predicted-only (n={len(predicted_only_props)})")

    # Difference heatmap
    ax2c = fig.add_subplot(gs[1, 2])
    if experimental_props and predicted_only_props:
        pivot_diff = pivot_combined - pivot_exp
        sns.heatmap(
            pivot_diff,
            annot=True,
            fmt="d",
            cmap="Reds",
            linewidths=0.5,
            ax=ax2c,
        )
        ax2c.set_title("Prediction Contribution")

    # Missing combinations
    ax3 = fig.add_subplot(gs[2, :])
    ax3.axis("off")

    all_combos = [(f, d) for f in all_functions for d in all_domains]
    observed = set((r["molecular_function"], r["domain"]) for r in combined_props)
    missing = [c for c in all_combos if c not in observed]

    text = "Missing Function-Domain Combinations:\n"
    if missing:
        by_func = {}
        for func, domain in missing:
            by_func.setdefault(func, []).append(domain)
        for func, domains in by_func.items():
            text += f"  {func}: {', '.join(domains)}\n"
        text += f"\nCoverage: {len(observed) / len(all_combos) * 100:.1f}%"
    else:
        text += "All combinations represented!"

    ax3.text(0.05, 0.95, text, transform=ax3.transAxes, fontsize=11, va="top", family="monospace")

    plt.suptitle("Microbial Opsin Property Analysis", fontsize=18)
    plt.tight_layout()

    fig_path = output_dir / "06_property_analysis.png"
    plt.savefig(fig_path, dpi=300, bbox_inches="tight")
    plt.close()
    print(f"[OK] Saved: {fig_path}")


def create_property_analysis(
    processed_structures: dict,
    property_df: pd.DataFrame,
    structure_mapping: dict,
    output_dir: Path,
):
    """Create property analysis visualizations."""
    print("\n[PROPERTY ANALYSIS] Generating property visualizations...")

    NATURAL_DOMAINS = ["Eukaryota", "Bacteria", "Archaea"]

    # Build property lookup
    props_by_id = {}
    for _, row in property_df.iterrows():
        row_dict = row.to_dict()
        props = {
            "domain": row_dict.get("Rhodopsin Type (Microbial)", "Unknown"),
            "molecular_function": row_dict.get("molecular_function", "Unknown"),
        }

        pdb_id = str(row_dict.get("PDB ID", "")).strip().lower()
        short_name = str(row_dict.get("short_name", "")).strip()

        if pdb_id:
            props_by_id[pdb_id] = props
        if short_name:
            props_by_id[short_name + "_model_0"] = props

    # Collect properties
    experimental_props = []
    predicted_only_props = []
    combined_props = []

    pred_with_exp = set(structure_mapping.values())

    for struct_id in processed_structures:
        props = props_by_id.get(struct_id, {})
        if not props or props.get("domain", "Unknown") not in NATURAL_DOMAINS:
            continue

        is_predicted = "_model_0" in struct_id or "_pred" in struct_id

        if not is_predicted:
            experimental_props.append(props)
            combined_props.append(props)
        else:
            if struct_id not in pred_with_exp:
                predicted_only_props.append(props)
                combined_props.append(props)

    print(f"  Experimental: {len(experimental_props)}")
    print(f"  Predicted-only: {len(predicted_only_props)}")
    print(f"  Combined: {len(combined_props)}")

    if combined_props:
        create_property_distribution_figure(
            experimental_props,
            predicted_only_props,
            combined_props,
            NATURAL_DOMAINS,
            output_dir,
        )


# =============================================================================
# Main Visualization Pipeline
# =============================================================================

def main(args=None):
    """Main visualization script."""
    parser = argparse.ArgumentParser(description="Generate MOGRN visualizations")
    parser.add_argument(
        "--input-dir", type=str, default="opsin_output",
        help="Input directory with workflow results"
    )
    parser.add_argument(
        "--output-dir", type=str, default="opsin_output/paper_figures",
        help="Output directory for figures"
    )
    parser.add_argument(
        "--chain-id", type=str, default="A",
        help="Chain ID (default: A)"
    )
    parser.add_argument(
        "--skip-property", action="store_true",
        help="Skip property analysis"
    )
    parser.add_argument(
        "--skip-interactive", action="store_true",
        help="Skip interactive visualizations"
    )
    parser.add_argument(
        "--only", type=str,
        choices=["overview", "errors", "rmsd", "distance", "logo", "property", "interactive"],
        help="Generate only specific visualization"
    )

    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)

    print("\n" + "=" * 60)
    print("MOGRN VISUALIZATION SCRIPT")
    print("=" * 60)

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = input_dir / "cache"

    print(f"Input: {input_dir}")
    print(f"Output: {output_dir}")

    # Load workflow data
    print("\n[LOADING] Loading workflow data...")
    data = load_workflow_data(input_dir, args.chain_id)

    processed_structures = data.get("processed_structures", {})
    structure_mapping = data.get("structure_mapping", {})
    rmsd_df = data.get("rmsd_matrix") or data.get("rmsd_df")

    if not processed_structures:
        print("[ERROR] No processed structures found!")
        return

    print(f"[INFO] Loaded {len(processed_structures)} structures")
    print(f"[INFO] Structure mapping pairs: {len(structure_mapping)}")

    # Load property data
    property_file = PROPERTY_DIR / "mo_exp.csv"
    property_df = load_property_data(property_file)

    # Create property lookup
    property_data = {}
    for _, row in property_df.iterrows():
        pdb_id = str(row.get("PDB ID", "")).strip().lower()
        short_name = str(row.get("short_name", "")).strip()
        props = {
            "domain": row.get("Rhodopsin Type (Microbial)", "Unknown"),
            "molecular_function": row.get("molecular_function", "Unknown"),
        }
        if pdb_id:
            property_data[pdb_id] = props
        if short_name:
            property_data[short_name + "_model_0"] = props

    # Create group/domain dicts for visualizations
    group_dict = {}
    domain_dict = {}
    for sid in processed_structures:
        props = property_data.get(sid, {})
        group_dict[sid] = props.get("molecular_function", "Unknown")
        domain_dict[sid] = {"domain": props.get("domain", "Unknown")}

    # =========================================================================
    # Visualization 1: Overview Plot
    # =========================================================================
    if args.only is None or args.only == "overview":
        print("\n[VIZ 1] Opsin Overview Plot...")
        try:
            overview_list = []
            for sid in processed_structures:
                props = property_data.get(sid, {})
                is_pred = "_model_0" in sid or "_pred" in sid
                overview_list.append({
                    "id": sid,
                    "short_name": sid.replace("_model_0", ""),
                    "molecular_function": props.get("molecular_function", "Unknown"),
                    "domain": props.get("domain", "Unknown"),
                    "experimentally_determined": not is_pred,
                })

            overview_df = pd.DataFrame(overview_list)
            if not overview_df.empty:
                fig = create_opsin_overview_plot(overview_df)
                fig_path = output_dir / "01_opsin_overview.png"
                fig.savefig(fig_path, dpi=300, bbox_inches="tight")
                plt.close(fig)
                print(f"[OK] Saved: {fig_path}")
        except Exception as e:
            print(f"[ERROR] Overview plot failed: {e}")
            traceback.print_exc()

    # =========================================================================
    # Visualization 2: Error Distribution
    # =========================================================================
    if args.only is None or args.only == "errors":
        print("\n[VIZ 2] Error Distribution...")
        try:
            set_a_path = input_dir / "set_a_errors.csv"
            set_b_path = input_dir / "set_b_errors.csv"

            df_a = pd.read_csv(set_a_path) if set_a_path.exists() else pd.DataFrame()
            df_b = pd.read_csv(set_b_path) if set_b_path.exists() else pd.DataFrame()

            if not df_a.empty or not df_b.empty:
                if "dataset_split" not in df_a.columns and not df_a.empty:
                    df_a["dataset_split"] = "A"
                if "dataset_split" not in df_b.columns and not df_b.empty:
                    df_b["dataset_split"] = "B"

                fig_path = output_dir / "02c_error_distribution.png"
                fig, summary = plot_error_box_comparison(
                    df_a, df_b,
                    metrics=["backbone_rmsd", "pocket_rmsd", "retinal_rmsd"],
                    dataset_labels=("Benchmark set", "Blind test set"),
                    output_path=fig_path,
                )
                plt.close(fig)
                print(f"[OK] Saved: {fig_path}")
                if not summary.empty:
                    print(summary.round(3))
            else:
                print("[WARN] No error CSV files found")
        except Exception as e:
            print(f"[ERROR] Error plot failed: {e}")
            traceback.print_exc()

    # =========================================================================
    # Visualization 3: RMSD Clustermap
    # =========================================================================
    if args.only is None or args.only == "rmsd":
        print("\n[VIZ 3] RMSD Clustermap...")
        if isinstance(rmsd_df, pd.DataFrame) and not rmsd_df.empty:
            try:
                # Prepare matrix for linkage
                rmsd_clean = rmsd_df.copy()
                if rmsd_clean.isnull().values.any():
                    fill_val = np.nanmax(rmsd_clean.values[np.isfinite(rmsd_clean.values)])
                    rmsd_clean.fillna(fill_val, inplace=True)
                np.fill_diagonal(rmsd_clean.values, 0.0)

                condensed = squareform(rmsd_clean.values, checks=False)
                Z = linkage(condensed, method="weighted")

                # Compute metrics
                metrics = compute_rmsd_metrics(
                    rmsd_df=rmsd_df,
                    linkage_matrix=Z,
                    thresholds=(2.0, 2.5, 3.0),
                    n_clusters=2,
                    outdir=output_dir / "rmsd_metrics",
                )

                # Continuous colormap version
                fig_path = output_dir / "02a_rmsd_clustermap.png"
                fig = visualize_rmsd_matrix_improved(
                    rmsd_df=rmsd_df,
                    group_dict=group_dict,
                    domain_dict=domain_dict,
                    linkage_matrix=Z,
                    figsize=(18, 15),
                    output_file=fig_path,
                )
                if fig:
                    if hasattr(fig, "fig"):
                        fig.fig.savefig(fig_path, dpi=300)
                        plt.close(fig.fig)
                    else:
                        fig.savefig(fig_path, dpi=300)
                        plt.close(fig)
                    print(f"[OK] Saved: {fig_path}")

                # Step colormap version
                fig_step_path = output_dir / "02b_rmsd_clustermap_step.png"
                fig_step = visualize_rmsd_matrix_improved(
                    rmsd_df=rmsd_df,
                    group_dict=group_dict,
                    domain_dict=domain_dict,
                    linkage_matrix=Z,
                    figsize=(18, 15),
                    output_file=fig_step_path,
                    color_mode="step",
                    step_cutoffs=[0.5, 1.5, 2.5],
                )
                if fig_step:
                    if hasattr(fig_step, "fig"):
                        fig_step.fig.savefig(fig_step_path, dpi=300)
                        plt.close(fig_step.fig)
                    else:
                        fig_step.savefig(fig_step_path, dpi=300)
                        plt.close(fig_step)
                    print(f"[OK] Saved: {fig_step_path}")

            except Exception as e:
                print(f"[ERROR] RMSD clustermap failed: {e}")
                traceback.print_exc()
        else:
            print("[WARN] RMSD matrix not found")

    # =========================================================================
    # Visualization 4-5: Distance Plots
    # =========================================================================
    if args.only is None or args.only == "distance":
        print("\n[VIZ 4-5] Distance Plots...")

        distance_table = data.get("distance_table")
        ca_distance_table = data.get("ca_distance_table")

        # Load residue table for filtering
        residue_table = data.get("msa_df") or data.get("msa_table") or data.get("residue_table")
        if residue_table is None or (isinstance(residue_table, pd.DataFrame) and residue_table.empty):
            curated_grn_path = input_dir / "curated_grn.csv"
            if curated_grn_path.exists():
                residue_table = pd.read_csv(curated_grn_path, index_col=0)

        # Filter by 7.50 position (Lysine required for retinal binding)
        df_filtered = pd.DataFrame()
        if isinstance(residue_table, pd.DataFrame) and not residue_table.empty:
            filter_col = None
            for candidate in ["7.50", "7.5"]:
                if candidate in residue_table.columns:
                    filter_col = candidate
                    break

            if filter_col:
                mask = residue_table[filter_col].astype(str).str.startswith("K", na=False)
                df_filtered = residue_table[mask]
                print(f"[INFO] Filtered to {len(df_filtered)} entries with K at 7.50")

        # All-atom distance plot
        if isinstance(distance_table, pd.DataFrame) and not distance_table.empty:
            try:
                if not df_filtered.empty:
                    common_idx = distance_table.index.intersection(df_filtered.index)
                    plot_table = distance_table.loc[common_idx]
                else:
                    plot_table = distance_table

                grns = get_tm_residues(sort_grns_str(plot_table.columns.astype(str).tolist()))

                fig = plot_distances_with_std(
                    plot_table[grns] if grns else plot_table,
                    title="All-Atom Distance to Retinal",
                    use_ca=False,
                    figsize=(14, 8),
                )
                fig_path = output_dir / "03_all_atom_distance.png"
                fig.savefig(fig_path, dpi=300)
                plt.close(fig)
                print(f"[OK] Saved: {fig_path}")
            except Exception as e:
                print(f"[ERROR] All-atom distance plot failed: {e}")
                traceback.print_exc()

        # CA distance plot
        if isinstance(ca_distance_table, pd.DataFrame) and not ca_distance_table.empty:
            try:
                if not df_filtered.empty:
                    common_idx = ca_distance_table.index.intersection(df_filtered.index)
                    plot_table = ca_distance_table.loc[common_idx]
                else:
                    plot_table = ca_distance_table

                fig = plot_distances_with_std(
                    plot_table,
                    title="CA-Atom Distance to Retinal",
                    use_ca=True,
                    figsize=(14, 8),
                )
                fig_path = output_dir / "04_ca_distance.png"
                fig.savefig(fig_path, dpi=300)
                plt.close(fig)
                print(f"[OK] Saved: {fig_path}")
            except Exception as e:
                print(f"[ERROR] CA distance plot failed: {e}")
                traceback.print_exc()

    # =========================================================================
    # Visualization 6: Helix Logo Plots
    # =========================================================================
    if args.only is None or args.only == "logo":
        print("\n[VIZ 6] Helix Logo Plots...")

        residue_table = data.get("msa_df") or data.get("msa_table") or data.get("residue_table")
        if residue_table is None or (isinstance(residue_table, pd.DataFrame) and residue_table.empty):
            curated_grn_path = input_dir / "curated_grn.csv"
            if curated_grn_path.exists():
                residue_table = pd.read_csv(curated_grn_path, index_col=0)

        if isinstance(residue_table, pd.DataFrame) and not residue_table.empty:
            # Filter for K at 7.50
            filter_col = None
            for candidate in ["7.50", "7.5"]:
                if candidate in residue_table.columns:
                    filter_col = candidate
                    break

            if filter_col:
                mask = residue_table[filter_col].astype(str).str.startswith("K", na=False)
                df_filtered = residue_table[mask]

                try:
                    fig = plot_helix_logo_plots(df_filtered, frequency_threshold=0.07)
                    fig_path = output_dir / "05_helix_logos.png"
                    fig.savefig(fig_path, dpi=300)
                    plt.close(fig)
                    print(f"[OK] Saved: {fig_path}")
                except Exception as e:
                    print(f"[ERROR] Helix logo plot failed: {e}")
                    traceback.print_exc()
        else:
            print("[WARN] No residue table for logo plots")

    # =========================================================================
    # Visualization 7: Property Analysis
    # =========================================================================
    if not args.skip_property and (args.only is None or args.only == "property"):
        print("\n[VIZ 7] Property Analysis...")
        if processed_structures and not property_df.empty:
            try:
                create_property_analysis(
                    processed_structures,
                    property_df,
                    structure_mapping,
                    output_dir,
                )
            except Exception as e:
                print(f"[ERROR] Property analysis failed: {e}")
                traceback.print_exc()

    # =========================================================================
    # Visualization 8: Interactive GRN Alignment
    # =========================================================================
    if not args.skip_interactive and HAS_INTERACTIVE and (args.only is None or args.only == "interactive"):
        print("\n[VIZ 8] Interactive GRN Visualization...")
        try:
            interactive_path = output_dir / "interactive_grn_alignment.html"
            fig = create_opsin_visualization_from_workflow(
                cache_dir=str(cache_dir),
                property_file=str(property_file),
                output_file=str(interactive_path),
                max_structures=100,
                show_membrane=True,
                membrane_opacity=0.05,
            )
            if fig:
                print(f"[OK] Saved: {interactive_path}")

            # Enhanced version
            interactive_path_b = output_dir / "interactive_grn_alignment_b.html"
            fig_b = create_opsin_visualization_from_workflow_b(
                cache_dir=str(cache_dir),
                property_file=str(property_file),
                output_file=str(interactive_path_b),
                max_structures=100,
                show_membrane=True,
                membrane_opacity=0.05,
            )
            if fig_b:
                print(f"[OK] Saved: {interactive_path_b}")

        except Exception as e:
            print(f"[ERROR] Interactive visualization failed: {e}")
            traceback.print_exc()

    # =========================================================================
    # Summary
    # =========================================================================
    print("\n" + "=" * 60)
    print("VISUALIZATION COMPLETE")
    print("=" * 60)
    print(f"Figures saved to: {output_dir}")


if __name__ == "__main__":
    main()
