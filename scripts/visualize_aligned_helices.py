#!/usr/bin/env python3
"""
Visualize aligned structures colored by helix number.

1. Load structures and annotate with helix numbers from helices_extended.json
2. Align structures based on TM bundle (helices 1-7) CA atoms
3. Visualize all aligned structures with helix coloring
   - Helix 1-7: distinct colors
   - Loops/tails: grey
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple
import numpy as np
import pandas as pd
import plotly.graph_objects as go

# Project setup
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT / "protos" / "src"))

import protos
protos.set_data_path(str(PROJECT_ROOT / "data"))

from protos.processing.structure import StructureProcessor
from protos.analysis.structure.alignment import kabsch_alignment


# Helix colors (rainbow-ish for visibility)
HELIX_COLORS = {
    "1": "#E41A1C",  # Red
    "2": "#FF7F00",  # Orange
    "3": "#FFFF33",  # Yellow
    "4": "#4DAF4A",  # Green
    "5": "#377EB8",  # Blue
    "6": "#984EA3",  # Purple
    "7": "#F781BF",  # Pink
}
LOOP_COLOR = "#AAAAAA"  # Grey for loops/tails


def load_helices() -> Dict[str, Dict[str, List[int]]]:
    """Load extended helix definitions."""
    helix_file = PROJECT_ROOT / "property" / "helices_extended.json"
    with open(helix_file) as f:
        return json.load(f)


def annotate_helix(df: pd.DataFrame, helices: Dict[str, List[int]]) -> pd.DataFrame:
    """
    Annotate structure with helix numbers.

    Args:
        df: Structure dataframe with auth_seq_id
        helices: Dict of {helix_num: [start, end]}

    Returns:
        DataFrame with 'helix' column added
    """
    df = df.copy()
    df['helix'] = None

    for helix_num, (start, end) in helices.items():
        mask = (df['auth_seq_id'] >= start) & (df['auth_seq_id'] <= end)
        df.loc[mask, 'helix'] = helix_num

    return df


def get_tm_ca_coords(df: pd.DataFrame, helices: Dict[str, List[int]]) -> pd.DataFrame:
    """
    Get CA coordinates for TM helices only.

    Args:
        df: Structure dataframe
        helices: Helix definitions

    Returns:
        DataFrame with CA coordinates from TM helices
    """
    atom_col = 'res_atom_name' if 'res_atom_name' in df.columns else 'atom_name'

    # Get CA atoms
    ca_df = df[df[atom_col] == 'CA'].copy()

    # Filter to TM helix residues
    tm_mask = pd.Series(False, index=ca_df.index)
    for helix_num, (start, end) in helices.items():
        tm_mask |= (ca_df['auth_seq_id'] >= start) & (ca_df['auth_seq_id'] <= end)

    return ca_df[tm_mask].sort_values('auth_seq_id')


def align_structure(ref_coords: np.ndarray,
                    struct_df: pd.DataFrame,
                    struct_helices: Dict[str, List[int]]) -> Tuple[pd.DataFrame, float]:
    """
    Align structure to reference based on TM bundle.

    Args:
        ref_coords: Reference CA coordinates (N x 3)
        struct_df: Structure to align
        struct_helices: Helix definitions for structure

    Returns:
        Tuple of (aligned structure DataFrame, RMSD)
    """
    # Get TM CA coordinates
    tm_ca = get_tm_ca_coords(struct_df, struct_helices)

    if len(tm_ca) == 0:
        return struct_df, float('nan')

    struct_coords = tm_ca[['x', 'y', 'z']].values

    # Match lengths (use minimum)
    min_len = min(len(ref_coords), len(struct_coords))
    ref_subset = ref_coords[:min_len]
    struct_subset = struct_coords[:min_len]

    # Calculate alignment
    try:
        rot, trans, rmsd = kabsch_alignment(ref_subset, struct_subset)
    except Exception as e:
        print(f"    Alignment failed: {e}")
        return struct_df, float('nan')

    # Apply transformation to entire structure
    aligned_df = struct_df.copy()
    coords = aligned_df[['x', 'y', 'z']].values
    aligned_coords = (rot @ coords.T).T + trans
    aligned_df[['x', 'y', 'z']] = aligned_coords

    return aligned_df, rmsd


def create_structure_trace(df: pd.DataFrame,
                           struct_id: str,
                           show_backbone: bool = True,
                           opacity: float = 0.7) -> List[go.Scatter3d]:
    """
    Create plotly traces for a structure colored by helix.

    Args:
        df: Annotated structure dataframe
        struct_id: Structure identifier
        show_backbone: Whether to show backbone trace
        opacity: Marker opacity

    Returns:
        List of plotly traces
    """
    atom_col = 'res_atom_name' if 'res_atom_name' in df.columns else 'atom_name'

    # Get CA atoms only for visualization
    ca_df = df[df[atom_col] == 'CA'].copy().sort_values('auth_seq_id')

    traces = []

    # Create trace for each helix + loops
    for helix_num in list(HELIX_COLORS.keys()) + [None]:
        if helix_num is None:
            mask = ca_df['helix'].isna()
            color = LOOP_COLOR
            name = f"{struct_id} - loops"
        else:
            mask = ca_df['helix'] == helix_num
            color = HELIX_COLORS[helix_num]
            name = f"{struct_id} - H{helix_num}"

        subset = ca_df[mask]
        if len(subset) == 0:
            continue

        trace = go.Scatter3d(
            x=subset['x'],
            y=subset['y'],
            z=subset['z'],
            mode='markers',
            marker=dict(size=3, color=color, opacity=opacity),
            name=name,
            legendgroup=struct_id,
            showlegend=False,
            hovertext=[f"{struct_id}<br>Res {r}<br>H{h if h else 'loop'}"
                      for r, h in zip(subset['auth_seq_id'], subset['helix'])],
            hoverinfo='text',
        )
        traces.append(trace)

    # Add backbone trace if requested
    if show_backbone:
        backbone_trace = go.Scatter3d(
            x=ca_df['x'],
            y=ca_df['y'],
            z=ca_df['z'],
            mode='lines',
            line=dict(color='grey', width=1),
            opacity=0.3,
            name=f"{struct_id} backbone",
            legendgroup=struct_id,
            showlegend=True,
            hoverinfo='skip',
        )
        traces.insert(0, backbone_trace)

    return traces


def visualize_aligned_structures(structures: Dict[str, pd.DataFrame],
                                  title: str = "Aligned Opsin Structures",
                                  max_structures: int = 20) -> go.Figure:
    """
    Create visualization of aligned structures.

    Args:
        structures: Dict of {struct_id: aligned_dataframe}
        title: Plot title
        max_structures: Maximum structures to show

    Returns:
        Plotly figure
    """
    fig = go.Figure()

    struct_ids = list(structures.keys())[:max_structures]

    for struct_id in struct_ids:
        df = structures[struct_id]
        traces = create_structure_trace(df, struct_id)
        for trace in traces:
            fig.add_trace(trace)

    # Add helix color legend
    for helix_num, color in HELIX_COLORS.items():
        fig.add_trace(go.Scatter3d(
            x=[None], y=[None], z=[None],
            mode='markers',
            marker=dict(size=10, color=color),
            name=f"Helix {helix_num}",
            showlegend=True,
        ))

    # Add loop legend
    fig.add_trace(go.Scatter3d(
        x=[None], y=[None], z=[None],
        mode='markers',
        marker=dict(size=10, color=LOOP_COLOR),
        name="Loops/Tails",
        showlegend=True,
    ))

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis=dict(visible=False, showgrid=False, zeroline=False),
            yaxis=dict(visible=False, showgrid=False, zeroline=False),
            zaxis=dict(visible=False, showgrid=False, zeroline=False),
            bgcolor='white',
            aspectmode='data',
        ),
        width=1200,
        height=900,
        paper_bgcolor='white',
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
        ),
    )

    return fig


def main():
    """Main function."""
    print("=" * 60)
    print("VISUALIZING ALIGNED STRUCTURES BY HELIX")
    print("=" * 60)

    # Load helix definitions
    helices = load_helices()
    print(f"[INFO] Loaded helix definitions for {len(helices)} structures")

    # Initialize processor
    processor = StructureProcessor("vis_align")

    # Select reference structure (use 4fbz or CnChR2)
    reference_id = "CnChR2_J230_refine9"
    if reference_id not in helices:
        reference_id = "4fbz"

    print(f"[INFO] Using reference: {reference_id}")

    # Load reference
    ref_df = processor.load_entity(reference_id)
    if ref_df is None:
        print(f"[ERROR] Could not load reference {reference_id}")
        return

    ref_df = ref_df.reset_index()
    if 'auth_chain_id' in ref_df.columns:
        ref_df = ref_df[ref_df['auth_chain_id'] == 'A']
    if 'atom_name' in ref_df.columns and 'res_atom_name' not in ref_df.columns:
        ref_df['res_atom_name'] = ref_df['atom_name']

    # Annotate reference
    ref_df = annotate_helix(ref_df, helices[reference_id])

    # Get reference TM CA coords
    ref_tm_ca = get_tm_ca_coords(ref_df, helices[reference_id])
    ref_coords = ref_tm_ca[['x', 'y', 'z']].values
    print(f"[INFO] Reference TM CA atoms: {len(ref_coords)}")

    # Select structures to visualize (sample from each dataset)
    datasets = ["mo_exp_A", "mo_exp_B", "mo_pred_exp", "mo_pred_novel"]
    sample_per_dataset = 5

    structures_to_load = [reference_id]
    for dataset in datasets:
        if processor.dataset_manager.dataset_exists(dataset):
            entities = processor.get_dataset_entities(dataset)
            # Filter to those with helix definitions
            entities = [e for e in entities if e in helices and e != reference_id]
            structures_to_load.extend(entities[:sample_per_dataset])

    print(f"\n[INFO] Processing {len(structures_to_load)} structures...")

    # Process and align structures
    aligned_structures = {}
    aligned_structures[reference_id] = ref_df  # Reference is already "aligned"

    rmsds = {}

    for struct_id in structures_to_load:
        if struct_id == reference_id:
            continue

        if struct_id not in helices:
            print(f"  {struct_id}: No helix definition, skipping")
            continue

        print(f"  {struct_id}...", end=" ", flush=True)

        try:
            df = processor.load_entity(struct_id)
            if df is None:
                print("LOAD FAILED")
                continue

            df = df.reset_index()

            # For predicted structures, set all chains to A
            if '_model_0' in struct_id:
                df['auth_chain_id'] = 'A'

            if 'auth_chain_id' in df.columns:
                df = df[df['auth_chain_id'] == 'A']

            if 'atom_name' in df.columns and 'res_atom_name' not in df.columns:
                df['res_atom_name'] = df['atom_name']

            # Annotate with helices
            df = annotate_helix(df, helices[struct_id])

            # Align to reference
            aligned_df, rmsd = align_structure(ref_coords, df, helices[struct_id])

            if not np.isnan(rmsd):
                aligned_structures[struct_id] = aligned_df
                rmsds[struct_id] = rmsd
                print(f"OK (RMSD: {rmsd:.2f} Å)")
            else:
                print("ALIGNMENT FAILED")

        except Exception as e:
            print(f"ERROR: {e}")

    print(f"\n[INFO] Successfully aligned {len(aligned_structures)} structures")

    # Show RMSD statistics
    if rmsds:
        rmsd_values = list(rmsds.values())
        print(f"[INFO] RMSD range: {min(rmsd_values):.2f} - {max(rmsd_values):.2f} Å")
        print(f"[INFO] Mean RMSD: {np.mean(rmsd_values):.2f} Å")

    # Create visualization
    print("\n[INFO] Creating visualization...")
    fig = visualize_aligned_structures(
        aligned_structures,
        title=f"Aligned Opsin Structures (n={len(aligned_structures)}, ref={reference_id})",
        max_structures=25
    )

    # Save to HTML
    output_file = PROJECT_ROOT / "opsin_output" / "aligned_helices_visualization.html"
    output_file.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_file))
    print(f"[INFO] Saved visualization to {output_file}")

    # Also save a static image if possible
    try:
        png_file = PROJECT_ROOT / "opsin_output" / "aligned_helices_visualization.png"
        fig.write_image(str(png_file), width=1200, height=900)
        print(f"[INFO] Saved PNG to {png_file}")
    except Exception as e:
        print(f"[INFO] Could not save PNG (install kaleido): {e}")

    return fig, aligned_structures


if __name__ == "__main__":
    main()
