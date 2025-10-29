#!/usr/bin/env python3
"""
Visualize the alignments between clustered MO structures and their best matches.
Uses protos visualization capabilities.
"""

import os
import numpy as np
import pandas as pd
import json
from pathlib import Path
import sys
import pickle
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Add project root to path
project_root = Path(__file__).resolve().parent
sys.path.append(str(project_root))

# Add protos to path if needed
protos_path = project_root / "protos" / "src"
if protos_path.exists():
    sys.path.insert(0, str(protos_path))

# Import PROTOS modules
from protos.processing.structure.struct_base_processor import CifBaseProcessor
from protos.visualization.structure_vis import (
    plot_structures,
    plot_structure_alignment,
    structure_vis,
    visualize_rmsd_heatmap,
    create_and_visualize_similarity_tree
)
from Bio.PDB.qcprot import QCPSuperimposer

# Import workflow functions
from src.data_processing import filter_structures_by_chain_and_retinal


def get_ca_atoms(df: pd.DataFrame) -> pd.DataFrame:
    """Get CA atoms from structure dataframe, handling different column names."""
    if 'res_atom_name' in df.columns:
        return df[df['res_atom_name'] == 'CA'].copy()
    elif 'atom_name' in df.columns:
        return df[df['atom_name'] == 'CA'].copy()
    else:
        return pd.DataFrame()


def load_alignment_results():
    """Load the alignment results and RMSD matrix."""
    # Load alignment results
    with open(project_root / "clustered_mo_alignment_results_v2.json", 'r') as f:
        alignment_results = json.load(f)
    
    # Load RMSD matrix
    rmsd_matrix_df = pd.read_csv(project_root / "clustered_mo_rmsd_matrix_v2.csv", index_col=0)
    
    return alignment_results, rmsd_matrix_df


def load_structures(cache_dir: Path, user_data_root: Path):
    """Load both clustered and processed structures."""
    
    # Load processed structures from cache
    processed_cache_path = cache_dir / "processed_structures_A.pkl"
    with open(processed_cache_path, 'rb') as f:
        processed_data = pickle.load(f)
    processed_structures = processed_data.get('processed_structures', {})
    
    # Load clustered structures - use protos default location
    cp = CifBaseProcessor(
        name="clustered_mo_processor",
        preload=False
    )
    cp.load_dataset('clustered_mo')
    
    clustered_structures = {}
    for pdb_id in cp.pdb_ids:
        df = cp.data[cp.data['pdb_id'] == pdb_id].copy()
        df_chain_a = df[df['auth_chain_id'] == 'A'].copy()
        clustered_structures[pdb_id] = df_chain_a
    
    return clustered_structures, processed_structures


def align_structures_for_vis(struct1_df, struct2_df):
    """Align two structures and return aligned coordinates."""
    # Filter structures by GRN if available
    struct1_filtered = struct1_df.copy()
    struct2_filtered = struct2_df.copy()
    
    # Check if struct1 (processed structure) has GRN column and filter
    if 'grn' in struct1_filtered.columns:
        # Filter out rows where grn is null, empty string, or NaN
        struct1_filtered = struct1_filtered[
            (struct1_filtered['grn'].notna()) & 
            (struct1_filtered['grn'] != '') &
            (struct1_filtered['grn'] != 'nan')
        ].copy()
    
    # Get CA atoms
    ca1 = get_ca_atoms(struct1_filtered)
    ca2 = get_ca_atoms(struct2_filtered)
    
    if len(ca1) == 0 or len(ca2) == 0:
        return struct1_df, struct2_df, None
    
    # Sort by residue number
    ca1 = ca1.sort_values('auth_seq_id')
    ca2 = ca2.sort_values('auth_seq_id')
    
    # Get coordinates
    coords1 = ca1[['x', 'y', 'z']].values.astype(float)
    coords2 = ca2[['x', 'y', 'z']].values.astype(float)
    
    # Truncate to same length for alignment
    min_len = min(len(coords1), len(coords2))
    coords1_trunc = coords1[:min_len]
    coords2_trunc = coords2[:min_len]
    
    if len(coords1_trunc) < 10:
        return struct1_df, struct2_df, None
    
    # Calculate transformation using QCPSuperimposer
    # Note: QCPSuperimposer expects (reference, mobile) order
    sup = QCPSuperimposer()
    sup.set(coords1_trunc, coords2_trunc)  # Reference, mobile
    sup.run()
    
    # Get rotation matrix and translation vector
    rot_matrix = sup.rot
    tran_vector = sup.tran
    rmsd = sup.get_rms()
    
    # Apply transformation to ALL atoms of struct2 (not just CA or filtered)
    struct2_aligned = struct2_df.copy()
    coords_all = struct2_df[['x', 'y', 'z']].values.astype(float)
    
    # Apply rotation first, then translation
    # QCPSuperimposer gives rotation matrix that should be applied as: rotated = original @ rotation + translation
    coords_transformed = np.dot(coords_all, rot_matrix) + tran_vector
    
    # Update coordinates in the dataframe
    struct2_aligned[['x', 'y', 'z']] = coords_transformed
    
    # Also return the filtered reference structure for visualization
    return struct1_filtered, struct2_aligned, rmsd


def visualize_single_alignment(clustered_id, best_match_id, clustered_struct, processed_struct, output_dir):
    """Create visualization for a single alignment."""
    
    # Align structures
    ref_struct, aligned_struct, rmsd = align_structures_for_vis(processed_struct, clustered_struct)
    
    if rmsd is None:
        print(f"  Skipping {clustered_id} - alignment failed")
        return
    
    # Create figure showing alignment
    fig = go.Figure()
    
    # Add reference structure (best match) - use filtered structure
    ca_ref = get_ca_atoms(ref_struct)
    if len(ca_ref) > 0:
        fig.add_trace(go.Scatter3d(
            x=ca_ref['x'],
            y=ca_ref['y'],
            z=ca_ref['z'],
            mode='lines+markers',
            marker=dict(size=3, color='blue'),
            line=dict(color='blue', width=2),
            name=f'{best_match_id} (Reference, GRN-filtered)',
            opacity=0.7
        ))
    
    # Add aligned clustered structure
    ca_aligned = get_ca_atoms(aligned_struct)
    if len(ca_aligned) > 0:
        fig.add_trace(go.Scatter3d(
            x=ca_aligned['x'],
            y=ca_aligned['y'],
            z=ca_aligned['z'],
            mode='lines+markers',
            marker=dict(size=3, color='red'),
            line=dict(color='red', width=2),
            name=f'{clustered_id} (Aligned)',
            opacity=0.7
        ))
    
    # Update layout
    fig.update_layout(
        title=f'Alignment: {clustered_id} → {best_match_id}<br>RMSD: {rmsd:.2f} Å',
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode='data'
        ),
        width=800,
        height=800
    )
    
    # Save as interactive HTML
    html_file = output_dir / f"alignment_{clustered_id}_to_{best_match_id}.docs"
    fig.write_html(str(html_file))
    
    # Save as static image
    png_file = output_dir / f"alignment_{clustered_id}_to_{best_match_id}.png"
    fig.write_image(str(png_file))
    
    print(f"  Saved visualization for {clustered_id} → {best_match_id}")


def create_overview_visualization(alignment_results, rmsd_matrix_df, output_dir):
    """Create overview visualizations of all alignments."""
    
    # 1. RMSD Heatmap
    plt.figure(figsize=(20, 6))
    
    # Extract RMSD values for best matches
    clustered_ids = list(alignment_results['alignments'].keys())
    best_match_ids = [alignment_results['alignments'][cid]['best_match_id'] for cid in clustered_ids]
    best_rmsds = [alignment_results['alignments'][cid]['best_rmsd'] for cid in clustered_ids]
    
    # Create bar plot of best RMSDs
    plt.subplot(1, 2, 1)
    bars = plt.bar(range(len(clustered_ids)), best_rmsds)
    plt.xlabel('Clustered Structure')
    plt.ylabel('Best RMSD (Å)')
    plt.title('Best RMSD for Each Clustered Structure')
    plt.xticks(range(len(clustered_ids)), [cid.replace('_model_0', '') for cid in clustered_ids], rotation=45, ha='right')
    
    # Color bars by RMSD value
    for i, (bar, rmsd) in enumerate(zip(bars, best_rmsds)):
        if rmsd < 5:
            bar.set_color('green')
        elif rmsd < 10:
            bar.set_color('yellow')
        elif rmsd < 15:
            bar.set_color('orange')
        else:
            bar.set_color('red')
    
    # 2. Distribution of RMSDs
    plt.subplot(1, 2, 2)
    plt.hist([r for r in best_rmsds if r is not None and r != float('inf')], bins=20, edgecolor='black')
    plt.xlabel('RMSD (Å)')
    plt.ylabel('Count')
    plt.title('Distribution of Best RMSDs')
    
    plt.tight_layout()
    plt.savefig(output_dir / "rmsd_overview.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 3. Create a smaller heatmap of clustered vs top processed structures
    # Find most common best matches
    from collections import Counter
    match_counts = Counter(best_match_ids)
    top_matches = [match for match, _ in match_counts.most_common(10)]
    
    # Extract relevant part of RMSD matrix
    if all(col in rmsd_matrix_df.columns for col in top_matches):
        subset_matrix = rmsd_matrix_df[top_matches]
        
        plt.figure(figsize=(12, 10))
        im = plt.imshow(subset_matrix.values, cmap='viridis_r', aspect='auto')
        plt.colorbar(im, label='RMSD (Å)')
        
        # Set ticks
        plt.yticks(range(len(subset_matrix.index)), [idx.replace('_model_0', '') for idx in subset_matrix.index])
        plt.xticks(range(len(top_matches)), [m.replace('_model_0', '') for m in top_matches], rotation=45, ha='right')
        
        plt.xlabel('Processed Structures (Top Matches)')
        plt.ylabel('Clustered Structures')
        plt.title('RMSD Matrix: Clustered vs Most Matched Structures')
        plt.tight_layout()
        plt.savefig(output_dir / "rmsd_heatmap_subset.png", dpi=300, bbox_inches='tight')
        plt.close()


def main():
    """Main function to create visualizations."""
    
    # Setup
    cache_dir = project_root / "opsin_output" / "cache"
    user_data_root = project_root / "data_clustered_mo"
    output_dir = project_root / "opsin_output" / "clustered"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Loading alignment results...")
    alignment_results, rmsd_matrix_df = load_alignment_results()
    
    print("Loading structures...")
    clustered_structures, processed_structures = load_structures(cache_dir, user_data_root)
    
    # Create overview visualizations
    print("\nCreating overview visualizations...")
    create_overview_visualization(alignment_results, rmsd_matrix_df, output_dir)
    
    # Visualize top 5 best alignments
    print("\nCreating individual alignment visualizations...")
    
    # Sort by RMSD to get best alignments
    alignments = alignment_results['alignments']
    sorted_alignments = sorted(
        [(cid, info) for cid, info in alignments.items() if info['best_match_id'] is not None],
        key=lambda x: x[1]['best_rmsd']
    )
    
    # Visualize top 5
    for i, (clustered_id, info) in enumerate(sorted_alignments[:5]):
        best_match_id = info['best_match_id']
        rmsd = info['best_rmsd']
        
        print(f"\n{i+1}. {clustered_id} → {best_match_id} (RMSD: {rmsd:.2f} Å)")
        
        if clustered_id in clustered_structures and best_match_id in processed_structures:
            clustered_struct = clustered_structures[clustered_id]
            
            # Get processed structure
            if isinstance(processed_structures[best_match_id], dict):
                processed_struct = processed_structures[best_match_id]['df']
            else:
                processed_struct = processed_structures[best_match_id]
            
            # Create visualization
            visualize_single_alignment(
                clustered_id, best_match_id,
                clustered_struct, processed_struct,
                output_dir
            )
    
    # Create summary figure showing all best matches
    print("\nCreating summary figure...")
    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=[f"{cid.replace('_model_0', '')}" for cid, _ in sorted_alignments[:6]],
        specs=[[{"type": "scatter3d"} for _ in range(3)] for _ in range(2)]
    )
    
    for idx, (clustered_id, info) in enumerate(sorted_alignments[:6]):
        row = idx // 3 + 1
        col = idx % 3 + 1
        
        if clustered_id in clustered_structures and info['best_match_id'] in processed_structures:
            # Get structures
            clustered_struct = clustered_structures[clustered_id]
            if isinstance(processed_structures[info['best_match_id']], dict):
                processed_struct = processed_structures[info['best_match_id']]['df']
            else:
                processed_struct = processed_structures[info['best_match_id']]
            
            # Align
            ref_struct, aligned_struct, _ = align_structures_for_vis(processed_struct, clustered_struct)
            
            # Add traces
            ca_ref = get_ca_atoms(ref_struct)
            ca_aligned = get_ca_atoms(aligned_struct)
            
            if len(ca_ref) > 0:
                fig.add_trace(
                    go.Scatter3d(
                        x=ca_ref['x'], y=ca_ref['y'], z=ca_ref['z'],
                        mode='lines',
                        line=dict(color='blue', width=2),
                        showlegend=False,
                        opacity=0.5
                    ),
                    row=row, col=col
                )
            
            if len(ca_aligned) > 0:
                fig.add_trace(
                    go.Scatter3d(
                        x=ca_aligned['x'], y=ca_aligned['y'], z=ca_aligned['z'],
                        mode='lines',
                        line=dict(color='red', width=2),
                        showlegend=False,
                        opacity=0.5
                    ),
                    row=row, col=col
                )
    
    fig.update_layout(
        title='Top 6 Best Alignments (Blue: Reference, Red: Clustered)',
        showlegend=False,
        height=800,
        width=1200
    )
    
    # Update all scenes to hide axes
    for i in range(1, 7):
        fig.update_scenes(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            aspectmode='data',
            row=(i-1)//3 + 1,
            col=(i-1)%3 + 1
        )
    
    fig.write_html(str(output_dir / "top_alignments_grid.docs"))
    fig.write_image(str(output_dir / "top_alignments_grid.png"))
    
    print(f"\nVisualizations saved to: {output_dir}")
    print("Done!")


if __name__ == "__main__":
    main()