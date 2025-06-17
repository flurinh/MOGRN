#!/usr/bin/env python3
"""
Simple rotation animation using matplotlib with robust error handling.
Creates 1080 frames of a rotating opsin structure with highlighted residues.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path
from tqdm import tqdm
import subprocess

from animate_grn_highlights import (
    load_rmsd_cache, load_processed_structures, load_grn_table,
    extract_ca_coordinates_with_grn, apply_alignment_transformations,
    apply_membrane_orientation
)
from src.opsin_color_scheme import HELIX_NUMBER_COLORS


def create_simple_frame(
    all_coords,
    highlight_coords,
    highlight_colors,
    camera_angle,
    target_grn,
    figsize=(19.2, 10.8),
    dpi=100
):
    """
    Create a simple 3D frame with basic error handling.
    
    Args:
        all_coords: All structure coordinates as numpy array
        highlight_coords: Highlighted residue coordinates
        highlight_colors: Colors for highlighted residues
        camera_angle: Camera rotation angle
        target_grn: GRN position being highlighted
        figsize: Figure size
        dpi: Resolution
    """
    
    # Create figure
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor='white')
    ax = fig.add_subplot(111, projection='3d')
    ax.set_facecolor('white')
    
    # Plot all structures as scatter points (simpler than lines)
    if len(all_coords) > 0:
        # Sample points for performance (every 3rd point)
        sample_coords = all_coords[::3]
        ax.scatter(sample_coords[:, 0], sample_coords[:, 1], sample_coords[:, 2],
                  c='lightgray', s=0.5, alpha=0.3)
    
    # Plot highlights
    if len(highlight_coords) > 0:
        coords_array = np.array(highlight_coords)
        ax.scatter(coords_array[:, 0], coords_array[:, 1], coords_array[:, 2],
                  c=highlight_colors, s=100, alpha=1.0, edgecolors='white', linewidth=2)
    
    # Set camera angle
    ax.view_init(elev=20, azim=camera_angle)
    
    # Set title
    helix_num = target_grn.split('.')[0]
    ax.set_title(f'MOGRN Helix {helix_num} Position 50 | Rotation: {camera_angle:.1f}°',
                fontsize=20, pad=20, fontweight='bold')
    
    # Set axis limits manually to avoid empty bbox error
    if len(all_coords) > 0:
        x_center = np.mean(all_coords[:, 0])
        y_center = np.mean(all_coords[:, 1])
        z_center = np.mean(all_coords[:, 2])
        
        # Use fixed ranges around center
        range_size = 40
        ax.set_xlim([x_center - range_size, x_center + range_size])
        ax.set_ylim([y_center - range_size, y_center + range_size])
        ax.set_zlim([z_center - range_size/2, z_center + range_size/2])
    else:
        # Fallback fixed limits
        ax.set_xlim([-40, 40])
        ax.set_ylim([-40, 40])
        ax.set_zlim([-20, 20])
    
    # Clean appearance
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    
    # Hide axis elements
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.line.set_visible(False)
    ax.yaxis.line.set_visible(False)
    ax.zaxis.line.set_visible(False)
    
    plt.tight_layout()
    
    return fig


def prepare_animation_data(aligned_structures, grn_df):
    """
    Prepare all data needed for animation to avoid processing during frame creation.
    """
    print("=== Preparing animation data ===")
    
    # Collect all coordinates
    all_coords = []
    structure_coords = {}
    
    # Color scheme
    helix_colors = {
        1: HELIX_NUMBER_COLORS[1], 2: HELIX_NUMBER_COLORS[2], 
        3: HELIX_NUMBER_COLORS[3], 4: HELIX_NUMBER_COLORS[4],
        5: HELIX_NUMBER_COLORS[5], 6: HELIX_NUMBER_COLORS[6],
        7: HELIX_NUMBER_COLORS[7], 0: '#D3D3D3'
    }
    
    # Process each structure
    for struct_id, data in aligned_structures.items():
        coords = data['coords']
        df = data.get('dataframe')
        
        if df is not None and len(coords) > 0:
            # Add to all coordinates
            all_coords.extend(coords)
            
            # Store structure data
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = coords
            structure_coords[struct_id] = df_aligned
    
    all_coords = np.array(all_coords)
    print(f"Total coordinates: {len(all_coords)}")
    
    # Pre-calculate highlight data for each GRN position
    grn_highlights = {}
    helix_grns = ["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"]
    
    for grn_pos in helix_grns:
        highlight_coords = []
        highlight_colors = []
        
        for struct_id, df_aligned in structure_coords.items():
            # Find residues with this GRN
            grn_matches = df_aligned[df_aligned['grn'] == grn_pos]
            
            for _, row in grn_matches.iterrows():
                highlight_coords.append([row['x'], row['y'], row['z']])
                helix_num = row['helix_num']
                highlight_colors.append(helix_colors.get(helix_num, '#D3D3D3'))
        
        grn_highlights[grn_pos] = {
            'coords': highlight_coords,
            'colors': highlight_colors
        }
        print(f"GRN {grn_pos}: {len(highlight_coords)} highlights")
    
    return all_coords, grn_highlights


def create_rotation_animation():
    """Create the full rotation animation with robust error handling."""
    
    print("=== MOGRN Simple Rotation Animation ===\n")
    
    # Parameters
    total_frames = 1080
    fps = 60
    output_dir = Path("rotation_simple_1080")
    width = 1920
    height = 1080
    dpi = 100
    
    fig_width = width / dpi
    fig_height = height / dpi
    
    print("=== Loading and Processing Data ===")
    try:
        cache_data = load_rmsd_cache()
        alignment_paths = cache_data.get('alignment_paths', {})
        processed_structures = load_processed_structures()
        grn_df = load_grn_table()
        
        structures = extract_ca_coordinates_with_grn(processed_structures, grn_df, chain_id='A', use_helix_only=True)
        aligned_structures = apply_alignment_transformations(structures, alignment_paths, 'MerMAID1_model_0')
        oriented_structures = apply_membrane_orientation(aligned_structures, 'MerMAID1_model_0')
        
        print(f"Processed {len(oriented_structures)} structures")
        
    except Exception as e:
        print(f"❌ Data loading failed: {e}")
        return
    
    # Prepare all animation data upfront
    try:
        all_coords, grn_highlights = prepare_animation_data(oriented_structures, grn_df)
        if len(all_coords) == 0:
            print("❌ No coordinate data found!")
            return
    except Exception as e:
        print(f"❌ Data preparation failed: {e}")
        return
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    
    print(f"\n=== Generating {total_frames} frames ===")
    print(f"Resolution: {width}x{height}")
    print(f"Duration: {total_frames/fps:.1f} seconds")
    
    # Define helix rotation pattern
    helix_grns = ["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"]
    
    successful_frames = 0
    
    # Generate frames
    for frame_num in tqdm(range(total_frames), desc="Creating frames"):
        angle = (frame_num / total_frames) * 360
        
        # Cycle through helices based on angle
        helix_idx = int((angle / 360) * 7) % 7
        target_grn = helix_grns[helix_idx]
        
        try:
            # Get highlight data for this GRN
            highlight_data = grn_highlights.get(target_grn, {'coords': [], 'colors': []})
            
            # Create frame
            fig = create_simple_frame(
                all_coords=all_coords,
                highlight_coords=highlight_data['coords'],
                highlight_colors=highlight_data['colors'],
                camera_angle=angle,
                target_grn=target_grn,
                figsize=(fig_width, fig_height),
                dpi=dpi
            )
            
            # Save frame
            frame_file = output_dir / f"frame_{frame_num:05d}.png"
            fig.savefig(frame_file, dpi=dpi, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            plt.close(fig)
            
            successful_frames += 1
            
        except Exception as e:
            if frame_num < 10:  # Only print first few errors
                print(f"\nError on frame {frame_num}: {e}")
            plt.close('all')  # Clean up
            continue
    
    print(f"\n✅ Successfully created {successful_frames}/{total_frames} frames")
    
    if successful_frames == 0:
        print("❌ No frames were created successfully!")
        return
    
    # Create video
    print("\n=== Creating video ===")
    video_file = output_dir / "mogrn_rotation.mp4"
    
    cmd = [
        'ffmpeg', '-y',
        '-framerate', str(fps),
        '-i', str(output_dir / 'frame_%05d.png'),
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p',
        '-crf', '18',
        str(video_file)
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"\n✅ Animation complete!")
        print(f"Video: {video_file}")
        print(f"Frames: {successful_frames}")
        print(f"Duration: {total_frames/fps:.1f} seconds")
    except Exception as e:
        print(f"\n❌ Video creation failed: {e}")
        print("PNG frames are available in:", output_dir)


if __name__ == "__main__":
    create_rotation_animation()