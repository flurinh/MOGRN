#!/usr/bin/env python3
"""
Create 2D front-view projections of rotating 3D structures.
This avoids matplotlib 3D issues by projecting to 2D coordinates.
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from tqdm import tqdm
import subprocess

from animate_grn_highlights import (
    load_rmsd_cache, load_processed_structures, load_grn_table,
    extract_ca_coordinates_with_grn, apply_alignment_transformations,
    apply_membrane_orientation
)
from src.opsin_color_scheme import HELIX_NUMBER_COLORS


def project_to_2d_frontview(coords_3d, camera_angle):
    """
    Project 3D coordinates to 2D front view at given camera angle.
    
    Args:
        coords_3d: Array of 3D coordinates [N, 3]
        camera_angle: Rotation angle in degrees around Z-axis
        
    Returns:
        Array of 2D coordinates [N, 2] for front view
    """
    angle_rad = np.radians(camera_angle)
    
    # Rotation matrix around Z-axis
    cos_a = np.cos(angle_rad)
    sin_a = np.sin(angle_rad)
    rotation_matrix = np.array([
        [cos_a, -sin_a, 0],
        [sin_a,  cos_a, 0],
        [0,      0,     1]
    ])
    
    # Rotate the coordinates
    rotated_coords = np.dot(coords_3d, rotation_matrix.T)
    
    # Project to front view (X-Z plane from rotated view)
    # This gives us a side view of the membrane protein
    front_view_2d = rotated_coords[:, [0, 2]]  # X and Z coordinates
    
    return front_view_2d


def create_frontview_frame(
    all_coords_3d,
    highlight_coords_3d,
    highlight_colors,
    camera_angle,
    target_grn,
    figsize=(19.2, 10.8),
    dpi=100
):
    """Create a 2D front-view frame from 3D coordinates."""
    
    # Create 2D plot
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi, facecolor='white')
    ax.set_facecolor('white')
    
    # Project background structures to 2D
    if len(all_coords_3d) > 0:
        # Sample for performance (every 4th point)
        sample_coords = all_coords_3d[::4]
        projected_bg = project_to_2d_frontview(sample_coords, camera_angle)
        
        # Plot background as small gray dots
        ax.scatter(projected_bg[:, 0], projected_bg[:, 1], 
                  c='lightgray', s=1, alpha=0.4, rasterized=True)
    
    # Project and plot highlights
    if len(highlight_coords_3d) > 0:
        highlight_array = np.array(highlight_coords_3d)
        projected_highlights = project_to_2d_frontview(highlight_array, camera_angle)
        
        # Plot highlighted residues as larger colored dots
        ax.scatter(projected_highlights[:, 0], projected_highlights[:, 1],
                  c=highlight_colors, s=120, alpha=1.0, 
                  edgecolors='white', linewidth=2, rasterized=True, zorder=10)
    
    # Set axis properties
    ax.set_xlim([-60, 60])
    ax.set_ylim([-35, 35])
    ax.set_aspect('equal')
    
    # Remove ticks and spines for clean look
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(False)
    
    # Add title
    helix_num = target_grn.split('.')[0]
    ax.set_title(f'MOGRN Helix {helix_num} Position 50 | Rotation: {camera_angle:.1f}°',
                fontsize=22, pad=20, fontweight='bold')
    
    # Add angle indicator in corner
    ax.text(0.98, 0.02, f'{camera_angle:.0f}°', transform=ax.transAxes,
           fontsize=18, ha='right', va='bottom',
           bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    # Add membrane indicator (horizontal line at Z=0)
    ax.axhline(y=0, color='lightblue', linestyle='--', alpha=0.5, linewidth=2)
    ax.text(0.02, 0.02, 'Membrane', transform=ax.transAxes,
           fontsize=12, ha='left', va='bottom', color='blue', alpha=0.7)
    
    plt.tight_layout()
    return fig


def create_frontview_animation():
    """Create front-view rotation animation."""
    
    print("=== MOGRN Front-View Rotation Animation ===")
    
    # Parameters
    total_frames = 1080
    fps = 60
    output_dir = Path("rotation_frontview_1080")
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
    
    print("=== Preparing Data ===")
    
    # Collect all 3D coordinates
    all_coords_3d = []
    for struct_id, data in oriented_structures.items():
        coords = data['coords']
        if len(coords) > 0:
            all_coords_3d.extend(coords)
    
    all_coords_3d = np.array(all_coords_3d)
    print(f"Total 3D coordinates: {len(all_coords_3d)}")
    
    # Prepare highlight data for each GRN
    helix_colors = {
        1: HELIX_NUMBER_COLORS[1], 2: HELIX_NUMBER_COLORS[2], 
        3: HELIX_NUMBER_COLORS[3], 4: HELIX_NUMBER_COLORS[4],
        5: HELIX_NUMBER_COLORS[5], 6: HELIX_NUMBER_COLORS[6],
        7: HELIX_NUMBER_COLORS[7], 0: '#D3D3D3'
    }
    
    grn_highlights = {}
    helix_grns = ["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"]
    
    for grn_pos in helix_grns:
        highlight_coords = []
        highlight_colors = []
        
        for struct_id, data in oriented_structures.items():
            df = data.get('dataframe')
            if df is not None:
                coords = data['coords']
                df_aligned = df.copy()
                df_aligned[['x', 'y', 'z']] = coords
                
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
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    
    print(f"\n=== Generating {total_frames} front-view frames ===")
    print(f"Resolution: {width}x{height}")
    print(f"Duration: {total_frames/fps:.1f} seconds")
    
    successful_frames = 0
    
    # Generate frames
    for frame_num in tqdm(range(total_frames), desc="Creating frames"):
        angle = (frame_num / total_frames) * 360
        
        # Cycle through helices
        helix_idx = int((angle / 360) * 7) % 7
        target_grn = helix_grns[helix_idx]
        
        try:
            # Get highlight data
            highlight_data = grn_highlights.get(target_grn, {'coords': [], 'colors': []})
            
            # Create front-view frame
            fig = create_frontview_frame(
                all_coords_3d=all_coords_3d,
                highlight_coords_3d=highlight_data['coords'],
                highlight_colors=highlight_data['colors'],
                camera_angle=angle,
                target_grn=target_grn,
                figsize=(fig_width, fig_height),
                dpi=dpi
            )
            
            # Save as PNG
            frame_file = output_dir / f"frame_{frame_num:05d}.png"
            fig.savefig(frame_file, dpi=dpi, bbox_inches='tight',
                       facecolor='white', edgecolor='none')
            plt.close(fig)
            
            successful_frames += 1
            
        except Exception as e:
            if frame_num < 5:
                print(f"\nError on frame {frame_num}: {e}")
            plt.close('all')
            continue
    
    print(f"\n✅ Successfully created {successful_frames}/{total_frames} frames")
    
    if successful_frames == 0:
        print("❌ No frames created!")
        return
    
    # Create video
    print("\n=== Creating video with ffmpeg ===")
    video_file = output_dir / "mogrn_frontview_rotation.mp4"
    
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
        print(f"\n✅ Front-view animation complete!")
        print(f"📹 Video: {video_file}")
        print(f"🎞️ Frames: {successful_frames}")
        print(f"⏱️ Duration: {total_frames/fps:.1f} seconds")
        print(f"👁️ View: 2D front projection showing rotation around Z-axis")
    except Exception as e:
        print(f"\n❌ Video creation failed: {e}")
        print(f"📁 PNG frames available in: {output_dir}")


if __name__ == "__main__":
    create_frontview_animation()