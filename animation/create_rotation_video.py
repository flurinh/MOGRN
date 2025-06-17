#!/usr/bin/env python3
"""
Create a 1080-frame rotation animation using matplotlib.
This bypasses Plotly's kaleido issues and directly generates PNG files.
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


def create_matplotlib_frame(
    aligned_structures,
    grn_df,
    target_grn,
    camera_angle=0,
    figsize=(19.2, 10.8),
    dpi=100
):
    """
    Create a single frame using matplotlib.
    
    Args:
        aligned_structures: Aligned structure data
        grn_df: GRN table
        target_grn: GRN position to highlight (e.g., "5.50")
        camera_angle: Camera rotation angle in degrees
        figsize: Figure size in inches (width, height)
        dpi: Dots per inch
        
    Returns:
        matplotlib figure
    """
    
    # Create figure
    fig = plt.figure(figsize=figsize, dpi=dpi, facecolor='white')
    ax = fig.add_subplot(111, projection='3d', facecolor='white')
    
    # Helix colors
    helix_colors = {
        1: HELIX_NUMBER_COLORS[1], 2: HELIX_NUMBER_COLORS[2], 
        3: HELIX_NUMBER_COLORS[3], 4: HELIX_NUMBER_COLORS[4],
        5: HELIX_NUMBER_COLORS[5], 6: HELIX_NUMBER_COLORS[6],
        7: HELIX_NUMBER_COLORS[7], 0: '#D3D3D3'
    }
    
    # Plot all structures with low opacity
    for struct_id, data in aligned_structures.items():
        coords = data['coords']
        df = data.get('dataframe')
        
        if df is not None:
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = coords
            
            # Plot each helix as a line
            for helix_num in range(8):
                helix_data = df_aligned[df_aligned['helix_num'] == helix_num]
                if len(helix_data) > 0:
                    # Sort by residue number for proper line connectivity
                    helix_data = helix_data.sort_values('auth_seq_id')
                    helix_coords = helix_data[['x', 'y', 'z']].values
                    
                    if len(helix_coords) > 1:
                        color = helix_colors.get(helix_num, '#D3D3D3')
                        ax.plot(helix_coords[:, 0], helix_coords[:, 1], helix_coords[:, 2],
                               color=color, alpha=0.1, linewidth=0.5)
    
    # Collect and plot highlighted residues
    highlight_coords = []
    highlight_colors = []
    
    for struct_id, data in aligned_structures.items():
        df = data.get('dataframe')
        if df is not None:
            coords = data['coords']
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = coords
            
            # Find residues with target GRN
            grn_matches = df_aligned[df_aligned['grn'] == target_grn]
            
            for _, row in grn_matches.iterrows():
                highlight_coords.append([row['x'], row['y'], row['z']])
                helix_num = row['helix_num']
                highlight_colors.append(helix_colors.get(helix_num, '#D3D3D3'))
    
    # Plot highlights
    if highlight_coords:
        coords_array = np.array(highlight_coords)
        ax.scatter(coords_array[:, 0], coords_array[:, 1], coords_array[:, 2], 
                  c=highlight_colors, s=60, alpha=1.0, edgecolors='white', linewidth=1.5)
    
    # Set camera angle
    ax.view_init(elev=20, azim=camera_angle)
    
    # Clean appearance
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_zlabel('')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    
    # Hide axis panes
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.line.set_visible(False)
    ax.yaxis.line.set_visible(False)
    ax.zaxis.line.set_visible(False)
    
    # Set title
    helix_num = int(target_grn.split('.')[0])
    ax.set_title(f'MOGRN Helix {helix_num} Position 50 | Rotation: {camera_angle:.1f}°',
                fontsize=20, pad=20, fontweight='bold')
    
    # Set consistent axis limits
    margin = 10
    ax.set_xlim([-40-margin, 40+margin])
    ax.set_ylim([-40-margin, 40+margin])
    ax.set_zlim([-20-margin, 20+margin])
    
    # Make equal aspect ratio
    ax.set_box_aspect([1,1,0.5])
    
    plt.tight_layout()
    
    return fig


def calculate_helix_centers(aligned_structures):
    """Calculate the center position of each helix."""
    helix_centers = {}
    
    for helix_num in range(1, 8):  # Helices 1-7
        helix_coords = []
        
        for struct_id, data in aligned_structures.items():
            df = data.get('dataframe')
            if df is not None:
                coords = data['coords']
                df_aligned = df.copy()
                df_aligned[['x', 'y', 'z']] = coords
                
                # Get coordinates for this helix
                helix_data = df_aligned[df_aligned['helix_num'] == helix_num]
                if len(helix_data) > 0:
                    helix_coords.extend(helix_data[['x', 'y', 'z']].values.tolist())
        
        if helix_coords:
            helix_centers[helix_num] = np.mean(np.array(helix_coords), axis=0)
    
    return helix_centers


def determine_visible_helix(camera_angle, helix_centers):
    """Determine which helix is most visible from the camera angle."""
    # Camera position (rotating around Z-axis)
    angle_rad = np.radians(camera_angle)
    camera_x = np.cos(angle_rad)
    camera_y = np.sin(angle_rad)
    camera_dir = np.array([camera_x, camera_y, 0])
    
    # Find helix with center most aligned with camera direction
    best_helix = 1
    max_dot = -1
    
    for helix_num, center in helix_centers.items():
        # Project center to XY plane and normalize
        center_xy = np.array([center[0], center[1], 0])
        if np.linalg.norm(center_xy) > 0:
            center_xy = center_xy / np.linalg.norm(center_xy)
            
            # Dot product gives alignment with camera
            dot = np.dot(camera_dir, center_xy)
            if dot > max_dot:
                max_dot = dot
                best_helix = helix_num
    
    return f"{best_helix}.50"


def create_full_rotation_animation():
    """Create the full 1080-frame rotation animation."""
    
    print("=== MOGRN Full Rotation Animation (1080 frames) ===\n")
    
    # Parameters
    total_frames = 1080
    fps = 60  # 60 fps = 18 second video
    output_dir = Path("rotation_1080_frames")
    width = 1920
    height = 1080
    dpi = 100
    
    # Calculate figure size in inches
    fig_width = width / dpi
    fig_height = height / dpi
    
    print("=== Loading Data ===")
    cache_data = load_rmsd_cache()
    alignment_paths = cache_data.get('alignment_paths', {})
    processed_structures = load_processed_structures()
    grn_df = load_grn_table()
    
    print("=== Processing Structures ===")
    structures = extract_ca_coordinates_with_grn(processed_structures, grn_df, chain_id='A', use_helix_only=True)
    aligned_structures = apply_alignment_transformations(structures, alignment_paths, 'MerMAID1_model_0')
    oriented_structures = apply_membrane_orientation(aligned_structures, 'MerMAID1_model_0')
    
    print("=== Calculating Helix Centers ===")
    helix_centers = calculate_helix_centers(oriented_structures)
    print(f"Found {len(helix_centers)} helix centers")
    
    # Create output directory
    output_dir.mkdir(exist_ok=True)
    
    print(f"\n=== Generating {total_frames} frames ===")
    print(f"Resolution: {width}x{height}")
    print(f"Duration: {total_frames/fps:.1f} seconds at {fps} fps")
    print(f"Output directory: {output_dir}\n")
    
    # Generate frames
    for frame_num in tqdm(range(total_frames), desc="Creating frames"):
        # Calculate rotation angle (0-360 degrees over all frames)
        angle = (frame_num / total_frames) * 360
        
        # Determine which helix to highlight based on visibility
        target_grn = determine_visible_helix(angle, helix_centers)
        
        try:
            # Create frame
            fig = create_matplotlib_frame(
                aligned_structures=oriented_structures,
                grn_df=grn_df,
                target_grn=target_grn,
                camera_angle=angle,
                figsize=(fig_width, fig_height),
                dpi=dpi
            )
            
            # Save as PNG
            frame_file = output_dir / f"frame_{frame_num:05d}.png"
            fig.savefig(frame_file, dpi=dpi, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            plt.close(fig)  # Free memory
            
        except Exception as e:
            print(f"\nError on frame {frame_num}: {e}")
            plt.close('all')  # Clean up any open figures
            continue
    
    print(f"\n✅ Frame generation complete!")
    
    # Create video with ffmpeg
    print("\n=== Creating video with ffmpeg ===")
    video_file = output_dir / "mogrn_rotation_1080.mp4"
    
    cmd = [
        'ffmpeg', '-y',
        '-framerate', str(fps),
        '-i', str(output_dir / 'frame_%05d.png'),
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p',
        '-crf', '18',  # High quality
        '-preset', 'slow',  # Better compression
        str(video_file)
    ]
    
    try:
        print(f"Running: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)
        print(f"\n✅ Video created successfully!")
        print(f"   File: {video_file}")
        print(f"   Duration: {total_frames/fps:.1f} seconds")
        print(f"   Frame rate: {fps} fps")
        print(f"   Resolution: {width}x{height}")
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"\n❌ ffmpeg error: {e}")
        print("\nMake sure ffmpeg is installed and in your PATH")
        print("You can still create the video manually with:")
        print(f"ffmpeg -framerate {fps} -i {output_dir}/frame_%05d.png -c:v libx264 -pix_fmt yuv420p output.mp4")


if __name__ == "__main__":
    create_full_rotation_animation()