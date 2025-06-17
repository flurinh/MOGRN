#!/usr/bin/env python3
"""
Create rotation animation using matplotlib for rendering.
This bypasses Plotly's image export issues.
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
    width=19.2,
    height=10.8,
    dpi=100,
    show_membrane=False
):
    """
    Create a frame using matplotlib instead of plotly.
    
    Args:
        aligned_structures: Aligned structure data
        grn_df: GRN table
        target_grn: GRN position to highlight
        camera_angle: Camera rotation angle
        width: Figure width in inches
        height: Figure height in inches
        dpi: Dots per inch
        
    Returns:
        matplotlib figure
    """
    
    # Create figure and 3D axis
    fig = plt.figure(figsize=(width, height), dpi=dpi)
    ax = fig.add_subplot(111, projection='3d')
    
    # Set background to white
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    
    # Color scheme
    helix_colors = {
        1: HELIX_NUMBER_COLORS[1], 2: HELIX_NUMBER_COLORS[2], 
        3: HELIX_NUMBER_COLORS[3], 4: HELIX_NUMBER_COLORS[4],
        5: HELIX_NUMBER_COLORS[5], 6: HELIX_NUMBER_COLORS[6],
        7: HELIX_NUMBER_COLORS[7], 0: '#D3D3D3'
    }
    
    # Plot structures with low opacity
    for struct_id, data in aligned_structures.items():
        coords = data['coords']
        df = data.get('dataframe')
        
        if df is not None:
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = coords
            
            # Plot each helix
            for helix_num in range(8):
                helix_data = df_aligned[df_aligned['helix_num'] == helix_num]
                if len(helix_data) > 0:
                    helix_coords = helix_data[['x', 'y', 'z']].values
                    color = helix_colors.get(helix_num, '#D3D3D3')
                    
                    # Plot as lines with low alpha
                    ax.plot(helix_coords[:, 0], helix_coords[:, 1], helix_coords[:, 2],
                           color=color, alpha=0.1, linewidth=0.5)
    
    # Highlight target GRN positions
    highlight_x, highlight_y, highlight_z = [], [], []
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
                highlight_x.append(row['x'])
                highlight_y.append(row['y'])
                highlight_z.append(row['z'])
                helix_num = row['helix_num']
                highlight_colors.append(helix_colors.get(helix_num, '#D3D3D3'))
    
    # Plot highlights
    if highlight_x:
        ax.scatter(highlight_x, highlight_y, highlight_z, 
                  c=highlight_colors, s=50, alpha=1.0, edgecolors='white', linewidth=1)
    
    # Set camera position
    angle_rad = np.radians(camera_angle)
    distance = 100
    ax.view_init(elev=20, azim=camera_angle)
    
    # Remove axes for clean look
    ax.set_xlabel('')
    ax.set_ylabel('')
    ax.set_zlabel('')
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    
    # Remove axis lines and panes
    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False
    ax.xaxis.line.set_visible(False)
    ax.yaxis.line.set_visible(False)
    ax.zaxis.line.set_visible(False)
    
    # Set title
    ax.set_title(f'MOGRN: Helix {target_grn[0]} Position 50 Highlighted | Rotation: {camera_angle:.0f}°',
                fontsize=16, pad=20)
    
    # Set axis limits
    all_coords = []
    for data in aligned_structures.values():
        all_coords.extend(data['coords'])
    all_coords = np.array(all_coords[:1000])  # Sample for speed
    
    margin = 5
    ax.set_xlim([all_coords[:, 0].min() - margin, all_coords[:, 0].max() + margin])
    ax.set_ylim([all_coords[:, 1].min() - margin, all_coords[:, 1].max() + margin])
    ax.set_zlim([all_coords[:, 2].min() - margin, all_coords[:, 2].max() + margin])
    
    plt.tight_layout()
    
    return fig


def create_rotation_animation_matplotlib(
    total_frames=360,  # 360 frames = 1 degree per frame
    fps=30,
    output_dir="matplotlib_rotation_frames",
    width=1920,
    height=1080,
    reference_id='MerMAID1_model_0'
):
    """
    Create rotation animation using matplotlib backend.
    """
    
    print("=== Loading Data ===")
    cache_data = load_rmsd_cache()
    alignment_paths = cache_data.get('alignment_paths', {})
    processed_structures = load_processed_structures()
    grn_df = load_grn_table()
    
    print("=== Processing Structures ===")
    structures = extract_ca_coordinates_with_grn(processed_structures, grn_df, chain_id='A', use_helix_only=True)
    aligned_structures = apply_alignment_transformations(structures, alignment_paths, reference_id)
    oriented_structures = apply_membrane_orientation(aligned_structures, reference_id)
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Calculate figure size and DPI to match desired resolution
    dpi = 100
    fig_width = width / dpi
    fig_height = height / dpi
    
    print(f"\n=== Generating {total_frames} frames ===")
    print(f"Resolution: {width}x{height} pixels")
    print(f"Duration: {total_frames/fps:.1f} seconds at {fps} fps")
    
    # Define helix positions to cycle through
    helix_grns = ["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"]
    
    frame_files = []
    
    for frame_num in tqdm(range(total_frames), desc="Generating frames"):
        angle = (frame_num / total_frames) * 360
        
        # Determine which helix to highlight based on angle
        # Each helix gets ~51 degrees of the rotation
        helix_idx = int((angle / 360) * 7) % 7
        target_grn = helix_grns[helix_idx]
        
        try:
            # Create frame using matplotlib
            fig = create_matplotlib_frame(
                aligned_structures=oriented_structures,
                grn_df=grn_df,
                target_grn=target_grn,
                camera_angle=angle,
                width=fig_width,
                height=fig_height,
                dpi=dpi,
                show_membrane=False
            )
            
            # Save frame as PNG
            frame_file = output_path / f"frame_{frame_num:05d}.png"
            fig.savefig(frame_file, dpi=dpi, bbox_inches='tight', 
                       facecolor='white', edgecolor='none')
            plt.close(fig)  # Free memory
            
            frame_files.append(str(frame_file))
            
        except Exception as e:
            print(f"\nError on frame {frame_num}: {e}")
            continue
    
    print(f"\n✅ Generated {len(frame_files)} PNG frames")
    
    # Create video with ffmpeg
    if frame_files:
        print("\n=== Creating video with ffmpeg ===")
        video_file = output_path / "opsin_rotation_matplotlib.mp4"
        
        cmd = [
            'ffmpeg', '-y',
            '-framerate', str(fps),
            '-i', str(output_path / 'frame_%05d.png'),
            '-c:v', 'libx264',
            '-pix_fmt', 'yuv420p',
            '-crf', '18',  # High quality
            str(video_file)
        ]
        
        try:
            subprocess.run(cmd, check=True)
            print(f"✅ Video created: {video_file}")
            print(f"Duration: {total_frames/fps:.1f} seconds")
            return video_file
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            print(f"❌ ffmpeg error: {e}")
            print("Make sure ffmpeg is installed and in your PATH")
            print("You can still manually create the video with:")
            print(f"ffmpeg -framerate {fps} -i {output_path}/frame_%05d.png -c:v libx264 -pix_fmt yuv420p output.mp4")
    
    return None


def main():
    """Create the rotation animation using matplotlib."""
    
    print("Creating rotation animation with matplotlib backend...")
    print("This bypasses Plotly's kaleido issues.\n")
    
    # Test with fewer frames first
    video = create_rotation_animation_matplotlib(
        total_frames=72,  # 72 frames = 5-degree increments
        fps=24,  # 3 seconds total
        output_dir="matplotlib_test_rotation",
        width=1920,
        height=1080
    )
    
    if video:
        print(f"\n✅ Test animation complete: {video}")
        print("\nTo create the full 1080-frame animation:")
        print("create_rotation_animation_matplotlib(total_frames=1080, fps=60)")
    else:
        print("\n⚠️ Video creation failed, but PNG frames should be available")


if __name__ == "__main__":
    main()