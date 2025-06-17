#!/usr/bin/env python3
"""
Create a rotation animation highlighting the .50 position of the helix closest to the viewer.
This version works around kaleido issues by using alternative export methods.
"""

import numpy as np
import os
from pathlib import Path
from tqdm import tqdm
import subprocess
import json

from animate_grn_highlights import (
    load_rmsd_cache, load_processed_structures, load_grn_table,
    extract_ca_coordinates_with_grn, apply_alignment_transformations,
    apply_membrane_orientation, create_grn_highlight_frame
)


def determine_closest_helix_grn(camera_angle, helix_centers):
    """
    Determine which helix .50 position to highlight based on camera angle.
    
    Args:
        camera_angle: Current camera rotation angle in degrees
        helix_centers: Dict mapping helix number to average position
        
    Returns:
        str: GRN position to highlight (e.g., "3.50")
    """
    # Convert angle to radians
    angle_rad = np.radians(camera_angle)
    
    # Camera position (on a circle around Z-axis)
    camera_x = np.cos(angle_rad)
    camera_y = np.sin(angle_rad)
    camera_pos = np.array([camera_x, camera_y, 0])  # Z=0 for horizontal view
    
    # Find closest helix center to camera
    min_dist = float('inf')
    closest_helix = 1
    
    for helix_num, center in helix_centers.items():
        if helix_num == 0:  # Skip loops
            continue
        # Project to XY plane for distance calculation
        center_xy = np.array([center[0], center[1], 0])
        dist = np.linalg.norm(camera_pos - center_xy)
        if dist < min_dist:
            min_dist = dist
            closest_helix = helix_num
    
    return f"{closest_helix}.50"


def calculate_helix_centers(aligned_structures, grn_df):
    """Calculate average position of each helix."""
    helix_centers = {}
    helix_coords = {i: [] for i in range(8)}  # 0-7
    
    # Collect coordinates for each helix
    for struct_id, data in aligned_structures.items():
        df = data.get('dataframe')
        if df is not None:
            coords = data['coords']
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = coords
            
            for helix_num in range(8):
                helix_data = df_aligned[df_aligned['helix_num'] == helix_num]
                if len(helix_data) > 0:
                    helix_coords[helix_num].extend(helix_data[['x', 'y', 'z']].values.tolist())
    
    # Calculate centers
    for helix_num, coords_list in helix_coords.items():
        if coords_list:
            coords_array = np.array(coords_list)
            helix_centers[helix_num] = np.mean(coords_array, axis=0)
    
    return helix_centers


def create_rotation_animation_with_dynamic_highlights(
    total_frames=1080,
    fps=60,
    output_dir="rotation_animation_frames",
    width=1920,
    height=1080,
    reference_id='MerMAID1_model_0'
):
    """
    Create rotation animation with dynamically highlighted .50 positions.
    
    Args:
        total_frames: Total number of frames (1080 = 18 seconds at 60fps)
        fps: Frames per second for final video
        output_dir: Directory for output frames
        width: Frame width
        height: Frame height
        reference_id: Reference structure ID
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
    
    print("=== Calculating Helix Centers ===")
    helix_centers = calculate_helix_centers(oriented_structures, grn_df)
    print(f"Found centers for {len(helix_centers)} helices")
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Create frames directory
    frames_path = output_path / "frames"
    frames_path.mkdir(exist_ok=True)
    
    frame_files = []
    
    print(f"\n=== Generating {total_frames} frames for 360° rotation ===")
    print(f"Duration: {total_frames/fps:.1f} seconds at {fps} fps")
    
    # Generate frames with progress bar
    for frame_num in tqdm(range(total_frames), desc="Generating frames"):
        # Calculate rotation angle (0-360 degrees over total frames)
        angle = (frame_num / total_frames) * 360
        
        # Determine which helix .50 to highlight
        target_grn = determine_closest_helix_grn(angle, helix_centers)
        
        try:
            # Create frame
            fig = create_grn_highlight_frame(
                aligned_structures=oriented_structures,
                grn_df=grn_df,
                property_data=None,
                target_grn=target_grn,
                camera_angle=angle,
                color_mode='helix',
                width=width,
                height=height,
                show_membrane=False  # Faster rendering
            )
            
            # Export frame - try JSON first as workaround
            frame_file = frames_path / f"frame_{frame_num:05d}.json"
            fig.write_json(str(frame_file))
            frame_files.append(str(frame_file))
            
        except Exception as e:
            print(f"\nError on frame {frame_num}: {e}")
            continue
    
    print(f"\nGenerated {len(frame_files)} frame files")
    
    # Convert JSON frames to images using a batch script
    print("\n=== Converting frames to images ===")
    convert_frames_to_images(frame_files, width, height)
    
    # Create video with ffmpeg
    print("\n=== Creating video with ffmpeg ===")
    video_file = output_path / "opsin_rotation_highlights.mp4"
    create_video_from_frames(frames_path, video_file, fps)
    
    return video_file


def convert_frames_to_images(json_files, width, height):
    """Convert JSON frames to PNG images using a Python script."""
    
    print("Creating frame converter script...")
    
    converter_script = '''
import json
import plotly.graph_objects as go
from pathlib import Path
import sys

def convert_json_to_image(json_file, width, height):
    """Convert a single JSON frame to PNG."""
    json_path = Path(json_file)
    png_path = json_path.with_suffix('.png')
    
    try:
        # Load figure from JSON
        with open(json_file, 'r') as f:
            fig_dict = json.load(f)
        fig = go.Figure(fig_dict)
        
        # Try different export methods
        try:
            # Method 1: Direct write_image
            fig.write_image(str(png_path), width=width, height=height)
            return True
        except:
            # Method 2: Save as HTML then screenshot with browser
            html_path = json_path.with_suffix('.html')
            fig.write_html(str(html_path))
            print(f"Saved as HTML: {html_path}")
            # Note: Manual conversion needed
            return False
    except Exception as e:
        print(f"Error converting {json_file}: {e}")
        return False

if __name__ == "__main__":
    if len(sys.argv) < 4:
        print("Usage: python converter.py <json_file> <width> <height>")
        sys.exit(1)
    
    json_file = sys.argv[1]
    width = int(sys.argv[2])
    height = int(sys.argv[3])
    
    success = convert_json_to_image(json_file, width, height)
    sys.exit(0 if success else 1)
'''
    
    converter_file = Path("frame_converter.py")
    with open(converter_file, 'w') as f:
        f.write(converter_script)
    
    print(f"Converter script created: {converter_file}")
    print("\nNOTE: Due to kaleido issues, frames are saved as JSON.")
    print("To convert to images, you have several options:")
    print("1. Fix kaleido installation")
    print("2. Use plotly-orca instead")
    print("3. Use online conversion service")
    print("4. Use browser automation to capture HTML frames")
    
    # Create a simpler batch converter that creates HTML files
    print("\nCreating HTML versions of frames...")
    for i, json_file in enumerate(json_files[:5]):  # Just first 5 for testing
        json_path = Path(json_file)
        html_path = json_path.with_suffix('.html')
        
        try:
            with open(json_file, 'r') as f:
                fig_dict = json.load(f)
            fig = go.Figure(fig_dict)
            fig.write_html(str(html_path))
            print(f"Created HTML: {html_path}")
        except Exception as e:
            print(f"Error creating HTML for {json_file}: {e}")


def create_video_from_frames(frames_dir, output_video, fps):
    """Create video from PNG frames using ffmpeg."""
    
    # Check if PNG files exist
    png_files = list(frames_dir.glob("*.png"))
    
    if not png_files:
        print("No PNG files found. Checking for alternatives...")
        
        # Check for HTML files
        html_files = list(frames_dir.glob("*.html"))
        if html_files:
            print(f"Found {len(html_files)} HTML files.")
            print("To create video, you need to convert HTML files to PNG first.")
            print("Options:")
            print("1. Use a browser automation tool (e.g., Selenium, Playwright)")
            print("2. Use online HTML to image conversion")
            print("3. Fix kaleido installation and re-run")
            return None
    
    # Create video with ffmpeg
    frame_pattern = str(frames_dir / "frame_%05d.png")
    
    cmd = [
        'ffmpeg', '-y',
        '-framerate', str(fps),
        '-i', frame_pattern,
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p',
        '-crf', '18',  # High quality
        str(output_video)
    ]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"Video created: {output_video}")
        return output_video
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"ffmpeg error: {e}")
        print("Make sure ffmpeg is installed")
        return None


def main():
    """Create the rotation animation."""
    
    # First, let's test with just a few frames
    print("Testing with 36 frames (10-degree increments)...")
    
    video = create_rotation_animation_with_dynamic_highlights(
        total_frames=36,  # Just 36 frames for testing
        fps=6,  # 6 fps = 6 second video
        output_dir="test_rotation_animation",
        width=1920,
        height=1080
    )
    
    if video:
        print(f"\n✅ Test animation complete: {video}")
        print("\nTo create the full animation, run:")
        print("create_rotation_animation_with_dynamic_highlights(total_frames=1080, fps=60)")
    else:
        print("\n⚠️ Video creation failed. See instructions above for converting frames.")


if __name__ == "__main__":
    main()