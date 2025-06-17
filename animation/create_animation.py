"""
Animation generation for MOGRN opsin structure visualization.

This module creates rotating animations of aligned opsin structures with 
GRN position highlights, outputting frames for ffmpeg video generation.
"""

import os
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import tempfile
import subprocess
from pathlib import Path
from tqdm import tqdm
import time

# Import our visualization functions
from visualize_alignment_grn import (
    load_rmsd_cache, 
    load_processed_structures, 
    load_grn_table,
    extract_ca_coordinates_with_grn,
    apply_alignment_transformations,
    apply_membrane_orientation,
    create_interactive_opsin_visualization
)


def create_static_frame(
    aligned_structures, 
    grn_df, 
    property_data,
    grn_position,
    camera_angle=0,
    color_mode='helix',
    width=1920,
    height=1080,
    show_membrane=True,
    membrane_opacity=0.03,
    title_suffix=""
):
    """
    Create a single static frame for animation.
    
    Args:
        aligned_structures: Dictionary of aligned structure data
        grn_df: GRN position table
        property_data: Property data for coloring
        grn_position: GRN position to highlight (e.g., 1.50, 2.50)
        camera_angle: Rotation angle around Z-axis in degrees
        color_mode: 'helix' or 'property'
        width: Frame width in pixels
        height: Frame height in pixels
        show_membrane: Whether to show membrane volume
        membrane_opacity: Membrane transparency
        title_suffix: Additional text for title
        
    Returns:
        plotly.graph_objects.Figure: Static frame figure
    """
    
    # Create the base visualization
    fig = create_interactive_opsin_visualization(
        aligned_structures=aligned_structures,
        grn_df=grn_df,
        property_data=property_data,
        title=f"MOGRN: GRN {grn_position} Highlighted{title_suffix}",
        width=width,
        height=height,
        membrane_opacity=membrane_opacity,
        show_membrane=show_membrane
    )
    
    # Calculate camera position for rotation around Z-axis
    # Convert angle to radians
    angle_rad = np.radians(camera_angle)
    
    # Camera distance from center
    camera_distance = 2.5
    
    # Calculate camera position (rotating around Z-axis)
    camera_x = camera_distance * np.cos(angle_rad)
    camera_y = camera_distance * np.sin(angle_rad)
    camera_z = 0.8  # Keep some elevation
    
    # Update camera for this frame
    fig.update_layout(
        scene_camera=dict(
            eye=dict(x=camera_x, y=camera_y, z=camera_z),
            center=dict(x=0, y=0, z=0),
            up=dict(x=0, y=0, z=1)
        ),
        # Remove interactive elements for static export
        sliders=[],
        updatemenus=[],
        showlegend=False,  # Hide legend for cleaner animation
        margin=dict(l=0, r=0, t=60, b=0)
    )
    
    # Manually set visibility for specific GRN position and color mode
    # This is a simplified version - you'd need to implement the logic
    # to show only the desired GRN position and color mode
    
    # For now, we'll create a simplified version focusing on the rotation
    # You can extend this to properly handle GRN highlighting
    
    return fig


def generate_animation_frames(
    aligned_structures,
    grn_df,
    property_data=None,
    grn_positions=None,
    color_mode='helix',
    rotation_steps=36,  # 10-degree increments for full rotation
    output_dir="animation_frames",
    width=1920,
    height=1080,
    show_membrane=True,
    membrane_opacity=0.03
):
    """
    Generate all frames for the animation.
    
    Args:
        aligned_structures: Dictionary of aligned structure data
        grn_df: GRN position table
        property_data: Property data for coloring
        grn_positions: List of GRN positions to highlight (e.g., [1.50, 2.50, 3.50])
        color_mode: 'helix' or 'property'
        rotation_steps: Number of rotation steps for 360-degree rotation
        output_dir: Directory to save frame images
        width: Frame width in pixels
        height: Frame height in pixels
        show_membrane: Whether to show membrane volume
        membrane_opacity: Membrane transparency
        
    Returns:
        List of frame file paths
    """
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Default GRN positions if not provided
    if grn_positions is None:
        # Select some key GRN positions
        grn_positions = [1.50, 2.50, 3.50, 4.50, 5.50, 6.50, 7.50]
    
    frame_files = []
    frame_number = 0
    
    print(f"Generating animation frames...")
    print(f"GRN positions: {grn_positions}")
    print(f"Rotation steps: {rotation_steps}")
    print(f"Total frames: {len(grn_positions) * rotation_steps}")
    
    # Generate frames for each GRN position
    for grn_idx, grn_pos in enumerate(grn_positions):
        print(f"\nGenerating frames for GRN {grn_pos} ({grn_idx + 1}/{len(grn_positions)})")
        
        # Generate rotation frames for this GRN position
        for step in tqdm(range(rotation_steps), desc=f"GRN {grn_pos}"):
            # Calculate rotation angle (360 degrees over rotation_steps)
            angle = (step / rotation_steps) * 360
            
            # Create frame
            fig = create_static_frame(
                aligned_structures=aligned_structures,
                grn_df=grn_df,
                property_data=property_data,
                grn_position=grn_pos,
                camera_angle=angle,
                color_mode=color_mode,
                width=width,
                height=height,
                show_membrane=show_membrane,
                membrane_opacity=membrane_opacity,
                title_suffix=f" | {color_mode.title()} Mode | Angle: {angle:.0f}°"
            )
            
            # Save frame
            frame_file = output_path / f"frame_{frame_number:04d}.png"
            
            try:
                fig.write_image(
                    str(frame_file),
                    width=width,
                    height=height,
                    scale=1,
                    format="png"
                )
                frame_files.append(str(frame_file))
                frame_number += 1
                
            except Exception as e:
                print(f"Error saving frame {frame_number}: {e}")
                continue
    
    print(f"\nGenerated {len(frame_files)} frames in {output_dir}")
    return frame_files


def create_video_with_ffmpeg(
    frame_files,
    output_video="opsin_animation.mp4",
    fps=24,
    quality="high",
    cleanup_frames=False
):
    """
    Create video from frame images using ffmpeg.
    
    Args:
        frame_files: List of frame file paths
        output_video: Output video filename
        fps: Frames per second
        quality: Video quality ('low', 'medium', 'high', 'ultra')
        cleanup_frames: Whether to delete frame files after video creation
        
    Returns:
        Path to output video file
    """
    
    if not frame_files:
        print("No frames to process!")
        return None
    
    # Quality settings
    quality_settings = {
        'low': ['-crf', '28', '-preset', 'fast'],
        'medium': ['-crf', '23', '-preset', 'medium'],
        'high': ['-crf', '18', '-preset', 'slow'],
        'ultra': ['-crf', '15', '-preset', 'slower']
    }
    
    # Get the directory containing frames
    frame_dir = Path(frame_files[0]).parent
    frame_pattern = str(frame_dir / "frame_%04d.png")
    
    # Build ffmpeg command
    cmd = [
        'ffmpeg',
        '-y',  # Overwrite output file
        '-framerate', str(fps),
        '-i', frame_pattern,
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p',  # For compatibility
    ]
    
    # Add quality settings
    cmd.extend(quality_settings.get(quality, quality_settings['high']))
    
    # Add output file
    cmd.append(output_video)
    
    print(f"Creating video with ffmpeg...")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        # Run ffmpeg
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        print(f"Video created successfully: {output_video}")
        
        # Print video info
        if os.path.exists(output_video):
            file_size = os.path.getsize(output_video) / (1024 * 1024)  # MB
            print(f"Video size: {file_size:.1f} MB")
        
        # Cleanup frames if requested
        if cleanup_frames:
            print("Cleaning up frame files...")
            for frame_file in frame_files:
                try:
                    os.remove(frame_file)
                except OSError:
                    pass
            
            # Remove frame directory if empty
            try:
                os.rmdir(frame_dir)
            except OSError:
                pass
        
        return output_video
        
    except subprocess.CalledProcessError as e:
        print(f"ffmpeg error: {e}")
        print(f"stdout: {e.stdout}")
        print(f"stderr: {e.stderr}")
        return None
    except FileNotFoundError:
        print("ffmpeg not found! Please install ffmpeg:")
        print("  Ubuntu/Debian: sudo apt install ffmpeg")
        print("  macOS: brew install ffmpeg")
        print("  Windows: Download from https://ffmpeg.org/download.html")
        return None


def create_opsin_animation(
    grn_positions=None,
    color_mode='helix',
    rotation_steps=36,
    fps=24,
    output_video="opsin_animation.mp4",
    quality="high",
    width=1920,
    height=1080,
    cleanup_frames=True,
    reference_id='MerMAID1_model_0'
):
    """
    Complete pipeline to create opsin structure animation.
    
    Args:
        grn_positions: List of GRN positions to highlight
        color_mode: 'helix' or 'property'
        rotation_steps: Number of rotation steps (36 = 10° increments)
        fps: Video framerate
        output_video: Output video filename
        quality: Video quality setting
        width: Frame width in pixels
        height: Frame height in pixels
        cleanup_frames: Delete frame files after video creation
        reference_id: Reference structure for alignment
        
    Returns:
        Path to created video file
    """
    
    print("=== MOGRN Animation Generator ===")
    print(f"Creating {color_mode} mode animation with {rotation_steps} rotation steps")
    
    # Load all required data
    print("\n=== Loading Data ===")
    cache_data = load_rmsd_cache()
    alignment_paths = cache_data.get('alignment_paths', {})
    processed_structures = load_processed_structures()
    grn_df = load_grn_table()
    
    # Load property data if needed
    property_data = None
    if color_mode == 'property':
        from src.data_processing import load_opsin_property_data
        from pathlib import Path
        property_file = Path('property/mo_exp.csv')
        
        if property_file.exists():
            try:
                property_result = load_opsin_property_data(property_file, processed_structures)
                if property_result and 'properties' in property_result:
                    property_data = property_result['properties']
                    print(f"Loaded property data for {len(property_data)} structures")
            except Exception as e:
                print(f"Failed to load property data: {e}")
    
    # Process structures
    print("\n=== Processing Structures ===")
    structures = extract_ca_coordinates_with_grn(
        processed_structures, grn_df, chain_id='A', use_helix_only=True
    )
    
    aligned_structures = apply_alignment_transformations(
        structures, alignment_paths, reference_id
    )
    
    oriented_structures = apply_membrane_orientation(aligned_structures, reference_id)
    
    print(f"Processed {len(oriented_structures)} structures")
    
    # Generate animation frames
    print("\n=== Generating Animation Frames ===")
    frame_files = generate_animation_frames(
        aligned_structures=oriented_structures,
        grn_df=grn_df,
        property_data=property_data,
        grn_positions=grn_positions,
        color_mode=color_mode,
        rotation_steps=rotation_steps,
        width=width,
        height=height
    )
    
    if not frame_files:
        print("No frames generated!")
        return None
    
    # Create video
    print("\n=== Creating Video ===")
    video_file = create_video_with_ffmpeg(
        frame_files=frame_files,
        output_video=output_video,
        fps=fps,
        quality=quality,
        cleanup_frames=cleanup_frames
    )
    
    if video_file:
        print(f"\n=== Animation Complete ===")
        print(f"Video saved to: {video_file}")
        print(f"Duration: ~{len(frame_files)/fps:.1f} seconds")
        print(f"Resolution: {width}x{height}")
        print(f"Framerate: {fps} fps")
    
    return video_file


def create_comparison_animation(
    grn_positions=None,
    rotation_steps=36,
    fps=24,
    output_video="opsin_comparison_animation.mp4",
    quality="high",
    width=1920,
    height=1080,
    cleanup_frames=True
):
    """
    Create side-by-side comparison animation of helix vs property coloring.
    
    Args:
        grn_positions: List of GRN positions to highlight
        rotation_steps: Number of rotation steps
        fps: Video framerate
        output_video: Output video filename
        quality: Video quality setting
        width: Frame width in pixels
        height: Frame height in pixels
        cleanup_frames: Delete frame files after video creation
        
    Returns:
        Path to created video file
    """
    
    print("=== MOGRN Comparison Animation Generator ===")
    print("Creating side-by-side helix vs property comparison")
    
    # This would create frames with subplot showing both modes
    # Implementation would be similar but create subplots
    # Left: helix mode, Right: property mode
    
    # For now, create separate animations
    helix_video = create_opsin_animation(
        grn_positions=grn_positions,
        color_mode='helix',
        rotation_steps=rotation_steps,
        fps=fps,
        output_video="helix_animation.mp4",
        quality=quality,
        width=width//2,
        height=height,
        cleanup_frames=cleanup_frames
    )
    
    property_video = create_opsin_animation(
        grn_positions=grn_positions,
        color_mode='property',
        rotation_steps=rotation_steps,
        fps=fps,
        output_video="property_animation.mp4",
        quality=quality,
        width=width//2,
        height=height,
        cleanup_frames=cleanup_frames
    )
    
    # Combine videos side by side using ffmpeg
    if helix_video and property_video:
        cmd = [
            'ffmpeg', '-y',
            '-i', helix_video,
            '-i', property_video,
            '-filter_complex', '[0:v][1:v]hstack=inputs=2[v]',
            '-map', '[v]',
            '-c:v', 'libx264',
            '-crf', '18',
            '-pix_fmt', 'yuv420p',
            output_video
        ]
        
        try:
            subprocess.run(cmd, check=True)
            print(f"Comparison video created: {output_video}")
            
            # Cleanup individual videos
            if cleanup_frames:
                os.remove(helix_video)
                os.remove(property_video)
            
            return output_video
            
        except subprocess.CalledProcessError as e:
            print(f"Error creating comparison video: {e}")
            return None
    
    return None


if __name__ == "__main__":
    # Example usage
    
    # Create standard animation
    video = create_opsin_animation(
        grn_positions=[1.50, 2.50, 3.50, 4.50, 5.50, 6.50, 7.50],
        color_mode='helix',
        rotation_steps=36,  # 10-degree increments
        fps=24,
        output_video="opsin_helix_animation.mp4",
        quality="high"
    )
    
    # Create property mode animation
    property_video = create_opsin_animation(
        grn_positions=[1.50, 2.50, 3.50, 4.50, 5.50, 6.50, 7.50],
        color_mode='property',
        rotation_steps=36,
        fps=24,
        output_video="opsin_property_animation.mp4",
        quality="high"
    )
    
    print("Animation generation complete!")