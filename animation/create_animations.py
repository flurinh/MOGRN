#!/usr/bin/env python3
"""
Quick script to generate MOGRN animations.

This script provides simple functions to create different types of animations
from your MOGRN opsin structure data.
"""

from animate_grn_highlights import generate_grn_animation, create_video
from pathlib import Path


def create_quick_animation(
    animation_type="helix_tour",
    output_name=None,
    quality="high",
    fps=24
):
    """
    Create predefined animation types.
    
    Args:
        animation_type: Type of animation to create
            - "helix_tour": Tour through all 7 helices (GRN 1.50-7.50)
            - "binding_pocket": Focus on binding pocket region (GRN 3.50-6.50)  
            - "single_rotation": Single 360° rotation at GRN 5.50
            - "property_comparison": Property-based coloring tour
        output_name: Custom output filename (optional)
        quality: Video quality ("low", "medium", "high", "ultra")
        fps: Frames per second
        
    Returns:
        Path to created video file
    """
    
    animations = {
        "helix_tour": {
            "grn_positions": [1.50, 2.50, 3.50, 4.50, 5.50, 6.50, 7.50],
            "color_mode": "helix",
            "rotation_steps": 24,  # 15° increments
            "output": "opsin_helix_tour.mp4"
        },
        "binding_pocket": {
            "grn_positions": [3.35, 4.50, 5.50, 6.44, 6.55],  # Key binding pocket positions
            "color_mode": "helix", 
            "rotation_steps": 36,  # 10° increments for smoother rotation
            "output": "opsin_binding_pocket.mp4"
        },
        "single_rotation": {
            "grn_positions": [5.50],  # Center of retinal binding pocket
            "color_mode": "helix",
            "rotation_steps": 36,  # Full 360° rotation
            "output": "opsin_single_rotation.mp4"
        },
        "property_comparison": {
            "grn_positions": [1.50, 2.50, 3.50, 4.50, 5.50, 6.50, 7.50],
            "color_mode": "property",
            "rotation_steps": 24,
            "output": "opsin_property_tour.mp4"
        }
    }
    
    if animation_type not in animations:
        print(f"Unknown animation type: {animation_type}")
        print(f"Available types: {list(animations.keys())}")
        return None
    
    config = animations[animation_type]
    output_file = output_name or config["output"]
    
    print(f"Creating {animation_type} animation...")
    print(f"GRN positions: {config['grn_positions']}")
    print(f"Color mode: {config['color_mode']}")
    print(f"Total frames: {len(config['grn_positions']) * config['rotation_steps']}")
    
    # Generate frames
    frames = generate_grn_animation(
        grn_positions=config["grn_positions"],
        color_mode=config["color_mode"],
        rotation_steps=config["rotation_steps"],
        output_dir=f"{animation_type}_frames"
    )
    
    # Create video
    if frames:
        video = create_video(
            frames,
            output_video=output_file,
            fps=fps,
            quality=quality
        )
        
        if video:
            print(f"✅ Animation created: {video}")
            duration = len(frames) / fps
            print(f"Duration: {duration:.1f} seconds")
            print(f"Resolution: 1920x1080")
            print(f"Quality: {quality}")
            return video
    
    return None


def create_custom_animation(
    grn_positions,
    color_mode="helix", 
    rotation_steps=24,
    output_file="custom_animation.mp4",
    fps=24,
    quality="high"
):
    """
    Create custom animation with specific parameters.
    
    Args:
        grn_positions: List of GRN positions (e.g., [1.50, 2.50, 3.50])
        color_mode: "helix" or "property"
        rotation_steps: Number of rotation steps per GRN (24 = 15° increments)
        output_file: Output video filename
        fps: Frames per second
        quality: Video quality
        
    Returns:
        Path to created video file
    """
    
    print(f"Creating custom animation...")
    print(f"GRN positions: {grn_positions}")
    print(f"Rotation steps: {rotation_steps}")
    print(f"Color mode: {color_mode}")
    
    frames = generate_grn_animation(
        grn_positions=grn_positions,
        color_mode=color_mode,
        rotation_steps=rotation_steps,
        output_dir="custom_frames"
    )
    
    if frames:
        video = create_video(frames, output_file, fps, quality)
        if video:
            print(f"✅ Custom animation created: {video}")
            return video
    
    return None


def create_presentation_set():
    """Create a set of animations suitable for presentations."""
    
    print("=== Creating Presentation Animation Set ===")
    
    animations_to_create = [
        ("single_rotation", "Quick overview (1 rotation at binding pocket center)"),
        ("binding_pocket", "Binding pocket tour (key GRN positions)"),
        ("helix_tour", "Complete helix tour (all 7 helices)"),
        ("property_comparison", "Functional classification view")
    ]
    
    created_videos = []
    
    for anim_type, description in animations_to_create:
        print(f"\n📹 {description}")
        video = create_quick_animation(anim_type, quality="high", fps=24)
        if video:
            created_videos.append(video)
    
    print(f"\n🎬 Presentation set complete!")
    print(f"Created {len(created_videos)} videos:")
    for video in created_videos:
        print(f"  - {video}")
    
    return created_videos


if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        animation_type = sys.argv[1]
        video = create_quick_animation(animation_type)
    else:
        print("MOGRN Animation Creator")
        print("=" * 40)
        print()
        print("Quick animations:")
        print("  python create_animations.py helix_tour")
        print("  python create_animations.py binding_pocket") 
        print("  python create_animations.py single_rotation")
        print("  python create_animations.py property_comparison")
        print()
        print("Create full presentation set:")
        print("  python create_animations.py presentation")
        print()
        
        # Default: create a quick overview
        print("Creating default single rotation animation...")
        video = create_quick_animation("single_rotation")
        
        if video:
            print(f"\n✅ Default animation created: {video}")
            print("\nTo create other animations, run:")
            print("  python create_animations.py <animation_type>")
        else:
            print("\n❌ Animation failed. Check error messages above.")
            print("\nCommon issues:")
            print("  - ffmpeg not installed")
            print("  - Missing data files (run opsin_analysis_workflow.py first)")
            print("  - Insufficient disk space for frames")