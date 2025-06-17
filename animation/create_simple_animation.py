#!/usr/bin/env python3
"""
Simplified animation without membrane volume to avoid rendering issues.
"""

from animate_grn_highlights import generate_grn_animation, create_video

def create_simple_animation():
    """Create animation without membrane volume for better compatibility."""
    
    print("Creating simplified animation (no membrane volume)...")
    
    # Generate frames with membrane disabled
    frames = generate_grn_animation(
        grn_positions=[5.5],
        color_mode='helix',
        rotation_steps=24,  # Fewer steps for faster testing
        output_dir="simple_frames",
        width=1280,  # Smaller resolution
        height=720,
        reference_id='MerMAID1_model_0'
    )
    
    if frames:
        print(f"Generated {len(frames)} frames")
        
        # Create video
        video = create_video(
            frames,
            output_video="simple_opsin_animation.mp4",
            fps=24,
            quality="medium"
        )
        
        if video:
            print(f"✅ Simple animation created: {video}")
            return video
    
    print("❌ Animation failed")
    return None

if __name__ == "__main__":
    create_simple_animation()