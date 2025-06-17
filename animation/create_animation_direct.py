#!/usr/bin/env python3
"""
Direct PNG export solution using plotly with offline mode.
This script creates frames and exports them directly as images.
"""

import os
import sys
import time
from pathlib import Path
import subprocess
import json
import numpy as np
from tqdm import tqdm

# Add a workaround for kaleido
os.environ['PLOTLY_RENDERER'] = 'browser'

from animate_grn_highlights import (
    load_rmsd_cache, load_processed_structures, load_grn_table,
    extract_ca_coordinates_with_grn, apply_alignment_transformations,
    apply_membrane_orientation, create_grn_highlight_frame
)


def test_kaleido_fix():
    """Test if kaleido can be fixed with environment settings."""
    try:
        # Try different kaleido settings
        import plotly.graph_objects as go
        
        # Create a simple test figure
        fig = go.Figure(data=[go.Scatter3d(x=[1,2,3], y=[1,2,3], z=[1,2,3])])
        
        # Try with specific kaleido settings
        import plotly.io as pio
        pio.kaleido.scope.mathjax = None
        pio.kaleido.scope.chromium_args = ['--disable-gpu', '--no-sandbox']
        
        # Test export
        print("Testing kaleido export with fixes...")
        fig.write_image("test_kaleido.png", width=100, height=100)
        print("✅ Kaleido export successful!")
        os.remove("test_kaleido.png")
        return True
        
    except Exception as e:
        print(f"❌ Kaleido still not working: {e}")
        return False


def export_json_to_static_image(json_file, output_png, width=1920, height=1080):
    """
    Convert JSON figure to static PNG using alternative method.
    """
    import plotly.graph_objects as go
    import plotly.io as pio
    
    # Load figure from JSON
    with open(json_file, 'r') as f:
        fig_dict = json.load(f)
    
    fig = go.Figure(fig_dict)
    
    # Try to export with timeout
    try:
        # Configure kaleido
        pio.kaleido.scope.mathjax = None
        pio.kaleido.scope.chromium_args = ['--disable-gpu', '--no-sandbox', '--single-process']
        
        # Export with short timeout
        fig.write_image(output_png, width=width, height=height, scale=1)
        return True
    except Exception as e:
        print(f"Direct export failed: {e}")
        
        # Fallback: Save as static HTML with embedded image
        try:
            static_html = fig.to_html(
                include_plotlyjs='inline',
                config={'staticPlot': True, 'toImageButtonOptions': {
                    'format': 'png',
                    'width': width,
                    'height': height,
                    'scale': 1
                }}
            )
            
            html_file = str(output_png).replace('.png', '_static.html')
            with open(html_file, 'w') as f:
                f.write(static_html)
            print(f"Saved static HTML: {html_file}")
            return False
        except:
            return False


def create_simple_rotation_animation_fixed(
    num_frames=36,
    output_dir="rotation_frames_png",
    width=1920,
    height=1080
):
    """
    Create rotation animation with fixed export process.
    """
    
    print("=== Testing Kaleido Fix ===")
    kaleido_works = test_kaleido_fix()
    
    if not kaleido_works:
        print("\n⚠️ Kaleido is not working. Will save frames as JSON and HTML.")
        print("You'll need to convert them manually or use a different system.\n")
    
    print("=== Loading Data ===")
    cache_data = load_rmsd_cache()
    alignment_paths = cache_data.get('alignment_paths', {})
    processed_structures = load_processed_structures()
    grn_df = load_grn_table()
    
    print("=== Processing Structures ===") 
    structures = extract_ca_coordinates_with_grn(processed_structures, grn_df, chain_id='A', use_helix_only=True)
    aligned_structures = apply_alignment_transformations(structures, alignment_paths, 'MerMAID1_model_0')
    oriented_structures = apply_membrane_orientation(aligned_structures, 'MerMAID1_model_0')
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    # Define GRN positions for each helix
    helix_grns = ["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"]
    
    print(f"\n=== Generating {num_frames} frames ===")
    
    successful_pngs = []
    all_files = []
    
    for frame_num in tqdm(range(num_frames), desc="Creating frames"):
        angle = (frame_num / num_frames) * 360
        
        # Simple rotation through helices
        helix_idx = int((angle / 360) * 7) % 7
        target_grn = helix_grns[helix_idx]
        
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
                show_membrane=False
            )
            
            # Save files
            base_name = output_path / f"frame_{frame_num:04d}"
            json_file = f"{base_name}.json"
            png_file = f"{base_name}.png"
            
            # Always save JSON as backup
            fig.write_json(json_file)
            all_files.append(json_file)
            
            # Try PNG export if kaleido works
            if kaleido_works:
                try:
                    import plotly.io as pio
                    pio.kaleido.scope.mathjax = None
                    pio.kaleido.scope.chromium_args = ['--disable-gpu', '--no-sandbox']
                    
                    fig.write_image(png_file, width=width, height=height)
                    successful_pngs.append(png_file)
                except Exception as e:
                    print(f"\nFrame {frame_num} PNG export failed: {e}")
            
        except Exception as e:
            print(f"\nError on frame {frame_num}: {e}")
            continue
    
    print(f"\n✅ Created {len(all_files)} frame files")
    print(f"   - {len(successful_pngs)} PNG files")
    print(f"   - {len(all_files)} JSON backup files")
    
    # Try to create video if we have PNGs
    if successful_pngs:
        print("\n=== Creating video with ffmpeg ===")
        video_file = output_path / "rotation_animation.mp4"
        
        cmd = [
            'ffmpeg', '-y',
            '-framerate', '10',
            '-i', str(output_path / 'frame_%04d.png'),
            '-c:v', 'libx264',
            '-pix_fmt', 'yuv420p',
            '-crf', '18',
            str(video_file)
        ]
        
        try:
            subprocess.run(cmd, check=True)
            print(f"✅ Video created: {video_file}")
        except Exception as e:
            print(f"❌ Video creation failed: {e}")
    else:
        print("\n⚠️ No PNG files created. Manual conversion needed.")
        print("\nNext steps:")
        print("1. Use the JSON files with a working Plotly installation")
        print("2. Or use the frame converter script")
        print("3. Or try on a different system where kaleido works")
        
        # Create a converter script
        create_json_to_png_converter(output_path)
    
    return all_files


def create_json_to_png_converter(output_path):
    """Create a converter script for the JSON files."""
    
    converter_content = '''#!/usr/bin/env python3
"""
Convert JSON frames to PNG images.
Run this on a system where kaleido works properly.
"""

import json
import glob
from pathlib import Path
import plotly.graph_objects as go
from tqdm import tqdm

def convert_all_json_to_png(input_dir=".", width=1920, height=1080):
    """Convert all JSON files in directory to PNG."""
    
    json_files = sorted(glob.glob(f"{input_dir}/frame_*.json"))
    print(f"Found {len(json_files)} JSON files to convert")
    
    successful = 0
    
    for json_file in tqdm(json_files, desc="Converting"):
        try:
            # Load figure
            with open(json_file, 'r') as f:
                fig_dict = json.load(f)
            fig = go.Figure(fig_dict)
            
            # Export as PNG
            png_file = json_file.replace('.json', '.png')
            fig.write_image(png_file, width=width, height=height)
            successful += 1
            
        except Exception as e:
            print(f"\\nError converting {json_file}: {e}")
    
    print(f"\\n✅ Successfully converted {successful}/{len(json_files)} files")
    
    if successful > 0:
        print("\\nCreate video with:")
        print(f"ffmpeg -framerate 10 -i frame_%04d.png -c:v libx264 -pix_fmt yuv420p output.mp4")

if __name__ == "__main__":
    convert_all_json_to_png()
'''
    
    converter_file = output_path / "convert_json_to_png.py"
    with open(converter_file, 'w') as f:
        f.write(converter_content)
    
    print(f"\n📝 Created converter script: {converter_file}")
    print("Run this script on a system where kaleido works to convert JSON to PNG")


def main():
    """Create the animation with workarounds."""
    
    print("MOGRN Rotation Animation Generator")
    print("=" * 50)
    print("This script will attempt to create a rotation animation.")
    print("If kaleido fails, it will save JSON files for later conversion.\n")
    
    # Create animation
    files = create_simple_rotation_animation_fixed(
        num_frames=36,  # 36 frames for testing
        output_dir="rotation_frames_final",
        width=1920,
        height=1080
    )
    
    print("\n" + "="*50)
    print("Process complete!")
    print("="*50)
    
    # Check what we have
    output_dir = Path("rotation_frames_final")
    png_files = list(output_dir.glob("*.png"))
    json_files = list(output_dir.glob("*.json"))
    
    if png_files:
        print(f"\n✅ Success! Found {len(png_files)} PNG files")
        print("You can now create a video with:")
        print(f"ffmpeg -framerate 10 -i {output_dir}/frame_%04d.png -c:v libx264 -pix_fmt yuv420p output.mp4")
    else:
        print(f"\n⚠️ No PNG files created, but {len(json_files)} JSON files are available")
        print("\nOptions:")
        print(f"1. Run the converter script: python {output_dir}/convert_json_to_png.py")
        print("2. Transfer files to a system where kaleido works")
        print("3. Use online Plotly chart studio to convert")


if __name__ == "__main__":
    main()