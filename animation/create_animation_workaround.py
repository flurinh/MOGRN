#!/usr/bin/env python3
"""
Workaround script for creating animations when kaleido has issues.
Uses plotly's built-in browser rendering to create frames.
"""

import numpy as np
import webbrowser
import time
import os
from pathlib import Path
from tqdm import tqdm
import subprocess

from animate_grn_highlights import (
    load_rmsd_cache, load_processed_structures, load_grn_table,
    extract_ca_coordinates_with_grn, apply_alignment_transformations,
    apply_membrane_orientation, create_grn_highlight_frame
)


def export_frame_via_browser(fig, output_file, width=1920, height=1080):
    """
    Export frame using Plotly's built-in browser screenshot capability.
    """
    # This uses plotly's show method with custom config
    config = {
        'toImageButtonOptions': {
            'format': 'png',
            'filename': output_file.stem,
            'height': height,
            'width': width,
            'scale': 1
        }
    }
    
    # Save as temporary HTML with auto-download config
    temp_html = output_file.with_suffix('.html')
    
    # Add JavaScript to automatically download the image
    html_str = fig.to_html(include_plotlyjs='cdn', config=config)
    
    # Inject auto-download script
    download_script = f"""
    <script>
    window.onload = function() {{
        setTimeout(function() {{
            // Get the plot
            var plot = document.getElementsByClassName('plotly-graph-div')[0];
            
            // Download as image
            Plotly.downloadImage(plot, {{
                format: 'png',
                width: {width},
                height: {height},
                filename: '{output_file.stem}'
            }});
            
            // Close window after download
            setTimeout(function() {{
                window.close();
            }}, 2000);
        }}, 1000);
    }};
    </script>
    """
    
    html_str = html_str.replace('</body>', download_script + '</body>')
    
    with open(temp_html, 'w') as f:
        f.write(html_str)
    
    return temp_html


def create_simple_rotation_animation(
    num_frames=36,  # 36 frames = 10-degree increments
    output_dir="rotation_frames_simple",
    width=1920,
    height=1080
):
    """
    Create a simple rotation animation highlighting different helices.
    This version creates files that can be manually processed.
    """
    
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
    
    # Define which GRN to highlight at each angle
    # For simplicity, cycle through helices 1-7
    helix_grns = ["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"]
    
    print(f"\n=== Generating {num_frames} frames ===")
    
    html_files = []
    instructions = []
    
    for frame_num in tqdm(range(num_frames), desc="Creating frames"):
        angle = (frame_num / num_frames) * 360
        
        # Determine which helix to highlight based on angle
        # Simple approach: divide 360 degrees into 7 segments
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
            
            # Save as static HTML with embedded image
            output_file = output_path / f"frame_{frame_num:04d}"
            
            # Create HTML file with embedded plot
            html_file = output_file.with_suffix('.html')
            
            # Create a static image HTML
            static_html = f"""
<!DOCTYPE html>
<html>
<head>
    <title>Frame {frame_num}</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body {{ margin: 0; padding: 0; }}
        #myDiv {{ width: {width}px; height: {height}px; }}
    </style>
</head>
<body>
    <div id="myDiv"></div>
    <script>
        var data = {fig.to_json()};
        var figure = JSON.parse(data);
        var layout = figure.layout;
        layout.width = {width};
        layout.height = {height};
        Plotly.newPlot('myDiv', figure.data, layout, {{staticPlot: true}});
    </script>
</body>
</html>
"""
            
            with open(html_file, 'w') as f:
                f.write(static_html)
            
            html_files.append(str(html_file))
            
            # Also save the raw figure data
            fig.write_json(str(output_file.with_suffix('.json')))
            
        except Exception as e:
            print(f"\nError on frame {frame_num}: {e}")
            continue
    
    # Create batch conversion script
    create_conversion_script(output_path, html_files)
    
    print(f"\n✅ Created {len(html_files)} frame files")
    print(f"\n📁 Output directory: {output_path}")
    print("\n🎬 Next steps to create video:")
    print("1. Open each HTML file in a browser and save as PNG")
    print("2. Or use the browser automation script: batch_convert.py")
    print("3. Once you have PNG files, run:")
    print(f"   ffmpeg -framerate 10 -i {output_path}/frame_%04d.png -c:v libx264 -pix_fmt yuv420p output.mp4")
    
    return html_files


def create_conversion_script(output_path, html_files):
    """Create a script to help with batch conversion."""
    
    script_content = f"""#!/usr/bin/env python3
'''
Batch conversion helper script.
This script provides options for converting HTML frames to PNG.
'''

import os
from pathlib import Path

output_dir = Path(r"{output_path}")
html_files = {html_files}

print("Frame Conversion Helper")
print("=" * 50)
print(f"Found {{len(html_files)}} HTML files in {{output_dir}}")
print()
print("Options for converting to PNG:")
print()
print("1. MANUAL METHOD (Most reliable):")
print("   - Open each HTML file in Chrome/Firefox")
print("   - Right-click > Save image as... > Save as PNG")
print("   - Name them frame_0000.png, frame_0001.png, etc.")
print()
print("2. BROWSER AUTOMATION (Requires Selenium):")
print("   pip install selenium pillow")
print("   Then run the selenium_convert() function below")
print()
print("3. ONLINE CONVERTER:")
print("   - Use online HTML to PNG converters")
print("   - Upload the HTML files in batches")
print()
print("Once you have PNG files, create video with:")
print(f"ffmpeg -framerate 10 -i {{output_dir}}/frame_%04d.png -c:v libx264 -pix_fmt yuv420p output.mp4")

def selenium_convert():
    '''Convert using Selenium (requires Chrome/Firefox driver)'''
    try:
        from selenium import webdriver
        from selenium.webdriver.chrome.options import Options
        from PIL import Image
        import time
        
        options = Options()
        options.add_argument('--headless')
        options.add_argument(f'--window-size={{1920}},{{1080}}')
        
        driver = webdriver.Chrome(options=options)
        
        for i, html_file in enumerate(html_files):
            print(f"Converting {{i+1}}/{{len(html_files)}}: {{Path(html_file).name}}")
            
            # Load HTML
            driver.get(f"file://{{os.path.abspath(html_file)}}")
            time.sleep(2)  # Wait for plot to render
            
            # Take screenshot
            png_file = Path(html_file).with_suffix('.png')
            driver.save_screenshot(str(png_file))
            
            # Crop if needed (remove browser chrome)
            img = Image.open(png_file)
            # Adjust crop values as needed
            # img = img.crop((0, 0, 1920, 1080))
            img.save(png_file)
        
        driver.quit()
        print("Conversion complete!")
        
    except ImportError:
        print("Selenium not installed. Run: pip install selenium")
    except Exception as e:
        print(f"Error: {{e}}")

if __name__ == "__main__":
    # Uncomment to run selenium conversion
    # selenium_convert()
    pass
"""
    
    script_file = output_path / "batch_convert.py"
    with open(script_file, 'w') as f:
        f.write(script_content)
    
    print(f"\n📝 Created helper script: {script_file}")


def main():
    """Create the animation frames."""
    
    # Create frames
    html_files = create_simple_rotation_animation(
        num_frames=36,  # 36 frames for testing
        output_dir="rotation_frames_test",
        width=1920,
        height=1080
    )
    
    print("\n" + "="*60)
    print("IMPORTANT: Due to kaleido issues, frames are saved as HTML")
    print("="*60)
    print("\nTo create the video:")
    print("1. Install working image export solution:")
    print("   - Fix kaleido: pip install -U kaleido")
    print("   - Or install plotly-orca")
    print("   - Or use Selenium for browser automation")
    print("\n2. Convert HTML files to PNG")
    print("\n3. Create video with ffmpeg")
    print("\nFor the full 1080-frame animation at 60fps, change:")
    print("  num_frames=1080")
    print("in the create_simple_rotation_animation() call")


if __name__ == "__main__":
    main()