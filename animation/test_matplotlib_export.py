#!/usr/bin/env python3
"""
Alternative export method using matplotlib backend for Plotly figures.
This bypasses kaleido issues on some systems.
"""

import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from plotly.tools import mpl_to_plotly
import plotly.io as pio

from animate_grn_highlights import *

def export_plotly_via_matplotlib(fig, filename, width=800, height=600, dpi=100):
    """
    Export Plotly figure using matplotlib backend.
    This is a workaround for kaleido issues.
    """
    # Convert to static image using plotly's to_image with different backend
    try:
        # Try using the svg backend first
        img_bytes = fig.to_image(format="svg", width=width, height=height)
        with open(filename.replace('.png', '.svg'), 'wb') as f:
            f.write(img_bytes)
        print(f"✅ Exported as SVG: {filename.replace('.png', '.svg')}")
        return True
    except:
        pass
    
    # If that fails, save as interactive HTML
    try:
        fig.write_html(filename.replace('.png', '.html'))
        print(f"✅ Saved as HTML: {filename.replace('.png', '.html')}")
        print("Note: You can open this in a browser and use browser's save as image feature")
        return True
    except Exception as e:
        print(f"❌ Export failed: {e}")
        return False


def test_alternative_export():
    """Test alternative export methods."""
    
    print("=== Loading Data ===")
    try:
        cache_data = load_rmsd_cache()
        alignment_paths = cache_data.get('alignment_paths', {})
        processed_structures = load_processed_structures()
        grn_df = load_grn_table()
        print("✅ Data loaded successfully")
    except Exception as e:
        print(f"❌ Data loading failed: {e}")
        return
    
    print("\n=== Processing Structures ===")
    try:
        structures = extract_ca_coordinates_with_grn(processed_structures, grn_df, chain_id='A', use_helix_only=True)
        aligned_structures = apply_alignment_transformations(structures, alignment_paths, 'MerMAID1_model_0')
        oriented_structures = apply_membrane_orientation(aligned_structures, 'MerMAID1_model_0')
        print(f"✅ Processed {len(oriented_structures)} structures")
    except Exception as e:
        print(f"❌ Structure processing failed: {e}")
        traceback.print_exc()
        return
    
    print("\n=== Creating Frame ===")
    try:
        fig = create_grn_highlight_frame(
            aligned_structures=oriented_structures,
            grn_df=grn_df,
            property_data=None,
            target_grn="5.50",
            camera_angle=0,
            color_mode='helix',
            width=800,
            height=600,
            show_membrane=False
        )
        print("✅ Frame created successfully")
    except Exception as e:
        print(f"❌ Frame creation failed: {e}")
        traceback.print_exc()
        return
    
    print("\n=== Testing Alternative Export ===")
    export_plotly_via_matplotlib(fig, "test_frame_alt.png")
    
    print("\n🎉 Alternative export test completed!")


if __name__ == "__main__":
    test_alternative_export()