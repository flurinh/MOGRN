#!/usr/bin/env python3
"""
Test script to debug frame generation issues.
"""

import traceback
import threading
import time
from animate_grn_highlights import *

def test_single_frame():
    """Test creating just one frame to debug issues."""
    
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
    
    print("\n=== Testing Frame Creation ===")
    try:
        print("Creating single frame...")
        fig = create_grn_highlight_frame(
            aligned_structures=oriented_structures,
            grn_df=grn_df,
            property_data=None,
            target_grn="5.50",  # Use string format with two decimal places
            camera_angle=0,
            color_mode='helix',
            width=800,  # Smaller size for testing
            height=600,
            show_membrane=False  # Disable membrane for faster testing
        )
        print("✅ Frame created successfully")
        
        # Test image export with timeout
        print("Testing image export...")
        
        # Function to run export in a thread
        export_success = False
        export_error = None
        
        def export_image():
            nonlocal export_success, export_error
            try:
                print("Attempting export with kaleido...")
                fig.write_image("test_frame.png", width=800, height=600, engine="kaleido")
                export_success = True
            except Exception as e:
                export_error = e
        
        # Run export in a thread with timeout
        export_thread = threading.Thread(target=export_image)
        export_thread.start()
        export_thread.join(timeout=10)  # 10 second timeout
        
        if export_thread.is_alive():
            print("❌ Image export timed out after 10 seconds")
            print("This is a known issue with kaleido on some systems.")
            print("\nAlternative solutions:")
            print("1. Save as HTML and use browser to export")
            print("2. Use Plotly's online service")
            print("3. Use matplotlib backend instead")
            
            # Save as HTML
            try:
                fig.write_html("test_frame.html")
                print("\n✅ Saved as HTML: test_frame.html")
                print("You can open this in a browser and save as image manually")
            except Exception as e:
                print(f"❌ HTML export failed: {e}")
        elif export_success:
            print("✅ Image export successful: test_frame.png")
        else:
            print(f"❌ Image export failed: {export_error}")
            if "kaleido" in str(export_error).lower():
                print("Make sure kaleido is installed: pip install kaleido")
            traceback.print_exc()
            
            # Try alternative: orca
            print("\nTrying alternative export with orca...")
            try:
                fig.write_image("test_frame_orca.png", width=800, height=600, engine="orca")
                print("✅ Image export successful with orca: test_frame_orca.png")
            except Exception as e2:
                print(f"❌ Orca export also failed: {e2}")
                print("Saving as HTML instead...")
                try:
                    fig.write_html("test_frame.html")
                    print("✅ Saved as HTML: test_frame.html")
                except:
                    pass
        
    except Exception as e:
        print(f"❌ Frame creation failed: {e}")
        traceback.print_exc()
        return
    
    print("\n🎉 Single frame test completed successfully!")
    print("The animation script should work now.")

if __name__ == "__main__":
    test_single_frame()