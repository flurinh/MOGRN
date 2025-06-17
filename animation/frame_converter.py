
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
