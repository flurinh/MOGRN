# MOGRN Animation System

Create rotating animations of opsin structures with GRN position highlights using ffmpeg.

## Quick Start

### Prerequisites
```bash
# Install ffmpeg
sudo apt install ffmpeg  # Ubuntu/Debian
brew install ffmpeg      # macOS
# Windows: Download from https://ffmpeg.org/

# Install Python dependencies
pip install plotly kaleido tqdm
```

### Simple Usage
```bash
# Create a quick overview animation
python create_animations.py single_rotation

# Create binding pocket tour
python create_animations.py binding_pocket

# Create complete helix tour
python create_animations.py helix_tour

# Create all presentation videos
python create_animations.py presentation
```

## Animation Types

### 1. Single Rotation (`single_rotation`)
- **GRN Position**: 5.50 (center of retinal binding pocket)
- **Duration**: ~1.5 seconds
- **Use Case**: Quick overview for presentations
- **Output**: `opsin_single_rotation.mp4`

### 2. Binding Pocket Tour (`binding_pocket`)
- **GRN Positions**: 3.35, 4.50, 5.50, 6.44, 6.55
- **Duration**: ~7.5 seconds
- **Use Case**: Focus on retinal binding region
- **Output**: `opsin_binding_pocket.mp4`

### 3. Helix Tour (`helix_tour`)
- **GRN Positions**: 1.50, 2.50, 3.50, 4.50, 5.50, 6.50, 7.50
- **Duration**: ~7 seconds
- **Use Case**: Complete structural overview
- **Output**: `opsin_helix_tour.mp4`

### 4. Property Comparison (`property_comparison`)
- **GRN Positions**: 1.50-7.50 (all helices)
- **Coloring**: Functional classification
- **Duration**: ~7 seconds
- **Use Case**: Show functional diversity
- **Output**: `opsin_property_tour.mp4`

## Custom Animations

### Python API
```python
from create_animations import create_custom_animation

# Custom GRN positions and rotation
video = create_custom_animation(
    grn_positions=[2.50, 3.50, 5.50],  # Specific positions
    color_mode="helix",                 # or "property"
    rotation_steps=36,                  # 10° increments
    output_file="my_animation.mp4",
    fps=24,
    quality="high"
)
```

### Advanced Control
```python
from animate_grn_highlights import generate_grn_animation, create_video

# Generate frames with full control
frames = generate_grn_animation(
    grn_positions=[1.50, 3.35, 5.50, 7.50],
    color_mode='helix',
    rotation_steps=24,           # 15° per step
    output_dir="my_frames",
    width=1920,                  # 4K: 3840
    height=1080,                 # 4K: 2160
    reference_id='MerMAID1_model_0'
)

# Create video with custom settings
video = create_video(
    frames,
    output_video="high_quality.mp4",
    fps=30,                      # Smooth motion
    quality="ultra"              # Highest quality
)
```

## Animation Structure

Each animation follows this pattern:
1. **GRN Sequence**: Cycles through specified GRN positions
2. **Rotation**: 360° rotation around Z-axis for each GRN
3. **Highlighting**: Target GRN residues shown with large, bright markers
4. **Background**: All structures shown at low opacity for context

### Camera Movement
- **Rotation Axis**: Z-axis (membrane normal)
- **Camera Distance**: 2.5 units from center
- **Elevation**: Slight upward angle (0.8 Z-offset)
- **Rotation Speed**: Consistent angular velocity

### Visual Elements
- **Background Structures**: Low opacity (0.1), thin lines
- **Highlighted Residues**: High opacity (1.0), large markers, white outlines
- **Membrane**: Translucent volume block (Z: -10 to +10)
- **Legend**: Shows helix or functional classifications

## File Organization

```
MOGRN/
├── create_animations.py          # Quick animation creator
├── animate_grn_highlights.py      # Core animation functions
├── create_animation.py           # Full-featured animation system
├── visualize_alignment_grn.py    # Base visualization functions
│
├── *_frames/                     # Frame directories (temporary)
│   ├── frame_0000.png
│   ├── frame_0001.png
│   └── ...
│
└── *.mp4                         # Output videos
    ├── opsin_single_rotation.mp4
    ├── opsin_binding_pocket.mp4
    ├── opsin_helix_tour.mp4
    └── opsin_property_tour.mp4
```

## Technical Details

### Frame Generation
- **Resolution**: 1920×1080 (Full HD) - configurable to 4K
- **Format**: PNG frames → MP4 video
- **Frame Rate**: 24 fps (configurable)
- **Quality**: H.264 encoding with CRF quality control

### Performance
- **Frame Generation**: ~2-5 seconds per frame
- **Total Time**: ~10-20 minutes for full helix tour
- **Disk Usage**: ~50MB per 100 frames (temporary)
- **Final Video**: 5-50MB depending on length and quality

### Memory Requirements
- **Frame Generation**: 4-8GB RAM (for 125 structures)
- **ffmpeg Encoding**: Minimal additional memory
- **Storage**: Plan for 1-2GB temporary frame storage

## Quality Settings

### Video Quality Options
```python
quality_settings = {
    'low': '-crf 28',      # ~2MB/minute, fast encoding
    'medium': '-crf 23',   # ~5MB/minute, balanced
    'high': '-crf 18',     # ~10MB/minute, high quality
    'ultra': '-crf 15'     # ~20MB/minute, maximum quality
}
```

### Resolution Options
- **1920×1080**: Standard Full HD (recommended)
- **3840×2160**: 4K Ultra HD (for high-end displays)
- **1280×720**: HD (for smaller file sizes)

## Usage Examples

### For Publications
```bash
# High-quality, single rotation for figure
python create_animations.py single_rotation
# Use in presentations, supplement videos

# Binding pocket focus for detailed analysis
python create_animations.py binding_pocket
# Show retinal binding site dynamics
```

### For Presentations
```bash
# Quick overview set
python create_animations.py presentation
# Creates 4 complementary videos

# Custom functional comparison
python -c "
from create_animations import create_custom_animation
create_custom_animation(
    grn_positions=[3.35, 5.50, 6.44],
    color_mode='property',
    rotation_steps=48,  # Smoother rotation
    output_file='functional_highlight.mp4'
)"
```

### For Web/Social Media
```bash
# Shorter, faster animations
python -c "
from animate_grn_highlights import generate_grn_animation, create_video
frames = generate_grn_animation(
    grn_positions=[5.50], 
    rotation_steps=20,  # Faster rotation
    width=1280, height=720  # Smaller size
)
create_video(frames, 'social_media.mp4', fps=30, quality='medium')
"
```

## Troubleshooting

### Common Issues

**1. ffmpeg not found**
```bash
# Install ffmpeg first
sudo apt install ffmpeg  # Linux
brew install ffmpeg      # macOS
```

**2. Out of memory**
```python
# Reduce structure count or resolution
generate_grn_animation(width=1280, height=720)  # Smaller frames
```

**3. Slow frame generation**
```python
# Reduce rotation steps for faster generation
generate_grn_animation(rotation_steps=12)  # 30° increments
```

**4. Large file sizes**
```python
# Use lower quality setting
create_video(frames, quality='medium')  # Instead of 'high'
```

### Performance Tips

1. **Start Small**: Test with `single_rotation` first
2. **Monitor Disk Space**: Frame generation can use 1-2GB temporarily
3. **Use SSD**: Faster disk I/O improves frame generation speed
4. **Close Other Apps**: Frame generation is memory-intensive

## Integration with GitHub

### Adding to README
```markdown
## 🎬 Animations

View our opsin structure animations:

- [**Single Rotation**](https://github.com/user/repo/raw/main/opsin_single_rotation.mp4) - Quick overview
- [**Binding Pocket Tour**](https://github.com/user/repo/raw/main/opsin_binding_pocket.mp4) - Retinal binding analysis  
- [**Complete Helix Tour**](https://github.com/user/repo/raw/main/opsin_helix_tour.mp4) - Full structural overview

*Download and open in video player for best quality*
```

### Git LFS for Large Files
```bash
# For videos >100MB, use Git LFS
git lfs install
git lfs track "*.mp4"
git add .gitattributes
git add *.mp4
git commit -m "Add animation videos"
```

## Future Enhancements

- **Side-by-side comparison**: Helix vs property coloring
- **Sequence-based animations**: Follow evolutionary relationships
- **Interactive web videos**: HTML5 with playback controls
- **Stereo 3D**: For VR/3D display systems
- **Custom camera paths**: Non-circular movements