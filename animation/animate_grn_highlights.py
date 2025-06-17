"""
Simplified animation script for GRN position highlights with rotation.

This script creates animations showing specific GRN positions highlighted
while rotating around the Z-axis to show different helix orientations.
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
from pathlib import Path
from tqdm import tqdm

# Import visualization components
from visualize_alignment_grn import (
    load_rmsd_cache, load_processed_structures, load_grn_table,
    extract_ca_coordinates_with_grn, apply_alignment_transformations,
    apply_membrane_orientation
)
from src.opsin_color_scheme import HELIX_NUMBER_COLORS, get_categorical_colors


def create_grn_highlight_frame(
    aligned_structures,
    grn_df,
    property_data,
    target_grn,
    camera_angle=0,
    color_mode='helix',
    width=1920,
    height=1080,
    show_membrane=True
):
    """
    Create a single frame with specific GRN position highlighted.
    
    Args:
        aligned_structures: Dictionary of aligned structure data
        grn_df: GRN position table
        property_data: Property data for structures
        target_grn: GRN position to highlight (e.g., 1.50)
        camera_angle: Rotation angle around Z-axis in degrees
        color_mode: 'helix' or 'property'
        width: Frame width
        height: Frame height
        show_membrane: Whether to show membrane volume
        
    Returns:
        plotly.graph_objects.Figure
    """
    
    print(f"    [DEBUG] Creating figure for GRN {target_grn} at angle {camera_angle}")
    
    # Create figure with single 3D plot
    fig = go.Figure()
    
    # Color schemes
    helix_colors = {
        1: HELIX_NUMBER_COLORS[1], 2: HELIX_NUMBER_COLORS[2], 3: HELIX_NUMBER_COLORS[3],
        4: HELIX_NUMBER_COLORS[4], 5: HELIX_NUMBER_COLORS[5], 6: HELIX_NUMBER_COLORS[6],
        7: HELIX_NUMBER_COLORS[7], 0: '#D3D3D3'
    }
    
    print(f"    [DEBUG] Color schemes loaded")
    
    # Property colors
    property_colors = {}
    if property_data and color_mode == 'property':
        molecular_functions = set()
        for struct_id, props in property_data.items():
            if 'molecular_function' in props:
                molecular_functions.add(props['molecular_function'])
        property_colors = get_categorical_colors(list(molecular_functions), property_type='property1')
    
    def get_color(struct_id, helix_num, mode):
        if mode == 'property' and property_data and struct_id in property_data:
            mol_func = property_data[struct_id].get('molecular_function', 'Unknown')
            return property_colors.get(mol_func, '#D3D3D3')
        return helix_colors.get(helix_num, '#D3D3D3')
    
    # Track what we've added to legend
    legend_added = set()
    
    print(f"    [DEBUG] Starting background structures processing...")
    
    # Add background structures (low opacity)
    structure_count = 0
    for struct_id, data in aligned_structures.items():
        if structure_count >= 125:
            break
            
        coords = data['coords']
        df = data.get('dataframe')
        
        if df is not None:
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = coords
            
            # Group by helix for line connectivity
            for helix_num, group in df_aligned.groupby('helix_num'):
                if len(group) == 0:
                    continue
                
                helix_coords = group[['x', 'y', 'z']].values
                color = get_color(struct_id, helix_num, color_mode)
                
                # Legend logic
                if color_mode == 'helix':
                    legend_key = f"Helix_{helix_num}"
                    legend_name = f"Helix {int(helix_num)}" if helix_num != 0 else "Loops"
                else:
                    mol_func = 'Unknown'
                    if property_data and struct_id in property_data:
                        mol_func = property_data[struct_id].get('molecular_function', 'Unknown')
                    legend_key = f"Property_{mol_func}"
                    legend_name = mol_func
                
                show_legend = legend_key not in legend_added
                if show_legend:
                    legend_added.add(legend_key)
                
                # Add background trace
                fig.add_trace(go.Scatter3d(
                    x=helix_coords[:, 0],
                    y=helix_coords[:, 1],
                    z=helix_coords[:, 2],
                    mode='markers+lines',
                    marker=dict(size=2, color=color, opacity=0.1),
                    line=dict(color=color, width=0.5),
                    name=legend_name,
                    legendgroup=legend_key,
                    showlegend=show_legend,
                    hoverinfo='skip',
                    opacity=0.1  # Move opacity to trace level
                ))
        
        structure_count += 1
    
    print(f"    [DEBUG] Processed {structure_count} background structures")
    
    # Add highlighted residues for target GRN
    print(f"    [DEBUG] Starting highlight residue collection for GRN {target_grn}...")
    
    highlight_coords = []
    highlight_colors = []
    highlight_text = []
    
    for struct_id, data in aligned_structures.items():
        df = data.get('dataframe')
        if df is not None:
            coords = data['coords']
            df_aligned = df.copy()
            df_aligned[['x', 'y', 'z']] = coords
            
            # Find residues with target GRN
            # Debug: Check GRN values
            grn_values = df_aligned['grn'].dropna().unique()
            if len(grn_values) < 10:
                print(f"    [DEBUG] GRN values in {struct_id}: {grn_values}")
            
            grn_matches = df_aligned[df_aligned['grn'] == target_grn]
            
            for _, row in grn_matches.iterrows():
                coord = [row['x'], row['y'], row['z']]
                helix_num = row['helix_num']
                color = get_color(struct_id, helix_num, color_mode)
                
                highlight_coords.append(coord)
                highlight_colors.append(color)
                
                if color_mode == 'helix':
                    highlight_text.append(f'{struct_id}<br>Helix: {int(helix_num) if helix_num != 0 else "Loop"}<br>GRN: {target_grn}')
                else:
                    mol_func = 'Unknown'
                    if property_data and struct_id in property_data:
                        mol_func = property_data[struct_id].get('molecular_function', 'Unknown')
                    highlight_text.append(f'{struct_id}<br>Function: {mol_func}<br>GRN: {target_grn}')
    
    print(f"    [DEBUG] Found {len(highlight_coords)} highlight coordinates")
    
    # Debug: If no highlights found, check why
    if len(highlight_coords) == 0:
        print(f"    [DEBUG] No highlights found for GRN {target_grn}")
        print(f"    [DEBUG] Checking a sample structure...")
        sample_struct = list(aligned_structures.keys())[0]
        sample_df = aligned_structures[sample_struct].get('dataframe')
        if sample_df is not None:
            grn_sample = sample_df['grn'].dropna().head(10)
            print(f"    [DEBUG] Sample GRN values from {sample_struct}: {grn_sample.tolist()}")
    
    # Add highlight trace
    if highlight_coords:
        coords_array = np.array(highlight_coords)
        fig.add_trace(go.Scatter3d(
            x=coords_array[:, 0],
            y=coords_array[:, 1],
            z=coords_array[:, 2],
            mode='markers',
            marker=dict(
                size=6,
                color=highlight_colors,
                opacity=1.0,
                line=dict(width=2, color='white')
            ),
            name=f'GRN {target_grn}',
            showlegend=False,
            text=highlight_text,
            hovertemplate='%{text}<extra></extra>'
        ))
    
    print(f"    [DEBUG] Highlight trace added")
    
    # Add membrane volume
    if show_membrane:
        print(f"    [DEBUG] Adding membrane volume...")
        # Simple membrane representation
        x_range = [-20, 20]
        y_range = [-20, 20]
        z_vol = np.linspace(-10, 10, 8)
        x_vol = np.linspace(x_range[0], x_range[1], 10)
        y_vol = np.linspace(y_range[0], y_range[1], 10)
        
        X_vol, Y_vol, Z_vol = np.meshgrid(x_vol, y_vol, z_vol, indexing='ij')
        membrane_values = np.ones_like(X_vol) * 0.5
        
        fig.add_trace(go.Volume(
            x=X_vol.flatten(),
            y=Y_vol.flatten(),
            z=Z_vol.flatten(),
            value=membrane_values.flatten(),
            isomin=0.3,
            isomax=0.7,
            opacity=0.02,
            surface_count=2,
            colorscale=[[0, 'lightgray'], [1, 'lightgray']],
            showscale=False,
            name='Membrane',
            showlegend=True,
            hoverinfo='skip'
        ))
        print(f"    [DEBUG] Membrane volume added")
    else:
        print(f"    [DEBUG] Membrane volume skipped")
    
    # Calculate camera position
    print(f"    [DEBUG] Setting up camera and layout...")
    angle_rad = np.radians(camera_angle)
    camera_distance = 2.5
    camera_x = camera_distance * np.cos(angle_rad)
    camera_y = camera_distance * np.sin(angle_rad)
    camera_z = 0.8
    
    # Update layout
    fig.update_layout(
        title=f'MOGRN: GRN {target_grn} Highlighted | {color_mode.title()} Mode | Angle: {camera_angle:.0f}°',
        width=width,
        height=height,
        scene=dict(
            camera=dict(
                eye=dict(x=camera_x, y=camera_y, z=camera_z),
                center=dict(x=0, y=0, z=0),
                up=dict(x=0, y=0, z=1)
            ),
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False),
            bgcolor='white',
            aspectmode='cube'
        ),
        paper_bgcolor='white',
        plot_bgcolor='white',
        showlegend=True,
        legend=dict(
            yanchor="top", y=0.99,
            xanchor="left", x=0.01,
            bgcolor='rgba(255,255,255,0.8)'
        ),
        margin=dict(l=0, r=0, t=60, b=0)
    )
    
    print(f"    [DEBUG] Figure creation complete")
    
    return fig


def generate_grn_animation(
    grn_positions=["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"],
    color_mode='helix',
    rotation_steps=24,  # 15-degree increments
    output_dir="grn_animation_frames",
    width=1920,
    height=1080,
    reference_id='MerMAID1_model_0',
    show_membrane=False  # Disable by default to avoid rendering issues
):
    """
    Generate animation frames for GRN highlights with rotation.
    
    Args:
        grn_positions: List of GRN positions to highlight
        color_mode: 'helix' or 'property'
        rotation_steps: Number of rotation steps per GRN position
        output_dir: Directory to save frames
        width: Frame width
        height: Frame height
        reference_id: Reference structure ID
        
    Returns:
        List of frame file paths
    """
    
    print("=== Loading Data ===")
    cache_data = load_rmsd_cache()
    alignment_paths = cache_data.get('alignment_paths', {})
    processed_structures = load_processed_structures()
    grn_df = load_grn_table()
    
    # Load property data
    property_data = None
    if color_mode == 'property':
        from src.data_processing import load_opsin_property_data
        property_file = Path('property/mo_exp.csv')
        if property_file.exists():
            try:
                property_result = load_opsin_property_data(property_file, processed_structures)
                if property_result and 'properties' in property_result:
                    property_data = property_result['properties']
            except Exception as e:
                print(f"Could not load property data: {e}")
    
    print("=== Processing Structures ===")
    structures = extract_ca_coordinates_with_grn(processed_structures, grn_df, chain_id='A', use_helix_only=True)
    aligned_structures = apply_alignment_transformations(structures, alignment_paths, reference_id)
    oriented_structures = apply_membrane_orientation(aligned_structures, reference_id)
    
    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)
    
    frame_files = []
    frame_number = 0
    
    print(f"=== Generating {len(grn_positions)} GRN positions × {rotation_steps} angles = {len(grn_positions) * rotation_steps} frames ===")
    
    for grn_idx, grn_pos in enumerate(grn_positions):
        print(f"\nGRN {grn_pos} ({grn_idx + 1}/{len(grn_positions)})")
        
        for step in tqdm(range(rotation_steps), desc=f"GRN {grn_pos}"):
            angle = (step / rotation_steps) * 360
            
            print(f"\n[DEBUG] Frame {frame_number}: Creating figure for angle {angle:.1f}°")
            
            try:
                fig = create_grn_highlight_frame(
                    aligned_structures=oriented_structures,
                    grn_df=grn_df,
                    property_data=property_data,
                    target_grn=grn_pos,
                    camera_angle=angle,
                    color_mode=color_mode,
                    width=width,
                    height=height,
                    show_membrane=show_membrane
                )
                print(f"[DEBUG] Frame {frame_number}: Figure created successfully")
                
            except Exception as e:
                print(f"[ERROR] Frame {frame_number}: Failed to create figure: {e}")
                import traceback
                traceback.print_exc()
                continue
            
            frame_file = output_path / f"frame_{frame_number:04d}.png"
            print(f"[DEBUG] Frame {frame_number}: Saving to {frame_file}")
            
            try:
                fig.write_image(str(frame_file), width=width, height=height, scale=1)
                print(f"[DEBUG] Frame {frame_number}: Image saved successfully")
                frame_files.append(str(frame_file))
                frame_number += 1
            except Exception as e:
                print(f"[ERROR] Frame {frame_number}: Failed to save image: {e}")
                import traceback
                traceback.print_exc()
                continue
    
    print(f"\nGenerated {len(frame_files)} frames")
    return frame_files


def create_video(frame_files, output_video="grn_highlights.mp4", fps=24, quality="high"):
    """Create video from frames using ffmpeg."""
    
    if not frame_files:
        print("No frames to process!")
        return None
    
    import subprocess
    
    frame_dir = Path(frame_files[0]).parent
    frame_pattern = str(frame_dir / "frame_%04d.png")
    
    quality_settings = {
        'low': ['-crf', '28'],
        'medium': ['-crf', '23'],
        'high': ['-crf', '18'],
        'ultra': ['-crf', '15']
    }
    
    cmd = [
        'ffmpeg', '-y',
        '-framerate', str(fps),
        '-i', frame_pattern,
        '-c:v', 'libx264',
        '-pix_fmt', 'yuv420p'
    ] + quality_settings.get(quality, quality_settings['high']) + [output_video]
    
    try:
        subprocess.run(cmd, check=True)
        print(f"Video created: {output_video}")
        return output_video
    except (subprocess.CalledProcessError, FileNotFoundError) as e:
        print(f"ffmpeg error: {e}")
        print("Make sure ffmpeg is installed and in your PATH")
        return None


def main():
    """Create GRN highlight animation."""
    
    # Generate frames
    frames = generate_grn_animation(
        grn_positions=["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"],
        color_mode='helix',
        rotation_steps=24,  # 15° increments
        width=1920,
        height=1080
    )
    
    # Create video
    if frames:
        video = create_video(
            frames,
            output_video="opsin_grn_helix_animation.mp4",
            fps=24,
            quality="high"
        )
        
        if video:
            print(f"\n✅ Animation complete: {video}")
            print(f"Duration: ~{len(frames)/24:.1f} seconds")
            
            # Create property mode version
            property_frames = generate_grn_animation(
                grn_positions=["1.50", "2.50", "3.50", "4.50", "5.50", "6.50", "7.50"],
                color_mode='property',
                rotation_steps=24,
                output_dir="grn_property_frames",
                width=1920,
                height=1080
            )
            
            if property_frames:
                property_video = create_video(
                    property_frames,
                    output_video="opsin_grn_property_animation.mp4",
                    fps=24,
                    quality="high"
                )
                print(f"✅ Property animation: {property_video}")


if __name__ == "__main__":
    main()