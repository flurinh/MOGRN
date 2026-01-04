"""Generate a reduced interactive opsin alignment figure with enriched hover text."""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict
import sys

PROJECT_DIR = Path(__file__).resolve().parents[1]
if str(PROJECT_DIR) not in sys.path:
    sys.path.append(str(PROJECT_DIR))

from src.visualize_alignment_grn import (
    load_rmsd_cache,
    load_processed_structures,
    load_grn_table,
    extract_ca_coordinates_with_grn,
    apply_alignment_transformations,
    apply_membrane_orientation,
    create_interactive_opsin_visualization_extended,
)
from src.data_processing import load_opsin_property_data


def _select_structures(structures: Dict[str, dict], max_structures: int, reference_id: str) -> Dict[str, dict]:
    """Return a deterministic subset of structures that always keeps the reference."""
    if max_structures <= 0 or len(structures) <= max_structures:
        return structures

    ordered_ids = list(structures.keys())
    selected_ids = []

    if reference_id in structures:
        selected_ids.append(reference_id)

    for struct_id in ordered_ids:
        if struct_id == reference_id:
            continue
        selected_ids.append(struct_id)
        if len(selected_ids) >= max_structures:
            break

    # If reference was missing, fall back to the leading IDs
    if not selected_ids:
        selected_ids = ordered_ids[:max_structures]

    return {sid: structures[sid] for sid in selected_ids}


def create_subset_alignment_visualization(
    cache_dir: str = "opsin_output/cache",
    property_file: str = "property/mo_exp.csv",
    grn_file: str = "opsin_output/curated_grn.csv",
    output_file: str = "opsin_output/interactive_grn_alignment_subset.html",
    reference_id: str = "MerMAID1_model_0",
    max_structures: int = 30,
    show_membrane: bool = True,
    membrane_opacity: float = 0.05,
) -> None:
    """Create an HTML alignment figure limited to the requested number of structures."""
    cache_data = load_rmsd_cache(cache_dir)
    alignment_paths = cache_data.get('alignment_paths', {})
    print(f"Found {len(alignment_paths)} alignment paths")

    processed_structures = load_processed_structures(cache_dir)
    grn_df = load_grn_table(grn_file)

    property_data = None
    property_path = Path(property_file)
    if property_path.exists():
        try:
            property_result = load_opsin_property_data(property_path, processed_structures)
            if property_result and 'properties' in property_result:
                property_data = property_result['properties']
                print(f"Loaded property data for {len(property_data)} structures")
        except Exception as exc:  # pragma: no cover - just logging
            print(f"Failed to load property data: {exc}")

    print("\n=== Extracting CA Coordinates with GRN Mapping ===")
    structures = extract_ca_coordinates_with_grn(
        processed_structures,
        grn_df,
        chain_id='A',
        use_helix_only=True
    )

    if not structures:
        raise RuntimeError("No structures available after GRN extraction")

    filtered_structures = _select_structures(structures, max_structures, reference_id)
    print(f"Selected {len(filtered_structures)} structures for subset visualization")

    print("\n=== Applying Alignment Transformations ===")
    aligned_structures = apply_alignment_transformations(filtered_structures, alignment_paths, reference_id)

    print("\n=== Applying Membrane Orientation ===")
    oriented_structures = apply_membrane_orientation(aligned_structures, reference_id)

    print("\n=== Creating Subset Interactive Visualization ===")
    fig = create_interactive_opsin_visualization_extended(
        oriented_structures,
        grn_df,
        property_data=property_data,
        max_structures=max_structures,
        show_membrane=show_membrane,
        membrane_opacity=membrane_opacity,
        hover_show_residue_name=True,
        include_retinal=True,
        retinal_reference_id=reference_id,
        title="Opsin Structure Alignment (30-structure subset)"
    )

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(output_path)
    print(f"Subset interactive visualization saved to: {output_path}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Create a reduced HTML structure alignment figure")
    parser.add_argument('--cache-dir', default='opsin_output/cache', help='Workflow cache directory')
    parser.add_argument('--property-file', default='property/mo_exp.csv', help='Path to property CSV file')
    parser.add_argument('--grn-file', default='opsin_output/curated_grn.csv', help='Path to curated GRN table')
    parser.add_argument('--output-file', default='opsin_output/interactive_grn_alignment_subset.html',
                        help='Destination HTML file')
    parser.add_argument('--reference-id', default='MerMAID1_model_0', help='Reference structure for alignment')
    parser.add_argument('--max-structures', type=int, default=30, help='Maximum structures to include')
    parser.add_argument('--hide-membrane', action='store_true', help='Disable the membrane volume overlay')
    parser.add_argument('--membrane-opacity', type=float, default=0.05, help='Membrane opacity in the figure')
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    create_subset_alignment_visualization(
        cache_dir=args.cache_dir,
        property_file=args.property_file,
        grn_file=args.grn_file,
        output_file=args.output_file,
        reference_id=args.reference_id,
        max_structures=args.max_structures,
        show_membrane=not args.hide_membrane,
        membrane_opacity=args.membrane_opacity,
    )


if __name__ == '__main__':
    main()
