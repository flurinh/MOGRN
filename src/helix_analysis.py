"""
Functions for loading and applying helix definitions to protein structures.

Helix definitions are loaded from property/helices_grn.json which is generated
by scripts/generate_helices_grn.py based on the postprocessed GRN table.
"""

import json
import os
from pathlib import Path


def load_helix_definitions(helix_file: str = None) -> dict:
    """
    Load helix definitions from JSON file.

    Args:
        helix_file: Path to helices JSON file. If None, uses property/helices_grn.json

    Returns:
        Dictionary mapping structure_id -> {helix_num: [start, end]}
    """
    if helix_file is None:
        helix_file = Path(__file__).parent.parent / 'property' / 'helices_grn.json'

    helix_file = Path(helix_file)

    if not helix_file.exists():
        print(f"[ERROR] Helix definitions file not found: {helix_file}")
        print("[INFO] Run scripts/generate_helices_grn.py to generate helix definitions")
        return {}

    with open(helix_file, 'r') as f:
        helix_data = json.load(f)

    print(f"[INFO] Loaded helix definitions for {len(helix_data)} structures from {helix_file}")
    return helix_data


def apply_helix_annotations(processed_structures: dict, helix_definitions: dict) -> dict:
    """
    Apply helix annotations to processed structures.

    Args:
        processed_structures: Dictionary of processed structures
        helix_definitions: Dictionary mapping structure_id -> {helix_num: [start, end]}

    Returns:
        Updated processed_structures with helix annotations
    """
    print(f"[INFO] Applying helix annotations to {len(processed_structures)} structures...")

    annotated_count = 0

    for struct_id, structure in processed_structures.items():
        # Try to find helix definitions (check lowercase)
        struct_id_lower = struct_id.lower()

        if struct_id in helix_definitions:
            struct_helix_defs = helix_definitions[struct_id]
        elif struct_id_lower in helix_definitions:
            struct_helix_defs = helix_definitions[struct_id_lower]
        else:
            # No helix definitions for this structure
            structure['helix_definitions'] = {}
            structure['residue_to_helix'] = {}
            structure['tm_helices'] = []
            continue

        # Convert to internal format and build residue mapping
        formatted_helix_defs = {}
        residue_to_helix = {}

        for helix_num, bounds in struct_helix_defs.items():
            if isinstance(bounds, list) and len(bounds) == 2:
                start, end = bounds[0], bounds[1]
                formatted_helix_defs[helix_num] = {'start': start, 'end': end}

                # Map each residue to its helix
                for res_id in range(start, end + 1):
                    residue_to_helix[res_id] = int(helix_num)

        # Store in structure
        structure['helix_definitions'] = formatted_helix_defs
        structure['residue_to_helix'] = residue_to_helix
        structure['tm_helices'] = list(range(1, 8))

        # Add helix_num column to dataframes
        for df_key in ['df', 'df_norm', 'df_ca', 'df_ca_norm']:
            if df_key in structure and structure[df_key] is not None and not structure[df_key].empty:
                df = structure[df_key]

                # Initialize helix_num column
                df['helix_num'] = 0

                # Apply helix numbers based on residue mapping
                if 'auth_seq_id' in df.columns:
                    for res_id, helix_num in residue_to_helix.items():
                        mask = df['auth_seq_id'] == res_id
                        df.loc[mask, 'helix_num'] = helix_num

                structure[df_key] = df

        annotated_count += 1

    print(f"[INFO] Applied helix annotations to {annotated_count} structures")
    return processed_structures


def load_and_apply_helix_annotations(data_dict: dict, helix_file: str = None) -> dict:
    """
    Load helix definitions and apply them to processed structures.

    This is the main entry point for helix annotation in the workflow.
    Helix definitions come from property/helices_grn.json which is generated
    by scripts/generate_helices_grn.py based on the postprocessed GRN table.

    Args:
        data_dict: Dictionary with 'processed_structures' key
        helix_file: Optional path to helices JSON file

    Returns:
        Dictionary with updated processed_structures and helix metadata
    """
    print("\n=== Loading and Applying Helix Annotations ===")

    processed_structures = data_dict.get('processed_structures', {})

    if not processed_structures:
        print("[ERROR] No structures available for annotation.")
        return {'processed_structures': processed_structures}

    # Load helix definitions
    helix_definitions = load_helix_definitions(helix_file)

    if not helix_definitions:
        print("[WARNING] No helix definitions loaded. Structures will not have helix annotations.")
        return {'processed_structures': processed_structures}

    # Apply annotations
    processed_structures = apply_helix_annotations(processed_structures, helix_definitions)

    return {
        'processed_structures': processed_structures,
        'helix_definitions': helix_definitions,
        'helix_annotations_file': str(helix_file) if helix_file else 'property/helices_grn.json'
    }
