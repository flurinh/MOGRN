"""
Functions for loading and applying helix definitions to protein structures.

Helix definitions are loaded from property/helices_grn.json which is generated
by scripts/generate_helices_grn.py based on the postprocessed GRN table.
"""

import json
import os
from pathlib import Path

import pandas as pd
from Bio.Align import PairwiseAligner


PROJECT_ROOT = Path(__file__).resolve().parent.parent
HELIX_STRUCTURE_ALIASES = (
    PROJECT_ROOT / "src" / "resources" / "helix_structure_aliases.json"
)
STRUCTURE_CACHE = PROJECT_ROOT / "data" / "structure" / "cache"


def _load_helix_structure_aliases(path: Path = None) -> dict:
    """Load structural replacements used to transfer baseline helix labels."""

    path = Path(path or HELIX_STRUCTURE_ALIASES)
    if not path.is_file():
        return {}
    payload = json.loads(path.read_text(encoding="utf-8"))
    aliases = payload.get("aliases", {})
    if not isinstance(aliases, dict):
        raise ValueError("Helix structure aliases must be a JSON object")
    return aliases


def _residue_sequence(frame: pd.DataFrame, chain_id: str = "A"):
    """Return author positions and sequence for one structure chain."""

    residues = frame.reset_index(drop=True)
    if "auth_chain_id" in residues.columns:
        residues = residues.loc[
            residues["auth_chain_id"].astype(str).eq(str(chain_id))
        ]
    if "atom_name" in residues.columns:
        ca = residues.loc[residues["atom_name"].astype(str).str.upper().eq("CA")]
        if not ca.empty:
            residues = ca
    residues = residues.copy()
    residues["_auth_seq_numeric"] = pd.to_numeric(
        residues["auth_seq_id"], errors="coerce"
    )
    residues = (
        residues.dropna(subset=["_auth_seq_numeric"])
        .sort_values("_auth_seq_numeric")
        .drop_duplicates("_auth_seq_numeric")
    )
    positions = residues["_auth_seq_numeric"].astype(int).tolist()
    sequence = "".join(
        residues["res_name1l"].fillna("X").astype(str).str.strip().str.upper()
    )
    return positions, sequence


def _sequence_position_map(
    source_frame: pd.DataFrame,
    target_frame: pd.DataFrame,
) -> dict[int, int]:
    """Map source author positions to target positions by global sequence."""

    source_positions, source_sequence = _residue_sequence(source_frame)
    target_positions, target_sequence = _residue_sequence(target_frame)
    if not source_sequence or not target_sequence:
        return {}

    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -3.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(source_sequence, target_sequence)[0]

    mapping = {}
    for (source_start, source_end), (target_start, target_end) in zip(
        alignment.aligned[0], alignment.aligned[1]
    ):
        length = min(source_end - source_start, target_end - target_start)
        for offset in range(length):
            mapping[source_positions[source_start + offset]] = target_positions[
                target_start + offset
            ]
    return mapping


def _add_replacement_helix_definitions(
    processed_structures: dict,
    helix_definitions: dict,
) -> dict:
    """Transfer raw helix ranges to replacement structures by sequence.

    This is baseline workflow metadata only. It does not read or apply the
    final curated ProtOS GRN reference.
    """

    aliases = _load_helix_structure_aliases()
    lookup = {str(key).lower(): key for key in processed_structures}
    definition_lookup = {str(key).lower(): key for key in helix_definitions}

    for target_id, spec in aliases.items():
        target_key = lookup.get(str(target_id).lower())
        source_id = str(spec.get("source_structure_id", "")).strip()
        source_definition_key = definition_lookup.get(source_id.lower())
        if target_key is None or source_definition_key is None:
            continue

        source_path = STRUCTURE_CACHE / f"{source_id}.pkl"
        if not source_path.is_file():
            continue
        source_frame = pd.read_pickle(source_path)
        target_structure = processed_structures[target_key]
        target_frame = target_structure.get("df")
        if target_frame is None or target_frame.empty:
            continue

        position_map = _sequence_position_map(source_frame, target_frame)
        translated = {}
        for helix_num, bounds in helix_definitions[source_definition_key].items():
            if not isinstance(bounds, list) or len(bounds) != 2:
                continue
            positions = [
                position_map[position]
                for position in range(int(bounds[0]), int(bounds[1]) + 1)
                if position in position_map
            ]
            if positions:
                translated[str(helix_num)] = [min(positions), max(positions)]
        if len(translated) == 7:
            helix_definitions[target_key] = translated
            print(
                f"[INFO] Transferred helix annotations "
                f"{source_definition_key} -> {target_key}"
            )
    return helix_definitions

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

    helix_definitions = _add_replacement_helix_definitions(
        processed_structures, helix_definitions
    )

    # Apply annotations
    processed_structures = apply_helix_annotations(processed_structures, helix_definitions)

    return {
        'processed_structures': processed_structures,
        'helix_definitions': helix_definitions,
        'helix_annotations_file': str(helix_file) if helix_file else 'property/helices_grn.json'
    }
