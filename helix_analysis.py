"""
Functions for identifying, annotating, and analyzing helices in protein structures.
These functions handle helix identification, annotation, and definitions for membrane proteins.
"""

import os
import pandas as pd
import numpy as np
import json
from pathlib import Path

from projects.opsin_analysis.data_processing import  ensure_structure_dtypes


def define_reference_helices(reference_structure, helix_ref_file=None):
    """
    Define the helix boundaries for the reference structure.
    This can either load from a file or use hardcoded values.
    
    Args:
        reference_structure: The reference structure data
        helix_ref_file: Optional path to a JSON file containing helix definitions
        
    Returns:
        Dictionary of helix definitions {helix_num: {'start': pos, 'end': pos}}
    """
    import json
    from pathlib import Path
    
    # First, try to load from file if provided
    if helix_ref_file and Path(helix_ref_file).exists():
        try:
            with open(helix_ref_file, 'r') as f:
                helix_data = json.load(f)
                
            # Check if the reference structure ID is in the helix data
            ref_id = next(iter(helix_data.keys()))
            
            # Convert to standard format
            helix_defs = {}
            for helix_num, bounds in helix_data[ref_id].items():
                if isinstance(bounds, list) and len(bounds) == 2:
                    helix_defs[helix_num] = {'start': bounds[0], 'end': bounds[1]}
            
            if helix_defs:
                print(f"[INFO] Loaded helix definitions from {helix_ref_file}: {helix_defs}")
                return helix_defs
        except Exception as e:
            print(f"[WARNING] Error loading helix definitions from file: {e}")
    
    # Fallback: Use hardcoded definitions for CnChR2_J230_refine9
    # These values are taken from the JSON file we fixed earlier
    ref_id = reference_structure.get('pdb_id', 'unknown')
    print(f"[INFO] Using hardcoded helix definitions for {ref_id}")
    
    if ref_id == 'CnChR2_J230_refine9':
        return {
            '1': {'start': 88, 'end': 111},
            '2': {'start': 120, 'end': 136},
            '3': {'start': 157, 'end': 174},
            '4': {'start': 188, 'end': 206},
            '5': {'start': 214, 'end': 232},
            '6': {'start': 250, 'end': 269},
            '7': {'start': 285, 'end': 304}
        }
    
    # If no valid definitions, return empty dict
    print(f"[WARNING] No helix definitions available for {ref_id}")
    return {}


def annotate_helices_from_alignments(processed_structures, reference_id, helix_definitions, 
                                     alignment_paths, chain_id='A'):
    """
    Annotate helices in all structures based on alignments to reference structure.
    
    Args:
        processed_structures: Dictionary of processed structures
        reference_id: ID of reference structure with defined helices
        helix_definitions: Dictionary of helix boundaries in reference
        alignment_paths: Paths from structure_comparison.get_structure_alignment
        chain_id: Chain ID to use
        
    Returns:
        Updated processed_structures with helix annotations
    """
    print(f"[INFO] Annotating helices based on reference structure: {reference_id}")
    
    if reference_id not in processed_structures:
        print(f"[ERROR] Reference structure {reference_id} not found in processed structures")
        return processed_structures
    
    # Get reference structure data
    ref_structure = processed_structures[reference_id]
    
    # Extract mapping from auth_seq_id to position for reference structure
    ref_df = None
    if 'df_norm' in ref_structure:
        ref_df = ref_structure['df_norm']
    elif 'df' in ref_structure:
        ref_df = ref_structure['df']
    
    if ref_df is None or ref_df.empty:
        print(f"[ERROR] No dataframe found for reference structure {reference_id}")
        return processed_structures
    
    # Filter for the specified chain
    ref_df = ref_df[ref_df['auth_chain_id'] == chain_id]
    
    # Filter for CA atoms
    ref_ca_df = ref_df[ref_df['res_atom_name'] == 'CA']
    if ref_ca_df.empty:
        try:
            # Try string comparison
            ref_ca_df = ref_df[ref_df['res_atom_name'].astype(str) == 'CA']
        except:
            pass
    
    if ref_ca_df.empty:
        print(f"[ERROR] No CA atoms found in reference structure {reference_id}")
        return processed_structures
    
    # Create mapping from auth_seq_id to indices for reference
    ref_seq_ids = ref_ca_df['auth_seq_id'].values
    ref_seq_map = {seq_id: i for i, seq_id in enumerate(ref_seq_ids)}
    
    # Create helix assignments for the reference structure
    helices_by_residue = {}
    for helix_num, helix_info in helix_definitions.items():
        start_pos = helix_info['start']
        end_pos = helix_info['end']
        
        # Assign all residues in range to this helix
        for res_id in range(start_pos, end_pos + 1):
            helices_by_residue[res_id] = int(helix_num)
    
    # Store helix assignments in the reference structure
    ref_structure['helix_assignments'] = helices_by_residue
    ref_structure['helix_definitions'] = helix_definitions
    
    # Now annotate all other structures using alignment paths
    for struct_id, structure in processed_structures.items():
        if struct_id == reference_id:
            continue  # Skip reference, already annotated
        
        # Get alignment path between reference and this structure
        if (reference_id, struct_id) in alignment_paths:
            path_info = alignment_paths[(reference_id, struct_id)]
        elif (struct_id, reference_id) in alignment_paths:
            # If found in reverse order, we need to flip the mapping
            path_info = alignment_paths[(struct_id, reference_id)]
            # Flip the residue mapping (B->A becomes A->B)
            if 'residue_mapping' in path_info:
                path_info['residue_mapping'] = [(b, a) for a, b in path_info['residue_mapping']]
        else:
            print(f"[WARNING] No alignment path found between {reference_id} and {struct_id}")
            continue
        
        # Extract residue mapping
        if 'residue_mapping' not in path_info:
            print(f"[WARNING] No residue mapping in alignment path for {struct_id}")
            continue
            
        residue_mapping = path_info['residue_mapping']
        
        # Create reverse mapping (reference -> structure)
        target_helix_assignments = {}
        
        for ref_res_id, target_res_id in residue_mapping:
            # Check if the reference residue has a helix assignment
            if ref_res_id in helices_by_residue:
                helix_num = helices_by_residue[ref_res_id]
                target_helix_assignments[target_res_id] = helix_num
        
        # Store the helix assignments in the structure
        structure['helix_assignments'] = target_helix_assignments
        
        # Create helix definitions for the target structure
        target_helix_defs = {}
        for helix_num in range(1, 8):  # Assume 7TM helices
            helix_residues = [res_id for res_id, h_num in target_helix_assignments.items() 
                             if h_num == helix_num]
            if helix_residues:
                target_helix_defs[str(helix_num)] = {
                    'start': min(helix_residues),
                    'end': max(helix_residues)
                }
        
        structure['helix_definitions'] = target_helix_defs
        print(f"[INFO] Annotated {len(target_helix_assignments)} residues with helix assignments in {struct_id}")
        
        # Add a tm_helices list for compatibility with visualization
        structure['tm_helices'] = list(range(1, 8))
    
    return processed_structures


def orient_and_annotate_structures(data_dict, output_dir='output', visualize=True):
    """
    Step 3 & 4: Orient structures and annotate helices using alignments

    Instead of using PCA-based orientation and algorithmic helix finding,
    this function uses structure alignments to transfer helix annotations
    from a reference structure to all other structures.

    Args:
        data_dict: Dictionary with data from previous steps
        output_dir: Directory to save output files
        visualize: Whether to generate visualizations

    Returns:
        Dictionary with oriented and annotated structures
    """
    print("\n=== Step 3 & 4: Annotating Helices Using Alignments ===")
    processed_structures = data_dict['processed_structures']
    alignment_paths = data_dict.get('alignment_paths', {})

    # Check if we have structures and alignment paths
    if not processed_structures:
        print("[ERROR] No structures available for annotation.")
        return {
            'processed_structures': processed_structures
        }

    if not alignment_paths:
        print("[ERROR] No alignment paths available. Make sure to run structure comparison first.")
        print("[ERROR] Falling back to old approach.")

        # Fallback to old approach (not recommended)
        # This is just a placeholder, the old code is removed as requested
        return {
            'processed_structures': processed_structures
        }

    # Ensure all structures have correct data types
    print("[INFO] Ensuring correct data types for all structures...")
    processed_structures = ensure_structure_dtypes(processed_structures)

    # Extract CA atoms for each structure (needed for helix annotation)
    for pdb_id, data in processed_structures.items():
        if 'df_norm' not in data or data['df_norm'].empty:
            if 'df' in data and not data['df'].empty:
                df_norm = data['df'].copy()
                data['df_norm'] = df_norm
            else:
                print(f"[WARNING] {pdb_id}: No structure data available.")
                continue

        df_norm = data['df_norm']
        try:
            # Try direct comparison first
            df_ca_norm = df_norm[df_norm['res_atom_name'] == 'CA'].copy()
            # If we didn't get any CA atoms, try string conversion
            if df_ca_norm.empty:
                print(f"[DEBUG] No CA atoms found with direct comparison for {pdb_id}, trying string conversion")
                df_ca_norm = df_norm[df_norm['res_atom_name'].astype(str) == 'CA'].copy()
        except TypeError:
            # If TypeError occurs, use string conversion
            print(f"[DEBUG] TypeError when extracting CA atoms for {pdb_id}, using string conversion")
            df_ca_norm = df_norm[df_norm['res_atom_name'].astype(str) == 'CA'].copy()

        if df_ca_norm.empty:
            print(f"[WARNING] No CA atoms found for {pdb_id} even after string conversion!")

        data['df_ca_norm'] = df_ca_norm

    # Choose reference structure for helix definitions
    # 1. Try to use the global reference if available
    reference_id = data_dict.get('global_ref')

    # 2. If not available, use CnChR2_J230_refine9 if present
    if not reference_id or reference_id not in processed_structures:
        if 'CnChR2_J230_refine9' in processed_structures:
            reference_id = 'CnChR2_J230_refine9'
        else:
            # 3. Otherwise use the first structure
            reference_id = next(iter(processed_structures.keys()))

    print(f"[INFO] Using {reference_id} as reference structure for helix annotation")

    # Get helix definitions for reference structure
    helix_ref_file = os.path.join(os.path.dirname(__file__), 'property', 'helix_ref_CnChR2_J230_refine9.json')
    helix_definitions = define_reference_helices(
        processed_structures[reference_id],
        helix_ref_file=helix_ref_file
    )

    if not helix_definitions:
        print("[ERROR] No helix definitions available. Cannot proceed with annotation.")
        return {
            'processed_structures': processed_structures,
            'reference_structure': reference_id
        }

    # Annotate helices in all structures based on reference
    processed_structures = annotate_helices_from_alignments(
        processed_structures,
        reference_id,
        helix_definitions,
        alignment_paths
    )

    return {
        'processed_structures': processed_structures,
        'reference_structure': reference_id,
        'helix_definitions': helix_definitions
    }


def align_to_reference_and_annotate_helices(data_dict, output_dir='output', visualize=True):
    """
    Step 4: Custom step that aligns all structures to a reference structure
    (from helix_ref_CnChR2_J230_refine9.json) and annotates all structures
    with helix numbers based on the reference helices.

    If helices.json already exists in property directory, it will load helix
    definitions from there instead of recalculating alignments.

    Args:
        data_dict: Dictionary with data from previous steps
        output_dir: Directory to save output files
        visualize: Whether to generate visualizations

    Returns:
        Dictionary with aligned structures and helix annotations
    """
    print("\n=== Step 4: Aligning to Reference and Annotating Helices ===")

    import json
    import os
    import numpy as np
    import pandas as pd
    from pathlib import Path
    from protos.processing.structure.struct_alignment import get_structure_alignment

    processed_structures = data_dict['processed_structures']

    # Check if we have structures to align
    if not processed_structures:
        print("[ERROR] No structures available for alignment.")
        return {
            'processed_structures': processed_structures
        }

    # Define paths for helix files
    property_dir = os.path.join(os.path.dirname(__file__), 'property')
    helix_ref_file = os.path.join(property_dir, 'helix_ref_CnChR2_J230_refine9.json')
    helix_cache_file = os.path.join(property_dir, 'helices.json')

    # Initialize variables
    global_helix_annotations = {}
    alignment_paths = {}
    formatted_helix_defs = {}
    ref_id = None

    # First, check if we already have helices.json with annotations for existing structures
    if os.path.exists(helix_cache_file):
        try:
            print(f"[INFO] Found existing helix definitions at {helix_cache_file}")

            # Check file size to detect potential truncation
            file_size = os.path.getsize(helix_cache_file)
            if file_size < 100:  # Suspiciously small for a JSON with structure definitions
                print(f"[WARNING] Helix cache file is suspiciously small ({file_size} bytes)")
                print(f"[WARNING] This may indicate a truncated or corrupted file")
                raise ValueError("Suspicious file size indicates potential corruption")

            # Try to read the first few lines to check format
            with open(helix_cache_file, 'r') as f:
                # Read first 1000 characters to check structure
                file_start = f.read(1000)
                if not file_start.strip().startswith('{') or '}' not in file_start:
                    print(f"[WARNING] Helix cache file does not appear to be valid JSON")
                    raise ValueError("Invalid JSON format detected")
                # Reset file pointer
                f.seek(0)

                # Try to load the full JSON
                try:
                    global_helix_annotations = json.load(f)
                except json.JSONDecodeError as json_err:
                    print(f"[WARNING] JSON parsing error: {json_err}")

                    # Try to read the whole file to get more info
                    f.seek(0)
                    file_content = f.read()
                    last_complete_brace = file_content.rfind('}')

                    if last_complete_brace > 0:
                        print("[INFO] Attempting to recover partial JSON data...")
                        try:
                            # Try to repair by closing the JSON properly
                            repaired_json = file_content[:last_complete_brace + 1]
                            global_helix_annotations = json.loads(repaired_json)
                            print(
                                f"[INFO] Successfully recovered {len(global_helix_annotations)} structure definitions")
                        except:
                            print("[WARNING] Failed to recover JSON data")
                            raise ValueError("Could not repair corrupted JSON file")
                    else:
                        raise ValueError("JSON file appears severely corrupted")

            # Extract reference structure ID - assume it's the first one
            existing_structures = list(global_helix_annotations.keys())
            if existing_structures:
                # Try to find the reference structure (CnChR2_J230_refine9 or similar)
                ref_candidates = [s for s in existing_structures if 'CnChR2_J230_refine9' in s]
                if ref_candidates:
                    ref_id = ref_candidates[0]
                else:
                    # Just use the first one
                    ref_id = existing_structures[0]

                print(f"[INFO] Using {ref_id} as reference structure from cached helix definitions")

                # Format helix definitions
                helix_definitions = global_helix_annotations[ref_id]
                for helix_num, bounds in helix_definitions.items():
                    if isinstance(bounds, list) and len(bounds) == 2:
                        formatted_helix_defs[helix_num] = {'start': bounds[0], 'end': bounds[1]}

                # Validate that we got all 7 helices for the reference structure
                if len(formatted_helix_defs) < 7:
                    print(
                        f"[WARNING] Reference structure only has {len(formatted_helix_defs)} helices defined, expected 7")
                    print(f"[WARNING] This may indicate incomplete or corrupted data")

                    # If we have fewer than 4 helices, consider the data unreliable
                    if len(formatted_helix_defs) < 4:
                        print(f"[WARNING] Too few helices defined. Will recalculate all helix definitions")
                        raise ValueError("Insufficient helix definitions in reference structure")
            else:
                print(f"[WARNING] No structures found in helix cache file")
                raise ValueError("Empty structure list in helix cache")

            # Check if we have definitions for all current structures
            current_struct_ids = set(processed_structures.keys())
            cached_struct_ids = set(global_helix_annotations.keys())
            missing_structs = current_struct_ids - cached_struct_ids

            # Check if any cached structures have incomplete definitions
            incomplete_structs = []
            for sid, struct_def in global_helix_annotations.items():
                if sid not in processed_structures:
                    continue  # Skip structures we don't need

                # Check if this structure has a reasonable number of helices
                if len(struct_def) < 4:  # Expect at least 4 of 7 helices
                    incomplete_structs.append(sid)

            if incomplete_structs:
                print(f"[WARNING] Found {len(incomplete_structs)} structures with incomplete helix definitions")
                print(f"[WARNING] These will be recalculated: {', '.join(incomplete_structs[:5])}")
                if len(incomplete_structs) > 5:
                    print(f"[WARNING] ... and {len(incomplete_structs) - 5} more")
                missing_structs.update(incomplete_structs)

            if missing_structs:
                print(f"[INFO] Found {len(missing_structs)} structures that need helix definitions")
                print(f"[INFO] Will calculate helix definitions for: {', '.join(list(missing_structs)[:5])}")
                if len(missing_structs) > 5:
                    print(f"[INFO] ... and {len(missing_structs) - 5} more")

                # We'll need to load the reference helix definitions and calculate alignments
                # for the missing structures
                need_alignment = True
            else:
                print(f"[INFO] Found helix definitions for all {len(current_struct_ids)} structures")
                need_alignment = False

        except Exception as e:
            print(f"[WARNING] Error loading cached helix definitions: {e}")
            print(f"[INFO] Will recalculate helix definitions for all structures")

            # Make a backup of the problematic file
            if os.path.exists(helix_cache_file):
                backup_file = f"{helix_cache_file}.bak"
                try:
                    import shutil
                    shutil.copy2(helix_cache_file, backup_file)
                    print(f"[INFO] Created backup of problematic helix cache at {backup_file}")
                except Exception as backup_err:
                    print(f"[WARNING] Failed to create backup: {backup_err}")

            # Reset to empty dictionary and recalculate
            global_helix_annotations = {}
            need_alignment = True
    else:
        print(f"[INFO] No cached helix definitions found at {helix_cache_file}")
        print(f"[INFO] Will calculate helix definitions for all structures")
        need_alignment = True

    # If we need to calculate alignments for all or some structures
    if need_alignment:
        # Load reference helix definitions if not already loaded
        if not ref_id or not formatted_helix_defs:
            if not os.path.exists(helix_ref_file):
                print(f"[ERROR] Reference helix file not found: {helix_ref_file}")
                return {
                    'processed_structures': processed_structures
                }

            try:
                with open(helix_ref_file, 'r') as f:
                    helix_data = json.load(f)

                # Get the reference structure ID
                ref_id = next(iter(helix_data.keys()))
                helix_definitions = helix_data[ref_id]

                print(f"[INFO] Loaded helix definitions for reference structure {ref_id}")

                # Format helix definitions as dictionary of {helix_num: {'start': pos, 'end': pos}}
                formatted_helix_defs = {}
                for helix_num, bounds in helix_definitions.items():
                    if isinstance(bounds, list) and len(bounds) == 2:
                        formatted_helix_defs[helix_num] = {'start': bounds[0], 'end': bounds[1]}

                if not formatted_helix_defs:
                    print("[ERROR] Invalid helix definition format.")
                    return {
                        'processed_structures': processed_structures
                    }

                # Add reference structure helices to global annotations
                global_helix_annotations[ref_id] = helix_definitions

            except Exception as e:
                print(f"[ERROR] Failed to load helix reference file: {e}")
                return {
                    'processed_structures': processed_structures
                }

        # Check if reference structure exists in our dataset
        if ref_id not in processed_structures:
            print(f"[WARNING] Reference structure {ref_id} not found in processed structures.")
            print(f"[WARNING] Will try to find a structure with similar name.")

            # Try to find a structure with a similar name
            potential_ref_ids = [sid for sid in processed_structures.keys() if ref_id in sid]
            if potential_ref_ids:
                # Use the first match
                ref_id = potential_ref_ids[0]
                print(f"[INFO] Using {ref_id} as reference structure.")
            else:
                print(f"[ERROR] Cannot find suitable reference structure. Aborting helix annotation.")
                return {
                    'processed_structures': processed_structures
                }

        # Ensure all structures have correct data types
        print("[INFO] Ensuring correct data types for all structures...")
        for pdb_id, data in processed_structures.items():
            if 'df' in data:
                df = data['df']
                for col in ['x', 'y', 'z']:
                    if col in df.columns:
                        try:
                            df[col] = pd.to_numeric(df[col], errors='coerce')
                        except Exception as e:
                            print(f"[WARNING] Error converting {col} in {pdb_id}: {e}")

                # Create df_norm if it doesn't exist
                if 'df_norm' not in data or data['df_norm'].empty:
                    data['df_norm'] = df.copy()

                # Extract CA atoms for alignment
                try:
                    df_ca = df[df['res_atom_name'] == 'CA'].copy()
                    if df_ca.empty:
                        # Try string conversion
                        df_ca = df[df['res_atom_name'].astype(str) == 'CA'].copy()
                    data['df_ca'] = df_ca
                except Exception as e:
                    print(f"[WARNING] Error extracting CA atoms in {pdb_id}: {e}")

        # Get reference structure data
        ref_structure = processed_structures[ref_id]

        # Extract reference CA atoms
        if 'df_ca' not in ref_structure or ref_structure['df_ca'].empty:
            print(f"[ERROR] No CA atoms found in reference structure {ref_id}")
            return {
                'processed_structures': processed_structures
            }

        ref_ca_df = ref_structure['df_ca']
        ref_ca_coords = ref_ca_df[['x', 'y', 'z']].astype(float).values
        ref_seq_ids = ref_ca_df['auth_seq_id'].values if 'auth_seq_id' in ref_ca_df.columns else list(
            range(len(ref_ca_coords)))

        # Determine which structures need alignment
        structures_to_align = []
        for struct_id in processed_structures.keys():
            if struct_id == ref_id:
                continue  # Skip reference structure

            # Check if structure already has helix definitions in the cache
            if struct_id in global_helix_annotations:
                # Verify the definitions are complete (all 7 helices)
                helix_defs = global_helix_annotations[struct_id]
                if len(helix_defs) < 7:
                    print(
                        f"[INFO] Structure {struct_id} has incomplete helix definitions ({len(helix_defs)}/7 helices)")
                    structures_to_align.append(struct_id)
                else:
                    print(f"[INFO] Using cached helix definitions for {struct_id}")
            else:
                # Structure not in cache, needs alignment
                structures_to_align.append(struct_id)

        # Align structures that don't have helix definitions yet
        if structures_to_align:
            print(f"[INFO] Aligning {len(structures_to_align)} structures to reference {ref_id} and mapping helices...")

            for struct_id in structures_to_align:
                structure = processed_structures[struct_id]

                # Skip structures without CA atoms
                if 'df_ca' not in structure or structure['df_ca'].empty:
                    print(f"[WARNING] Skipping {struct_id} - No CA atoms found.")
                    continue

                # Get target structure CA atoms
                target_ca_df = structure['df_ca']
                target_ca_coords = target_ca_df[['x', 'y', 'z']].astype(float).values
                target_seq_ids = target_ca_df['auth_seq_id'].values if 'auth_seq_id' in target_ca_df.columns else list(
                    range(len(target_ca_coords)))

                try:
                    # Perform structure alignment
                    rotation, translation, best_path, rmsd = get_structure_alignment(ref_ca_coords, target_ca_coords)

                    # Extract indices from alignment path
                    ref_indices, target_indices = best_path

                    # Map indices to auth_seq_id values
                    ref_res_ids = [ref_seq_ids[idx] for idx in ref_indices]
                    target_res_ids = [target_seq_ids[idx] for idx in target_indices]

                    # Create mapping between reference and target residue IDs
                    residue_mapping = list(zip(ref_res_ids, target_res_ids))

                    # Store alignment path
                    alignment_paths[(ref_id, struct_id)] = {
                        'rotation': rotation.tolist(),
                        'translation': translation.tolist(),
                        'coord_indices': best_path,
                        'residue_mapping': residue_mapping,
                        'rmsd': rmsd
                    }

                    # Now map helix definitions from reference to target
                    target_helix_defs = {}
                    target_residue_helices = {}

                    # For each helix in reference, find corresponding residues in target
                    for helix_num, bounds in formatted_helix_defs.items():
                        ref_start = bounds['start']
                        ref_end = bounds['end']

                        # Find all aligned residues in this helix
                        helix_mappings = []
                        for ref_res, target_res in residue_mapping:
                            if ref_start <= ref_res <= ref_end:
                                helix_mappings.append((ref_res, target_res))
                                target_residue_helices[target_res] = int(helix_num)

                        # If we found mappings for this helix, create target helix definition
                        if helix_mappings:
                            target_res_in_helix = [t for _, t in helix_mappings]
                            target_helix_defs[helix_num] = {
                                'start': min(target_res_in_helix),
                                'end': max(target_res_in_helix)
                            }

                    # Add to global annotations
                    global_helix_annotations[struct_id] = {
                        str(h_num): [bounds['start'], bounds['end']]
                        for h_num, bounds in target_helix_defs.items()
                    }

                    print(f"[INFO] Successfully aligned and annotated {struct_id} (RMSD: {rmsd:.3f}Å)")

                except Exception as e:
                    print(f"[ERROR] Failed to align and annotate {struct_id}: {e}")

        # Save global helix annotations to property directory
        try:
            # Make sure all helix definitions are serializable
            clean_helix_annotations = {}
            for struct_id, helix_defs in global_helix_annotations.items():
                # Create a clean copy with only strings and simple types
                clean_struct_defs = {}
                for helix_id, bounds in helix_defs.items():
                    # Convert bounds to simple list with integers if needed
                    if isinstance(bounds, list) and len(bounds) == 2:
                        # Ensure bounds are integers
                        clean_bounds = [int(bounds[0]), int(bounds[1])]
                        clean_struct_defs[str(helix_id)] = clean_bounds
                # Only add structures with at least one helix definition
                if clean_struct_defs:
                    clean_helix_annotations[struct_id] = clean_struct_defs

            # First, write to a temporary file to avoid corruption
            temp_file = f"{helix_cache_file}.tmp"
            with open(temp_file, 'w') as f:
                json.dump(clean_helix_annotations, f, indent=2)

            # Now rename the temporary file to the actual file (atomic operation)
            import os
            if os.path.exists(helix_cache_file):
                # Create a backup of the current file first
                backup_file = f"{helix_cache_file}.bak"
                try:
                    import shutil
                    shutil.copy2(helix_cache_file, backup_file)
                except Exception as backup_err:
                    print(f"[WARNING] Failed to create backup: {backup_err}")

            # Rename temp file to actual file name
            import shutil
            shutil.move(temp_file, helix_cache_file)

            print(f"[INFO] Saved helix annotations for {len(clean_helix_annotations)} structures to {helix_cache_file}")
        except Exception as e:
            print(f"[WARNING] Failed to save helix annotations: {e}")
            import traceback
            traceback.print_exc()

    # Now apply helix annotations to all structures
    print("[INFO] Applying helix annotations to all structures...")

    for struct_id, structure in processed_structures.items():
        # Skip structures without helix definitions
        if struct_id not in global_helix_annotations:
            print(f"[WARNING] No helix definition found for {struct_id}")

            # Instead of skipping, initialize with empty helix data
            # This ensures all structures have the expected fields
            structure['helix_definitions'] = {}
            structure['residue_to_helix'] = {}
            structure['tm_helices'] = []

            # Initialize dataframes with helix_num = 0
            for df_key in ['df', 'df_norm', 'df_ca', 'df_ca_norm']:
                if df_key in structure and not structure[df_key].empty:
                    if 'helix_num' not in structure[df_key].columns:
                        structure[df_key]['helix_num'] = 0

            continue

        # Get helix definitions for this structure
        struct_helix_defs = global_helix_annotations[struct_id]

        # Convert to internal format
        formatted_struct_helix_defs = {}
        target_residue_helices = {}

        for helix_num, bounds in struct_helix_defs.items():
            if isinstance(bounds, list) and len(bounds) == 2:
                formatted_struct_helix_defs[helix_num] = {'start': bounds[0], 'end': bounds[1]}

                # Create mapping from residue ID to helix number
                start_pos = bounds[0]
                end_pos = bounds[1]
                for res_id in range(start_pos, end_pos + 1):
                    target_residue_helices[res_id] = int(helix_num)

        # Store helix definitions and mapping in structure data
        structure['helix_definitions'] = formatted_struct_helix_defs
        structure['residue_to_helix'] = target_residue_helices

        # Also store as list for compatibility with visualization functions
        structure['tm_helices'] = list(range(1, 8))

        # Process every dataframe in the structure to ensure consistent helix annotations
        for df_key in ['df', 'df_norm', 'df_ca', 'df_ca_norm']:
            if df_key in structure and not structure[df_key].empty:
                df = structure[df_key]

                # Initialize helix_num column to 0 (no helix)
                if 'helix_num' in df.columns:
                    # Reset existing column
                    df['helix_num'] = 0
                else:
                    # Create new column
                    df['helix_num'] = 0

                # Update helix_num based on mapped residues
                if 'auth_seq_id' in df.columns:
                    # Use vectorized operations when possible for efficiency
                    if len(target_residue_helices) > 0:
                        # Create a mapping series for faster lookup
                        import pandas as pd
                        import numpy as np

                        # Extract all residue IDs and convert to numeric if needed
                        residue_ids = df['auth_seq_id'].unique()

                        # Update each residue's helix number
                        for res_id, helix_num in target_residue_helices.items():
                            # Apply the helix number to all atoms in this residue
                            mask = df['auth_seq_id'] == res_id
                            df.loc[mask, 'helix_num'] = int(helix_num)

                # Store updated dataframe back in the structure
                structure[df_key] = df

        # Also update df_ret if it exists to maintain consistency
        if 'df_ret' in structure and not structure['df_ret'].empty:
            # Retinal is not part of helices, but should have helix_num column for consistency
            if 'helix_num' not in structure['df_ret'].columns:
                structure['df_ret']['helix_num'] = 0

        # Debugging print statement showing helix details for this structure
        unique_helix_numbers = sorted(list(set([int(h) for h in formatted_struct_helix_defs.keys()])))
        helix_info = f"{struct_id}: {len(formatted_struct_helix_defs)} helices => {unique_helix_numbers}"

        # Add details for each helix
        for h_num in range(1, 8):
            h_str = str(h_num)
            if h_str in formatted_struct_helix_defs:
                start = formatted_struct_helix_defs[h_str]['start']
                end = formatted_struct_helix_defs[h_str]['end']
                helix_info += f" | Helix {h_num}: {start}-{end}"
            else:
                helix_info += f" | Helix {h_num}: missing"

        print(f"[DEBUG] {helix_info}")

    return {
        'processed_structures': processed_structures,
        'reference_structure': ref_id,
        'helix_definitions': formatted_helix_defs,
        'alignment_paths': alignment_paths,
        'helix_annotations_file': helix_cache_file
    }
