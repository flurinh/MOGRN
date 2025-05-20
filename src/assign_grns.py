import pandas as pd
import os
import matplotlib.pyplot as plt

from src.reference_alignment import (
    find_type_references,
    find_global_reference,
    create_seq_alignment_dicts_from_paths
)

from src.tree_based_alignment import (
    align_and_assign_grn_tree_based,
    build_similarity_tree,
    create_guide_tree,
    generate_transitive_alignment_paths,
    create_tree_based_seq_alignment_dicts
)

from src.msa_grn import (
    analyze_residue_composition,
    generate_grn_msa_tables
)


def align_and_assign_grn(data_dict, output_dir='output', visualize=True):
    """
    Step 6: Structure alignment and GRN assignment

    This function assigns Generic Residue Numbers (GRN) using cached alignment paths
    without recalculating alignments, and also runs the tree-based alignment method
    as an additional approach.

    Args:
        data_dict: Dictionary with data from previous steps
        output_dir: Directory to save output files
        visualize: Whether to generate visualizations

    Returns:
        Dictionary with alignment data
    """
    print("\n=== Step 6: Structure Alignment & GRN Assignment ===")

    processed_structures_complete = data_dict['processed_structures']
    alignment_paths = data_dict['alignment_paths']

    rmsd_df = data_dict.get('rmsd_df', pd.DataFrame())
    pdb_list = data_dict.get('pdb_list', [])
    group_dict = data_dict.get('group_dict', {})

    # Check if we have structures to align
    if not processed_structures_complete or len(processed_structures_complete) < 2:
        print(f"Need at least 2 structures for alignment. Found {len(processed_structures_complete)}.")
        # Return empty data structures
        return {
            'processed_structures': processed_structures_complete,
            'seq_alignment_dicts': {},
            'msa_df': pd.DataFrame(),
            'distance_table': pd.DataFrame(),
            'ca_msa_df': pd.DataFrame(),
            'ca_distance_table': pd.DataFrame(),
            'global_ref': None,
            'type_reference_dict': {}
        }

    # Check if RMSD matrix is valid
    if rmsd_df.empty or len(pdb_list) < 2:
        print("No valid RMSD matrix available for alignment.")
        return {
            'processed_structures': processed_structures_complete,
            'seq_alignment_dicts': {},
            'msa_df': pd.DataFrame(),
            'distance_table': pd.DataFrame(),
            'ca_msa_df': pd.DataFrame(),
            'ca_distance_table': pd.DataFrame(),
            'global_ref': None,
            'type_reference_dict': {}
        }
    
    # Filter out structures with suspiciously low RMSD values (<0.1)
    filtered_pdb_list = []
    low_rmsd_structures = []
    low_rmsd_values = {}
    
    for pdb_id in pdb_list:
        if pdb_id in rmsd_df.index:
            avg_rmsd = rmsd_df.loc[pdb_id].mean()
            if avg_rmsd < 0.1:
                print(f"[WARNING] Excluding {pdb_id} from GRN assignment: suspiciously low RMSD ({avg_rmsd:.4f})")
                low_rmsd_structures.append(pdb_id)
                low_rmsd_values[pdb_id] = avg_rmsd
                continue
        filtered_pdb_list.append(pdb_id)
    
    if len(filtered_pdb_list) < len(pdb_list):
        print(f"[INFO] Filtered out {len(pdb_list) - len(filtered_pdb_list)} structures with RMSD < 0.1")
        # Debug: Print all low RMSD structures with their values
        if low_rmsd_structures:
            print(f"[DEBUG] Low RMSD structures: " + ", ".join([f"{sid} ({low_rmsd_values[sid]:.4f}Å)" for sid in low_rmsd_structures]))
        pdb_list = filtered_pdb_list

    # Find type references and global reference
    type_reference_dict = find_type_references(rmsd_df.loc[pdb_list, pdb_list], group_dict)

    print(type_reference_dict)

    global_ref = find_global_reference(rmsd_df.loc[pdb_list, pdb_list], type_reference_dict)

    # extremely important, the reference below if hardcoded must be updated in msa_grn.py!!!!
    # global_ref = 'CnChR2_J230_refine9'

    print(f"Global reference structure: {global_ref}")

    # Use cached alignment paths instead of recalculating alignments
    print("Using cached alignment paths from step 5 instead of recalculating alignments...")
    seq_alignment_dicts = create_seq_alignment_dicts_from_paths(
        alignment_paths=alignment_paths,
        structure_ids=pdb_list,
        global_ref=global_ref,
        type_reference_dict=type_reference_dict
    )

    # Generate all MSA tables with GRN labeling
    print("Generating MSA and distance tables with GRN labeling...")

    try:
        # Get structure mapping from data_dict if available
        structure_mapping = data_dict.get('structure_mapping', {})
        
        # Pass the RMSD matrix and structure mapping to filter structures properly
        tables = generate_grn_msa_tables(
            seq_alignment_dicts,
            processed_structures_complete,
            global_ref,
            rmsd_df=rmsd_df,  # Pass RMSD matrix for filtering
            max_rmsd_threshold=3.0,  # Filter structures with RMSD > 3.0 to reference
            structure_mapping=structure_mapping  # Pass structure mapping to prioritize experimental structures
        )

        # Extract tables
        msa_df = tables["residue_table"]
        distance_table = tables["distance_table"]
        ca_msa_df = tables["ca_residue_table"]
        ca_distance_table = tables["ca_distance_table"]

        # Report on any excluded structures
        if "excluded_structures" in tables and tables["excluded_structures"]:
            excluded_count = len(tables["excluded_structures"])
            print(f"\n[INFO] {excluded_count} structures were excluded from MSA due to high RMSD (>3.0Å)")
            print(f"[INFO] Final MSA includes {len(msa_df)} structures and {len(msa_df.columns)} positions")

        # Save the tables using direct file operations
        msa_df.to_csv(os.path.join(output_dir, "msa_table_grn.csv"))
        distance_table.to_csv(os.path.join(output_dir, "distance_table_grn.csv"))
        ca_msa_df.to_csv(os.path.join(output_dir, "ca_msa_table_grn.csv"))
        ca_distance_table.to_csv(os.path.join(output_dir, "ca_distance_table_grn.csv"))

        _50_positions = []

        print("MSA and distance tables generated and saved.")

        # Display statistics about the tables
        print("\nMSA and Distance Table Statistics:")
        print(f"Number of structures: {len(msa_df)}")
        print(f"Number of aligned positions: {len(msa_df.columns)}")

        # Count TM helices in GRN positions
        tm_positions = [col for col in msa_df.columns if '.' in col and not col.startswith('L.')]
        print(f"TM residue positions: {len(tm_positions)}")

        # Count positions by helix
        for helix in range(1, 8):
            helix_positions = [col for col in msa_df.columns if col.startswith(f"{helix}.")]
            print(f"  Helix {helix}: {len(helix_positions)} positions")

        # Distance statistics - handle NaN values properly
        print(f"\nDistance statistics:")
        if not distance_table.empty and not ca_distance_table.empty:
            avg_sidechain = distance_table.mean(skipna=True).mean(skipna=True)
            avg_backbone = ca_distance_table.mean(skipna=True).mean(skipna=True)
            print(f"  Average distance to RET (sidechain): {avg_sidechain:.2f}Å")
            print(f"  Average distance to RET (backbone): {avg_backbone:.2f}Å")

            # Find closest residues to RET (handling NaN values)
            closest_residues = []
            for col in distance_table.columns:
                if '.' in col and not col.startswith('L.'):
                    avg_dist = distance_table[col].mean(skipna=True)
                    if not pd.isna(avg_dist):  # Skip columns with all NaN
                        closest_residues.append((col, avg_dist))

            closest_residues.sort(key=lambda x: x[1])
            print("\nTop 10 closest residues to RET (across all structures):")
            for pos, dist in closest_residues[:10]:
                print(f"  {pos}: {dist:.2f}Å")

        # Now run the tree-based alignment method as an additional approach
        print("\n=== Running Tree-Based Alignment as Additional Method ===")
        tree_based_results = run_tree_based_alignment(data_dict, output_dir, visualize)

        # Combine results from both methods
        result_dict = {
            'processed_structures': processed_structures_complete,
            'seq_alignment_dicts': seq_alignment_dicts,
            'msa_df': msa_df,
            'distance_table': distance_table,
            'ca_msa_df': ca_msa_df,
            'ca_distance_table': ca_distance_table,
            'global_ref': global_ref,
            'type_reference_dict': type_reference_dict
        }
        
        # Add tree-based results with distinct keys
        result_dict.update({
            'tree_seq_alignment_dicts': tree_based_results.get('seq_alignment_dicts', {}),
            'tree_msa_df': tree_based_results.get('msa_df', pd.DataFrame()),
            'tree_distance_table': tree_based_results.get('distance_table', pd.DataFrame()),
            'tree_ca_msa_df': tree_based_results.get('ca_msa_df', pd.DataFrame()),
            'tree_ca_distance_table': tree_based_results.get('ca_distance_table', pd.DataFrame()),
            'tree_central_ref': tree_based_results.get('central_ref', None),
            'tree_enhanced_paths': tree_based_results.get('enhanced_paths', {})
        })
        
        return result_dict

    except Exception as e:
        print(f"Error during GRN assignment: {e}")
        import traceback
        traceback.print_exc()

        return {
            'processed_structures': processed_structures_complete,
            'seq_alignment_dicts': {},
            'msa_df': pd.DataFrame(),
            'distance_table': pd.DataFrame(),
            'ca_msa_df': pd.DataFrame(),
            'ca_distance_table': pd.DataFrame(),
            'global_ref': global_ref if 'global_ref' in locals() else None,
            'type_reference_dict': type_reference_dict if 'type_reference_dict' in locals() else {}
        }


def run_tree_based_alignment(data_dict, output_dir='output', visualize=True, method='average'):
    """
    Run the tree-based alignment approach to generate alternative GRN assignments.
    
    This function uses hierarchical clustering to build a guide tree for progressive alignment,
    which can handle structures that don't align well to a single reference.
    
    Args:
        data_dict: Dictionary with data from previous steps
        output_dir: Directory to save output files
        visualize: Whether to generate visualizations
        method: Linkage method for hierarchical clustering ('average', 'single', 'complete', etc.)
        
    Returns:
        Dictionary with tree-based alignment results
    """
    print("\nRunning tree-based alignment for GRN assignment...")
    
    # Create a subdirectory for tree-based outputs
    tree_output_dir = os.path.join(output_dir, 'opsin_grn_based')
    os.makedirs(tree_output_dir, exist_ok=True)
    print(f"Saving tree-based outputs to: {tree_output_dir}")
    
    processed_structures = data_dict.get('processed_structures', {})
    alignment_paths = data_dict.get('alignment_paths', {})
    rmsd_df = data_dict.get('rmsd_df', pd.DataFrame())
    pdb_list = data_dict.get('pdb_list', [])
    
    # Generate tree-based alignment dictionaries
    try:
        # Create tree-based sequence alignment dictionaries
        seq_alignment_dicts, central_ref, enhanced_paths = create_tree_based_seq_alignment_dicts(
            rmsd_df.loc[pdb_list, pdb_list], 
            alignment_paths, 
            method=method
        )
        
        print(f"Tree-based central reference structure: {central_ref}")
        print(f"Generated {len(enhanced_paths) - len(alignment_paths)} additional transitive alignments")
        
        # Get structure mapping from data_dict if available
        structure_mapping = data_dict.get('structure_mapping', {})
        
        # Generate GRN tables using the tree-based alignments
        tables = generate_grn_msa_tables(
            seq_alignment_dicts,
            processed_structures,
            central_ref,
            rmsd_df=rmsd_df,
            max_rmsd_threshold=3.0,
            structure_mapping=structure_mapping
        )
        
        # Extract tables
        msa_df = tables["residue_table"]
        distance_table = tables["distance_table"]
        ca_msa_df = tables["ca_residue_table"]
        ca_distance_table = tables["ca_distance_table"]
        
        # Report on any excluded structures
        if "excluded_structures" in tables and tables["excluded_structures"]:
            excluded_count = len(tables["excluded_structures"])
            print(f"\n[INFO] {excluded_count} structures were excluded from tree-based MSA due to high RMSD (>3.0Å)")
            print(f"[INFO] Final tree-based MSA includes {len(msa_df)} structures and {len(msa_df.columns)} positions")
        
        # Save the tree-based tables
        msa_df.to_csv(os.path.join(tree_output_dir, "msa_table_grn_tree.csv"))
        distance_table.to_csv(os.path.join(tree_output_dir, "distance_table_grn_tree.csv"))
        ca_msa_df.to_csv(os.path.join(tree_output_dir, "ca_msa_table_grn_tree.csv"))
        ca_distance_table.to_csv(os.path.join(tree_output_dir, "ca_distance_table_grn_tree.csv"))
        
        print("Tree-based MSA and distance tables generated and saved.")
        
        # Display statistics about the tables
        print("\nTree-Based MSA Statistics:")
        print(f"Number of structures: {len(msa_df)}")
        print(f"Number of aligned positions: {len(msa_df.columns)}")
        
        # Count TM helices in GRN positions
        tm_positions = [col for col in msa_df.columns if '.' in col and not col.startswith('L.')]
        print(f"TM residue positions: {len(tm_positions)}")
        
        # Count positions by helix
        for helix in range(1, 8):
            helix_positions = [col for col in msa_df.columns if col.startswith(f"{helix}.")]
            print(f"  Helix {helix}: {len(helix_positions)} positions")
        
        return {
            'seq_alignment_dicts': seq_alignment_dicts,
            'msa_df': msa_df,
            'distance_table': distance_table,
            'ca_msa_df': ca_msa_df,
            'ca_distance_table': ca_distance_table,
            'central_ref': central_ref,
            'enhanced_paths': enhanced_paths
        }
        
    except Exception as e:
        print(f"Error during tree-based alignment: {e}")
        import traceback
        traceback.print_exc()
        
        return {
            'seq_alignment_dicts': {},
            'msa_df': pd.DataFrame(),
            'distance_table': pd.DataFrame(),
            'ca_msa_df': pd.DataFrame(),
            'ca_distance_table': pd.DataFrame(),
            'central_ref': None,
            'enhanced_paths': {}
        }
