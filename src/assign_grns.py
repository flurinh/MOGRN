import pandas as pd
import os
import matplotlib.pyplot as plt

from src.reference_alignment import (
    find_type_references,
    find_global_reference,
    create_seq_alignment_dicts_from_paths
)

from src.msa_grn import (
    analyze_residue_composition,
    generate_grn_msa_tables
)

from src.visualization_functions import (
    plot_average_distances_by_helix,
    plot_distance_heatmap,
    print_residue_composition,
    create_residue_conservation_plot,
    plot_helix_logo_plots
)


def align_and_assign_grn(data_dict, output_dir='output', visualize=True):
    """
    Step 6: Structure alignment and GRN assignment

    This function assigns Generic Residue Numbers (GRN) using cached alignment paths
    without recalculating alignments.

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

    # extremely important, the reference we use is hardcoded here!!!!
    global_ref = 'CnChR2_J230_refine9'

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

        # Visualizations
        if visualize:
            try:
                # Plot distance heatmap and line plots
                if not distance_table.empty:
                    fig1 = plot_average_distances_by_helix(distance_table)
                    fig1.savefig(os.path.join(output_dir, 'distances_by_helix.png'), dpi=300, bbox_inches='tight')
                    plt.close(fig1)
                    print(f"Saved distance by helix plot to {os.path.join(output_dir, 'distances_by_helix.png')}")

                    fig2 = plot_distance_heatmap(distance_table)
                    fig2.savefig(os.path.join(output_dir, 'distance_heatmap.png'), dpi=300, bbox_inches='tight')
                    plt.close(fig2)
                    print(f"Saved distance heatmap to {os.path.join(output_dir, 'distance_heatmap.png')}")

                # Analyze residue composition at key positions
                if not msa_df.empty:
                    positions_to_analyze = [f"{helix}.50" for helix in range(1, 8)]
                    positions_in_df = [pos for pos in positions_to_analyze if pos in msa_df.columns]

                    if positions_in_df:
                        residue_composition = analyze_residue_composition(msa_df, positions_in_df)
                        print_residue_composition(residue_composition)

                    # Also analyze binding pocket positions
                    binding_pocket_positions = ["3.37", "3.40", "6.48", "6.51", "7.43", "7.46", "7.51"]
                    positions_in_df = [pos for pos in binding_pocket_positions if pos in msa_df.columns]

                    if positions_in_df:
                        binding_pocket_composition = analyze_residue_composition(msa_df, positions_in_df)
                        print_residue_composition(binding_pocket_composition)

                    # Create helix logo plots instead of sequence logo
                    try:
                        fig = plot_helix_logo_plots(msa_df)
                        fig.savefig(os.path.join(output_dir, "helix_logo_plots.png"), dpi=300, bbox_inches='tight')
                        plt.close(fig)
                        print(f"Saved helix logo plots to {os.path.join(output_dir, 'helix_logo_plots.png')}")
                    except Exception as e:
                        print(f"[WARNING] Error creating helix logo plots: {str(e)}")

                    # Create residue conservation plot
                    fig = create_residue_conservation_plot(msa_df)
                    fig.savefig(os.path.join(output_dir, "residue_conservation.png"), dpi=300, bbox_inches='tight')
                    plt.close(fig)
                    print(f"Saved residue conservation plot to {os.path.join(output_dir, 'residue_conservation.png')}")

            except Exception as e:
                print(f"Warning: Error during visualization: {e}")

        return {
            'processed_structures': processed_structures_complete,
            'seq_alignment_dicts': seq_alignment_dicts,
            'msa_df': msa_df,
            'distance_table': distance_table,
            'ca_msa_df': ca_msa_df,
            'ca_distance_table': ca_distance_table,
            'global_ref': global_ref,
            'type_reference_dict': type_reference_dict
        }

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
