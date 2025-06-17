# plot.py

import os
import sys
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import traceback
from pathlib import Path
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage

# --- Setup Project Paths ---
try:
    current_file_path = Path(__file__).resolve()
    PROJECT_DIR = current_file_path.parent
except NameError:
    print("Running in interactive mode. Using Path.cwd() for project_dir.")
    PROJECT_DIR = Path.cwd()

print(f"[INFO] Project directory set to: {PROJECT_DIR}")

SRC_DIR = PROJECT_DIR / 'src'
if str(SRC_DIR) not in sys.path:
    sys.path.append(str(SRC_DIR))
    print(f"[INFO] Added '{SRC_DIR}' to sys.path")

FIGURES_OUTPUT_DIR = PROJECT_DIR / 'opsin_output' / 'figures'
os.makedirs(FIGURES_OUTPUT_DIR, exist_ok=True)
print(f"[INFO] Figures will be saved to: {FIGURES_OUTPUT_DIR}")

CACHE_DIR = PROJECT_DIR / 'opsin_output' / 'cache'
GRN_TABLES_DIR = PROJECT_DIR / 'opsin_grn_tables'

try:
    from src.data_processing import load_opsin_property_data
    from src.visualization_functions import (
        create_opsin_overview_plot,
        visualize_rmsd_matrix_improved,
        plot_distances_with_std,
        plot_helix_logo_plots
    )
    from protos.processing.grn.grn_utils import sort_grns_str, get_tm_residues

    print("[INFO] Successfully imported custom modules.")
except ImportError as e:
    print(f"[ERROR] Failed to import custom modules: {e}")
    traceback.print_exc()
    sys.exit(1)


# --- Helper Functions for Data Loading ---
def load_cached_data(cache_path, description="data"):
    if os.path.exists(cache_path):
        print(f"[INFO] Loading {description} from cache: {cache_path}")
        try:
            with open(cache_path, 'rb') as f:
                result = pickle.load(f)
            print(f"[INFO] Successfully loaded {description}")
            return result
        except Exception as e:
            print(f"[ERROR] Error loading {description} from cache '{cache_path}': {e}")
            traceback.print_exc()
    else:
        print(f"[WARN] Cache file not found: {cache_path}")
    return None


def load_grn_tables_data():  # Matching notebook name
    grn_tables_pkl = GRN_TABLES_DIR / 'grn_tables_data.pkl'  # Corrected path
    return load_cached_data(grn_tables_pkl, "GRN tables data")


# --- Main Script Logic ---
def main():
    print("\n" + "=" * 30 + " STARTING VISUALIZATION SCRIPT " + "=" * 30)

    data = {}  # Matching notebook variable name
    chain_id = 'A'

    cache_files = {
        "raw_structures": (f"raw_structures_{chain_id}.pkl", "raw structures"),
        "processed_structures": (f"processed_structures_{chain_id}.pkl", "processed structures"),
        "structure_errors": (f"structure_errors_{chain_id}.pkl", "structure errors"),
        "helix_annotations": (f"helix_annotations_{chain_id}.pkl", "helix annotations"),
        "structure_comparison": (f"structure_comparison_{chain_id}.pkl", "structure comparison"),
        "grn_assignment": (f"grn_assignment_{chain_id}.pkl", "GRN assignment")
    }

    print("\n[PHASE] Loading cached workflow data...")
    for data_type, (file_name, description) in cache_files.items():
        cache_path = CACHE_DIR / file_name  # Corrected path
        component_data = load_cached_data(cache_path, description)
        if component_data is not None:
            if isinstance(component_data, dict):
                data.update(component_data)
            else:
                data[data_type] = component_data

    grn_data = load_grn_tables_data()
    if grn_data: data.update(grn_data)

    print(f"\n[INFO] Loaded data with {len(data.keys())} top-level keys")
    if 'processed_structures' in data:
        print(f"[INFO] Found {len(data['processed_structures'])} processed structures")
    else:
        print("[ERROR] 'processed_structures' not found. This is crucial.")
        # return # Exiting if crucial data is missing

    print("\n[PHASE] Loading opsin property data...")
    property_csv_path = PROJECT_DIR / 'property' / 'mo_exp.csv'  # Corrected path

    processed_structures_for_prop = data.get('processed_structures', {})
    if not processed_structures_for_prop:
        print("[WARN] No processed structures available for property loading.")

    property_data_loaded = load_opsin_property_data(property_csv_path, processed_structures_for_prop)  # Renamed

    if property_data_loaded and 'properties' in property_data_loaded:
        data['property_data'] = {k.replace('_smile_', '_'): v for k, v in property_data_loaded['properties'].items()}
        data['structure_mapping'] = property_data_loaded['structure_mapping']
        print(f"[INFO] Loaded property data for {len(data['property_data'])} entries.")
    else:
        print("[WARN] Failed to load or parse property data.")
        data['property_data'] = {}
        data['structure_mapping'] = {}

    # --- Visualization 1: Opsins Overview Plot ---
    print("\n[VISUALIZATION 1] Generating Opsin Overview Plot...")
    try:
        overview_data_list = []
        processed_structures_map = data.get('processed_structures', {})
        property_data_map = data.get('property_data', {})

        for sid, struct_info in processed_structures_map.items():
            mf, domain, is_experimental, display_name = "Unknown", "Unknown", False, sid
            prop_entry = property_data_map.get(sid)
            if prop_entry and isinstance(prop_entry, dict):
                mf_val = prop_entry.get('molecular_function')  # As per notebook
                if pd.notna(mf_val) and str(mf_val).strip() and str(mf_val).lower() != "unknown": mf = str(mf_val)
                domain_val = prop_entry.get('domain')
                if pd.notna(domain_val) and str(domain_val).strip() and str(
                    domain_val).lower() != "unknown": domain = str(domain_val)
                is_experimental_prop = prop_entry.get('experimentally_determined', False)
                is_experimental = str(is_experimental_prop).lower() == 'true' if isinstance(is_experimental_prop,
                                                                                            str) else bool(
                    is_experimental_prop)
                display_name_prop = prop_entry.get('short_name', prop_entry.get('display_name', sid))
                if pd.notna(display_name_prop) and str(display_name_prop).strip(): display_name = str(display_name_prop)
            elif isinstance(struct_info, dict):  # Fallback to structure_info as per notebook
                struct_props = struct_info.get('properties', {})
                if isinstance(struct_props, dict):
                    if mf == "Unknown":
                        mf_val_struct = struct_props.get('molecular_function')
                        if pd.notna(mf_val_struct) and str(mf_val_struct).strip() and str(
                            mf_val_struct).lower() != "unknown": mf = str(mf_val_struct)
                    if domain == "Unknown":
                        domain_val_struct = struct_props.get('domain')
                        if pd.notna(domain_val_struct) and str(domain_val_struct).strip() and str(
                            domain_val_struct).lower() != "unknown": domain = str(domain_val_struct)

            overview_data_list.append({
                'id': sid, 'short_name': display_name,
                'molecular_function_normalized': mf, 'domain': domain,
                'experimentally_determined': is_experimental
            })
        overview_df = pd.DataFrame(overview_data_list)

        if not overview_df.empty:
            # Update 'experimentally_determined' based on structure_mapping as per notebook
            flat_list = set()  # Use a set for faster lookups
            structure_map = data.get('structure_mapping', {})
            for k, v_list_or_str in structure_map.items():
                flat_list.add(k)
                if isinstance(v_list_or_str, list):
                    flat_list.update(v_list_or_str)
                elif isinstance(v_list_or_str, str):
                    flat_list.add(v_list_or_str)
            overview_df['experimentally_determined'] = overview_df['id'].isin(flat_list)

            # Call create_opsin_overview_plot - assuming it now handles its own figsize or uses a default
            fig1 = create_opsin_overview_plot(overview_df)
            # If you need to set title explicitly after creation:
            # fig1.suptitle("Opsin Overview", fontsize=16) # Use suptitle for figure-level title

            fig1_path = FIGURES_OUTPUT_DIR / "01_opsin_overview.png"
            fig1.savefig(fig1_path, dpi=300, bbox_inches='tight')
            print(f"[SUCCESS] Saved Opsin Overview Plot to: {fig1_path}")
            plt.close(fig1)
        else:
            print("[WARN] Overview DataFrame is empty. Skipping Opsin Overview Plot.")
    except Exception as e:
        print(f"[ERROR] Failed to generate Opsin Overview Plot: {e}")
        traceback.print_exc()

    # --- Prepare group_dict and domain_dict as per notebook ---
    group_dict = {}
    domain_dict = {}  # This will store {'domain': name, 'average_error': val}

    # Use 'overview_data_list' as it was constructed with latest property info
    for item in overview_data_list:  # Ensure overview_data_list is used
        group_dict[item['id']] = item['molecular_function_normalized']
        # For domain_dict, fetch average_error if available
        avg_error = data.get('structure_errors', {}).get(item['id'], {}).get('average_error')
        domain_dict[item['id']] = {'domain': item['domain'], 'average_error': avg_error}

    # --- Get RMSD matrix from data as per notebook ---
    rmsd_df = data.get('rmsd_matrix')
    if not isinstance(rmsd_df, pd.DataFrame) or rmsd_df.empty:
        rmsd_df = data.get('rmsd_df')  # Fallback key

    if not isinstance(rmsd_df, pd.DataFrame) or rmsd_df.empty:  # Check again
        print("[WARN] RMSD matrix (rmsd_df or rmsd_matrix) not found or empty in data. Skipping RMSD Heatmap.")
        rmsd_df = pd.DataFrame()  # Ensure it's an empty DataFrame if not found
    else:
        print(f"[INFO] Found RMSD matrix with shape: {rmsd_df.shape}")

    # --- Visualization 2: RMSD Heatmap with Improved Visualization (Clustermap) ---
    # Renumbered from Vis 3 in notebook to Vis 2 here
    print("\n[VISUALIZATION 2] Generating RMSD Clustermap...")
    if not rmsd_df.empty:
        try:
            # Prepare linkage matrix
            rmsd_matrix_for_linkage = rmsd_df.copy()
            # Fill NaNs robustly for linkage calculation
            if rmsd_matrix_for_linkage.isnull().values.any():
                finite_vals = rmsd_matrix_for_linkage.values[np.isfinite(rmsd_matrix_for_linkage.values)]
                fill_val = np.nanmax(finite_vals) if finite_vals.size > 0 else 10.0  # Default if all NaNs
                rmsd_matrix_for_linkage.fillna(fill_val, inplace=True)
            np.fill_diagonal(rmsd_matrix_for_linkage.values, 0.0)  # Ensure diagonal is zero

            condensed_matrix = squareform(rmsd_matrix_for_linkage.values, checks=False)
            Z_linkage = linkage(condensed_matrix, method='average')

            # Call the improved visualization function
            # Assuming visualize_rmsd_matrix_improved handles its own saving if output_file is passed
            fig2 = visualize_rmsd_matrix_improved(
                rmsd_df=rmsd_df,
                group_dict=group_dict,
                domain_dict=domain_dict,  # Pass the dict with 'domain' and 'average_error' keys
                linkage_matrix=Z_linkage,
                figsize=(18, 15),  # Adjusted figsize
                output_file=FIGURES_OUTPUT_DIR / "02_rmsd_clustermap.png"  # Pass output file directly
            )
            if fig2 and hasattr(fig2, 'fig'):  # If it returns ClusterGrid object
                # If the function didn't save, or to be sure
                # fig2.fig.savefig(FIGURES_OUTPUT_DIR / "02_rmsd_clustermap.png", dpi=300)
                print(f"[INFO] RMSD Clustermap generation initiated (saved by function).")
                plt.close(fig2.fig)
            elif fig2:  # If it returns a Figure object directly
                # fig2.savefig(FIGURES_OUTPUT_DIR / "02_rmsd_clustermap.png", dpi=300)
                print(f"[INFO] RMSD Clustermap generation initiated (saved by function).")
                plt.close(fig2)
            else:
                print("[WARN] RMSD Clustermap was not generated.")
        except Exception as e:
            print(f"[ERROR] Failed to generate RMSD Clustermap: {e}")
            traceback.print_exc()
    else:
        print("[WARN] RMSD matrix is empty. Skipping RMSD Clustermap.")

    # --- Prepare data for subsequent distance and logo plots ---
    print("\n[PHASE] Preparing data for distance and logo plots...")
    distance_table = data.get('distance_table')
    ca_distance_table = data.get('ca_distance_table')

    # Determine residue_table_for_logo as per notebook
    residue_table_for_logo = data.get('msa_df')
    msa_table_notebook = data.get('msa_table')  # Check if 'msa_table' was loaded
    if not isinstance(residue_table_for_logo, pd.DataFrame) or residue_table_for_logo.empty:
        if isinstance(msa_table_notebook, pd.DataFrame) and not msa_table_notebook.empty:
            residue_table_for_logo = msa_table_notebook
            print("[INFO] Using 'msa_table' for logo plots.")
        elif 'residue_table' in data and isinstance(data['residue_table'], pd.DataFrame):
            residue_table_for_logo = data['residue_table']  # Fallback if 'msa_table' also not good
            print("[INFO] Using 'residue_table' as fallback for logo plots.")
        else:
            print("[WARN] No suitable 'msa_df', 'msa_table', or 'residue_table' found for logo plots.")
            residue_table_for_logo = pd.DataFrame()

    # Initialize filtered versions
    df_filtered_residue = pd.DataFrame()
    df_filtered_distance = pd.DataFrame()
    df_filtered_ca_distance = pd.DataFrame()
    grns = []

    if not residue_table_for_logo.empty:
        col_name_filter = '7.50'  # Renamed from col_name
        if col_name_filter in residue_table_for_logo.columns:
            print(f"[INFO] Filtering by '{col_name_filter}' starting with 'K'...")
            boolean_mask = residue_table_for_logo[col_name_filter].astype(str).str.startswith('K', na=False)
            df_filtered_residue = residue_table_for_logo[boolean_mask]
            print(f"[INFO] {len(df_filtered_residue)} entries after filtering residues.")

            if not df_filtered_residue.empty:
                if isinstance(distance_table, pd.DataFrame) and not distance_table.empty:
                    common_idx = distance_table.index.intersection(df_filtered_residue.index)
                    df_filtered_distance = distance_table.loc[common_idx]
                else:
                    print("[WARN] Original distance_table is empty or not a DataFrame.")

                if isinstance(ca_distance_table, pd.DataFrame) and not ca_distance_table.empty:
                    common_idx_ca = ca_distance_table.index.intersection(df_filtered_residue.index)
                    df_filtered_ca_distance = ca_distance_table.loc[common_idx_ca]
                else:
                    print("[WARN] Original ca_distance_table is empty or not a DataFrame.")

                if not df_filtered_distance.empty:
                    # Get GRNs from columns of df_filtered_distance as per notebook
                    grns = get_tm_residues(sort_grns_str(df_filtered_distance.columns.astype(str).tolist()))
                    print(f"[INFO] Derived {len(grns)} TM GRNs for plotting from filtered distance table.")
                else:
                    print("[WARN] Filtered distance table is empty, cannot derive GRNs for plotting.")
            else:
                print(
                    "[WARN] No residues remained after filtering. Distance tables will not be filtered based on this.")
        else:
            print(f"[WARN] Filter column '{col_name_filter}' not found in residue_table_for_logo.")
    else:
        print("[WARN] residue_table_for_logo is empty. Cannot perform filtering for distance/logo plots.")

    # --- Visualization 3: All-Atom Distance to Retinal Plot (was Vis 4) ---
    print("\n[VISUALIZATION 3] Generating All-Atom Distance to Retinal Plot...")
    if not df_filtered_distance.empty and grns:  # Ensure grns list is also populated
        try:
            fig3 = plot_distances_with_std(
                df_filtered_distance[grns],  # Select only the sorted TM GRN columns
                title="All-Atom Distance to Retinal by Position",
                use_ca=False, figsize=(14, 8)  # Adjusted figsize
            )
            fig3_path = FIGURES_OUTPUT_DIR / "03_all_atom_distance_std.png"
            fig3.savefig(fig3_path, dpi=300);
            plt.close(fig3)
            print(f"[SUCCESS] Saved All-Atom Distance Plot to: {fig3_path}")
        except KeyError as e:
            print(f"[ERROR] KeyError during All-Atom Distance Plot (some GRNs might be missing after filtering): {e}")
            traceback.print_exc()
        except Exception as e:
            print(f"[ERROR] Failed to generate All-Atom Distance Plot: {e}")
            traceback.print_exc()
    else:
        print("[WARN] Filtered distance data or GRNs list is empty. Skipping All-Atom Distance Plot.")

    # --- Visualization 4: CA-Atom Distance to Retinal Plot (was Vis 5) ---
    print("\n[VISUALIZATION 4] Generating CA-Atom Distance to Retinal Plot...")
    # For CA plot, we need GRNs relevant to df_filtered_ca_distance if its columns differ
    # However, plot_distances_with_std itself filters for TM columns.
    if not df_filtered_ca_distance.empty:
        try:
            # It's safer if grns for ca_distance are derived from its own columns if they differ
            # For now, assume plot_distances_with_std handles column selection.
            fig4 = plot_distances_with_std(
                df_filtered_ca_distance,  # Pass the whole filtered CA table
                title="CA-Atom Distance to Retinal by Position",
                use_ca=True, figsize=(14, 8)  # Adjusted figsize
            )
            fig4_path = FIGURES_OUTPUT_DIR / "04_ca_atom_distance_std.png"
            fig4.savefig(fig4_path, dpi=300);
            plt.close(fig4)
            print(f"[SUCCESS] Saved CA-Atom Distance Plot to: {fig4_path}")
        except Exception as e:
            print(f"[ERROR] Failed to generate CA-Atom Distance Plot: {e}")
            traceback.print_exc()
    else:
        print("[WARN] Filtered CA distance data is empty. Skipping CA-Atom Distance Plot.")

    # --- Visualization 5: Helix Logo Plots (was Vis 7) ---
    print("\n[VISUALIZATION 5] Generating Helix Logo Plots...")
    if not df_filtered_residue.empty:  # Use df_filtered_residue as per notebook for logos
        try:
            fig5 = plot_helix_logo_plots(
                df_filtered_residue,
                frequency_threshold=0.07  # As per notebook
            )
            fig5_path = FIGURES_OUTPUT_DIR / "05_helix_logos_x50.png"
            fig5.savefig(fig5_path, dpi=300);
            plt.close(fig5)
            print(f"[SUCCESS] Saved Helix Logo Plots to: {fig5_path}")
        except Exception as e:
            print(f"[ERROR] Failed to generate Helix Logo Plots: {e}")
            traceback.print_exc()
    else:
        print("[WARN] Filtered residue data for logos is empty. Skipping Helix Logo Plots.")

    print("\n" + "=" * 30 + " VISUALIZATION SCRIPT COMPLETE " + "=" * 30)


if __name__ == "__main__":
    main()