# plot.py
# Consolidated visualization script for MOGRN project
# Includes all visualization functionality:
# - Opsin overview plots
# - RMSD heatmaps and clustering
# - Distance plots (all-atom and CA)
# - Helix logo plots
# - Property analysis (natural domains, contributions, missing combinations)
# - Interactive GRN alignment visualization

import os
import sys
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import traceback
import argparse
from pathlib import Path
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage
from collections import Counter

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

# Default directories - will be updated based on command-line args
FIGURES_OUTPUT_DIR = PROJECT_DIR / 'opsin_output' / 'paper_figures'
CACHE_DIR = PROJECT_DIR / 'opsin_output' / 'cache'
GRN_TABLES_DIR = PROJECT_DIR / 'opsin_grn_tables'

from src.visualize_alignment_grn import create_opsin_visualization_from_workflow

try:
    from src.property_mapping import create_unified_property_mapper
except:
    print("cannot import unified property mapper")
from src.opsin_color_scheme import OPSIN_COLORS
from src.visualization_functions import (
        compute_rmsd_metrics,
        create_opsin_overview_plot,
        visualize_rmsd_matrix_improved,
        _annotate_metrics_on_clustergrid,
        plot_distances_with_std,
        plot_helix_logo_plots,
        plot_error_box_comparison
    )
try:
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


def load_grn_tables_data():
    # Try multiple possible locations for GRN tables
    possible_paths = [
        GRN_TABLES_DIR / 'grn_tables_data.pkl',
        PROJECT_DIR / 'opsin_grn_tables' / 'grn_tables_data.pkl',
        PROJECT_DIR / 'opsin_output' / 'grn_tables_data.pkl'
    ]
    
    for grn_path in possible_paths:
        if grn_path.exists():
            return load_cached_data(grn_path, "GRN tables data")
    
    print(f"[WARN] GRN tables data not found in any of the expected locations")
    return None


# --- Property Analysis Functions ---
def create_property_analysis_visualizations(mapper, processed_structures, output_dir):
    """Create property analysis visualizations."""
    print("\n[PROPERTY ANALYSIS] Starting property analysis visualizations...")
    
    # Natural domains to include
    NATURAL_DOMAINS = ['Eukaryota', 'Bacteria', 'Archaea']
    
    # Get all experimental-predicted pairs
    all_structure_ids = list(processed_structures.keys())
    all_pairs = mapper.get_all_experimental_predicted_pairs(all_structure_ids)
    pred_with_exp = set(all_pairs.values())
    
    # Collect properties for different datasets
    experimental_props = []
    predicted_only_props = []
    combined_props = []
    
    for struct_id in processed_structures:
        props = mapper.get_properties(struct_id)
        if props and props['domain'] in NATURAL_DOMAINS:
            struct_type, _ = mapper.identify_structure_type(struct_id)
            
            if struct_type in ['standard_exp', 'hideaki_exp']:
                experimental_props.append(props)
                combined_props.append(props)
            elif struct_type in ['standard_pred', 'hideaki_pred']:
                if struct_id not in pred_with_exp:
                    predicted_only_props.append(props)
                    combined_props.append(props)
    
    print(f"Dataset composition (natural domains only):")
    print(f"  - Experimental structures: {len(experimental_props)}")
    print(f"  - Predicted-only structures: {len(predicted_only_props)}")
    print(f"  - Combined dataset: {len(combined_props)}")
    
    # Create comprehensive property analysis figure
    create_property_distribution_figure(experimental_props, predicted_only_props, 
                                      combined_props, NATURAL_DOMAINS, output_dir)
    
    # Create prediction contribution bar chart
    create_prediction_contribution_chart(experimental_props, predicted_only_props, 
                                       NATURAL_DOMAINS, output_dir)


def create_property_distribution_figure(experimental_props, predicted_only_props, 
                                       combined_props, natural_domains, output_dir):
    """Create comprehensive property distribution analysis figure."""
    from src.opsin_color_scheme import RMSD_WHITE_TO_DARKGRAY_CMAP
    
    all_functions = sorted(list(set(p['molecular_function'] for p in combined_props)))
    all_domains = sorted(natural_domains)
    
    fig = plt.figure(figsize=(20, 14))
    gs = fig.add_gridspec(3, 3, height_ratios=[2, 2, 1], hspace=0.3, wspace=0.25)
    
    # 1. Combined Results Heatmap
    ax1 = fig.add_subplot(gs[0, :])
    df_combined = pd.DataFrame(combined_props)
    pivot_combined = pd.crosstab(df_combined['molecular_function'], df_combined['domain'])
    pivot_combined = pivot_combined.reindex(index=all_functions, columns=all_domains, fill_value=0)
    
    annot_combined = pivot_combined.astype(str)
    annot_combined[pivot_combined == 0] = ''
    
    sns.heatmap(pivot_combined, annot=annot_combined, fmt='', cmap=RMSD_WHITE_TO_DARKGRAY_CMAP,
                linewidths=0.5, linecolor='gray', cbar_kws={'label': 'Count'},
                vmin=0, ax=ax1)
    
    ax1.set_xlabel('Domain', fontsize=14)
    ax1.set_ylabel('Molecular Function', fontsize=14)
    ax1.set_title('Combined Dataset Results (Experimental + Predicted without Exp Counterpart)', fontsize=16)
    ax1.text(1.02, 0.5, f'Total: {len(combined_props)}\n({len(experimental_props)} experimental +\n{len(predicted_only_props)} predicted-only)',
             transform=ax1.transAxes, verticalalignment='center',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 2a. Experimental Only Heatmap
    ax2a = fig.add_subplot(gs[1, 0])
    if experimental_props:
        df_exp = pd.DataFrame(experimental_props)
        pivot_exp = pd.crosstab(df_exp['molecular_function'], df_exp['domain'])
        pivot_exp = pivot_exp.reindex(index=all_functions, columns=all_domains, fill_value=0)
        
        annot_exp = pivot_exp.astype(str)
        annot_exp[pivot_exp == 0] = ''
        
        sns.heatmap(pivot_exp, annot=annot_exp, fmt='', cmap=RMSD_WHITE_TO_DARKGRAY_CMAP,
                    linewidths=0.5, linecolor='gray', cbar_kws={'label': 'Count'},
                    vmin=0, ax=ax2a)
        ax2a.set_title(f'Experimental Only (n={len(experimental_props)})', fontsize=14)
    
    # 2b. Predicted Only Heatmap
    ax2b = fig.add_subplot(gs[1, 1])
    if predicted_only_props:
        df_pred_only = pd.DataFrame(predicted_only_props)
        pivot_pred_only = pd.crosstab(df_pred_only['molecular_function'], df_pred_only['domain'])
        pivot_pred_only = pivot_pred_only.reindex(index=all_functions, columns=all_domains, fill_value=0)
        
        annot_pred = pivot_pred_only.astype(str)
        annot_pred[pivot_pred_only == 0] = ''
        
        sns.heatmap(pivot_pred_only, annot=annot_pred, fmt='', cmap='YlOrRd',
                    linewidths=0.5, linecolor='gray', cbar_kws={'label': 'Count'},
                    vmin=0, ax=ax2b)
        ax2b.set_title(f'Predicted-Only Contribution (n={len(predicted_only_props)})', fontsize=14)
    
    # 2c. Difference Plot
    ax2c = fig.add_subplot(gs[1, 2])
    if experimental_props and predicted_only_props:
        pivot_diff = pivot_combined - pivot_exp
        
        colors_list = ['white', 'lightcoral', 'red', 'darkred']
        cmap = plt.cm.colors.LinearSegmentedColormap.from_list('diff_cmap', colors_list, N=10)
        
        sns.heatmap(pivot_diff, annot=True, fmt='d', cmap=cmap,
                    linewidths=0.5, linecolor='gray', cbar_kws={'label': 'Structures Added'},
                    vmin=0, ax=ax2c)
        ax2c.set_title('Contribution of Predictions\n(Combined - Experimental)', fontsize=14)
    
    # 3. Missing Combinations Analysis
    ax3 = fig.add_subplot(gs[2, :])
    ax3.axis('off')
    
    all_possible_combinations = [(f, d) for f in all_functions for d in all_domains]
    observed_combinations = set((row['molecular_function'], row['domain']) 
                               for _, row in df_combined.iterrows())
    missing_combinations = [combo for combo in all_possible_combinations 
                           if combo not in observed_combinations]
    
    missing_text = "Missing Function-Domain Combinations:\n\n"
    if missing_combinations:
        missing_by_function = {}
        for func, domain in missing_combinations:
            if func not in missing_by_function:
                missing_by_function[func] = []
            missing_by_function[func].append(domain)
        
        for func, domains in missing_by_function.items():
            missing_text += f"• {func}: {', '.join(domains)}\n"
        
        missing_text += f"\nTotal missing: {len(missing_combinations)} out of {len(all_possible_combinations)} possible combinations"
        missing_text += f"\nCoverage: {(len(observed_combinations)/len(all_possible_combinations)*100):.1f}%"
    else:
        missing_text += "All possible combinations are represented in the dataset!"
    
    ax3.text(0.05, 0.95, missing_text, transform=ax3.transAxes, fontsize=12,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightgray', alpha=0.8))
    
    plt.suptitle('Microbial Opsin Property Analysis: Natural Domains Only', fontsize=18, y=0.98)
    plt.tight_layout()
    
    fig_path = output_dir / '06_property_analysis_natural_domains.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[SUCCESS] Saved property analysis to: {fig_path}")


def create_prediction_contribution_chart(experimental_props, predicted_only_props, 
                                       natural_domains, output_dir):
    """Create bar chart showing prediction contributions."""
    from src.opsin_color_scheme import get_categorical_colors
    
    # Count by function-domain combinations
    experimental_counts = {}
    predicted_only_counts = {}
    
    for props in experimental_props:
        key = (props['molecular_function'], props['domain'])
        experimental_counts[key] = experimental_counts.get(key, 0) + 1
    
    for props in predicted_only_props:
        key = (props['molecular_function'], props['domain'])
        predicted_only_counts[key] = predicted_only_counts.get(key, 0) + 1
    
    # Get all combinations
    all_combos = sorted(set(list(experimental_counts.keys()) + list(predicted_only_counts.keys())))
    
    # Create figure
    fig, ax = plt.subplots(figsize=(14, 8))
    
    # Prepare data
    combo_labels = [f"{func}\n{domain}" for func, domain in all_combos]
    exp_values = [experimental_counts.get(combo, 0) for combo in all_combos]
    pred_values = [predicted_only_counts.get(combo, 0) for combo in all_combos]
    total_values = [e + p for e, p in zip(exp_values, pred_values)]
    
    x = np.arange(len(combo_labels))
    width = 0.35
    
    # Create bars
    bars1 = ax.bar(x - width/2, exp_values, width, label='Experimental', color='steelblue', alpha=0.8)
    bars2 = ax.bar(x + width/2, pred_values, width, label='Predicted-only', color='coral', alpha=0.8)
    
    # Add value labels
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax.text(bar.get_x() + bar.get_width()/2., height,
                       f'{int(height)}', ha='center', va='bottom', fontsize=9)
    
    # Add totals
    for i, (exp, pred) in enumerate(zip(exp_values, pred_values)):
        total = exp + pred
        max_height = max(exp, pred)
        if total > 0 and pred > 0:
            ax.text(i, max_height + 0.5, f'Total: {total}', ha='center', va='bottom', 
                    fontsize=10, fontweight='bold')
    
    # Customize plot
    ax.set_xlabel('Function-Domain Combination', fontsize=12)
    ax.set_ylabel('Number of Structures', fontsize=12)
    ax.set_title('Contribution of Predicted Structures to Each Function-Domain Combination\n(Natural Domains Only)', fontsize=14)
    ax.set_xticks(x)
    ax.set_xticklabels(combo_labels, rotation=45, ha='right')
    ax.legend(loc='upper left')
    ax.grid(axis='y', alpha=0.3)
    
    # Add percentage contribution text
    total_exp = sum(exp_values)
    total_pred = sum(pred_values)
    contribution_text = f"Overall Contribution:\n{total_pred/(total_exp+total_pred)*100:.1f}% from predictions"
    
    ax.text(0.98, 0.97, contribution_text, transform=ax.transAxes,
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            fontsize=11)
    
    plt.tight_layout()
    fig_path = output_dir / '07_prediction_contribution_bars.png'
    plt.savefig(fig_path, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"[SUCCESS] Saved prediction contribution chart to: {fig_path}")


# --- Main Script Logic ---
def main(args=None):
    """Main visualization script with command-line argument support."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Generate MOGRN visualizations')
    parser.add_argument('--input-dir', type=str, default='opsin_output',
                       help='Input directory with analysis results')
    parser.add_argument('--outputs-dir', type=str, default='opsin_output/paper_figures',
                       help='Output directory for figures')
    parser.add_argument('--skip-property-analysis', action='store_true',
                       help='Skip property analysis visualizations')
    parser.add_argument('--skip-interactive', action='store_true',
                       help='Skip interactive GRN visualization')
    parser.add_argument('--only', type=str, choices=['overview', 'errors', 'rmsd', 'distance', 'logo', 'property', 'interactive'],
                       help='Generate only specific visualization type')
    
    if args is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(args)
    
    print("\n" + "=" * 30 + " STARTING VISUALIZATION SCRIPT " + "=" * 30)
    print(f"Input directory: {args.input_dir}")
    print(f"Output directory: {args.output_dir}")

    data = {}  # Matching notebook variable name
    chain_id = 'A'

    # Update paths based on input directory
    global CACHE_DIR, GRN_TABLES_DIR, FIGURES_OUTPUT_DIR
    CACHE_DIR = Path(args.input_dir) / 'cache'
    GRN_TABLES_DIR = Path(args.input_dir) / 'tree_based_grn'
    FIGURES_OUTPUT_DIR = Path(args.output_dir)
    os.makedirs(FIGURES_OUTPUT_DIR, exist_ok=True)
    
    # Skip specific visualizations if requested
    if args.only:
        print(f"Generating only '{args.only}' visualization as requested")
    
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

    print("\n[PHASE] Loading opsin property data with unified mapper...")
    property_csv_path = PROJECT_DIR / 'property' / 'mo_exp.csv'

    processed_structures = data.get('processed_structures', {})
    if not processed_structures:
        print("[WARN] No processed structures available for property loading.")
        data['property_data'] = {}
        data['structure_mapping'] = {}
    else:
        # Create unified property mapper
        mapper = create_unified_property_mapper(property_csv_path)
        
        # Update processed structures with properties
        updated_structures = mapper.update_processed_structures(processed_structures)
        data['processed_structures'] = updated_structures
        
        # Get all properties for backward compatibility
        data['property_data'] = mapper.get_all_properties()
        
        # Get experimental-predicted mapping
        all_structure_ids = list(processed_structures.keys())
        data['structure_mapping'] = mapper.get_all_experimental_predicted_pairs(all_structure_ids)
        
        # Validate mapping
        stats = mapper.validate_mapping(updated_structures)
        print(f"[INFO] Property mapping statistics:")
        print(f"  - Total structures: {stats['total_structures']}")
        print(f"  - With properties: {stats['structures_with_properties']} ({stats['structures_with_properties']/stats['total_structures']*100:.1f}%)")
        print(f"  - Experimental-predicted pairs: {stats['total_pairs']}")
        
        # For backward compatibility, also create the filtered property_data without _smile
        data['property_data'] = {k.replace('_smile', ''): v for k, v in data['property_data'].items()}

    # --- Visualization 1: Opsins Overview Plot ---
    if args.only is None or args.only == 'overview':
        print("\n[VISUALIZATION 1] Generating Opsin Overview Plot...")
        try:
            overview_data_list = []
            processed_structures_map = data.get('processed_structures', {})
            property_data_map = data.get('property_data', {})

            for sid, struct_info in processed_structures_map.items():
                molecular_function = "Unknown"
                domain = "Unknown"
                is_experimental = False
                display_name = sid

                property_candidates = []
                if isinstance(struct_info, dict):
                    property_candidates.append(struct_info.get('properties', {}))
                property_candidates.append(property_data_map.get(sid, property_data_map.get(sid.replace('_model_0', ''), {})))
                if mapper is not None:
                    property_candidates.append(mapper.get_properties(sid) or {})

                for props in property_candidates:
                    if not isinstance(props, dict) or not props:
                        continue

                    mf_val = props.get('molecular_function') or props.get('molecular_function_normalized')
                    if pd.notna(mf_val) and str(mf_val).strip() and str(mf_val).lower() != "unknown":
                        molecular_function = str(mf_val)

                    domain_val = props.get('domain')
                    if pd.notna(domain_val) and str(domain_val).strip() and str(domain_val).lower() != "unknown":
                        domain = str(domain_val)

                    if 'experimentally_determined' in props:
                        is_experimental = bool(props.get('experimentally_determined'))

                    display_name_prop = props.get('short_name', props.get('opsin_name', display_name))
                    if pd.notna(display_name_prop) and str(display_name_prop).strip():
                        display_name = str(display_name_prop)

                    if molecular_function.lower() != 'unknown' and domain.lower() != 'unknown' and display_name:
                        break

                overview_data_list.append({
                    'id': sid,
                    'short_name': display_name,
                    'molecular_function': molecular_function,
                    'domain': domain,
                    'experimentally_determined': is_experimental,
                    'is_predicted': not is_experimental
                })
            overview_df = pd.DataFrame(overview_data_list)

            if not overview_df.empty:
                # Update 'experimentally_determined' based on structure type
                # Use the mapper to properly identify experimental vs predicted structures
                for idx, row in overview_df.iterrows():
                    struct_id = row['id']
                    struct_type, _ = mapper.identify_structure_type(struct_id)
                    # Experimental structures are standard_exp or hideaki_exp
                    is_exp = struct_type in ['standard_exp', 'hideaki_exp']
                    overview_df.at[idx, 'experimentally_determined'] = is_exp

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

    # --- Visualization 2: Error Distribution Violin Plot ---
    if args.only is None or args.only == 'errors':
        print("\n[VISUALIZATION 2] Generating Error Distribution Plot...")
        try:
            set_a_path = Path(args.input_dir) / 'set_a_errors.csv'
            set_b_path = Path(args.input_dir) / 'set_b_errors.csv'

            error_frames = []
            if set_a_path.exists():
                df_a = pd.read_csv(set_a_path)
                df_a['dataset_split'] = df_a.get('dataset_split', 'A')
            else:
                print(f"[WARN] Set A error file not found: {set_a_path}")
                df_a = pd.DataFrame()

            if set_b_path.exists():
                df_b = pd.read_csv(set_b_path)
                df_b['dataset_split'] = df_b.get('dataset_split', 'B')
            else:
                print(f"[WARN] Set B error file not found: {set_b_path}")
                df_b = pd.DataFrame()

            if df_a.empty and df_b.empty:
                print("[WARN] No error CSV files available. Skipping comparison plot.")
            else:
                comparison_path = FIGURES_OUTPUT_DIR / "02c_error_distribution_box.png"
                fig_errors, summary = plot_error_box_comparison(
                    df_a,
                    df_b,
                    metrics=['backbone_rmsd', 'pocket_rmsd', 'retinal_rmsd'],
                    dataset_labels=('Benchmark set', 'Blind test set'),
                    output_path=comparison_path
                )
                plt.close(fig_errors)
                print(f"[SUCCESS] Saved error comparison plot to: {comparison_path}")
                if not summary.empty:
                    print("Summary statistics (Å):")
                    print(summary.round(3))
        except Exception as e:
            print(f"[ERROR] Failed to generate error distribution plot: {e}")
            traceback.print_exc()

    # --- Prepare group_dict and domain_dict as per notebook ---
    group_dict = {}
    domain_dict = {}  # This will store {'domain': name, 'average_error': val}

    # Use 'overview_data_list' as it was constructed with latest property info
    if 'overview_data_list' not in locals():
        overview_data_list = []
        
    for item in overview_data_list:  # Ensure overview_data_list is used
        group_dict[item['id']] = item['molecular_function']
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

    # --- Visualization 3: RMSD Heatmap with Improved Visualization (Clustermap) ---
    if args.only is None or args.only == 'rmsd':
        print("\n[VISUALIZATION 3] Generating RMSD Clustermap...")
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
                Z_linkage = linkage(condensed_matrix, method='weighted')

                # Compute linkage metrics

                metrics = compute_rmsd_metrics(
                    rmsd_df=rmsd_df,
                    linkage_matrix=Z_linkage,
                    thresholds=(2.0, 2.5, 3.0),
                    n_clusters=2,
                    outdir=FIGURES_OUTPUT_DIR / "rmsd_metrics"
                )

                # Attach function/domain annotations if you have dicts
                names = rmsd_df.index.to_list()
                labels = metrics["labels"]

                def get_annotation(name):
                    func = group_dict.get(name, "NA")
                    dom = domain_dict.get(name, {}).get("domain", "NA")
                    return func, dom

                print("\nTop 3 reference-like opsins (lowest mean cross-cluster RMSD):")
                for name, val in metrics["ref_like_top10"][:3]:
                    idx = names.index(name)
                    cluster = labels[idx]
                    func, dom = get_annotation(name)
                    print(
                        f"  {name} | cross-cluster RMSD = {val:.2f} Å | cluster {cluster} | function={func} | domain={dom}")

                print("\nTop 3 outlier opsins (highest mean global RMSD):")
                for name, val in metrics["outliers_top10"][:3]:
                    idx = names.index(name)
                    cluster = labels[idx]
                    func, dom = get_annotation(name)
                    print(f"  {name} | global RMSD = {val:.2f} Å | cluster {cluster} | function={func} | domain={dom}")

                # Call the improved visualization function
                # Assuming visualize_rmsd_matrix_improved handles its own saving if output_file is passed
                fig2_path = FIGURES_OUTPUT_DIR / "02a_rmsd_clustermap.png"
                fig2 = visualize_rmsd_matrix_improved(
                    rmsd_df=rmsd_df,
                    group_dict=group_dict,
                    domain_dict=domain_dict,  # Pass the dict with 'domain' and 'average_error' keys
                    linkage_matrix=Z_linkage,
                    figsize=(18, 15),  # Adjusted figsize
                    output_file=fig2_path  # Pass outputs file directly
                )

                if fig2 and hasattr(fig2, "ax_heatmap"):
                    # Ensure saving (if your visualize function didn’t)
                    fig2.fig.savefig(fig2_path, dpi=300)
                    plt.close(fig2.fig)
                elif fig2:  # raw Figure
                    fig2.savefig(fig2_path, dpi=300)
                    plt.close(fig2)
                else:
                    print("[WARN] RMSD Clustermap was not generated.")

                fig2_step_path = FIGURES_OUTPUT_DIR / "02b_rmsd_clustermap_step.png"
                fig2_step = visualize_rmsd_matrix_improved(
                    rmsd_df=rmsd_df,
                    group_dict=group_dict,
                    domain_dict=domain_dict,
                    linkage_matrix=Z_linkage,
                    figsize=(18, 15),
                    output_file=fig2_step_path,
                    color_mode='step',
                    step_cutoffs=[0.5, 1.5, 2.5]
                )

                if fig2_step and hasattr(fig2_step, "ax_heatmap"):
                    fig2_step.fig.savefig(fig2_step_path, dpi=300)
                    plt.close(fig2_step.fig)
                elif fig2_step:
                    fig2_step.savefig(fig2_step_path, dpi=300)
                    plt.close(fig2_step)
                else:
                    print("[WARN] RMSD Clustermap (step) was not generated.")
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

    if residue_table_for_logo.empty:
        curated_grn_path = Path(args.input_dir) / 'curated_grn.csv'
        if curated_grn_path.exists():
            try:
                residue_table_for_logo = pd.read_csv(curated_grn_path, index_col=0)
                print(f"[INFO] Loaded residue table for logo plots from '{curated_grn_path}'.")
            except Exception as exc:
                print(f"[ERROR] Failed to load curated GRN table from '{curated_grn_path}': {exc}")

    # Initialize filtered versions
    df_filtered_residue = pd.DataFrame()
    df_filtered_distance = pd.DataFrame()
    df_filtered_ca_distance = pd.DataFrame()
    grns = []

    if not residue_table_for_logo.empty:
        col_name_filter = None
        filter_candidates = ['7.50', '7.5']
        for candidate in filter_candidates:
            if candidate in residue_table_for_logo.columns:
                col_name_filter = candidate
                break

        if col_name_filter is None:
            fallback_candidates = [
                col for col in residue_table_for_logo.columns
                if str(col).startswith('7.5')
            ]
            if fallback_candidates:
                col_name_filter = sorted(fallback_candidates, key=str)[0]

        if col_name_filter and col_name_filter in residue_table_for_logo.columns:
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
            print("[WARN] No suitable filter column (e.g., '7.5') found in residue_table_for_logo.")
    else:
        print("[WARN] residue_table_for_logo is empty. Cannot perform filtering for distance/logo plots.")

    # --- Visualization 4: All-Atom Distance to Retinal Plot (was Vis 4) ---
    if args.only is None or args.only == 'distance':
        print("\n[VISUALIZATION 4] Generating All-Atom Distance to Retinal Plot...")
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

    # --- Visualization 5: CA-Atom Distance to Retinal Plot (was Vis 5) ---
    if args.only is None or args.only == 'distance':
        print("\n[VISUALIZATION 5] Generating CA-Atom Distance to Retinal Plot...")
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

    # --- Visualization 6: Helix Logo Plots (was Vis 7) ---
    if args.only is None or args.only == 'logo':
        print("\n[VISUALIZATION 6] Generating Helix Logo Plots...")
        if not df_filtered_residue.empty:  # Use df_filtered_residue as per notebook for logos
            try:
                fig5 = plot_helix_logo_plots(
                    df_filtered_residue,
                    frequency_threshold=0.07  # As per notebook
                )
                fig5_path = FIGURES_OUTPUT_DIR / "05_helix_logos_x50.png"
                fig5.savefig(fig5_path, dpi=300);
                fig5_path_b = FIGURES_OUTPUT_DIR / "05b_helix_logos_x50.png"
                fig5.savefig(fig5_path_b, dpi=300);
                plt.close(fig5)
                print(f"[SUCCESS] Saved Helix Logo Plots to: {fig5_path}")
                print(f"[SUCCESS] Saved Helix Logo Plots (05b layout) to: {fig5_path_b}")
            except Exception as e:
                print(f"[ERROR] Failed to generate Helix Logo Plots: {e}")
                traceback.print_exc()
        else:
            print("[WARN] Filtered residue data for logos is empty. Skipping Helix Logo Plots.")

    # --- Additional Visualizations ---
    
    # Property Analysis Visualizations
    if not args.skip_property_analysis and (args.only is None or args.only == 'property'):
        if 'processed_structures' in data and property_csv_path.exists():
            try:
                print("\n[VISUALIZATION 7-8] Generating Property Analysis Visualizations...")
                create_property_analysis_visualizations(mapper, data['processed_structures'], FIGURES_OUTPUT_DIR)
            except Exception as e:
                print(f"[ERROR] Failed to generate property analysis visualizations: {e}")
                traceback.print_exc()
        else:
            print("[WARN] Skipping property analysis - missing required data")
    
    # Interactive GRN Visualization
    if not args.skip_interactive and (args.only is None or args.only == 'interactive'):
        try:
            print("\n[VISUALIZATION 9] Generating Interactive GRN Alignment Visualization...")
            # Import the interactive visualization module
            
            interactive_output = FIGURES_OUTPUT_DIR / 'interactive_grn_alignment.html'
            fig_interactive = create_opsin_visualization_from_workflow(
                cache_dir=str(CACHE_DIR),
                property_file=str(property_csv_path),
                output_file=str(interactive_output),
                max_structures=100,
                show_membrane=True,
                membrane_opacity=0.05
            )
            
            if fig_interactive:
                print(f"[SUCCESS] Saved interactive GRN visualization to: {interactive_output}")
            else:
                print("[WARN] Interactive visualization generation failed")
        except ImportError:
            print("[WARN] Interactive visualization module not available")
        except Exception as e:
            print(f"[ERROR] Failed to generate interactive visualization: {e}")
            traceback.print_exc()
    
    print("\n" + "=" * 30 + " VISUALIZATION SCRIPT COMPLETE " + "=" * 30)
    print(f"\nGenerated visualizations saved to: {FIGURES_OUTPUT_DIR}")
    print("\nVisualization summary:")
    print("  1. Opsin overview plot")
    print("  2. RMSD clustermap")
    print("  3. All-atom distance plot")
    print("  4. CA-atom distance plot")
    print("  5. Helix logo plots")
    if not args.skip_property_analysis:
        print("  6. Property analysis (natural domains)")
        print("  7. Prediction contribution chart")
    if not args.skip_interactive:
        print("  8. Interactive GRN alignment (HTML)")


if __name__ == "__main__":
    main()
