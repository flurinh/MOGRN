import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import json
from pathlib import Path
import pickle

# Import core workflow functions from modular structure
# Error Analysis
from src.error_analysis import calculate_structure_errors

# Structure Comparison
from src.structure_comparison import (
    compare_structures,
    create_unified_structure_mapping
)

# Data Processing Functions
from src.data_processing import (
    load_opsin_property_data,
    load_opsin_structures
)

# Helix Analysis Functions
from src.helix_analysis import align_to_reference_and_annotate_helices

# GRN Assignment
from src.assign_grns import align_and_assign_grn


def run_opsin_analysis_workflow(output_dir=None,
                                visualize=True,
                                use_foldmason=True,
                                use_cache=True,
                                cache_raw=True,
                                property_file=None,
                                chain_id='A',
                                retinal_name='RET',
                                retinal_cutoff=6.0,
                                global_ref_override=None,
                                helices_file='property/helices_curated.json'):
    """
    Run the full opsin analysis workflow with standardized dataset handling
    without using PropertyProcessor
    
    Args:
        output_dir: Directory to save output files. If None, uses projects/opsin_analysis/opsin_output.
        visualize: Whether to generate visualizations
        use_foldmason: Whether to use FoldMason for structure alignment
        use_cache: Whether to use cached structure data (default: True)
        cache_raw: Whether to cache raw unfiltered data (default: True)
        property_file: Path to CSV file with property data (optional)
        chain_id: Chain ID to use for analysis (default: 'A')
        retinal_name: Name of retinal residue (default: 'RET')
        retinal_cutoff: Distance cutoff in Angstroms for retinal selection (default: 6.0)
        global_ref_override: Optional override for global reference structure ID
        helices_file: Path to JSON file containing helix boundaries (default: 'property/helices_curated.json')
        
    Returns:
        Dictionary with analysis results
    """
    # Step 1: Load and process structures
    print("\n" + "=" * 80)
    print("RUNNING OPSIN ANALYSIS WORKFLOW")
    print("=" * 80 + "\n")
    
    # Set default directories
    if output_dir is None:
        output_dir = os.path.join(os.path.dirname(__file__), 'opsin_output')
    
    # Define data directory (where PROTOS data lives)
    data_dir = os.path.join(os.path.dirname(__file__), 'data')
    
    # Create output directory and cache directory (for analysis results and figures)
    os.makedirs(output_dir, exist_ok=True)
    cache_dir = os.path.join(output_dir, 'cache')
    os.makedirs(cache_dir, exist_ok=True)
    
    # Make sure data directory exists
    os.makedirs(data_dir, exist_ok=True)
    
    print(f"Using output directory for results and cache: {output_dir}")
    print(f"Using data directory for PROTOS data: {data_dir}")
    
    # Make sure we have the necessary structure directory in data path
    structure_dir = os.path.join(data_dir, 'structure')
    os.makedirs(structure_dir, exist_ok=True)
    
    # Create output structure directory for results if needed
    structure_output_dir = os.path.join(output_dir, 'structure')
    os.makedirs(structure_output_dir, exist_ok=True)
    
    # Initialize data dictionary to store all results
    data = {}
    chain_id = 'A'
    
    # Step 1: Load opsin structures from datasets (with two-stage caching)
    # This step has its own internal caching mechanism
    print(f"Calling load_opsin_structures with data_dir={data_dir}, output_dir={output_dir}")
    data = load_opsin_structures(data_dir=data_dir, output_dir=output_dir, chain_id=chain_id,
                                 visualize=visualize, use_cache=use_cache, cache_raw=cache_raw,
                                 retinal_name=retinal_name, retinal_cutoff=retinal_cutoff)
    
    # Step 2: Load property data if provided - cache the results
    property_data = None
    property_data = load_opsin_property_data("property/mo_exp.csv", data['processed_structures'])
    print("Property data loaded:", property_data)
    
    # Step 3: Create a unified structure mapping
    structure_mapping = None
    # If no cache or cache failed, create the mapping
    if structure_mapping is None:
        print("Creating unified structure mapping...")
        structure_mapping = create_unified_structure_mapping(data, property_data)
        data['structure_mapping'] = structure_mapping
    
    print(f"\nUnified structure mapping with {len(structure_mapping)} pairs")
    print(f"This mapping will be used throughout the workflow for all structure comparisons")
    
    # Step 4a: Calculate errors between experimental and predicted structures
    errors_cache_path = os.path.join(cache_dir, f"structure_errors_{chain_id}.pkl")
    errors_data = None
    
    if use_cache and os.path.exists(errors_cache_path):
        # Try to load from cache first
        try:
            print(f"Loading structure errors from cache: {errors_cache_path}")
            with open(errors_cache_path, 'rb') as f:
                errors_data = pickle.load(f)
                data.update(errors_data)
                print(f"Loaded structure errors from cache")
        except Exception as e:
            print(f"Warning: Failed to load structure errors from cache: {e}")
            errors_data = None
    
    # If no cache or cache failed, calculate errors
    if errors_data is None:
        print("Calculating structure errors...")
        errors_data = calculate_structure_errors(data, output_dir=output_dir, visualize=visualize)
        data.update(errors_data)
        
        # Cache the errors data
        if use_cache:
            try:
                with open(errors_cache_path, 'wb') as f:
                    pickle.dump(errors_data, f)
                print(f"Saved structure errors to cache: {errors_cache_path}")
            except Exception as e:
                print(f"Warning: Failed to cache structure errors: {e}")
    
    # Step 4b: Align structures to reference and annotate helices
    helix_cache_path = os.path.join(cache_dir, f"helix_annotations_{chain_id}.pkl")
    helix_data = None
    
    if use_cache and os.path.exists(helix_cache_path):
        # Try to load from cache first
        try:
            print(f"Loading helix annotations from cache: {helix_cache_path}")
            with open(helix_cache_path, 'rb') as f:
                helix_data = pickle.load(f)
                data.update(helix_data)
                print(f"Loaded helix annotations from cache")
        except Exception as e:
            print(f"Warning: Failed to load helix annotations from cache: {e}")
            helix_data = None
    
    # If no cache or cache failed, perform alignment and annotation
    if helix_data is None:
        print("Aligning structures and annotating helices...")
        helix_data = align_to_reference_and_annotate_helices(data, output_dir, visualize=visualize)
        data.update(helix_data)
        
        # Cache the helix data
        if use_cache:
            try:
                with open(helix_cache_path, 'wb') as f:
                    pickle.dump(helix_data, f)
                print(f"Saved helix annotations to cache: {helix_cache_path}")
            except Exception as e:
                print(f"Warning: Failed to cache helix annotations: {e}")
    
    # Step 5: Structure comparison
    comparison_cache_path = os.path.join(cache_dir, f"structure_comparison_{chain_id}.pkl")
    comparison_data = None
    
    if use_cache and os.path.exists(comparison_cache_path):
        # Try to load from cache first
        try:
            print(f"Loading structure comparison from cache: {comparison_cache_path}")
            with open(comparison_cache_path, 'rb') as f:
                comparison_data = pickle.load(f)
                data.update(comparison_data)
                print(f"Loaded structure comparison from cache")
        except Exception as e:
            print(f"Warning: Failed to load structure comparison from cache: {e}")
            comparison_data = None
    
    # If no cache or cache failed, perform comparison
    if comparison_data is None:
        print("Comparing structures...")
        comparison_data = compare_structures(data, output_dir, visualize=False)
        data.update(comparison_data)
        
        # Cache the comparison data
        if use_cache:
            try:
                with open(comparison_cache_path, 'wb') as f:
                    pickle.dump(comparison_data, f)
                print(f"Saved structure comparison to cache: {comparison_cache_path}")
            except Exception as e:
                print(f"Warning: Failed to cache structure comparison: {e}")
    
    # Step 6: Structure alignment and GRN assignment
    grn_cache_path = os.path.join(cache_dir, f"grn_assignment_{chain_id}.pkl")
    grn_data = None
    
    # If no cache or cache failed, perform GRN assignment
    if grn_data is None:
        print("Assigning GRNs...")
        grn_data = align_and_assign_grn(data, output_dir, visualize=False,
                                       global_ref_override=global_ref_override,
                                       helices_file=helices_file)
        data.update(grn_data)
        
        # Cache the GRN data
        if use_cache:
            try:
                with open(grn_cache_path, 'wb') as f:
                    pickle.dump(grn_data, f)
                print(f"Saved GRN assignment to cache: {grn_cache_path}")
            except Exception as e:
                print(f"Warning: Failed to cache GRN assignment: {e}")
    
    # Save a summary of the results
    try:
        summary = {
            "datasets": list(data.get('datasets', {}).keys()),
            "structures_count": len(data.get('processed_structures', {})),
            "exp_pred_pairs": len(data.get('structure_mapping', {})),
            "timestamp": pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S'),
            "cache_enabled": use_cache,
            "cache_dir": cache_dir
        }
        
        with open(os.path.join(output_dir, 'analysis_summary.json'), 'w') as f:
            json.dump(summary, f, indent=2)
    except Exception as e:
        print(f"Warning: Could not save summary: {e}")
    
    print("\n" + "=" * 80)
    print(f"ANALYSIS COMPLETE. Results saved to {output_dir}")
    print("=" * 80)
    
    return data


if __name__ == '__main__':
    import argparse
    
    # Set up command line arguments
    parser = argparse.ArgumentParser(description='Run opsin analysis workflow')
    parser.add_argument('--output-dir', type=str, default=None,
                        help='Directory to save output files. If not specified, uses projects/opsin_analysis/opsin_output')
    parser.add_argument('--no-visualize', action='store_false', dest='visualize',
                        help='Disable visualization generation')
    parser.add_argument('--no-foldmason', action='store_false', dest='use_foldmason',
                        help='Disable FoldMason alignment')
    parser.add_argument('--no-cache', action='store_false', dest='use_cache',
                        help='Disable structure caching')
    parser.add_argument('--no-raw-cache', action='store_false', dest='cache_raw',
                        help='Disable raw data caching')
    parser.add_argument('--property-file', type=str, default=None,
                        help='Path to CSV file with property data')
    parser.add_argument('--chain-id', type=str, default='A',
                        help='Chain ID to use for analysis (default: A)')
    parser.add_argument('--global-ref', type=str, default=None,
                        help='Override for global reference structure ID')
    parser.add_argument('--helices-file', type=str, default='property/helices_curated.json',
                        help='Path to JSON file containing helix boundaries')
    
    parser.set_defaults(visualize=False, use_foldmason=False, use_cache=True, cache_raw=True)
    
    # Parse arguments
    args = parser.parse_args()
    
    # Run the workflow
    results = run_opsin_analysis_workflow(
        output_dir=args.output_dir,
        visualize=args.visualize,
        use_foldmason=args.use_foldmason,
        use_cache=args.use_cache,
        cache_raw=args.cache_raw,
        property_file=args.property_file,
        chain_id=args.chain_id,
        global_ref_override=args.global_ref,
        helices_file=args.helices_file
    )