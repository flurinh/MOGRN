#!/usr/bin/env python3
"""
Test script for retinal carbon assignment with detailed reporting.
"""

import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from src.retinal_carbon_mapping import assign_retinal_carbons_correct
from src.data_processing import load_opsin_structures
import pandas as pd
import numpy as np


def format_table(df: pd.DataFrame, title: str = None) -> str:
    """Format a DataFrame as a text table."""
    if title:
        print(f"\n{title}")
        print("=" * len(title))

    # Convert DataFrame to string with formatting
    return df.to_string(index=False, max_rows=None, max_cols=None)


def get_retinal_atoms(structure_df: pd.DataFrame) -> pd.DataFrame:
    """Extract retinal atoms from structure DataFrame."""
    retinal_mask = structure_df['res_name3l'].isin(['RET', 'LYR'])
    retinal_df = structure_df[retinal_mask].copy()

    if 'LYR' in retinal_df['res_name3l'].values:
        lys_atoms = {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'}
        retinal_df = retinal_df[~retinal_df['res_atom_name'].isin(lys_atoms)]

    return retinal_df.reset_index(drop=True)


def analyze_carbon_assignment(struct_id: str, struct_data: dict) -> pd.DataFrame:
    """
    Analyze carbon assignment for a single structure and return detailed results.

    Returns:
        DataFrame with assignment details
    """
    # Extract retinal atoms
    df = struct_data['df']
    retinal_df = get_retinal_atoms(df)

    # Get carbon assignments from the algorithm
    try:
        assigned_carbons = assign_retinal_carbons_correct(retinal_df)
    except Exception as e:
        print(f"Error in {struct_id}: {e}")
        return pd.DataFrame()

    if not assigned_carbons:
        print(f"No assignments for {struct_id}")
        return pd.DataFrame()

    # Create a detailed mapping table
    results = []

    # Get all carbon atoms
    carbon_mask = retinal_df['atom_name'] == 'C'
    carbon_atoms = retinal_df[carbon_mask]

    for idx, atom in carbon_atoms.iterrows():
        original_idx = idx
        canonical_name = atom['res_atom_name']

        # Find what this carbon was assigned as
        assigned_as = None
        for carbon_label, assigned_idx in assigned_carbons.items():
            if carbon_label.startswith('C') and assigned_idx == original_idx:
                assigned_as = carbon_label
                break

        # Check if assignment matches canonical
        matches = (assigned_as == canonical_name) if assigned_as else False

        # Get position
        x, y, z = atom['x'], atom['y'], atom['z']

        results.append({
            'Structure': struct_id,
            'Index': original_idx,
            'Canonical': canonical_name,
            'Assigned': assigned_as or 'UNASSIGNED',
            'Match': '✓' if matches else '✗' if assigned_as else '-',
            'X': f"{x:.2f}",
            'Y': f"{y:.2f}",
            'Z': f"{z:.2f}"
        })

    # Sort by canonical name (C1, C2, ..., C20)
    results_df = pd.DataFrame(results)
    if not results_df.empty:
        # Custom sort function for carbon names
        def sort_carbon_name(name):
            if name.startswith('C') and name[1:].isdigit():
                return int(name[1:])
            return 999  # Put non-standard names at the end

        results_df['sort_key'] = results_df['Canonical'].apply(sort_carbon_name)
        results_df = results_df.sort_values('sort_key').drop('sort_key', axis=1)

    return results_df


def create_summary_report(all_results: pd.DataFrame) -> pd.DataFrame:
    """Create a summary report of assignment accuracy across structures."""
    summary = []

    for struct_id in all_results['Structure'].unique():
        struct_data = all_results[all_results['Structure'] == struct_id]

        total_carbons = len(struct_data)
        assigned = len(struct_data[struct_data['Assigned'] != 'UNASSIGNED'])
        correct = len(struct_data[struct_data['Match'] == '✓'])

        accuracy = (correct / assigned * 100) if assigned > 0 else 0
        completeness = (assigned / total_carbons * 100) if total_carbons > 0 else 0

        summary.append({
            'Structure': struct_id,
            'Total C': total_carbons,
            'Assigned': assigned,
            'Correct': correct,
            'Accuracy %': f"{accuracy:.1f}",
            'Complete %': f"{completeness:.1f}"
        })

    return pd.DataFrame(summary)


def analyze_mismatches(all_results: pd.DataFrame):
    """Analyze patterns in mismatched assignments."""
    mismatches = all_results[all_results['Match'] == '✗']

    if mismatches.empty:
        print("\nNo mismatches found!")
        return

    print("\n=== Mismatch Analysis ===")

    # Group by canonical carbon
    mismatch_patterns = {}
    for _, row in mismatches.iterrows():
        canonical = row['Canonical']
        assigned = row['Assigned']

        if canonical not in mismatch_patterns:
            mismatch_patterns[canonical] = {}

        if assigned not in mismatch_patterns[canonical]:
            mismatch_patterns[canonical][assigned] = 0

        mismatch_patterns[canonical][assigned] += 1

    # Print patterns
    for canonical, assignments in sorted(mismatch_patterns.items()):
        print(f"\n{canonical} was misassigned as:")
        for assigned, count in sorted(assignments.items()):
            print(f"  {assigned}: {count} times")


def main():
    """Run detailed carbon assignment analysis."""

    print("Loading opsin structures...")

    # Set up paths
    output_dir = os.path.join(os.path.dirname(__file__), 'opsin_output')
    data_dir = os.path.join(os.path.dirname(__file__), 'data')

    # Load structures
    data = load_opsin_structures(
        data_dir=data_dir,
        output_dir=output_dir,
        chain_id='A',
        use_cache=True,
        visualize=False
    )

    processed_structures = data.get('processed_structures', {})

    if not processed_structures:
        print("No structures loaded!")
        return

    print(f"\nLoaded {len(processed_structures)} structures")

    # Analyze all structures
    all_results = []

    # Process structures in batches for better output
    struct_ids = list(processed_structures.keys())
    batch_size = 5

    for i in range(0, min(20, len(struct_ids)), batch_size):  # Analyze first 20 structures
        batch = struct_ids[i:i + batch_size]

        print(f"\n{'=' * 80}")
        print(f"Processing structures {i + 1}-{min(i + batch_size, len(struct_ids))}")
        print(f"{'=' * 80}")

        for struct_id in batch:
            print(f"\nAnalyzing {struct_id}...")
            results_df = analyze_carbon_assignment(struct_id, processed_structures[struct_id])

            if not results_df.empty:
                all_results.append(results_df)

                # Print detailed table for this structure
                print(f"\nCarbon assignments for {struct_id}:")
                print(format_table(results_df))

                # Quick summary
                matches = len(results_df[results_df['Match'] == '✓'])
                total = len(results_df)
                print(f"\nSummary: {matches}/{total} correct assignments ({matches / total * 100:.1f}%)")

    # Combine all results
    if all_results:
        all_results_df = pd.concat(all_results, ignore_index=True)

        # Create and print summary report
        print(f"\n{'=' * 80}")
        print("OVERALL SUMMARY")
        print(f"{'=' * 80}")

        summary_df = create_summary_report(all_results_df)
        print("\n" + format_table(summary_df, "Structure Summary"))

        # Analyze mismatches
        analyze_mismatches(all_results_df)

        # Save detailed results
        output_file = os.path.join(output_dir, 'carbon_assignment_detailed.csv')
        all_results_df.to_csv(output_file, index=False)
        print(f"\nDetailed results saved to: {output_file}")

        # Save summary
        summary_file = os.path.join(output_dir, 'carbon_assignment_summary.csv')
        summary_df.to_csv(summary_file, index=False)
        print(f"Summary saved to: {summary_file}")

        # Create a pivot table showing assignment patterns
        print(f"\n{'=' * 80}")
        print("ASSIGNMENT PATTERN MATRIX")
        print(f"{'=' * 80}")

        # Create confusion matrix
        confusion_data = []
        for _, row in all_results_df.iterrows():
            if row['Assigned'] != 'UNASSIGNED':
                confusion_data.append({
                    'Canonical': row['Canonical'],
                    'Assigned': row['Assigned']
                })

        if confusion_data:
            confusion_df = pd.DataFrame(confusion_data)
            confusion_matrix = pd.crosstab(confusion_df['Canonical'], confusion_df['Assigned'])

            # Only show carbons with issues
            problem_carbons = confusion_matrix.index[confusion_matrix.apply(
                lambda x: x.argmax() != x.name[1:] if x.name.startswith('C') and x.name[1:].isdigit() else True,
                axis=1)]

            if len(problem_carbons) > 0:
                print("\nCarbons with assignment issues:")
                print(confusion_matrix.loc[problem_carbons])
            else:
                print("\nAll carbons correctly assigned!")


if __name__ == "__main__":
    main()