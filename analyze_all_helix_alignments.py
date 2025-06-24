#!/usr/bin/env python3
"""
Analyze helix alignments for all proteins to find systematic patterns
"""

import pandas as pd
import json
from pathlib import Path
from typing import Dict, List, Tuple
from collections import defaultdict
import numpy as np

def load_grn_entry(df: pd.DataFrame, protein_id: str) -> Dict[str, str]:
    """Load a single GRN table entry and extract sequence info"""
    if protein_id not in df.index:
        return None
    
    # Get the row as a dict
    row = df.loc[protein_id]
    sequence_dict = {}
    
    for col, value in row.items():
        if pd.notna(value) and value != '-':
            sequence_dict[col] = value
    
    return sequence_dict

def extract_residue_info(grn_value: str) -> Tuple[str, int]:
    """Extract residue type and position from GRN value like 'K296'"""
    if not grn_value or grn_value == '-':
        return None, None
    
    residue = grn_value[0]
    position = int(grn_value[1:])
    return residue, position

def analyze_protein_alignment(sequence_dict: Dict[str, str], helix_boundaries: Dict[str, List[int]]) -> Dict:
    """Analyze alignment for a single protein"""
    
    helix_stats = {}
    
    for helix_num in ['1', '2', '3', '4', '5', '6', '7']:
        if helix_num not in helix_boundaries:
            continue
            
        start, end = helix_boundaries[helix_num]
        
        # Find all positions for this helix
        helix_positions = []
        for grn_pos, grn_value in sequence_dict.items():
            if grn_pos.startswith(f"{helix_num}."):
                residue, position = extract_residue_info(grn_value)
                if position is not None:
                    helix_positions.append(position)
        
        if not helix_positions:
            continue
            
        actual_start = min(helix_positions)
        actual_end = max(helix_positions)
        
        helix_stats[helix_num] = {
            'defined_start': start,
            'defined_end': end,
            'actual_start': actual_start,
            'actual_end': actual_end,
            'start_diff': actual_start - start,
            'end_diff': actual_end - end,
            'n_positions': len(helix_positions)
        }
    
    return helix_stats

def main():
    # Load data
    grn_table_path = "opsin_output/opsin_grn_tables/residue_table_grn.csv"
    helices_path = "property/helices_curated.json"
    
    print("Loading GRN table...")
    grn_df = pd.read_csv(grn_table_path, index_col=0)
    
    print("Loading helix boundaries...")
    with open(helices_path, 'r') as f:
        all_helices = json.load(f)
    
    # Analyze all proteins
    all_results = {}
    missing_proteins = []
    
    for protein_id in grn_df.index:
        sequence_dict = load_grn_entry(grn_df, protein_id)
        
        if protein_id not in all_helices:
            missing_proteins.append(protein_id)
            continue
        
        helix_boundaries = all_helices[protein_id]
        helix_stats = analyze_protein_alignment(sequence_dict, helix_boundaries)
        all_results[protein_id] = helix_stats
    
    print(f"\nAnalyzed {len(all_results)} proteins")
    print(f"Missing helix definitions for {len(missing_proteins)} proteins")
    
    # Aggregate statistics
    print("\n" + "="*80)
    print("AGGREGATE STATISTICS FOR HELIX BOUNDARY MISALIGNMENTS")
    print("="*80)
    
    helix_diffs = defaultdict(lambda: {'start_diffs': [], 'end_diffs': []})
    
    for protein_id, helix_stats in all_results.items():
        for helix_num, stats in helix_stats.items():
            helix_diffs[helix_num]['start_diffs'].append(stats['start_diff'])
            helix_diffs[helix_num]['end_diffs'].append(stats['end_diff'])
    
    # Print statistics for each helix
    for helix_num in ['1', '2', '3', '4', '5', '6', '7']:
        if helix_num not in helix_diffs:
            continue
            
        start_diffs = np.array(helix_diffs[helix_num]['start_diffs'])
        end_diffs = np.array(helix_diffs[helix_num]['end_diffs'])
        
        print(f"\nHelix {helix_num} (n={len(start_diffs)} proteins):")
        print(f"  Start differences:")
        print(f"    Mean: {np.mean(start_diffs):.1f}")
        print(f"    Median: {np.median(start_diffs):.0f}")
        print(f"    Range: [{np.min(start_diffs)}, {np.max(start_diffs)}]")
        print(f"    Std: {np.std(start_diffs):.1f}")
        
        print(f"  End differences:")
        print(f"    Mean: {np.mean(end_diffs):.1f}")
        print(f"    Median: {np.median(end_diffs):.0f}")
        print(f"    Range: [{np.min(end_diffs)}, {np.max(end_diffs)}]")
        print(f"    Std: {np.std(end_diffs):.1f}")
    
    # Find proteins with large deviations
    print("\n" + "="*80)
    print("PROTEINS WITH LARGE DEVIATIONS (>5 residues)")
    print("="*80)
    
    outliers = []
    for protein_id, helix_stats in all_results.items():
        max_deviation = 0
        worst_helix = None
        
        for helix_num, stats in helix_stats.items():
            deviation = max(abs(stats['start_diff']), abs(stats['end_diff']))
            if deviation > max_deviation:
                max_deviation = deviation
                worst_helix = helix_num
        
        if max_deviation > 5:
            outliers.append({
                'protein': protein_id,
                'max_deviation': max_deviation,
                'worst_helix': worst_helix,
                'details': helix_stats[worst_helix]
            })
    
    # Sort by deviation
    outliers.sort(key=lambda x: x['max_deviation'], reverse=True)
    
    print(f"\nFound {len(outliers)} proteins with deviations > 5 residues:")
    for i, outlier in enumerate(outliers[:10]):  # Show top 10
        print(f"\n{i+1}. {outlier['protein']}:")
        print(f"   Worst helix: {outlier['worst_helix']}")
        print(f"   Max deviation: {outlier['max_deviation']}")
        details = outlier['details']
        print(f"   Start diff: {details['start_diff']:+d}, End diff: {details['end_diff']:+d}")
    
    # Save detailed results
    output_data = {
        'summary_statistics': {
            helix: {
                'start_diff_mean': float(np.mean(helix_diffs[helix]['start_diffs'])),
                'start_diff_median': float(np.median(helix_diffs[helix]['start_diffs'])),
                'start_diff_std': float(np.std(helix_diffs[helix]['start_diffs'])),
                'end_diff_mean': float(np.mean(helix_diffs[helix]['end_diffs'])),
                'end_diff_median': float(np.median(helix_diffs[helix]['end_diffs'])),
                'end_diff_std': float(np.std(helix_diffs[helix]['end_diffs'])),
                'n_proteins': len(helix_diffs[helix]['start_diffs'])
            }
            for helix in helix_diffs
        },
        'all_results': all_results,
        'missing_proteins': missing_proteins,
        'outliers': outliers
    }
    
    with open('helix_alignment_analysis_all_proteins.json', 'w') as f:
        json.dump(output_data, f, indent=2)
    
    print(f"\nDetailed results saved to: helix_alignment_analysis_all_proteins.json")
    
    # Suggest corrections
    print("\n" + "="*80)
    print("SUGGESTED SYSTEMATIC CORRECTIONS")
    print("="*80)
    
    for helix_num in ['1', '2', '3', '4', '5', '6', '7']:
        if helix_num not in helix_diffs:
            continue
            
        start_median = int(np.median(helix_diffs[helix_num]['start_diffs']))
        end_median = int(np.median(helix_diffs[helix_num]['end_diffs']))
        
        if abs(start_median) >= 1 or abs(end_median) >= 1:
            print(f"\nHelix {helix_num}:")
            print(f"  Suggested start adjustment: {start_median:+d}")
            print(f"  Suggested end adjustment: {end_median:+d}")

if __name__ == "__main__":
    main()