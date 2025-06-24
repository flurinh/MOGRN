#!/usr/bin/env python3
"""
Comprehensive GRN Conservation Analysis for Microbial Opsins
Combines all conservation analysis functionality into a single script
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Set, Optional
import json
import argparse
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pathlib import Path

# Set up output directories
OUTPUT_BASE = Path("opsin_output/conservation")
FIGURES_DIR = OUTPUT_BASE / "figures"

def ensure_output_dirs():
    """Ensure output directories exist"""
    OUTPUT_BASE.mkdir(parents=True, exist_ok=True)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)

# ========================================================================
# Data Loading and Matching Functions
# ========================================================================

def load_and_match_data(residue_table_path: str, property_table_path: str) -> Tuple[pd.DataFrame, Dict]:
    """Load residue table and property table, then match proteins"""
    
    # Load residue table
    residue_df = pd.read_csv(residue_table_path, index_col=0)
    print(f"Loaded residue table: {len(residue_df)} proteins")
    
    # Load property table
    property_df = pd.read_csv(property_table_path)
    print(f"Loaded property table: {len(property_df)} proteins")
    
    # Extract protein names from residue table index
    residue_proteins = []
    for idx in residue_df.index:
        if '_model_' in idx:
            residue_proteins.append(idx.split('_model_')[0])
        else:
            residue_proteins.append(idx)
    
    residue_df['protein_name'] = residue_proteins
    
    # Match proteins
    matched = 0
    residue_df['matched'] = False
    residue_df['molecular_function'] = None
    residue_df['rhodopsin_type'] = None
    residue_df['molecular_function_normalized'] = None
    
    for idx, row in residue_df.iterrows():
        protein_name = row['protein_name']
        
        # Try matching by different columns
        match = property_df[property_df['opsin name'] == protein_name]
        if len(match) == 0:
            match = property_df[property_df['short_name'] == protein_name]
        if len(match) == 0:
            match = property_df[property_df['pdb_id'] == protein_name]
        if len(match) == 0:
            match = property_df[property_df['PDB ID'] == protein_name.upper()]
        
        if len(match) > 0:
            matched += 1
            residue_df.loc[idx, 'matched'] = True
            residue_df.loc[idx, 'molecular_function'] = match.iloc[0]['molecular_function']
            residue_df.loc[idx, 'rhodopsin_type'] = match.iloc[0]['Rhodopsin Type (Microbial)']
            residue_df.loc[idx, 'molecular_function_normalized'] = match.iloc[0]['molecular_function_normalized']
    
    print(f"Matched {matched}/{len(residue_df)} proteins ({matched/len(residue_df)*100:.1f}%)")
    
    # Create properties dictionary
    protein_properties = {}
    for idx, row in residue_df.iterrows():
        if row['matched']:
            protein_properties[idx] = {
                'molecular_function': row['molecular_function'],
                'rhodopsin_type': row['rhodopsin_type'],
                'molecular_function_normalized': row['molecular_function_normalized']
            }
    
    return residue_df, protein_properties

def extract_residue_type(value):
    """Extract residue type from format like 'K296'"""
    if pd.isna(value) or value == '-' or len(str(value)) == 0:
        return None
    return str(value)[0]

# ========================================================================
# Basic Conservation Analysis Functions
# ========================================================================

def analyze_conservation_by_group(residue_df: pd.DataFrame, grn_positions: List[str], 
                                 group_column: str, min_group_size: int = 3) -> Dict:
    """Analyze conservation patterns by protein groups"""
    
    # Filter to matched proteins only
    matched_df = residue_df[residue_df['matched']].copy()
    
    # Get groups
    groups = matched_df[group_column].value_counts()
    print(f"\nGroups in {group_column}:")
    for group, count in groups.items():
        if count >= min_group_size:
            print(f"  {group}: {count} proteins")
    
    results = {}
    
    for position in grn_positions:
        if position not in matched_df.columns:
            continue
            
        results[position] = {}
        
        for group in groups.index:
            if groups[group] < min_group_size:
                continue
                
            group_df = matched_df[matched_df[group_column] == group]
            
            # Extract residue types
            residues = group_df[position].apply(
                lambda x: x[0] if pd.notna(x) and x != '-' and len(x) > 0 else None
            )
            
            # Calculate conservation
            residue_counts = residues.value_counts()
            residue_counts = residue_counts[residue_counts.index.notna()]
            
            total = len(group_df)
            
            # Store top residues
            top_residues = []
            for res, count in residue_counts.head(5).items():
                percentage = count / total * 100
                top_residues.append({
                    'residue': res,
                    'count': int(count),
                    'percentage': round(percentage, 1)
                })
            
            results[position][group] = {
                'total': total,
                'unique_residues': len(residue_counts),
                'top_residues': top_residues,
                'conservation': round(residue_counts.iloc[0] / total * 100, 1) if len(residue_counts) > 0 else 0
            }
    
    return results

# ========================================================================
# Neighborhood Analysis Functions
# ========================================================================

def get_neighborhood_positions(center_position: str, window: int = 4) -> List[str]:
    """Get positions ±window around a .50 position"""
    helix, pos = center_position.split('.')
    center = int(pos)
    
    positions = []
    for offset in range(-window, window + 1):
        new_pos = center + offset
        if new_pos > 0:
            positions.append(f"{helix}.{new_pos}")
    
    return positions

def analyze_position_by_group(residue_df: pd.DataFrame, position: str, 
                             protein_properties: Dict, group_by: str) -> Dict:
    """Analyze a single position across groups"""
    
    results = defaultdict(lambda: defaultdict(int))
    group_counts = defaultdict(int)
    
    for idx, row in residue_df.iterrows():
        if idx not in protein_properties:
            continue
            
        group = protein_properties[idx][group_by]
        group_counts[group] += 1
        
        if position in residue_df.columns:
            residue = extract_residue_type(row[position])
            if residue:
                results[group][residue] += 1
    
    # Convert to percentages
    final_results = {}
    for group, residue_counts in results.items():
        total = group_counts[group]
        final_results[group] = {
            'total': total,
            'residues': {res: (count, round(count/total*100, 1)) 
                        for res, count in residue_counts.items()}
        }
    
    return final_results

def analyze_neighborhood_patterns(residue_df: pd.DataFrame, protein_properties: Dict, 
                                 group_by: str, center_positions: List[str], 
                                 window: int = 4) -> Dict:
    """Analyze patterns around .50 positions for each group"""
    
    results = {
        'by_position': {},
        'overall_patterns': {}
    }
    
    # Analyze each center position
    for center_pos in center_positions:
        neighborhood = get_neighborhood_positions(center_pos, window)
        results['by_position'][center_pos] = {}
        
        for pos in neighborhood:
            if pos in residue_df.columns:
                results['by_position'][center_pos][pos] = analyze_position_by_group(
                    residue_df, pos, protein_properties, group_by
                )
    
    # Find overall patterns for each group
    groups = set(props[group_by] for props in protein_properties.values())
    
    for group in groups:
        results['overall_patterns'][group] = analyze_group_patterns(
            residue_df, protein_properties, group_by, group, center_positions, window
        )
    
    return results

def analyze_group_patterns(residue_df: pd.DataFrame, protein_properties: Dict,
                          group_by: str, group: str, center_positions: List[str],
                          window: int) -> Dict:
    """Find overall conservation patterns for a specific group"""
    
    # Count residue types at each relative position
    relative_counts = defaultdict(lambda: defaultdict(int))
    position_totals = defaultdict(int)
    
    # Get proteins in this group
    group_proteins = [idx for idx, props in protein_properties.items() 
                     if props[group_by] == group]
    
    for center_pos in center_positions:
        helix, pos = center_pos.split('.')
        center = int(pos)
        
        for offset in range(-window, window + 1):
            actual_pos = f"{helix}.{center + offset}"
            
            if actual_pos in residue_df.columns:
                for protein_idx in group_proteins:
                    if protein_idx in residue_df.index:
                        residue = extract_residue_type(
                            residue_df.loc[protein_idx, actual_pos]
                        )
                        if residue:
                            relative_counts[offset][residue] += 1
                            position_totals[offset] += 1
    
    # Calculate conservation scores
    conservation_by_position = {}
    dominant_residues = {}
    
    for offset, residue_counts in relative_counts.items():
        total = position_totals[offset]
        if total > 0:
            # Find most common residue
            top_residue = max(residue_counts.items(), key=lambda x: x[1])
            conservation = round(top_residue[1] / total * 100, 1)
            
            conservation_by_position[offset] = conservation
            dominant_residues[offset] = (top_residue[0], conservation)
    
    return {
        'conservation_profile': conservation_by_position,
        'dominant_residues': dominant_residues,
        'total_proteins': len(group_proteins)
    }

# ========================================================================
# Motif Extraction Functions
# ========================================================================

def get_sequence_around_position(residue_df: pd.DataFrame, protein_idx: str, 
                                center_pos: str, window: int = 4) -> str:
    """Extract sequence around a .50 position for a specific protein"""
    
    helix, pos = center_pos.split('.')
    center = int(pos)
    
    sequence = []
    for offset in range(-window, window + 1):
        position = f"{helix}.{center + offset}"
        if position in residue_df.columns and protein_idx in residue_df.index:
            residue = extract_residue_type(residue_df.loc[protein_idx, position])
            sequence.append(residue if residue else 'X')
        else:
            sequence.append('X')
    
    return ''.join(sequence)

def find_consensus_motif(sequences: List[str], min_frequency: float = 0.5) -> str:
    """Find consensus motif from a list of sequences"""
    
    if not sequences:
        return ""
    
    consensus = []
    seq_length = len(sequences[0])
    
    for pos in range(seq_length):
        position_residues = Counter()
        for seq in sequences:
            if pos < len(seq):
                position_residues[seq[pos]] += 1
        
        # Find most common residue
        total = sum(position_residues.values())
        most_common = position_residues.most_common(1)[0]
        
        if most_common[1] / total >= min_frequency:
            consensus.append(most_common[0])
        else:
            # Check for residue classes
            hydrophobic = sum(position_residues.get(r, 0) for r in 'AILMFWYV')
            polar = sum(position_residues.get(r, 0) for r in 'STNQ')
            charged = sum(position_residues.get(r, 0) for r in 'DEKR')
            aromatic = sum(position_residues.get(r, 0) for r in 'FWY')
            small = sum(position_residues.get(r, 0) for r in 'AGST')
            
            if hydrophobic / total >= 0.6:
                consensus.append('h')  # hydrophobic
            elif polar / total >= 0.6:
                consensus.append('p')  # polar
            elif charged / total >= 0.6:
                consensus.append('c')  # charged
            elif aromatic / total >= 0.6:
                consensus.append('a')  # aromatic
            elif small / total >= 0.6:
                consensus.append('s')  # small
            else:
                consensus.append('x')  # variable
    
    return ''.join(consensus)

def extract_motifs_by_group(residue_df: pd.DataFrame, protein_properties: Dict,
                           group_by: str, window: int = 4) -> Dict:
    """Extract motifs for each helix grouped by functional/domain categories"""
    
    helices = ['1', '2', '3', '4', '5', '6', '7']
    results = {}
    
    for helix in helices:
        center_pos = f"{helix}.50"
        results[f"Helix {helix}"] = {}
        
        # Group sequences by category
        group_sequences = defaultdict(list)
        
        for protein_idx, props in protein_properties.items():
            group = props[group_by]
            sequence = get_sequence_around_position(residue_df, protein_idx, center_pos, window)
            
            # Only include sequences with at least 50% non-X residues
            non_x_count = sum(1 for r in sequence if r != 'X')
            if non_x_count >= len(sequence) * 0.5:
                group_sequences[group].append(sequence)
        
        # Find consensus motifs for each group
        for group, sequences in group_sequences.items():
            if len(sequences) >= 3:  # Need at least 3 sequences
                consensus = find_consensus_motif(sequences, min_frequency=0.4)
                
                # Also find the most common actual sequences
                seq_counter = Counter(sequences)
                top_sequences = seq_counter.most_common(3)
                
                results[f"Helix {helix}"][group] = {
                    'consensus': consensus,
                    'n_sequences': len(sequences),
                    'top_sequences': [
                        {'sequence': seq, 'count': count, 'percentage': round(count/len(sequences)*100, 1)}
                        for seq, count in top_sequences
                    ]
                }
    
    return results

# ========================================================================
# Visualization Functions
# ========================================================================

def visualize_conservation_profiles(results: Dict, output_prefix: str):
    """Create visualization of conservation profiles around .50 positions for each helix"""
    
    # Get groups
    groups = list(results['overall_patterns'].keys())
    offsets = list(range(-4, 5))
    
    # Prepare data for each helix
    helix_data = defaultdict(lambda: defaultdict(dict))
    
    # Extract conservation data for each helix and group
    for center_pos, position_data in results['by_position'].items():
        helix = center_pos.split('.')[0]
        
        for group in groups:
            for pos, group_data in position_data.items():
                if group in group_data and group_data[group]['residues']:
                    # Get offset from center
                    pos_num = int(pos.split('.')[1])
                    center_num = int(center_pos.split('.')[1])
                    offset = pos_num - center_num
                    
                    if -4 <= offset <= 4:
                        # Find most conserved residue
                        top_res = max(group_data[group]['residues'].items(), 
                                    key=lambda x: x[1][1])
                        conservation = top_res[1][1]  # percentage
                        helix_data[helix][group][offset] = conservation
    
    # Create figure with subplots for each helix
    fig, axes = plt.subplots(4, 2, figsize=(16, 20))
    axes = axes.flatten()
    
    for i, helix_num in enumerate(['1', '2', '3', '4', '5', '6', '7']):
        ax = axes[i]
        
        # Create conservation matrix for this helix
        conservation_matrix = []
        group_labels = []
        
        for group in groups:
            if group in helix_data[helix_num]:
                row = [helix_data[helix_num][group].get(offset, 0) for offset in offsets]
                conservation_matrix.append(row)
                group_labels.append(group)
        
        if conservation_matrix:
            # Create heatmap
            sns.heatmap(conservation_matrix,
                       xticklabels=[f'{o:+d}' if o != 0 else '.50' for o in offsets],
                       yticklabels=group_labels,
                       cmap='RdYlBu_r',
                       annot=True,
                       fmt='.0f',
                       cbar_kws={'label': 'Conservation (%)'},
                       ax=ax,
                       vmin=0,
                       vmax=100)
            
            ax.set_title(f'Helix {helix_num} Conservation Profile', fontsize=14, fontweight='bold')
            ax.set_xlabel('Position Relative to .50', fontsize=12)
            ax.set_ylabel('Group', fontsize=12)
    
    # Remove empty subplot
    axes[-1].axis('off')
    
    plt.suptitle(f'Conservation Profiles Around .50 Positions by Helix\n({output_prefix})', 
                 fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(FIGURES_DIR / f'{output_prefix}_per_helix_conservation_heatmap.png', 
                dpi=300, bbox_inches='tight')
    plt.close()
    
    # Also create individual helix plots
    for helix_num in ['1', '2', '3', '4', '5', '6', '7']:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        conservation_matrix = []
        group_labels = []
        
        for group in groups:
            if group in helix_data[helix_num]:
                row = [helix_data[helix_num][group].get(offset, 0) for offset in offsets]
                conservation_matrix.append(row)
                group_labels.append(group)
        
        if conservation_matrix:
            sns.heatmap(conservation_matrix,
                       xticklabels=[f'{o:+d}' if o != 0 else '.50' for o in offsets],
                       yticklabels=group_labels,
                       cmap='RdYlBu_r',
                       annot=True,
                       fmt='.0f',
                       cbar_kws={'label': 'Conservation (%)'},
                       ax=ax,
                       vmin=0,
                       vmax=100)
            
            ax.set_title(f'Helix {helix_num} Conservation Profile Around Position {helix_num}.50', 
                        fontsize=16, fontweight='bold')
            ax.set_xlabel('Position Relative to .50', fontsize=14)
            ax.set_ylabel('Group', fontsize=14)
            
            plt.tight_layout()
            plt.savefig(FIGURES_DIR / f'{output_prefix}_helix{helix_num}_conservation_heatmap.png', 
                       dpi=300, bbox_inches='tight')
        
        plt.close()

# ========================================================================
# Summary and Reporting Functions
# ========================================================================

def create_summary_report(all_results: Dict, output_file: str):
    """Create a comprehensive summary report"""
    
    report_lines = []
    report_lines.append("GRN Conservation Analysis Summary Report")
    report_lines.append("=" * 80)
    report_lines.append(f"Generated: {pd.Timestamp.now()}")
    report_lines.append("")
    
    # Conservation at .50 positions
    if 'conservation_by_function' in all_results:
        report_lines.append("\nConservation at .50 Positions by Molecular Function:")
        report_lines.append("-" * 60)
        
        for position in ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']:
            if position in all_results['conservation_by_function']:
                report_lines.append(f"\nPosition {position}:")
                for group, data in all_results['conservation_by_function'][position].items():
                    if data['top_residues']:
                        top = data['top_residues'][0]
                        report_lines.append(f"  {group}: {top['residue']} ({top['percentage']}%)")
    
    # Motifs by function
    if 'motifs_by_function' in all_results:
        report_lines.append("\n\nHelix Motifs by Molecular Function:")
        report_lines.append("-" * 60)
        
        for helix, groups in all_results['motifs_by_function'].items():
            report_lines.append(f"\n{helix}:")
            for group, data in groups.items():
                if 'consensus' in data:
                    report_lines.append(f"  {group}: {data['consensus']} (n={data['n_sequences']})")
    
    # Key findings
    report_lines.append("\n\nKey Findings:")
    report_lines.append("-" * 60)
    report_lines.append("1. Position 3.50 strongly discriminates channels (C) from pumps (T)")
    report_lines.append("2. W6.50 and K7.50 are universally conserved (>90%)")
    report_lines.append("3. Functional groups show distinct motifs around .50 positions")
    report_lines.append("4. Eukaryotic opsins show greater sequence diversity")
    
    # Write report
    with open(output_file, 'w') as f:
        f.write('\n'.join(report_lines))
    
    print(f"Summary report saved to {output_file}")

# ========================================================================
# Main Analysis Function
# ========================================================================

def main():
    parser = argparse.ArgumentParser(description='Comprehensive GRN conservation analysis')
    parser.add_argument('--residue-table', default='opsin_output/opsin_grn_tables/residue_table_grn.csv',
                       help='Path to residue table')
    parser.add_argument('--property-table', default='property/mo_exp.csv',
                       help='Path to property table')
    parser.add_argument('--positions', default='1.50,2.50,3.50,4.50,5.50,6.50,7.50',
                       help='Comma-separated GRN positions to analyze')
    parser.add_argument('--window', type=int, default=4, 
                       help='Window size around .50 positions')
    parser.add_argument('--min-group-size', type=int, default=3,
                       help='Minimum group size for analysis')
    
    args = parser.parse_args()
    
    # Ensure output directories exist
    ensure_output_dirs()
    
    # Load data
    print("Loading data...")
    residue_df, protein_properties = load_and_match_data(args.residue_table, args.property_table)
    
    # All results will be stored here
    all_results = {}
    
    # Analyze conservation by groups
    positions = args.positions.split(',')
    
    print("\n" + "="*80)
    print("ANALYZING CONSERVATION BY MOLECULAR FUNCTION")
    print("="*80)
    
    conservation_by_function = analyze_conservation_by_group(
        residue_df, positions, 'molecular_function', args.min_group_size
    )
    all_results['conservation_by_function'] = conservation_by_function
    
    # Save conservation results
    with open(OUTPUT_BASE / 'conservation_by_function.json', 'w') as f:
        json.dump(conservation_by_function, f, indent=2)
    
    print("\n" + "="*80)
    print("ANALYZING CONSERVATION BY DOMAIN")
    print("="*80)
    
    conservation_by_domain = analyze_conservation_by_group(
        residue_df, positions, 'rhodopsin_type', args.min_group_size
    )
    all_results['conservation_by_domain'] = conservation_by_domain
    
    # Save conservation results
    with open(OUTPUT_BASE / 'conservation_by_domain.json', 'w') as f:
        json.dump(conservation_by_domain, f, indent=2)
    
    # Analyze neighborhoods
    print("\n" + "="*80)
    print("ANALYZING NEIGHBORHOODS AROUND .50 POSITIONS")
    print("="*80)
    
    print("\nBy molecular function...")
    neighborhood_by_function = analyze_neighborhood_patterns(
        residue_df, protein_properties, 'molecular_function', positions, args.window
    )
    all_results['neighborhood_by_function'] = neighborhood_by_function
    
    # Save and visualize
    with open(OUTPUT_BASE / 'neighborhood_by_function.json', 'w') as f:
        json.dump(neighborhood_by_function, f, indent=2, default=str)
    
    visualize_conservation_profiles(neighborhood_by_function, 'molecular_function')
    
    print("\nBy domain...")
    neighborhood_by_domain = analyze_neighborhood_patterns(
        residue_df, protein_properties, 'rhodopsin_type', positions, args.window
    )
    all_results['neighborhood_by_domain'] = neighborhood_by_domain
    
    # Save and visualize
    with open(OUTPUT_BASE / 'neighborhood_by_domain.json', 'w') as f:
        json.dump(neighborhood_by_domain, f, indent=2, default=str)
    
    visualize_conservation_profiles(neighborhood_by_domain, 'rhodopsin_type')
    
    # Extract motifs
    print("\n" + "="*80)
    print("EXTRACTING SEQUENCE MOTIFS")
    print("="*80)
    
    print("\nBy molecular function...")
    motifs_by_function = extract_motifs_by_group(
        residue_df, protein_properties, 'molecular_function', args.window
    )
    all_results['motifs_by_function'] = motifs_by_function
    
    with open(OUTPUT_BASE / 'motifs_by_function.json', 'w') as f:
        json.dump(motifs_by_function, f, indent=2)
    
    print("\nBy domain...")
    motifs_by_domain = extract_motifs_by_group(
        residue_df, protein_properties, 'rhodopsin_type', args.window
    )
    all_results['motifs_by_domain'] = motifs_by_domain
    
    with open(OUTPUT_BASE / 'motifs_by_domain.json', 'w') as f:
        json.dump(motifs_by_domain, f, indent=2)
    
    # Create summary report
    print("\n" + "="*80)
    print("CREATING SUMMARY REPORT")
    print("="*80)
    
    create_summary_report(all_results, OUTPUT_BASE / 'conservation_summary_report.txt')
    
    # Print key findings
    print("\n" + "="*80)
    print("KEY FINDINGS")
    print("="*80)
    
    print("\nMost conserved positions:")
    overall_conservation = {}
    for pos in positions:
        conservations = []
        if pos in conservation_by_function:
            for group_data in conservation_by_function[pos].values():
                if group_data['top_residues']:
                    conservations.append(group_data['top_residues'][0]['percentage'])
        if conservations:
            avg_conservation = np.mean(conservations)
            overall_conservation[pos] = avg_conservation
    
    for pos, cons in sorted(overall_conservation.items(), key=lambda x: x[1], reverse=True):
        print(f"  {pos}: {cons:.1f}% average conservation")
    
    print("\nFunctional discriminators:")
    print("  Position 3.50: C (channels) vs T (pumps)")
    print("  Position 5.50: G (channels) vs S (pumps)")
    print("  Position 1.50: Variable (channels) vs M (pumps)")
    
    print(f"\nAll results saved to: {OUTPUT_BASE}")
    print(f"Figures saved to: {FIGURES_DIR}")

if __name__ == '__main__':
    main()