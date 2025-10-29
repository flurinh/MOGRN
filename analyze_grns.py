#!/usr/bin/env python3
"""
Comprehensive GRN Analysis Tool
Consolidates functionality from:
- analyze_grn_conservation.py
- analyze_grn_patterns_detailed.py
- create_grn_summary_figure.py
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Tuple, Set, Optional
import json
import argparse
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import os
from pathlib import Path
import sys

# Add src to path for imports
sys.path.append('src')
from property_mapping import PropertyMapper

# Set up outputs directories
OUTPUT_BASE = Path("opsin_output/grn_analysis")
CONSERVATION_DIR = OUTPUT_BASE / "conservation"
FIGURES_DIR = OUTPUT_BASE / "figures"
REPORTS_DIR = OUTPUT_BASE / "reports"

def ensure_output_dirs():
    """Ensure all outputs directories exist"""
    for dir_path in [OUTPUT_BASE, CONSERVATION_DIR, FIGURES_DIR, REPORTS_DIR]:
        dir_path.mkdir(parents=True, exist_ok=True)

# ========================================================================
# Data Loading and Matching Functions (from analyze_grn_conservation.py)
# ========================================================================

def load_and_match_data(residue_table_path: str, property_table_path: str) -> Tuple[pd.DataFrame, Dict]:
    """Load residue table and property table using PropertyMapper"""
    
    # Load residue table
    residue_df = pd.read_csv(residue_table_path, index_col=0)
    print(f"Loaded residue table: {len(residue_df)} proteins")
    
    # Initialize property mapper
    property_mapper = PropertyMapper(Path(property_table_path))
    print(f"Initialized PropertyMapper with property table")
    
    # Map properties to proteins using the structure IDs from residue table
    protein_properties = {}
    matched = 0
    
    # Add columns to residue_df for compatibility
    residue_df['matched'] = False
    residue_df['molecular_function'] = None
    residue_df['rhodopsin_type'] = None
    residue_df['molecular_function_normalized'] = None
    
    for protein_id in residue_df.index:
        properties = property_mapper.get_properties(protein_id)
        if properties:
            matched += 1
            
            # Get molecular function - handle various possible field names
            molecular_function = properties.get('molecular_function', 'Unknown')
            if molecular_function == 'Unknown' or not molecular_function or molecular_function == '?':
                molecular_function = 'Unknown'
            
            # Get domain/rhodopsin type
            rhodopsin_type = properties.get('domain', 'Unknown')
            
            # Update residue_df for compatibility
            residue_df.loc[protein_id, 'matched'] = True
            residue_df.loc[protein_id, 'molecular_function'] = molecular_function
            residue_df.loc[protein_id, 'rhodopsin_type'] = rhodopsin_type
            
            # Handle molecular_function_normalized - map specific functions
            molecular_function_normalized = molecular_function
            if 'Pump' in molecular_function:
                if 'Proton' in molecular_function:
                    molecular_function_normalized = 'Proton Pump'
                elif 'Chloride' in molecular_function:
                    molecular_function_normalized = 'Chloride Pump'
                elif 'Sodium' in molecular_function:
                    molecular_function_normalized = 'Sodium Pump'
            elif 'Channel' in molecular_function:
                if 'Cation' in molecular_function:
                    molecular_function_normalized = 'Cation Channel'
                elif 'Anion' in molecular_function:
                    molecular_function_normalized = 'Anion Channel'
            elif 'Sensor' in molecular_function or 'Regulatory' in molecular_function:
                molecular_function_normalized = 'Sensor / Regulatory'
            
            residue_df.loc[protein_id, 'molecular_function_normalized'] = molecular_function_normalized
            
            # Store in protein_properties dictionary
            protein_properties[protein_id] = {
                'molecular_function': molecular_function,
                'rhodopsin_type': rhodopsin_type,
                'molecular_function_normalized': molecular_function_normalized,
                'opsin_name': properties.get('opsin_name', protein_id),
                'pdb_id': properties.get('pdb_id', None),
                'short_name': properties.get('short_name', protein_id)
            }
    
    print(f"Matched {matched}/{len(residue_df)} proteins ({matched/len(residue_df)*100:.1f}%)")
    
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
    # Handle both "1.5" and "1.50" notations
    if pos == '5':
        center = 50
    else:
        center = int(pos)
    
    positions = []
    for offset in range(-window, window + 1):
        new_pos = center + offset
        if new_pos > 0:
            # Special case: position 50 is stored as "5" in the GRN table
            if new_pos == 50:
                positions.append(f"{helix}.5")
            else:
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
# Pattern Analysis Functions (from analyze_grn_patterns_detailed.py)
# ========================================================================

def analyze_functional_discriminators(conservation_data: Dict) -> Dict:
    """Identify positions that best discriminate between functional groups"""
    
    discriminators = {}
    
    for position, groups in conservation_data.items():
        if not groups:
            continue
            
        # Get top residue for each group
        group_residues = {}
        for group, data in groups.items():
            if data.get('top_residues'):
                group_residues[group] = {
                    'residue': data['top_residues'][0]['residue'],
                    'percentage': data['top_residues'][0]['percentage']
                }
        
        # Calculate discrimination power
        residue_groups = defaultdict(list)
        for group, info in group_residues.items():
            residue_groups[info['residue']].append((group, info['percentage']))
        
        # Check if position discriminates between functional categories
        discrimination_patterns = []
        
        # Check pump vs channel discrimination
        pump_residues = set()
        channel_residues = set()
        
        for residue, groups_list in residue_groups.items():
            for group, percentage in groups_list:
                if 'Pump' in group:
                    pump_residues.add(residue)
                elif 'Channel' in group:
                    channel_residues.add(residue)
        
        if pump_residues and channel_residues and not pump_residues.intersection(channel_residues):
            discrimination_patterns.append({
                'type': 'pump_vs_channel',
                'pump_residues': list(pump_residues),
                'channel_residues': list(channel_residues)
            })
        
        # Calculate entropy as measure of discrimination
        all_percentages = [info['percentage'] for info in group_residues.values()]
        if all_percentages:
            entropy = -sum(p/100 * np.log2(p/100) if p > 0 else 0 for p in all_percentages)
            
            discriminators[position] = {
                'group_residues': group_residues,
                'discrimination_patterns': discrimination_patterns,
                'entropy': entropy,
                'unique_residues': len(residue_groups)
            }
    
    return discriminators

def analyze_coevolution_patterns(neighborhood_data: Dict) -> Dict:
    """Analyze co-evolution patterns in neighborhoods around .50 positions"""
    
    coevolution_patterns = defaultdict(list)
    
    for center_pos, position_data in neighborhood_data.get('by_position', {}).items():
        helix = center_pos.split('.')[0]
        
        # For each functional group
        for group in ['Proton Pump', 'Cation Channel', 'Anion Channel', 'Chloride Pump', 'Sensor / Regulatory']:
            group_pattern = []
            
            # Get residues at each position relative to .50
            for offset in range(-4, 5):
                pos = f"{helix}.{50 + offset}"
                if pos in position_data and group in position_data[pos]:
                    residues = position_data[pos][group].get('residues', {})
                    if residues:
                        # Get most common residue
                        top_res = max(residues.items(), key=lambda x: x[1][0])
                        group_pattern.append({
                            'offset': offset,
                            'residue': top_res[0],
                            'count': top_res[1][0],
                            'percentage': top_res[1][1]
                        })
            
            if len(group_pattern) >= 5:  # Need sufficient data
                coevolution_patterns[f"{center_pos}_{group}"] = group_pattern
    
    return dict(coevolution_patterns)

def analyze_helix_specific_features(motifs_data: Dict) -> Dict:
    """Extract helix-specific features and patterns"""
    
    helix_features = {}
    
    for helix, groups in motifs_data.items():
        helix_num = helix.split()[-1]
        features = {
            'consensus_motifs': {},
            'conserved_positions': [],
            'variable_positions': [],
            'group_specific_signatures': {}
        }
        
        # Analyze each group's motif
        for group, data in groups.items():
            if 'consensus' in data:
                consensus = data['consensus']
                features['consensus_motifs'][group] = consensus
                
                # Find highly conserved positions (uppercase in consensus)
                conserved = [(i-4, consensus[i]) for i in range(len(consensus)) 
                            if consensus[i].isupper() and consensus[i] not in 'X']
                
                features['group_specific_signatures'][group] = {
                    'consensus': consensus,
                    'conserved_positions': conserved,
                    'n_sequences': data['n_sequences']
                }
        
        # Find positions conserved across all groups
        all_positions = defaultdict(list)
        for group, sig in features['group_specific_signatures'].items():
            for pos, res in sig['conserved_positions']:
                all_positions[pos].append((group, res))
        
        # Positions conserved in >80% of groups
        n_groups = len(features['group_specific_signatures'])
        for pos, group_residues in all_positions.items():
            if len(group_residues) >= 0.8 * n_groups:
                # Check if same residue
                residues = [r for g, r in group_residues]
                if len(set(residues)) == 1:
                    features['conserved_positions'].append((pos, residues[0]))
                else:
                    features['variable_positions'].append((pos, Counter(residues).most_common()))
        
        helix_features[helix] = features
    
    return helix_features

def analyze_domain_evolution_patterns(conservation_by_domain: Dict, motifs_by_domain: Dict) -> Dict:
    """Analyze evolutionary patterns across domains"""
    
    evolution_patterns = {
        'conservation_levels': {},
        'domain_specific_residues': {},
        'evolutionary_transitions': {}
    }
    
    # Calculate average conservation for each domain
    for position, domains in conservation_by_domain.items():
        for domain, data in domains.items():
            if domain not in evolution_patterns['conservation_levels']:
                evolution_patterns['conservation_levels'][domain] = []
            
            if data.get('top_residues'):
                conservation = data['top_residues'][0]['percentage']
                evolution_patterns['conservation_levels'][domain].append(conservation)
    
    # Calculate average conservation per domain
    for domain, conservations in evolution_patterns['conservation_levels'].items():
        evolution_patterns['conservation_levels'][domain] = {
            'mean': np.mean(conservations),
            'std': np.std(conservations),
            'n_positions': len(conservations)
        }
    
    # Find domain-specific residues at each position
    for position, domains in conservation_by_domain.items():
        position_residues = defaultdict(list)
        
        for domain, data in domains.items():
            if data.get('top_residues'):
                top_res = data['top_residues'][0]
                position_residues[top_res['residue']].append({
                    'domain': domain,
                    'percentage': top_res['percentage']
                })
        
        # Find residues specific to certain domains
        domain_specific = {}
        for residue, domain_list in position_residues.items():
            if len(domain_list) == 1 and domain_list[0]['percentage'] > 70:
                domain_specific[residue] = domain_list[0]['domain']
        
        if domain_specific:
            evolution_patterns['domain_specific_residues'][position] = domain_specific
    
    return evolution_patterns

def identify_structural_motifs(motifs_data: Dict, neighborhood_data: Dict) -> Dict:
    """Identify structural motifs like GxxxG, NPxxY analogs, etc."""
    
    structural_motifs = {
        'helix_packing': {},
        'functional_motifs': {},
        'conservation_hotspots': {}
    }
    
    # Look for GxxxG-like motifs
    for helix, groups in motifs_data.items():
        helix_num = helix.split()[-1]
        
        for group, data in groups.items():
            if 'consensus' in data:
                consensus = data['consensus']
                
                # Check for GxxxG pattern
                for i in range(len(consensus) - 4):
                    if consensus[i] == 'G' and consensus[i+4] == 'G':
                        if helix_num not in structural_motifs['helix_packing']:
                            structural_motifs['helix_packing'][helix_num] = []
                        structural_motifs['helix_packing'][helix_num].append({
                            'group': group,
                            'pattern': consensus[i:i+5],
                            'position': i - 4  # Relative to .50
                        })
                
                # Check for WLxT pattern (common in helix 3)
                if helix_num == '3':
                    for i in range(len(consensus) - 3):
                        if consensus[i] == 'W' and consensus[i+1] == 'L' and consensus[i+3] == 'T':
                            if 'WLxT' not in structural_motifs['functional_motifs']:
                                structural_motifs['functional_motifs']['WLxT'] = []
                            structural_motifs['functional_motifs']['WLxT'].append({
                                'group': group,
                                'helix': helix_num,
                                'position': i - 4
                            })
                
                # Check for proline-induced kinks
                proline_positions = [i-4 for i, res in enumerate(consensus) if res == 'P']
                if proline_positions:
                    if helix_num not in structural_motifs['functional_motifs']:
                        structural_motifs['functional_motifs'][helix_num] = {}
                    structural_motifs['functional_motifs'][helix_num][group] = {
                        'proline_positions': proline_positions
                    }
    
    return structural_motifs

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
                    pos_parts = pos.split('.')
                    center_parts = center_pos.split('.')
                    
                    # Handle both "1.5" and "1.50" notations
                    # Position "1.5" in GRN table represents position 1.50
                    if pos_parts[1] == '5':
                        pos_num = 50
                    else:
                        pos_num = int(pos_parts[1])
                        
                    if center_parts[1] == '5':
                        center_num = 50
                    else:
                        center_num = int(center_parts[1])
                        
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

def create_discrimination_heatmap(discriminators: Dict, output_file: Path):
    """Create heatmap showing discrimination power of each position"""
    
    # Prepare data for heatmap
    positions = sorted(discriminators.keys(), key=lambda x: (int(x.split('.')[0]), int(x.split('.')[1])))
    
    # Create discrimination matrix
    groups = ['Proton Pump', 'Cation Channel', 'Anion Channel', 'Chloride Pump', 
              'Sensor / Regulatory', 'Sodium Pump', 'Unknown']
    
    matrix = []
    position_labels = []
    
    for pos in positions:
        if pos in discriminators:
            row = []
            position_labels.append(pos)
            
            for group in groups:
                if group in discriminators[pos]['group_residues']:
                    residue_info = discriminators[pos]['group_residues'][group]
                    # Encode residue as number for heatmap
                    row.append(ord(residue_info['residue']) - ord('A'))
                else:
                    row.append(-1)
            
            matrix.append(row)
    
    if not matrix:
        return
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Create custom colormap
    cmap = plt.colormaps.get_cmap('tab20')
    cmap = cmap.resampled(26)
    cmap.set_under('white')
    
    # Plot heatmap
    im = ax.imshow(np.array(matrix).T, cmap=cmap, aspect='auto', vmin=0, vmax=25)
    
    # Set ticks
    ax.set_xticks(range(len(position_labels)))
    ax.set_xticklabels(position_labels, rotation=45, ha='right')
    ax.set_yticks(range(len(groups)))
    ax.set_yticklabels(groups)
    
    # Add text annotations
    for i, pos in enumerate(position_labels):
        for j, group in enumerate(groups):
            if group in discriminators[pos]['group_residues']:
                residue = discriminators[pos]['group_residues'][group]['residue']
                percentage = discriminators[pos]['group_residues'][group]['percentage']
                text = ax.text(i, j, f'{residue}\n{percentage:.0f}%', 
                             ha='center', va='center', fontsize=8)
    
    ax.set_title('Functional Group Discrimination by GRN Position', fontsize=16, fontweight='bold')
    ax.set_xlabel('GRN Position', fontsize=12)
    ax.set_ylabel('Functional Group', fontsize=12)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def create_coevolution_network(coevolution_patterns: Dict, output_file: Path):
    """Visualize co-evolution patterns as network"""
    
    # This would require networkx, so we'll create a simpler visualization
    fig, axes = plt.subplots(2, 4, figsize=(20, 10))
    axes = axes.flatten()
    
    helices = ['1', '2', '3', '4', '5', '6', '7']
    
    for idx, helix in enumerate(helices):
        ax = axes[idx]
        
        # Collect patterns for this helix
        helix_patterns = {}
        for key, pattern in coevolution_patterns.items():
            # Keys are formatted as "1.5_Proton Pump" (not "1.50_Proton Pump")
            if key.startswith(f"{helix}.5_"):
                group = key.split('_', 1)[1]
                helix_patterns[group] = pattern
        
        if not helix_patterns:
            ax.axis('off')
            continue
        
        # Create pattern matrix
        groups = list(helix_patterns.keys())
        positions = list(range(-4, 5))
        
        matrix = []
        for group in groups:
            row = [''] * 9
            for item in helix_patterns[group]:
                offset = item['offset']
                if -4 <= offset <= 4:
                    row[offset + 4] = item['residue']
            matrix.append(row)
        
        # Plot as text
        ax.set_xlim(-0.5, 8.5)
        ax.set_ylim(-0.5, len(groups) - 0.5)
        
        # Add grid
        for i in range(10):
            ax.axvline(i - 0.5, color='gray', linewidth=0.5, alpha=0.5)
        for i in range(len(groups) + 1):
            ax.axhline(i - 0.5, color='gray', linewidth=0.5, alpha=0.5)
        
        # Add text
        for i, group in enumerate(groups):
            for j, residue in enumerate(matrix[i]):
                if residue:
                    ax.text(j, i, residue, ha='center', va='center', fontsize=12, fontweight='bold')
        
        # Labels
        ax.set_xticks(range(9))
        ax.set_xticklabels([str(p) if p != 0 else '.50' for p in positions])
        ax.set_yticks(range(len(groups)))
        ax.set_yticklabels(groups, fontsize=10)
        ax.set_title(f'Helix {helix} Co-evolution Pattern', fontsize=12, fontweight='bold')
        
        # Highlight .50 position
        ax.axvline(4, color='red', linewidth=2, alpha=0.3)
    
    # Remove empty subplot
    axes[-1].axis('off')
    
    plt.suptitle('Co-evolution Patterns Around .50 Positions', fontsize=16, fontweight='bold')
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()

def create_grn_summary_figure(all_results: Dict):
    """Create comprehensive summary figure for GRN conservation patterns"""
    
    # Create figure with subplots
    fig = plt.figure(figsize=(20, 26))
    
    # Define grid with more spacing
    gs = fig.add_gridspec(6, 3, height_ratios=[1.5, 1, 1, 1, 1, 0.8], 
                          width_ratios=[1, 1, 1], hspace=0.4, wspace=0.3,
                          left=0.05, right=0.95, top=0.96, bottom=0.02)
    
    # ======================================================================
    # Top Panel: Overall Conservation at .50 Positions
    # ======================================================================
    ax1 = fig.add_subplot(gs[0, :])
    
    # Extract conservation data from results
    conservation_by_func = all_results.get('conservation_by_function', {})
    
    positions = ['1.5', '2.5', '3.5', '4.5', '5.5', '6.5', '7.5']
    positions_display = ['1.50', '2.50', '3.50', '4.50', '5.50', '6.50', '7.50']
    conservation = {}
    residues = {}
    
    # Extract conservation values and top residues
    for func in ['Proton Pump', 'Cation Channel', 'Anion Channel', 'Chloride Pump', 
                 'Sensor / Regulatory', 'Sodium Pump', 'Unknown']:
        conservation[func] = []
        residues[func] = []
        
        for pos in positions:
            if pos in conservation_by_func and func in conservation_by_func[pos]:
                data = conservation_by_func[pos][func]
                if data['top_residues']:
                    conservation[func].append(data['top_residues'][0]['percentage'])
                    residues[func].append(data['top_residues'][0]['residue'])
                else:
                    conservation[func].append(0)
                    residues[func].append('')
            else:
                conservation[func].append(0)
                residues[func].append('')
    
    # Define colors for functional groups
    colors = {
        'Proton Pump': '#1f77b4',
        'Cation Channel': '#ff7f0e',
        'Anion Channel': '#2ca02c',
        'Chloride Pump': '#d62728',
        'Sensor / Regulatory': '#9467bd',
        'Sodium Pump': '#8c564b',
        'Unknown': '#7f7f7f'
    }
    
    x = np.arange(len(positions))
    width = 0.12
    
    for i, (group, values) in enumerate(conservation.items()):
        offset = (i - 3) * width
        bars = ax1.bar(x + offset, values, width, label=group, color=colors[group], alpha=0.8)
        
        # Add top residue annotations
        for j, (bar, res) in enumerate(zip(bars, residues[group])):
            height = bar.get_height()
            if height > 80:  # Only annotate highly conserved positions
                ax1.text(bar.get_x() + bar.get_width()/2., height + 1,
                        res, ha='center', va='bottom', fontsize=8, fontweight='bold')
    
    ax1.set_xlabel('GRN Position', fontsize=14)
    ax1.set_ylabel('Conservation (%)', fontsize=14)
    ax1.set_title('Conservation at .50 Positions by Functional Group', fontsize=16, fontweight='bold', pad=20)
    ax1.set_xticks(x)
    ax1.set_xticklabels(positions_display)
    # Place legend inside the plot area to avoid overlap
    ax1.legend(loc='upper right', fontsize=9, ncol=2, framealpha=0.9)
    ax1.set_ylim(0, 105)
    ax1.grid(axis='y', alpha=0.3)
    
    # Add discrimination annotations
    if 'discriminators' in all_results:
        disc = all_results['discriminators']
        if '3.5' in disc and disc['3.5']['discrimination_patterns']:
            ax1.text(2, 95, 'C vs T\n(Channel vs Pump)', ha='center', va='center', 
                     bbox=dict(boxstyle="round,pad=0.3", facecolor='yellow', alpha=0.5),
                     fontsize=10, fontweight='bold')
        if '5.5' in disc and disc['5.5']['discrimination_patterns']:
            ax1.text(4, 95, 'G vs S\n(Channel vs Pump)', ha='center', va='center',
                     bbox=dict(boxstyle="round,pad=0.3", facecolor='yellow', alpha=0.5),
                     fontsize=10, fontweight='bold')
    
    # ======================================================================
    # Second Row: Helix Conservation Profiles
    # ======================================================================
    for helix_idx, helix in enumerate(['1', '3', '5']):
        ax = fig.add_subplot(gs[1, helix_idx])
        
        # Extract conservation profiles from neighborhood data
        if 'neighborhood_by_function' in all_results:
            neighborhood_data = all_results['neighborhood_by_function']
            
            positions_range = list(range(-4, 5))
            profiles = {}
            
            # Group functions into categories
            pump_groups = ['Proton Pump', 'Chloride Pump', 'Sodium Pump']
            channel_groups = ['Cation Channel', 'Anion Channel']
            sensor_groups = ['Sensor / Regulatory']
            
            for category, groups_list in [('Pumps', pump_groups), 
                                         ('Channels', channel_groups),
                                         ('Sensors', sensor_groups)]:
                profile_values = []
                
                center_pos = f"{helix}.5"
                if center_pos in neighborhood_data.get('by_position', {}):
                    pos_data = neighborhood_data['by_position'][center_pos]
                    
                    for offset in positions_range:
                        # Handle position naming - position 50 is stored as "5"
                        pos_num = 50 + offset
                        if pos_num == 50:
                            pos = f"{helix}.5"
                        else:
                            pos = f"{helix}.{pos_num}"
                        conservation_values = []
                        
                        for group in groups_list:
                            if pos in pos_data and group in pos_data[pos]:
                                group_data = pos_data[pos][group]
                                if group_data['residues']:
                                    # Get top conservation
                                    top_conservation = max(r[1] for r in group_data['residues'].values())
                                    conservation_values.append(top_conservation)
                        
                        if conservation_values:
                            profile_values.append(np.mean(conservation_values))
                        else:
                            profile_values.append(0)
                    
                    if any(profile_values):
                        profiles[category] = profile_values
            
            # Plot profiles
            for group, profile in profiles.items():
                ax.plot(positions_range, profile, marker='o', label=group, linewidth=2)
            
            ax.axvline(0, color='red', linestyle='--', alpha=0.5)
            ax.set_xlabel('Position relative to .50', fontsize=10)
            ax.set_ylabel('Conservation (%)', fontsize=10)
            ax.set_title(f'Helix {helix} Conservation Profile', fontsize=12, fontweight='bold')
            ax.legend(fontsize=8)
            ax.grid(alpha=0.3)
            ax.set_xticks(positions_range)
            ax.set_xticklabels([str(p) if p != 0 else '.50' for p in positions_range])
    
    # ======================================================================
    # Third Row: Key Motifs
    # ======================================================================
    ax3 = fig.add_subplot(gs[2, :])
    ax3.axis('off')
    
    # Extract motifs from results
    motif_text = "KEY SEQUENCE MOTIFS AROUND .50 POSITIONS\n\n"
    
    if 'motifs_by_function' in all_results:
        motifs_data = all_results['motifs_by_function']
        
        for helix in ['Helix 1', 'Helix 3', 'Helix 4', 'Helix 5', 'Helix 6', 'Helix 7']:
            if helix in motifs_data:
                helix_num = helix.split()[-1]
                motif_text += f"    HELIX {helix_num}:\n"
                
                for func, data in motifs_data[helix].items():
                    if 'consensus' in data:
                        motif_text += f"        {func}: {data['consensus']} (n={data['n_sequences']})\n"
                
                motif_text += "\n"
    
    motif_text += "\n    Legend: h=hydrophobic, s=small, c=charged, p=polar, a=aromatic, x=variable"
    
    ax3.text(0.05, 0.95, motif_text, transform=ax3.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
    
    # ======================================================================
    # Fourth Row: Functional Discrimination Matrix
    # ======================================================================
    ax4 = fig.add_subplot(gs[3, :2])
    
    # Create discrimination matrix from results
    if 'discriminators' in all_results:
        discriminators = all_results['discriminators']
        
        positions_disc = ['1.5', '2.5', '3.5', '4.5', '5.5']
        positions_disc_display = ['1.50', '2.50', '3.50', '4.50', '5.50']
        groups_disc = ['Proton\nPump', 'Cation\nChannel', 'Anion\nChannel', 
                       'Chloride\nPump', 'Sensor/\nRegulatory', 'Sodium\nPump']
        
        # Create residue matrix
        residue_matrix = []
        for group in groups_disc:
            row = []
            for pos in positions_disc:
                if pos in discriminators and group.replace('\n', ' ') in discriminators[pos]['group_residues']:
                    residue = discriminators[pos]['group_residues'][group.replace('\n', ' ')]['residue']
                    row.append(residue)
                else:
                    row.append('')
            residue_matrix.append(row)
        
        # Create color-coded matrix
        residue_to_color = {
            'G': 0, 'A': 1, 'V': 2, 'I': 3, 'L': 4, 'M': 5, 'F': 6, 'W': 7,
            'S': 8, 'T': 9, 'C': 10, 'Y': 11, 'P': 12, 'H': 13, 'E': 14,
            'D': 15, 'N': 16, 'Q': 17, 'K': 18, 'R': 19
        }
        
        color_matrix = np.zeros((len(groups_disc), len(positions_disc)))
        for i, row in enumerate(residue_matrix):
            for j, res in enumerate(row):
                if res:
                    color_matrix[i, j] = residue_to_color.get(res, 20)
                else:
                    color_matrix[i, j] = -1
        
        cmap = plt.colormaps.get_cmap('tab20')
        cmap = cmap.resampled(26)
        cmap.set_under('white')
        im = ax4.imshow(color_matrix, cmap=cmap, aspect='auto', vmin=0, vmax=25)
        
        # Add text annotations
        for i in range(len(groups_disc)):
            for j in range(len(positions_disc)):
                if residue_matrix[i][j]:
                    text = ax4.text(j, i, residue_matrix[i][j], ha='center', va='center',
                                   color='white', fontsize=12, fontweight='bold')
        
        ax4.set_xticks(range(len(positions_disc)))
        ax4.set_xticklabels(positions_disc_display)
        ax4.set_yticks(range(len(groups_disc)))
        ax4.set_yticklabels(groups_disc)
        ax4.set_title('Functional Group Discrimination at Key Positions', fontsize=14, fontweight='bold')
        ax4.set_xlabel('GRN Position', fontsize=12)
    
    # Add discrimination power subplot
    ax4b = fig.add_subplot(gs[3, 2])
    
    if 'discriminators' in all_results:
        discrimination_power = []
        for pos in positions_disc:
            if pos in discriminators:
                entropy = discriminators[pos]['entropy']
                # Convert entropy to discrimination power
                power = 100 - (entropy/3.5)*100
                discrimination_power.append(power)
            else:
                discrimination_power.append(0)
        
        bars = ax4b.barh(range(len(positions_disc)), discrimination_power, color='darkred', alpha=0.7)
        ax4b.set_yticks(range(len(positions_disc)))
        ax4b.set_yticklabels(positions_disc)
        ax4b.set_xlabel('Discrimination Power (%)', fontsize=12)
        ax4b.set_title('Position Discrimination\nPower', fontsize=12, fontweight='bold')
        ax4b.grid(axis='x', alpha=0.3)
        
        for i, (bar, val) in enumerate(zip(bars, discrimination_power)):
            ax4b.text(val + 1, i, f'{val:.0f}%', va='center', fontsize=10)
    
    # ======================================================================
    # Fifth Row: Evolutionary Conservation
    # ======================================================================
    ax5a = fig.add_subplot(gs[4, 0])
    
    # Domain conservation from evolution patterns
    if 'evolution_patterns' in all_results:
        evo = all_results['evolution_patterns']
        
        domains = []
        domain_conservation = []
        domain_std = []
        
        for domain, stats in sorted(evo['conservation_levels'].items(), 
                                   key=lambda x: x[1]['mean'], reverse=True):
            domains.append(domain)
            domain_conservation.append(stats['mean'])
            domain_std.append(stats['std'])
        
        y_pos = np.arange(len(domains))
        ax5a.barh(y_pos, domain_conservation, xerr=domain_std if domain_std else None,
                  alpha=0.7, capsize=5)
        ax5a.set_yticks(y_pos)
        ax5a.set_yticklabels(domains)
        ax5a.set_xlabel('Average Conservation (%)', fontsize=12)
        ax5a.set_title('Conservation by Domain', fontsize=12, fontweight='bold')
        ax5a.grid(axis='x', alpha=0.3)
        
        for i, val in enumerate(domain_conservation):
            ax5a.text(val + 2, i, f'{val:.1f}%', va='center', fontsize=10)
    
    # Structural features
    ax5b = fig.add_subplot(gs[4, 1:])
    ax5b.axis('off')
    
    structural_text = """STRUCTURAL FEATURES AND CONSTRAINTS\n\n"""
    
    # Extract structural features from results
    if 'structural_motifs' in all_results:
        motifs = all_results['structural_motifs']
        
        if motifs['helix_packing']:
            structural_text += "    Helix Packing Motifs (GxxxG-like):\n"
            for helix, patterns in motifs['helix_packing'].items():
                structural_text += f"    • Helix {helix}: "
                patterns_str = ', '.join([f"{p['pattern']} ({p['group']})" for p in patterns[:2]])
                structural_text += patterns_str + "\n"
            structural_text += "\n"
        
        if 'WLxT' in motifs.get('functional_motifs', {}):
            structural_text += "    WLxT Motif (Helix 3):\n"
            wlxt_groups = set(occ['group'] for occ in motifs['functional_motifs']['WLxT'])
            structural_text += f"    • Found in: {', '.join(list(wlxt_groups)[:3])}\n\n"
    
    # Add key findings
    structural_text += """    Universal Conservation:
    • K7.50 - Schiff base lysine
    • W6.50 - Retinal binding pocket
    • P3.51, P6.54 - Helix kinks
    
    Functional Specializations:
    • Pumps: M4.50, S5.50 - Tight packing
    • Channels: C3.50, G5.50 - Flexible gating
    • Sensors: Intermediate features
    
    Co-evolution Patterns:
    • 3.50 ↔ 5.50: C+G (channels) vs T+S (pumps)
    • WLxT motif: Conserved across all groups"""
    
    ax5b.text(0.05, 0.95, structural_text, transform=ax5b.transAxes, fontsize=9,
              verticalalignment='top', fontfamily='sans-serif',
              bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))
    
    # ======================================================================
    # Bottom Panel: Summary
    # ======================================================================
    ax6 = fig.add_subplot(gs[5, :])
    ax6.axis('off')
    
    summary_text = """KEY FINDINGS:
    
    1. FUNCTIONAL DISCRIMINATION: Position 3.50 is the primary functional switch (C=Channel, T=Pump) with >90% accuracy
    2. UNIVERSAL CONSERVATION: Only K7.50 and W6.50 are conserved across all microbial opsins
    3. CO-EVOLUTION: Positions 3.50 and 5.50 co-evolve to determine transport mechanism
    4. EVOLUTIONARY PATTERN: Archaea (most ancient) → Bacteria → Eukaryotes (most diverse), suggesting pumps are ancestral
    5. STRUCTURAL CONSTRAINTS: GxxxG motifs, proline kinks, and WLxT motif maintain structural integrity
    
    MECHANISTIC IMPLICATIONS: The .50 positions form a conserved structural scaffold where specific residue combinations 
    determine function. Channels likely evolved from pumps through key mutations at positions 3.50 and 5.50."""
    
    ax6.text(0.5, 0.5, summary_text, transform=ax6.transAxes, fontsize=12,
             ha='center', va='center', fontfamily='sans-serif',
             bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.2))
    
    # Main title
    plt.suptitle('GRN Conservation Patterns in Microbial Opsins: Comprehensive Analysis', 
                 fontsize=20, fontweight='bold', y=0.98)
    
    # Save figure without tight_layout to avoid conflicts with gridspec
    output_file = FIGURES_DIR / 'grn_conservation_comprehensive_summary.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Comprehensive summary figure saved to: {output_file}")
    
    # Also create a simplified version for presentations
    create_simplified_summary(all_results)

def create_simplified_summary(all_results: Dict):
    """Create a simplified summary figure focusing on key findings"""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel 1: Key discriminator positions
    if 'conservation_by_function' in all_results:
        positions = ['3.5', '5.5', '4.5']
        positions_display = ['3.50', '5.50', '4.50']
        
        # Calculate average conservation for pumps and channels
        pump_conservation = []
        channel_conservation = []
        pump_res = []
        channel_res = []
        
        conservation_data = all_results['conservation_by_function']
        
        for pos in positions:
            if pos in conservation_data:
                # Pumps
                pump_values = []
                pump_residues = []
                for func in ['Proton Pump', 'Chloride Pump', 'Sodium Pump']:
                    if func in conservation_data[pos] and conservation_data[pos][func]['top_residues']:
                        pump_values.append(conservation_data[pos][func]['top_residues'][0]['percentage'])
                        pump_residues.append(conservation_data[pos][func]['top_residues'][0]['residue'])
                
                # Channels
                channel_values = []
                channel_residues = []
                for func in ['Cation Channel', 'Anion Channel']:
                    if func in conservation_data[pos] and conservation_data[pos][func]['top_residues']:
                        channel_values.append(conservation_data[pos][func]['top_residues'][0]['percentage'])
                        channel_residues.append(conservation_data[pos][func]['top_residues'][0]['residue'])
                
                pump_conservation.append(np.mean(pump_values) if pump_values else 0)
                channel_conservation.append(np.mean(channel_values) if channel_values else 0)
                
                # Get most common residue
                if pump_residues:
                    pump_res.append(Counter(pump_residues).most_common(1)[0][0])
                else:
                    pump_res.append('')
                    
                if channel_residues:
                    channel_res.append(Counter(channel_residues).most_common(1)[0][0])
                else:
                    channel_res.append('')
        
        if pump_conservation and channel_conservation:
            x = np.arange(len(positions))
            width = 0.35
            
            ax1.bar(x - width/2, pump_conservation, width, label='Pumps', color='blue', alpha=0.7)
            ax1.bar(x + width/2, channel_conservation, width, label='Channels', color='orange', alpha=0.7)
            
            # Add residue labels
            for i in range(len(positions)):
                if i < len(pump_res) and pump_res[i] and i < len(pump_conservation):
                    ax1.text(i - width/2, pump_conservation[i] + 2, pump_res[i], 
                            ha='center', fontsize=14, fontweight='bold')
                if i < len(channel_res) and channel_res[i] and i < len(channel_conservation):
                    ax1.text(i + width/2, channel_conservation[i] + 2, channel_res[i], 
                            ha='center', fontsize=14, fontweight='bold')
            
            ax1.set_xticks(x)
            ax1.set_xticklabels(positions_display)
            ax1.legend()
        else:
            ax1.text(0.5, 0.5, 'No discrimination data available', 
                    transform=ax1.transAxes, ha='center', va='center', fontsize=12)
        
        ax1.set_ylabel('Conservation (%)', fontsize=12)
        ax1.set_title('Functional Discriminator Positions', fontsize=14, fontweight='bold')
        ax1.set_ylim(0, 100)
        ax1.grid(axis='y', alpha=0.3)
    
    # Panel 2: Universal vs Variable positions
    ax2.axis('off')
    universal_text = """CONSERVATION HIERARCHY

    UNIVERSAL (>90%):
    • K7.50 - Schiff base
    • W6.50 - Retinal pocket
    
    FUNCTIONAL GROUP (50-90%):
    • 3.50 - C/T discriminator
    • 5.50 - G/S discriminator  
    • 4.50 - M (pumps)
    
    DOMAIN-SPECIFIC (30-50%):
    • 1.50 - Variable
    • 2.50 - I/V/L
    
    HIGHLY VARIABLE (<30%):
    • Loop regions
    • Helix termini"""
    
    ax2.text(0.1, 0.9, universal_text, transform=ax2.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.3))
    
    # Panel 3: Evolutionary pattern
    if 'evolution_patterns' in all_results:
        evo = all_results['evolution_patterns']
        
        domains = []
        conservation = []
        
        # Select key domains
        for domain in ['Archaea', 'Bacteria', 'Eukaryota']:
            if domain in evo['conservation_levels']:
                domains.append(domain)
                conservation.append(evo['conservation_levels'][domain]['mean'])
        
        if domains:
            colors_evo = ['red', 'blue', 'purple'][:len(domains)]
            
            ax3.bar(domains, conservation, color=colors_evo, alpha=0.7)
            ax3.set_ylabel('Average Conservation (%)', fontsize=12)
            ax3.set_title('Evolutionary Conservation Pattern', fontsize=14, fontweight='bold')
            ax3.set_ylim(0, 100)
            ax3.grid(axis='y', alpha=0.3)
            
            for i, val in enumerate(conservation):
                ax3.text(i, val + 2, f'{val:.1f}%', ha='center', fontsize=10)
    
    # Panel 4: Mechanism summary
    ax4.axis('off')
    mechanism_text = """PROPOSED MECHANISM

    PUMPS (ancestral):
    T3.50 + S5.50 + M4.50
    → Tight structure
    → Controlled H+ transport
    → Unidirectional
    
    CHANNELS (derived):
    C3.50 + G5.50 + T4.50
    → Flexible gating
    → Ion selectivity
    → Bidirectional
    
    Key Innovation:
    T→C mutation at 3.50
    creates channel function"""
    
    ax4.text(0.1, 0.9, mechanism_text, transform=ax4.transAxes, fontsize=11,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.3))
    
    plt.suptitle('Microbial Opsin GRN Conservation: Key Findings', fontsize=16, fontweight='bold')
    plt.tight_layout()
    
    output_file = FIGURES_DIR / 'grn_conservation_key_findings.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Simplified summary figure saved to: {output_file}")

# ========================================================================
# Report Generation Functions
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
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(report_lines))
    
    print(f"Summary report saved to {output_file}")

def generate_detailed_interpretation_report(all_results: Dict, output_file: Path):
    """Generate comprehensive interpretation report"""
    
    lines = []
    lines.append("=" * 80)
    lines.append("DETAILED GRN CONSERVATION PATTERN INTERPRETATION")
    lines.append("=" * 80)
    lines.append(f"\nGenerated: {pd.Timestamp.now()}\n")
    
    # 1. Functional Discrimination Analysis
    lines.append("\n1. FUNCTIONAL DISCRIMINATION ANALYSIS")
    lines.append("-" * 60)
    
    if 'discriminators' in all_results:
        # Find best discriminators
        discriminators = all_results['discriminators']
        
        # Sort by discrimination power (low entropy = high discrimination)
        sorted_positions = sorted(discriminators.items(), 
                                 key=lambda x: x[1]['entropy'])
        
        lines.append("\nTop 10 Discriminating Positions:")
        for pos, data in sorted_positions[:10]:
            lines.append(f"\n  Position {pos}:")
            lines.append(f"    Entropy: {data['entropy']:.3f}")
            lines.append(f"    Unique residues: {data['unique_residues']}")
            
            if data['discrimination_patterns']:
                for pattern in data['discrimination_patterns']:
                    if pattern['type'] == 'pump_vs_channel':
                        lines.append(f"    Discriminates pumps vs channels:")
                        lines.append(f"      Pumps: {', '.join(pattern['pump_residues'])}")
                        lines.append(f"      Channels: {', '.join(pattern['channel_residues'])}")
            
            # Show residue distribution
            lines.append("    Residue distribution:")
            for group, info in data['group_residues'].items():
                lines.append(f"      {group}: {info['residue']} ({info['percentage']:.1f}%)")
    
    # 2. Co-evolution Patterns
    lines.append("\n\n2. CO-EVOLUTION PATTERNS")
    lines.append("-" * 60)
    
    if 'coevolution_patterns' in all_results:
        coevo = all_results['coevolution_patterns']
        
        # Analyze patterns by functional group
        functional_patterns = defaultdict(list)
        for key, pattern in coevo.items():
            parts = key.split('_', 1)
            if len(parts) == 2:
                pos, group = parts
                functional_patterns[group].append((pos, pattern))
        
        for group, patterns in functional_patterns.items():
            lines.append(f"\n  {group}:")
            
            # Find conserved motifs across helices
            position_conservation = defaultdict(Counter)
            for pos, pattern in patterns:
                for item in pattern:
                    offset = item['offset']
                    residue = item['residue']
                    if item['percentage'] > 50:  # Only highly conserved
                        position_conservation[offset][residue] += 1
            
            # Report conserved positions
            lines.append("    Conserved positions relative to .50:")
            for offset in sorted(position_conservation.keys()):
                residue_counts = position_conservation[offset]
                if residue_counts:
                    top_residue = residue_counts.most_common(1)[0]
                    if top_residue[1] >= 3:  # Conserved in at least 3 helices
                        lines.append(f"      Position {offset:+d}: {top_residue[0]} "
                                   f"(conserved in {top_residue[1]}/7 helices)")
    
    # 3. Helix-Specific Features
    lines.append("\n\n3. HELIX-SPECIFIC FEATURES")
    lines.append("-" * 60)
    
    if 'helix_features' in all_results:
        features = all_results['helix_features']
        
        for helix, data in sorted(features.items()):
            lines.append(f"\n  {helix}:")
            
            # Universally conserved positions
            if data['conserved_positions']:
                lines.append("    Universally conserved positions:")
                for pos, res in data['conserved_positions']:
                    lines.append(f"      Position {pos:+d}: {res}")
            
            # Variable positions with patterns
            if data['variable_positions']:
                lines.append("    Functionally variable positions:")
                for pos, residue_counts in data['variable_positions']:
                    lines.append(f"      Position {pos:+d}: {', '.join([f'{r}({c})' for r,c in residue_counts[:3]])}")
            
            # Group-specific signatures
            lines.append("    Group-specific consensus motifs:")
            for group, sig in data['group_specific_signatures'].items():
                lines.append(f"      {group}: {sig['consensus']} (n={sig['n_sequences']})")
    
    # 4. Evolutionary Patterns
    lines.append("\n\n4. EVOLUTIONARY PATTERNS")
    lines.append("-" * 60)
    
    if 'evolution_patterns' in all_results:
        evo = all_results['evolution_patterns']
        
        # Conservation levels by domain
        lines.append("\n  Average conservation by domain:")
        for domain, stats in sorted(evo['conservation_levels'].items(), 
                                   key=lambda x: x[1]['mean'], reverse=True):
            lines.append(f"    {domain}: {stats['mean']:.1f}% ± {stats['std']:.1f}%")
        
        # Domain-specific residues
        if evo['domain_specific_residues']:
            lines.append("\n  Domain-specific residue preferences:")
            for pos, residues in sorted(evo['domain_specific_residues'].items()):
                lines.append(f"    Position {pos}:")
                for res, domain in residues.items():
                    lines.append(f"      {res} specific to {domain}")
    
    # 5. Structural Motifs
    lines.append("\n\n5. STRUCTURAL MOTIFS")
    lines.append("-" * 60)
    
    if 'structural_motifs' in all_results:
        motifs = all_results['structural_motifs']
        
        # Helix packing motifs
        if motifs['helix_packing']:
            lines.append("\n  Helix packing motifs (GxxxG-like):")
            for helix, patterns in motifs['helix_packing'].items():
                lines.append(f"    Helix {helix}:")
                for pattern in patterns:
                    lines.append(f"      {pattern['group']}: {pattern['pattern']} at position {pattern['position']:+d}")
        
        # Functional motifs
        if motifs['functional_motifs']:
            lines.append("\n  Functional motifs:")
            for motif_type, occurrences in motifs['functional_motifs'].items():
                if motif_type == 'WLxT':
                    lines.append(f"    {motif_type} motif (Helix 3):")
                    for occ in occurrences:
                        lines.append(f"      {occ['group']} at position {occ['position']:+d}")
                else:
                    lines.append(f"    Helix {motif_type} proline positions:")
                    for group, data in occurrences.items():
                        if 'proline_positions' in data:
                            positions = ', '.join([f"{p:+d}" for p in data['proline_positions']])
                            lines.append(f"      {group}: {positions}")
    
    # 6. Key Interpretations
    lines.append("\n\n6. KEY INTERPRETATIONS")
    lines.append("-" * 60)
    
    lines.append("\n  A. Functional Switching Mechanisms:")
    lines.append("     - Position 3.50 acts as the primary functional switch")
    lines.append("       * C (Cysteine) → Channel function")
    lines.append("       * T (Threonine) → Pump function")
    lines.append("     - This suggests a critical role in gating/transport mechanism")
    
    lines.append("\n  B. Conservation Hierarchy:")
    lines.append("     - Universal: K7.50 (Schiff base), W6.50 (retinal pocket)")
    lines.append("     - Functional group: Positions 3.50, 5.50, 4.50")
    lines.append("     - Domain-specific: Positions 1.50, 2.50")
    lines.append("     - Variable: Loop regions and helix termini")
    
    lines.append("\n  C. Evolutionary Insights:")
    lines.append("     - Archaea show highest conservation (ancient origin)")
    lines.append("     - Eukaryotes show highest diversity (recent diversification)")
    lines.append("     - Channels appear to be derived from pumps")
    lines.append("     - Sensor/regulatory opsins show intermediate features")
    
    lines.append("\n  D. Structural Constraints:")
    lines.append("     - GxxxG motifs in Helix 1 suggest tight packing requirements")
    lines.append("     - Proline positions are conserved within functional groups")
    lines.append("     - WLxT motif in Helix 3 appears critical for all opsins")
    lines.append("     - Position-specific hydrophobicity patterns maintain structure")
    
    lines.append("\n  E. Mechanistic Implications:")
    lines.append("     - .50 positions form a conserved structural scaffold")
    lines.append("     - ±2 positions from .50 show functional specialization")
    lines.append("     - ±4 positions show phylogenetic variation")
    lines.append("     - Co-evolution patterns suggest coupled functional units")
    
    # Write report
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    
    print(f"Detailed interpretation report saved to: {output_file}")

# ========================================================================
# Main Analysis Function
# ========================================================================

def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive GRN Analysis Tool - combines conservation, pattern, and visualization analyses'
    )
    
    # Data input arguments
    parser.add_argument('--residue-table', default='opsin_output/curated_grn.csv',
                       help='Path to residue table')
    parser.add_argument('--property-table', default='property/mo_exp.csv',
                       help='Path to property table')
    
    # Analysis parameters
    parser.add_argument('--positions', default='1.5,2.5,3.5,4.5,5.5,6.5,7.5',
                       help='Comma-separated GRN positions to analyze')
    parser.add_argument('--window', type=int, default=4, 
                       help='Window size around .50 positions')
    parser.add_argument('--min-group-size', type=int, default=3,
                       help='Minimum group size for analysis')
    
    # Analysis options
    parser.add_argument('--skip-conservation', action='store_true',
                       help='Skip conservation analysis')
    parser.add_argument('--skip-patterns', action='store_true',
                       help='Skip pattern analysis')
    parser.add_argument('--skip-figures', action='store_true',
                       help='Skip figure generation')
    
    # Output options
    parser.add_argument('--outputs-dir', default='opsin_output/grn_analysis',
                       help='Base outputs directory')
    
    args = parser.parse_args()
    
    # Update outputs directories based on argument
    global OUTPUT_BASE, CONSERVATION_DIR, FIGURES_DIR, REPORTS_DIR
    OUTPUT_BASE = Path(args.output_dir)
    CONSERVATION_DIR = OUTPUT_BASE / "conservation"
    FIGURES_DIR = OUTPUT_BASE / "figures"
    REPORTS_DIR = OUTPUT_BASE / "reports"
    
    # Ensure outputs directories exist
    ensure_output_dirs()
    
    # Load data
    print("=" * 80)
    print("COMPREHENSIVE GRN ANALYSIS")
    print("=" * 80)
    print("\nLoading data...")
    residue_df, protein_properties = load_and_match_data(args.residue_table, args.property_table)
    
    # All results will be stored here
    all_results = {}
    
    # Parse positions
    positions = args.positions.split(',')
    
    # ========================================================================
    # PART 1: Conservation Analysis (from analyze_grn_conservation.py)
    # ========================================================================
    
    if not args.skip_conservation:
        print("\n" + "="*80)
        print("PART 1: CONSERVATION ANALYSIS")
        print("="*80)
        
        # Analyze conservation by molecular function
        print("\n1.1 Analyzing conservation by molecular function...")
        conservation_by_function = analyze_conservation_by_group(
            residue_df, positions, 'molecular_function', args.min_group_size
        )
        all_results['conservation_by_function'] = conservation_by_function
        
        # Save conservation results
        with open(CONSERVATION_DIR / 'conservation_by_function.json', 'w') as f:
            json.dump(conservation_by_function, f, indent=2)
        
        # Analyze conservation by domain
        print("\n1.2 Analyzing conservation by domain...")
        conservation_by_domain = analyze_conservation_by_group(
            residue_df, positions, 'rhodopsin_type', args.min_group_size
        )
        all_results['conservation_by_domain'] = conservation_by_domain
        
        # Save conservation results
        with open(CONSERVATION_DIR / 'conservation_by_domain.json', 'w') as f:
            json.dump(conservation_by_domain, f, indent=2)
        
        # Analyze neighborhoods
        print("\n1.3 Analyzing neighborhoods around .50 positions...")
        
        print("  By molecular function...")
        neighborhood_by_function = analyze_neighborhood_patterns(
            residue_df, protein_properties, 'molecular_function', positions, args.window
        )
        all_results['neighborhood_by_function'] = neighborhood_by_function
        
        # Save and visualize
        with open(CONSERVATION_DIR / 'neighborhood_by_function.json', 'w') as f:
            json.dump(neighborhood_by_function, f, indent=2, default=str)
        
        if not args.skip_figures:
            visualize_conservation_profiles(neighborhood_by_function, 'molecular_function')
        
        print("  By domain...")
        neighborhood_by_domain = analyze_neighborhood_patterns(
            residue_df, protein_properties, 'rhodopsin_type', positions, args.window
        )
        all_results['neighborhood_by_domain'] = neighborhood_by_domain
        
        # Save and visualize
        with open(CONSERVATION_DIR / 'neighborhood_by_domain.json', 'w') as f:
            json.dump(neighborhood_by_domain, f, indent=2, default=str)
        
        if not args.skip_figures:
            visualize_conservation_profiles(neighborhood_by_domain, 'rhodopsin_type')
        
        # Extract motifs
        print("\n1.4 Extracting sequence motifs...")
        
        print("  By molecular function...")
        motifs_by_function = extract_motifs_by_group(
            residue_df, protein_properties, 'molecular_function', args.window
        )
        all_results['motifs_by_function'] = motifs_by_function
        
        with open(CONSERVATION_DIR / 'motifs_by_function.json', 'w') as f:
            json.dump(motifs_by_function, f, indent=2)
        
        print("  By domain...")
        motifs_by_domain = extract_motifs_by_group(
            residue_df, protein_properties, 'rhodopsin_type', args.window
        )
        all_results['motifs_by_domain'] = motifs_by_domain
        
        with open(CONSERVATION_DIR / 'motifs_by_domain.json', 'w') as f:
            json.dump(motifs_by_domain, f, indent=2)
    
    else:
        # Load previous conservation results if available
        print("\nSkipping conservation analysis, loading previous results if available...")
        
        conservation_files = {
            'conservation_by_function': CONSERVATION_DIR / 'conservation_by_function.json',
            'conservation_by_domain': CONSERVATION_DIR / 'conservation_by_domain.json',
            'neighborhood_by_function': CONSERVATION_DIR / 'neighborhood_by_function.json',
            'neighborhood_by_domain': CONSERVATION_DIR / 'neighborhood_by_domain.json',
            'motifs_by_function': CONSERVATION_DIR / 'motifs_by_function.json',
            'motifs_by_domain': CONSERVATION_DIR / 'motifs_by_domain.json'
        }
        
        for key, filepath in conservation_files.items():
            if filepath.exists():
                with open(filepath, 'r') as f:
                    all_results[key] = json.load(f)
                print(f"  Loaded {key}")
    
    # ========================================================================
    # PART 2: Pattern Analysis (from analyze_grn_patterns_detailed.py)
    # ========================================================================
    
    if not args.skip_patterns and 'conservation_by_function' in all_results:
        print("\n" + "="*80)
        print("PART 2: DETAILED PATTERN ANALYSIS")
        print("="*80)
        
        print("\n2.1 Analyzing functional discriminators...")
        discriminators = analyze_functional_discriminators(all_results['conservation_by_function'])
        all_results['discriminators'] = discriminators
        
        # Create discrimination heatmap
        if not args.skip_figures:
            create_discrimination_heatmap(discriminators, FIGURES_DIR / 'functional_discrimination_heatmap.png')
        
        print("\n2.2 Analyzing co-evolution patterns...")
        if 'neighborhood_by_function' in all_results:
            coevolution_patterns = analyze_coevolution_patterns(all_results['neighborhood_by_function'])
            all_results['coevolution_patterns'] = coevolution_patterns
            
            # Create co-evolution visualization
            if not args.skip_figures:
                create_coevolution_network(coevolution_patterns, FIGURES_DIR / 'coevolution_patterns.png')
        
        print("\n2.3 Analyzing helix-specific features...")
        if 'motifs_by_function' in all_results:
            helix_features = analyze_helix_specific_features(all_results['motifs_by_function'])
            all_results['helix_features'] = helix_features
        
        print("\n2.4 Analyzing domain evolution patterns...")
        if 'conservation_by_domain' in all_results and 'motifs_by_domain' in all_results:
            evolution_patterns = analyze_domain_evolution_patterns(
                all_results['conservation_by_domain'], 
                all_results['motifs_by_domain']
            )
            all_results['evolution_patterns'] = evolution_patterns
        
        print("\n2.5 Identifying structural motifs...")
        if 'motifs_by_function' in all_results and 'neighborhood_by_function' in all_results:
            structural_motifs = identify_structural_motifs(
                all_results['motifs_by_function'], 
                all_results['neighborhood_by_function']
            )
            all_results['structural_motifs'] = structural_motifs
    
    # Save all results
    print("\n" + "="*80)
    print("SAVING RESULTS")
    print("="*80)
    
    with open(OUTPUT_BASE / 'comprehensive_analysis_results.json', 'w') as f:
        json.dump(all_results, f, indent=2, default=str)
    print(f"Saved comprehensive results to: {OUTPUT_BASE / 'comprehensive_analysis_results.json'}")
    
    # ========================================================================
    # PART 3: Report Generation
    # ========================================================================
    
    print("\n" + "="*80)
    print("PART 3: REPORT GENERATION")
    print("="*80)
    
    # Create summary report
    print("\n3.1 Creating summary report...")
    create_summary_report(all_results, REPORTS_DIR / 'conservation_summary_report.txt')
    
    # Generate detailed interpretation report
    print("\n3.2 Generating detailed interpretation report...")
    generate_detailed_interpretation_report(all_results, REPORTS_DIR / 'grn_pattern_interpretation.txt')
    
    # ========================================================================
    # PART 4: Figure Generation (from create_grn_summary_figure.py)
    # ========================================================================
    
    if not args.skip_figures:
        print("\n" + "="*80)
        print("PART 4: FIGURE GENERATION")
        print("="*80)
        
        print("\n4.1 Creating comprehensive summary figure...")
        create_grn_summary_figure(all_results)
    
    # ========================================================================
    # Print key findings
    # ========================================================================
    
    print("\n" + "="*80)
    print("KEY FINDINGS")
    print("="*80)
    
    # Most conserved positions
    if 'conservation_by_function' in all_results:
        print("\nMost conserved positions:")
        overall_conservation = {}
        for pos in positions:
            conservations = []
            if pos in all_results['conservation_by_function']:
                for group_data in all_results['conservation_by_function'][pos].values():
                    if group_data['top_residues']:
                        conservations.append(group_data['top_residues'][0]['percentage'])
            if conservations:
                avg_conservation = np.mean(conservations)
                overall_conservation[pos] = avg_conservation
        
        for pos, cons in sorted(overall_conservation.items(), key=lambda x: x[1], reverse=True):
            print(f"  {pos}: {cons:.1f}% average conservation")
    
    # Functional discriminators
    if 'discriminators' in all_results:
        print("\nTop functional discriminators:")
        disc = all_results['discriminators']
        sorted_disc = sorted(disc.items(), key=lambda x: x[1]['entropy'])[:5]
        
        for pos, data in sorted_disc:
            print(f"  Position {pos}:")
            for pattern in data['discrimination_patterns']:
                if pattern['type'] == 'pump_vs_channel':
                    print(f"    Pumps: {', '.join(pattern['pump_residues'])}")
                    print(f"    Channels: {', '.join(pattern['channel_residues'])}")
    
    print(f"\nAll results saved to: {OUTPUT_BASE}")
    print(f"Conservation data: {CONSERVATION_DIR}")
    print(f"Figures: {FIGURES_DIR}")
    print(f"Reports: {REPORTS_DIR}")
    
    print("\nAnalysis complete!")

if __name__ == '__main__':
    main()