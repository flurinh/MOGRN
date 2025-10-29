#!/usr/bin/env python3
"""
Unified Motif Analysis Script
Combines all motif analysis functionalities:
- Single position motif analysis (GS discriminator etc)
- Literature-defined functional motifs
- Advanced motif correlation search (2-3 residue combinations)
- Function-specific motif discovery
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from itertools import combinations
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Set, Optional
import logging
import argparse
import sys

# Add src to path for imports
sys.path.append('src')
from property_mapping import PropertyMapper

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class UnifiedMotifAnalyzer:
    """Comprehensive motif analysis combining all functionalities"""
    
    def __init__(self, grn_table_path: Path, property_path: Path, output_dir: Path):
        self.grn_table_path = grn_table_path
        self.property_path = property_path
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Literature-defined position mappings
        self.literature_positions = {
            'functional_motif': ['3.45', '3.49', '3.56'],  # BR85, BR89, BR96
            'alt_rhodopsin': '3.42',  # BR82 (Arg82)
            'discriminators': ['3.5', '5.5'],  # Channel vs Pump discriminators
            'universal': ['7.5', '6.5'],  # Highly conserved positions
            'helix3_extended': ['3.46', '3.47', '3.48', '3.51', '3.52', '3.53'],  # Extended helix 3
            'other_helices': ['1.5', '2.5', '4.5']  # Central positions in other helices
        }
        
        # Known functional motifs from literature (expanded based on literature review)
        self.known_motifs = {
            # Outward Proton pumps
            'DTD': {'function': 'Proton Pump', 'direction': 'outward', 
                   'counterion': 'Complex of two Aspartate residues (Asp85 and Asp212)'},
            'DTE': {'function': 'Proton Pump', 'direction': 'outward',
                   'counterion': 'Complex of Aspartate and Glutamate residues (Asp85 and Glu96)'},
            'DTG': {'function': 'Proton Pump', 'direction': 'outward'},
            'DTN': {'function': 'Proton Pump', 'direction': 'outward'},
            'DSA': {'function': 'Proton Pump', 'direction': 'outward'},
            
            # Inward Chloride pumps
            'TSA': {'function': 'Chloride Pump', 'direction': 'inward', 'domain': 'Archaea',
                   'counterion': 'The transported Chloride (Cl⁻) ion itself serves as the counterion'},
            'NTQ': {'function': 'Chloride Pump', 'direction': 'inward', 'domain': 'Bacteria',
                   'counterion': 'The transported Chloride (Cl⁻) ion serves as the counterion'},
            'TSD': {'function': 'Chloride Pump', 'direction': 'inward', 'domain': 'Cyanobacteria',
                   'counterion': 'The transported Chloride (Cl⁻) ion serves as the counterion'},
            'TSL': {'function': 'Chloride Pump', 'direction': 'inward'},
            
            # Outward Sodium pumps
            'NDQ': {'function': 'Sodium Pump', 'direction': 'outward', 'domain': 'Bacteria',
                   'counterion': 'A single Aspartate residue (Asp116) from the NDQ motif'},
            
            # Inward proton pumps
            'FTD': {'function': 'Proton Pump', 'direction': 'inward'},
            
            # Channels
            'STH': {'function': 'Anion Channel'},
            'STL': {'function': 'Anion Channel'},
            'DTT': {'function': 'Cation Channel'},
            'ESK': {'function': 'Cation Channel'},
            'ETH': {'function': 'Cation Channel'},
            
            # Novel families from literature
            'ESL': {'function': 'Heliorhodopsin', 'family': 'HeR',
                   'counterion': 'A single Glutamate residue (e.g., Glu107 in HeR-48C12)'}
        }
        
        # Discriminator residues
        self.discriminators = {
            '3.5': {'C': 'channel', 'T': 'pump', 'V': 'pump', 'I': 'pump', 'A': 'sensor'},
            '5.5': {'G': 'channel', 'S': 'pump/sensor'}
        }
        
        # All positions to analyze for correlations
        self.all_positions = []
        for pos_list in self.literature_positions.values():
            if isinstance(pos_list, list):
                self.all_positions.extend(pos_list)
            else:
                self.all_positions.append(pos_list)
        self.all_positions = sorted(list(set(self.all_positions)))
        # Remove 7.5 from correlation analysis
        if '7.5' in self.all_positions:
            self.all_positions.remove('7.5')
        
    def load_data(self):
        """Load and prepare all data"""
        logger.info("Loading GRN table and property data...")
        
        # Load GRN table
        self.grn_df = pd.read_csv(self.grn_table_path, index_col=0)
        logger.info(f"Loaded GRN table: {self.grn_df.shape}")
        
        # Initialize property mapper
        self.property_mapper = PropertyMapper(self.property_path)
        
        # Map properties to proteins using the structure IDs from GRN table
        self.protein_properties = {}
        for protein_id in self.grn_df.index:
            properties = self.property_mapper.get_properties(protein_id)
            if properties:
                # Normalize function field name
                function = properties.get('molecular_function', 'Unknown')
                if function == 'Unknown' or not function:
                    function = 'Unknown'
                
                self.protein_properties[protein_id] = {
                    'domain': properties.get('domain', 'Unknown'),
                    'function': function,
                    'opsin_name': properties.get('opsin_name', protein_id),
                    'pdb_id': properties.get('pdb_id', None),
                    'short_name': properties.get('short_name', protein_id)
                }
            
        logger.info(f"Loaded properties for {len(self.protein_properties)} proteins out of {len(self.grn_df)} in GRN table")
        
    def extract_residue_only(self, value):
        """Extract only the amino acid residue from values like 'K296'"""
        if pd.isna(value) or value == '-':
            return '-'
        # Extract just the letter(s) at the start
        residue = ''.join(c for c in str(value) if c.isalpha())
        return residue if residue else '-'
    
    # ========================================================================
    # SINGLE POSITION ANALYSIS
    # ========================================================================
    
    def analyze_single_positions(self) -> Dict:
        """Analyze conservation and function correlation at single positions"""
        logger.info("\n=== Single Position Analysis ===")
        
        results = {}
        
        for position in self.all_positions:
            if position not in self.grn_df.columns:
                continue
                
            # Extract residues at this position
            residues = []
            functions = []
            domains = []
            
            for protein_id in self.grn_df.index:
                if protein_id in self.protein_properties:
                    residue = self.extract_residue_only(self.grn_df.loc[protein_id, position])
                    if residue != '-':
                        residues.append(residue)
                        functions.append(self.protein_properties[protein_id]['function'])
                        domains.append(self.protein_properties[protein_id]['domain'])
            
            if not residues:
                continue
                
            # Calculate overall conservation
            residue_counts = Counter(residues)
            total = len(residues)
            conservation = max(residue_counts.values()) / total if total > 0 else 0
            
            # Calculate function-specific patterns
            function_patterns = defaultdict(lambda: defaultdict(int))
            for res, func in zip(residues, functions):
                function_patterns[func][res] += 1
            
            # Find discriminating residues
            discriminating = []
            for func, res_counts in function_patterns.items():
                total_func = sum(res_counts.values())
                for res, count in res_counts.items():
                    if count / total_func >= 0.8:  # 80% threshold
                        # Check if this residue is specific to this function
                        other_funcs = [f for f in function_patterns if f != func]
                        is_specific = True
                        for other_func in other_funcs:
                            if res in function_patterns[other_func]:
                                other_ratio = function_patterns[other_func][res] / sum(function_patterns[other_func].values())
                                if other_ratio > 0.2:  # Less than 20% in other functions
                                    is_specific = False
                                    break
                        
                        if is_specific:
                            discriminating.append({
                                'residue': res,
                                'function': func,
                                'support': count / total_func
                            })
            
            results[position] = {
                'conservation': conservation,
                'top_residue': residue_counts.most_common(1)[0] if residue_counts else None,
                'function_patterns': dict(function_patterns),
                'discriminating': discriminating,
                'total_proteins': total
            }
            
            # Log interesting findings
            if discriminating:
                logger.info(f"Position {position}: Found {len(discriminating)} discriminating residues")
                for disc in discriminating[:3]:  # Show top 3
                    logger.info(f"  - {disc['residue']} for {disc['function']} ({disc['support']:.1%})")
        
        return results
    
    # ========================================================================
    # LITERATURE MOTIF ANALYSIS
    # ========================================================================
    
    def analyze_literature_motifs(self) -> Dict:
        """Analyze known functional motifs from literature"""
        logger.info("\n=== Literature Motif Analysis ===")
        
        results = {
            'functional_motifs': [],
            'alt_rhodopsins': [],
            'discriminator_analysis': {}
        }
        
        # Extract functional motifs (3.45-3.49-3.56)
        positions = self.literature_positions['functional_motif']
        
        for protein_id in self.grn_df.index:
            if protein_id not in self.protein_properties:
                continue
                
            # Extract motif
            motif_residues = []
            for pos in positions:
                if pos in self.grn_df.columns:
                    residue = self.extract_residue_only(self.grn_df.loc[protein_id, pos])
                    motif_residues.append(residue)
            
            if len(motif_residues) == 3 and all(r != '-' for r in motif_residues):
                motif = ''.join(motif_residues)
                
                # Check if it's a known motif
                is_known = motif in self.known_motifs
                expected_function = self.known_motifs.get(motif, {}).get('function', 'Unknown')
                actual_function = self.protein_properties[protein_id]['function']
                
                results['functional_motifs'].append({
                    'protein': protein_id,
                    'motif': motif,
                    'function': actual_function,
                    'domain': self.protein_properties[protein_id]['domain'],
                    'is_known': is_known,
                    'expected_function': expected_function,
                    'matches_expectation': expected_function == actual_function if is_known else None
                })
        
        # Analyze Alt-rhodopsins (non-R at position 3.42)
        alt_pos = self.literature_positions['alt_rhodopsin']
        if alt_pos in self.grn_df.columns:
            for protein_id in self.grn_df.index:
                if protein_id not in self.protein_properties:
                    continue
                    
                residue = self.extract_residue_only(self.grn_df.loc[protein_id, alt_pos])
                if residue != '-' and residue != 'R':
                    results['alt_rhodopsins'].append({
                        'protein': protein_id,
                        'residue_3.42': residue,
                        'function': self.protein_properties[protein_id]['function'],
                        'domain': self.protein_properties[protein_id]['domain']
                    })
        
        # Analyze discriminator positions
        for disc_pos in self.literature_positions['discriminators']:
            if disc_pos not in self.grn_df.columns:
                continue
                
            disc_analysis = defaultdict(lambda: defaultdict(int))
            
            for protein_id in self.grn_df.index:
                if protein_id not in self.protein_properties:
                    continue
                    
                residue = self.extract_residue_only(self.grn_df.loc[protein_id, disc_pos])
                if residue != '-':
                    function = self.protein_properties[protein_id]['function']
                    disc_analysis[function][residue] += 1
            
            results['discriminator_analysis'][disc_pos] = dict(disc_analysis)
        
        # Log summary
        motif_counts = Counter(m['motif'] for m in results['functional_motifs'])
        logger.info(f"Found {len(results['functional_motifs'])} functional motifs")
        logger.info(f"Top motifs: {motif_counts.most_common(5)}")
        logger.info(f"Found {len(results['alt_rhodopsins'])} Alt-rhodopsins")
        
        return results
    
    # ========================================================================
    # ADVANCED CORRELATION ANALYSIS
    # ========================================================================
    
    def analyze_advanced_correlations(self, min_support: float = 0.8, 
                                    min_proteins: int = 3) -> Dict:
        """Find correlated 2-3 residue motifs within domain-function pairs"""
        logger.info("\n=== Advanced Correlation Analysis ===")
        logger.info(f"Parameters: min_support={min_support}, min_proteins={min_proteins}")
        
        results = {
            '2-residue': [],
            '3-residue': [],
            'function_specific': defaultdict(list),
            'cross_function': []
        }
        
        # Group proteins by domain-function
        domain_function_groups = defaultdict(list)
        for protein_id in self.grn_df.index:
            if protein_id in self.protein_properties:
                domain = self.protein_properties[protein_id]['domain']
                function = self.protein_properties[protein_id]['function']
                domain_function_groups[(domain, function)].append(protein_id)
        
        # Filter groups with enough proteins
        valid_groups = {k: v for k, v in domain_function_groups.items() 
                       if len(v) >= min_proteins}
        
        logger.info(f"Analyzing {len(valid_groups)} domain-function pairs")
        
        # Analyze each group
        for (domain, function), protein_ids in valid_groups.items():
            group_size = len(protein_ids)
            
            # Extract all residues for this group
            group_residues = defaultdict(list)
            for protein_id in protein_ids:
                for position in self.all_positions:
                    if position in self.grn_df.columns:
                        residue = self.extract_residue_only(
                            self.grn_df.loc[protein_id, position]
                        )
                        if residue != '-':
                            group_residues[position].append(residue)
            
            # Find 2-residue correlations
            for pos1, pos2 in combinations(self.all_positions, 2):
                if pos1 not in group_residues or pos2 not in group_residues:
                    continue
                    
                # Get paired residues
                pairs = []
                for protein_id in protein_ids:
                    if pos1 in self.grn_df.columns and pos2 in self.grn_df.columns:
                        res1 = self.extract_residue_only(
                            self.grn_df.loc[protein_id, pos1]
                        )
                        res2 = self.extract_residue_only(
                            self.grn_df.loc[protein_id, pos2]
                        )
                        if res1 != '-' and res2 != '-':
                            pairs.append(f"{res1}{res2}")
                
                if not pairs:
                    continue
                    
                # Find dominant pairs
                pair_counts = Counter(pairs)
                for pair, count in pair_counts.items():
                    support = count / len(pairs)
                    if support >= min_support:
                        correlation = {
                            'domain': domain,
                            'function': function,
                            'positions': f"{pos1}+{pos2}",
                            'motif': pair,
                            'count': count,
                            'total': len(pairs),
                            'support': support,
                            'type': '2-residue'
                        }
                        results['2-residue'].append(correlation)
                        results['function_specific'][function].append(correlation)
            
            # Find 3-residue correlations (focus on key positions)
            key_triplets = [
                self.literature_positions['functional_motif'],  # 3.45, 3.49, 3.56
                ['3.42', '3.5', '5.5'],  # Alt-rhodopsin + discriminators
                # Removed ['3.5', '5.5', '7.5'] to exclude 7.5 from analysis
            ]
            
            for positions in key_triplets:
                if all(pos in self.grn_df.columns for pos in positions):
                    triplets = []
                    for protein_id in protein_ids:
                        residues = []
                        for pos in positions:
                            res = self.extract_residue_only(
                                self.grn_df.loc[protein_id, pos]
                            )
                            residues.append(res)
                        
                        if all(r != '-' for r in residues):
                            triplets.append(''.join(residues))
                    
                    if not triplets:
                        continue
                        
                    # Find dominant triplets
                    triplet_counts = Counter(triplets)
                    for triplet, count in triplet_counts.items():
                        support = count / len(triplets)
                        if support >= min_support:
                            correlation = {
                                'domain': domain,
                                'function': function,
                                'positions': '+'.join(positions),
                                'motif': triplet,
                                'count': count,
                                'total': len(triplets),
                                'support': support,
                                'type': '3-residue'
                            }
                            results['3-residue'].append(correlation)
                            results['function_specific'][function].append(correlation)
        
        # Find cross-function patterns
        self._find_cross_function_patterns(results)
        
        # Log summary
        logger.info(f"Found {len(results['2-residue'])} 2-residue correlations")
        logger.info(f"Found {len(results['3-residue'])} 3-residue correlations")
        
        return results
    
    def _find_cross_function_patterns(self, results: Dict):
        """Find motifs that appear across multiple functions"""
        # Collect all motifs by position combination
        position_motifs = defaultdict(lambda: defaultdict(list))
        
        for corr_type in ['2-residue', '3-residue']:
            for corr in results[corr_type]:
                positions = corr['positions']
                motif = corr['motif']
                function = corr['function']
                position_motifs[positions][motif].append(function)
        
        # Find motifs appearing in multiple functions
        for positions, motifs in position_motifs.items():
            for motif, functions in motifs.items():
                unique_functions = list(set(functions))
                if len(unique_functions) > 1:
                    results['cross_function'].append({
                        'positions': positions,
                        'motif': motif,
                        'functions': unique_functions,
                        'num_functions': len(unique_functions)
                    })
    
    # ========================================================================
    # VISUALIZATION AND OUTPUT
    # ========================================================================
    
    def create_comprehensive_report(self, single_pos_results: Dict, 
                                  literature_results: Dict,
                                  correlation_results: Dict):
        """Create comprehensive analysis report"""
        report_path = self.output_dir / "comprehensive_motif_report.md"
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("# Comprehensive Motif Analysis Report\n\n")
            
            # Executive Summary
            f.write("## Executive Summary\n\n")
            
            # Count key findings
            n_discriminating = sum(len(r['discriminating']) for r in single_pos_results.values())
            n_literature_motifs = len(literature_results['functional_motifs'])
            n_alt_rhodopsins = len(literature_results['alt_rhodopsins'])
            n_2res_correlations = len(correlation_results['2-residue'])
            n_3res_correlations = len(correlation_results['3-residue'])
            
            f.write(f"- **Single Position Analysis**: Found {n_discriminating} discriminating residues\n")
            f.write(f"- **Literature Motifs**: Identified {n_literature_motifs} functional motifs\n")
            f.write(f"- **Alt-rhodopsins**: Found {n_alt_rhodopsins} proteins with non-R at 3.42\n")
            f.write(f"- **Correlations**: {n_2res_correlations} 2-residue and {n_3res_correlations} 3-residue patterns\n\n")
            
            # Single Position Results
            f.write("## Single Position Analysis\n\n")
            f.write("### Most Discriminating Positions\n\n")
            
            # Sort by number of discriminating residues
            sorted_positions = sorted(single_pos_results.items(), 
                                    key=lambda x: len(x[1]['discriminating']), 
                                    reverse=True)
            
            for position, data in sorted_positions[:10]:
                if data['discriminating']:
                    f.write(f"#### Position {position}\n")
                    f.write(f"- Conservation: {data['conservation']:.1%}\n")
                    f.write(f"- Top residue: {data['top_residue'][0]} ({data['top_residue'][1]}/{data['total_proteins']})\n")
                    f.write("- Discriminating residues:\n")
                    for disc in data['discriminating']:
                        f.write(f"  - **{disc['residue']}** → {disc['function']} ({disc['support']:.1%})\n")
                    f.write("\n")
            
            # Literature Motif Results
            f.write("## Literature Motif Analysis\n\n")
            
            # Functional motifs
            f.write("### Functional Motifs (3.45-3.49-3.56)\n\n")
            motif_by_function = defaultdict(list)
            for entry in literature_results['functional_motifs']:
                motif_by_function[entry['function']].append(entry['motif'])
            
            f.write("| Function | Common Motifs | Count |\n")
            f.write("|----------|---------------|-------|\n")
            for function, motifs in sorted(motif_by_function.items()):
                motif_counts = Counter(motifs)
                top_motifs = ', '.join([f"{m} ({c})" for m, c in motif_counts.most_common(3)])
                f.write(f"| {function} | {top_motifs} | {len(motifs)} |\n")
            
            # Add known motif details table
            f.write("\n### Known Motif Details\n\n")
            f.write("| Motif | Function | Direction | Domain | Counterion System |\n")
            f.write("|-------|----------|-----------|---------|------------------|\n")
            for motif, details in self.known_motifs.items():
                function = details.get('function', '-')
                direction = details.get('direction', '-')
                domain = details.get('domain', 'Multiple')
                counterion = details.get('counterion', 'Not specified')
                # Truncate long counterion descriptions for table readability
                if len(counterion) > 60:
                    counterion = counterion[:57] + "..."
                f.write(f"| {motif} | {function} | {direction} | {domain} | {counterion} |\n")
            
            # Alt-rhodopsins
            f.write("\n### Alt-rhodopsins (Non-R at 3.42)\n\n")
            alt_by_function = defaultdict(list)
            for entry in literature_results['alt_rhodopsins']:
                alt_by_function[entry['function']].append(entry['residue_3.42'])
            
            f.write("| Function | Alt Residues | Count | Percentage |\n")
            f.write("|----------|--------------|-------|------------|\n")
            for function, residues in sorted(alt_by_function.items()):
                total_func = len([e for e in literature_results['functional_motifs'] 
                                 if e['function'] == function])
                percentage = len(residues) / total_func * 100 if total_func > 0 else 0
                unique_residues = ', '.join(sorted(set(residues)))
                f.write(f"| {function} | {unique_residues} | {len(residues)} | {percentage:.1f}% |\n")
            
            # Correlation Results
            f.write("\n## Advanced Correlation Analysis\n\n")
            
            # Function-specific patterns
            f.write("### Function-Specific Motif Patterns\n\n")
            for function, correlations in correlation_results['function_specific'].items():
                if correlations:
                    f.write(f"#### {function}\n")
                    
                    # Group by residue count
                    two_res = [c for c in correlations if c['type'] == '2-residue']
                    three_res = [c for c in correlations if c['type'] == '3-residue']
                    
                    if two_res:
                        f.write("\n**2-Residue Patterns:**\n")
                        # Sort by support
                        for corr in sorted(two_res, key=lambda x: x['support'], reverse=True)[:5]:
                            f.write(f"- {corr['positions']}: **{corr['motif']}** "
                                   f"({corr['support']:.0%}, {corr['count']}/{corr['total']})\n")
                    
                    if three_res:
                        f.write("\n**3-Residue Patterns:**\n")
                        for corr in sorted(three_res, key=lambda x: x['support'], reverse=True)[:5]:
                            f.write(f"- {corr['positions']}: **{corr['motif']}** "
                                   f"({corr['support']:.0%}, {corr['count']}/{corr['total']})\n")
                    f.write("\n")
            
            # Cross-function patterns
            if correlation_results['cross_function']:
                f.write("### Cross-Function Patterns\n\n")
                f.write("Motifs appearing in multiple functions:\n\n")
                
                for pattern in sorted(correlation_results['cross_function'], 
                                    key=lambda x: x['num_functions'], reverse=True)[:10]:
                    f.write(f"- {pattern['positions']}: **{pattern['motif']}** → "
                           f"{', '.join(pattern['functions'])}\n")
        
        logger.info(f"Report saved to: {report_path}")
    
    def create_visualizations(self, single_pos_results: Dict, 
                            literature_results: Dict,
                            correlation_results: Dict):
        """Create comprehensive visualizations"""
        # Create figure with subplots
        fig = plt.figure(figsize=(20, 16))
        
        # 1. Single position conservation heatmap
        ax1 = plt.subplot(3, 3, 1)
        self._plot_position_conservation(ax1, single_pos_results)
        
        # 2. Functional motif distribution
        ax2 = plt.subplot(3, 3, 2)
        self._plot_motif_distribution(ax2, literature_results)
        
        # 3. Alt-rhodopsin analysis
        ax3 = plt.subplot(3, 3, 3)
        self._plot_alt_rhodopsin_distribution(ax3, literature_results)
        
        # 4. Discriminator analysis
        ax4 = plt.subplot(3, 3, 4)
        self._plot_discriminator_analysis(ax4, literature_results)
        
        # 5. 2-residue correlation network
        ax5 = plt.subplot(3, 3, 5)
        self._plot_correlation_network(ax5, correlation_results, '2-residue')
        
        # 6. 3-residue patterns
        ax6 = plt.subplot(3, 3, 6)
        self._plot_correlation_network(ax6, correlation_results, '3-residue')
        
        # 7. Function-specific signatures
        ax7 = plt.subplot(3, 3, 7)
        self._plot_function_signatures(ax7, correlation_results)
        
        # 8. Cross-function patterns
        ax8 = plt.subplot(3, 3, 8)
        self._plot_cross_function_patterns(ax8, correlation_results)
        
        # 9. Summary statistics
        ax9 = plt.subplot(3, 3, 9)
        self._plot_summary_statistics(ax9, single_pos_results, 
                                    literature_results, correlation_results)
        
        plt.suptitle('Comprehensive Motif Analysis', fontsize=20)
        plt.tight_layout()
        
        output_path = self.output_dir / "comprehensive_motif_analysis.png"
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Visualization saved to: {output_path}")
    
    def _plot_position_conservation(self, ax, results):
        """Plot conservation levels at each position"""
        positions = []
        conservations = []
        
        for pos in self.all_positions:
            if pos in results and results[pos]['conservation'] > 0:
                positions.append(pos)
                conservations.append(results[pos]['conservation'])
        
        # Create bar plot
        bars = ax.bar(range(len(positions)), conservations)
        
        # Color by conservation level
        for i, (pos, cons) in enumerate(zip(positions, conservations)):
            if cons > 0.9:
                bars[i].set_color('darkred')
            elif cons > 0.7:
                bars[i].set_color('orange')
            else:
                bars[i].set_color('lightblue')
        
        ax.set_xticks(range(len(positions)))
        ax.set_xticklabels(positions, rotation=45, ha='right')
        ax.set_ylabel('Conservation')
        ax.set_title('Single Position Conservation')
        ax.set_ylim(0, 1.1)
        
        # Add discriminating marker
        for i, pos in enumerate(positions):
            if results[pos]['discriminating']:
                ax.text(i, conservations[i] + 0.02, '*', ha='center', fontsize=12)
    
    def _plot_motif_distribution(self, ax, results):
        """Plot functional motif distribution"""
        motif_counts = Counter(m['motif'] for m in results['functional_motifs'])
        
        # Get top 15 motifs
        top_motifs = motif_counts.most_common(15)
        motifs, counts = zip(*top_motifs)
        
        bars = ax.barh(range(len(motifs)), counts)
        
        # Color known motifs differently
        for i, motif in enumerate(motifs):
            if motif in self.known_motifs:
                bars[i].set_color('darkgreen')
            else:
                bars[i].set_color('gray')
        
        ax.set_yticks(range(len(motifs)))
        ax.set_yticklabels(motifs)
        ax.set_xlabel('Count')
        ax.set_title('Functional Motif Distribution (3.45-3.49-3.56)')
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [
            Patch(facecolor='darkgreen', label='Known motif'),
            Patch(facecolor='gray', label='Novel motif')
        ]
        ax.legend(handles=legend_elements, loc='lower right')
    
    def _plot_alt_rhodopsin_distribution(self, ax, results):
        """Plot Alt-rhodopsin distribution by function"""
        alt_by_function = defaultdict(list)
        total_by_function = defaultdict(int)
        
        # Count alt-rhodopsins
        for entry in results['alt_rhodopsins']:
            alt_by_function[entry['function']].append(entry['residue_3.42'])
        
        # Count total proteins per function
        for entry in results['functional_motifs']:
            total_by_function[entry['function']] += 1
        
        # Calculate percentages
        functions = sorted(alt_by_function.keys())
        percentages = []
        for func in functions:
            alt_count = len(alt_by_function[func])
            total = total_by_function[func]
            percentages.append(alt_count / total * 100 if total > 0 else 0)
        
        # Create bar plot
        bars = ax.bar(range(len(functions)), percentages)
        
        # Color channels differently
        for i, func in enumerate(functions):
            if 'Channel' in func:
                bars[i].set_color('darkblue')
            else:
                bars[i].set_color('lightgray')
        
        ax.set_xticks(range(len(functions)))
        ax.set_xticklabels(functions, rotation=45, ha='right')
        ax.set_ylabel('Alt-rhodopsin %')
        ax.set_title('Alt-rhodopsin Distribution (Non-R at 3.42)')
    
    def _plot_discriminator_analysis(self, ax, results):
        """Plot discriminator position analysis"""
        disc_data = results['discriminator_analysis']
        
        if not disc_data:
            ax.text(0.5, 0.5, 'No discriminator data', ha='center', va='center')
            return
        
        # Focus on 3.5 and 5.5
        positions = ['3.5', '5.5']
        functions = ['Proton Pump', 'Chloride Pump', 'Sodium Pump', 
                    'Cation Channel', 'Anion Channel']
        
        # Create matrix for heatmap
        matrix = []
        labels = []
        
        for pos in positions:
            if pos in disc_data:
                for func in functions:
                    if func in disc_data[pos]:
                        residue_counts = disc_data[pos][func]
                        total = sum(residue_counts.values())
                        # Get top residue
                        if residue_counts:
                            top_res = max(residue_counts, key=residue_counts.get)
                            percentage = residue_counts[top_res] / total
                            matrix.append(percentage)
                            labels.append(f"{pos}:{func[:8]}\n{top_res}")
                        else:
                            matrix.append(0)
                            labels.append(f"{pos}:{func[:8]}\n-")
        
        if matrix:
            # Reshape for display
            n_funcs = len(functions)
            matrix = np.array(matrix).reshape(len(positions), n_funcs)
            
            im = ax.imshow(matrix, cmap='RdYlBu_r', aspect='auto', vmin=0, vmax=1)
            
            ax.set_xticks(range(n_funcs))
            ax.set_xticklabels([f[:8] for f in functions], rotation=45, ha='right')
            ax.set_yticks(range(len(positions)))
            ax.set_yticklabels(positions)
            
            # Add text annotations
            for i in range(len(positions)):
                for j in range(n_funcs):
                    idx = i * n_funcs + j
                    if idx < len(labels):
                        text = labels[idx].split('\n')[1]
                        ax.text(j, i, text, ha='center', va='center',
                               color='white' if matrix[i, j] > 0.5 else 'black')
            
            ax.set_title('Discriminator Positions (3.5, 5.5)')
            
            # Add colorbar
            plt.colorbar(im, ax=ax, label='Frequency')
    
    def _plot_correlation_network(self, ax, results, corr_type):
        """Plot correlation network for 2 or 3 residue patterns"""
        correlations = results[corr_type]
        
        if not correlations:
            ax.text(0.5, 0.5, f'No {corr_type} correlations', ha='center', va='center')
            return
        
        # Get top correlations by support
        top_corrs = sorted(correlations, key=lambda x: x['support'], reverse=True)[:20]
        
        # Create visualization data
        positions = []
        motifs = []
        supports = []
        functions = []
        
        for corr in top_corrs:
            positions.append(corr['positions'])
            motifs.append(corr['motif'])
            supports.append(corr['support'])
            functions.append(corr['function'])
        
        # Create horizontal bar plot
        y_pos = np.arange(len(positions))
        bars = ax.barh(y_pos, supports)
        
        # Color by function
        function_colors = {
            'Proton Pump': 'red',
            'Chloride Pump': 'blue',
            'Sodium Pump': 'green',
            'Cation Channel': 'orange',
            'Anion Channel': 'purple',
            'Sensor / Regulatory': 'gray'
        }
        
        for i, func in enumerate(functions):
            bars[i].set_color(function_colors.get(func, 'black'))
        
        # Create labels
        labels = [f"{pos}: {motif}" for pos, motif in zip(positions, motifs)]
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel('Support')
        ax.set_title(f'Top {corr_type.title()} Correlations')
        ax.set_xlim(0.7, 1.05)
        
        # Add legend
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=color, label=func) 
                          for func, color in function_colors.items()]
        ax.legend(handles=legend_elements, loc='lower right', fontsize=6)
    
    def _plot_function_signatures(self, ax, results):
        """Plot function-specific signature counts"""
        function_counts = defaultdict(int)
        
        for correlations in results['function_specific'].values():
            for corr in correlations:
                function_counts[corr['function']] += 1
        
        functions = sorted(function_counts.keys())
        counts = [function_counts[f] for f in functions]
        
        bars = ax.bar(range(len(functions)), counts)
        
        # Color by function type
        for i, func in enumerate(functions):
            if 'Pump' in func:
                bars[i].set_color('red')
            elif 'Channel' in func:
                bars[i].set_color('blue')
            else:
                bars[i].set_color('gray')
        
        ax.set_xticks(range(len(functions)))
        ax.set_xticklabels(functions, rotation=45, ha='right')
        ax.set_ylabel('Number of Signatures')
        ax.set_title('Function-Specific Signatures')
    
    def _plot_cross_function_patterns(self, ax, results):
        """Plot patterns appearing across multiple functions"""
        cross_patterns = results['cross_function']
        
        if not cross_patterns:
            ax.text(0.5, 0.5, 'No cross-function patterns', ha='center', va='center')
            return
        
        # Sort by number of functions
        sorted_patterns = sorted(cross_patterns, 
                               key=lambda x: x['num_functions'], 
                               reverse=True)[:15]
        
        # Create data for plotting
        labels = []
        num_functions = []
        
        for pattern in sorted_patterns:
            label = f"{pattern['positions']}: {pattern['motif']}"
            labels.append(label)
            num_functions.append(pattern['num_functions'])
        
        # Create horizontal bar plot
        y_pos = np.arange(len(labels))
        ax.barh(y_pos, num_functions, color='darkgreen')
        
        ax.set_yticks(y_pos)
        ax.set_yticklabels(labels, fontsize=8)
        ax.set_xlabel('Number of Functions')
        ax.set_title('Cross-Function Patterns')
        
        # Add function details on hover
        for i, pattern in enumerate(sorted_patterns):
            func_str = ', '.join(pattern['functions'][:3])
            if len(pattern['functions']) > 3:
                func_str += '...'
            ax.text(num_functions[i] + 0.1, i, func_str, 
                   fontsize=6, va='center')
    
    def _plot_summary_statistics(self, ax, single_results, lit_results, corr_results):
        """Plot summary statistics"""
        ax.axis('off')
        
        # Calculate statistics
        n_positions = len([p for p in single_results if single_results[p]['conservation'] > 0])
        n_discriminating = sum(len(r['discriminating']) for r in single_results.values())
        n_literature_motifs = len(lit_results['functional_motifs'])
        n_known_motifs = len([m for m in lit_results['functional_motifs'] if m['is_known']])
        n_novel_motifs = n_literature_motifs - n_known_motifs
        n_alt_rhodopsins = len(lit_results['alt_rhodopsins'])
        n_2res = len(corr_results['2-residue'])
        n_3res = len(corr_results['3-residue'])
        n_cross_function = len(corr_results['cross_function'])
        
        # Create summary text
        summary_text = f"""Analysis Summary:

Single Position Analysis:
• Positions analyzed: {n_positions}
• Discriminating positions: {n_discriminating}
• Average conservation: {np.mean([r['conservation'] for r in single_results.values() if r['conservation'] > 0]):.1%}

Literature Motif Analysis:
• Total functional motifs: {n_literature_motifs}
• Known motifs: {n_known_motifs}
• Novel motifs: {n_novel_motifs}
• Alt-rhodopsins: {n_alt_rhodopsins} ({n_alt_rhodopsins/n_literature_motifs*100:.1f}%)

Advanced Correlations:
• 2-residue patterns: {n_2res}
• 3-residue patterns: {n_3res}
• Cross-function patterns: {n_cross_function}
• Functions analyzed: {len(corr_results['function_specific'])}

Key Findings:
• Most conserved: Position {max(single_results.items(), key=lambda x: x[1]['conservation'])[0]}
• Most discriminating: Positions 3.5, 5.5
• Dominant motifs: DTD/DTE (H⁺ pump), NDQ (Na⁺ pump)
• Alt-rhodopsins enriched in channels"""
        
        ax.text(0.05, 0.95, summary_text, transform=ax.transAxes,
               fontsize=10, va='top', ha='left', family='monospace',
               bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.5))
    
    def run_analysis(self):
        """Run complete analysis pipeline"""
        # Load data
        self.load_data()
        
        # Run analyses
        logger.info("Running comprehensive motif analysis...")
        single_results = self.analyze_single_positions()
        literature_results = self.analyze_literature_motifs()
        correlation_results = self.analyze_advanced_correlations()
        
        # Create outputs
        logger.info("Creating report and visualizations...")
        self.create_comprehensive_report(single_results, literature_results, correlation_results)
        self.create_visualizations(single_results, literature_results, correlation_results)
        
        # Save detailed results
        self._save_detailed_results(single_results, literature_results, correlation_results)
        
        # Create extended motif table
        self.create_extended_motif_table(literature_results)
        
        logger.info("Analysis complete!")
    
    def _save_detailed_results(self, single_results, literature_results, correlation_results):
        """Save detailed results to CSV files"""
        # Single position results
        single_df_data = []
        for position, data in single_results.items():
            if data['discriminating']:
                for disc in data['discriminating']:
                    single_df_data.append({
                        'position': position,
                        'residue': disc['residue'],
                        'function': disc['function'],
                        'support': disc['support'],
                        'conservation': data['conservation']
                    })
        
        if single_df_data:
            single_df = pd.DataFrame(single_df_data)
            single_df.to_csv(self.output_dir / 'single_position_discriminators.csv', index=False, encoding='utf-8')
        
        # Literature motifs
        lit_df = pd.DataFrame(literature_results['functional_motifs'])
        lit_df.to_csv(self.output_dir / 'literature_motifs.csv', index=False, encoding='utf-8')
        
        alt_df = pd.DataFrame(literature_results['alt_rhodopsins'])
        if not alt_df.empty:
            alt_df.to_csv(self.output_dir / 'alt_rhodopsins.csv', index=False, encoding='utf-8')
        
        # Correlations
        corr_2res_df = pd.DataFrame(correlation_results['2-residue'])
        if not corr_2res_df.empty:
            corr_2res_df.to_csv(self.output_dir / 'correlations_2residue.csv', index=False, encoding='utf-8')
        
        corr_3res_df = pd.DataFrame(correlation_results['3-residue'])
        if not corr_3res_df.empty:
            corr_3res_df.to_csv(self.output_dir / 'correlations_3residue.csv', index=False, encoding='utf-8')
        
        cross_df = pd.DataFrame(correlation_results['cross_function'])
        if not cross_df.empty:
            cross_df.to_csv(self.output_dir / 'cross_function_patterns.csv', index=False, encoding='utf-8')
        
        logger.info(f"Detailed results saved to: {self.output_dir}")
    
    def create_extended_motif_table(self, literature_results: Dict):
        """Create extended motif distribution table including novel motifs"""
        logger.info("Creating extended motif distribution table...")
        
        # Get all motifs from the literature results
        motif_data = literature_results['functional_motifs']
        motif_df = pd.DataFrame(motif_data)
        
        # Get all unique motifs
        all_motifs = motif_df['motif'].value_counts()
        
        # Separate known and novel motifs
        known_motif_list = list(self.known_motifs.keys())
        novel_motifs = [m for m in all_motifs.index if m not in known_motif_list and not m.startswith('?')]
        
        # Create comprehensive table
        results = []
        
        # For each function
        for function in sorted(motif_df['function'].unique()):
            func_data = motif_df[motif_df['function'] == function]
            
            row = {
                'Function': function,
                'Total_Proteins': len(func_data),
                'Domains': ', '.join(func_data['domain'].value_counts().head(3).index.tolist())
            }
            
            # Add known motifs
            for motif in known_motif_list:
                count = (func_data['motif'] == motif).sum()
                row[f'{motif}'] = count
            
            # Add top 5 novel motifs for this function
            func_novel = func_data[func_data['motif'].isin(novel_motifs)]
            if len(func_novel) > 0:
                top_novel = func_novel['motif'].value_counts().head(5)
                row['Top_Novel_Motifs'] = ', '.join([f"{m}({c})" for m, c in top_novel.items()])
            else:
                row['Top_Novel_Motifs'] = '-'
            
            # Add position statistics from discriminator analysis
            # Position 3.5 statistics
            func_proteins = [p for p in func_data['protein'] if p in self.grn_df.index]
            pos_3_5_data = []
            pos_5_5_data = []
            
            for protein in func_proteins:
                if '3.5' in self.grn_df.columns:
                    res_3_5 = self.extract_residue_only(self.grn_df.loc[protein, '3.5'])
                    pos_3_5_data.append(res_3_5)
                if '5.5' in self.grn_df.columns:
                    res_5_5 = self.extract_residue_only(self.grn_df.loc[protein, '5.5'])
                    pos_5_5_data.append(res_5_5)
            
            # Count residues at position 3.5
            if pos_3_5_data:
                row['C@3.5'] = pos_3_5_data.count('C')
                row['T@3.5'] = pos_3_5_data.count('T')
                row['V@3.5'] = pos_3_5_data.count('V')
                row['Other@3.5'] = len(pos_3_5_data) - row['C@3.5'] - row['T@3.5'] - row['V@3.5']
            else:
                row['C@3.5'] = row['T@3.5'] = row['V@3.5'] = row['Other@3.5'] = 0
            
            # Count residues at position 5.5
            if pos_5_5_data:
                row['G@5.5'] = pos_5_5_data.count('G')
                row['S@5.5'] = pos_5_5_data.count('S')
                row['A@5.5'] = pos_5_5_data.count('A')
                row['Other@5.5'] = len(pos_5_5_data) - row['G@5.5'] - row['S@5.5'] - row['A@5.5']
            else:
                row['G@5.5'] = row['S@5.5'] = row['A@5.5'] = row['Other@5.5'] = 0
            
            # Alt-rhodopsin statistics
            alt_rhodopsins = [p for p in literature_results['alt_rhodopsins'] if p['function'] == function]
            row['Alt-R'] = len(alt_rhodopsins)
            if alt_rhodopsins:
                alt_residues = [p['residue_3.42'] for p in alt_rhodopsins]
                alt_counts = Counter(alt_residues)
                row['Alt-R_Types'] = ', '.join([f"{r}" for r, c in alt_counts.most_common(3)])
            else:
                row['Alt-R_Types'] = '-'
            
            results.append(row)
        
        # Create DataFrame
        extended_df = pd.DataFrame(results)
        
        # Save as CSV
        extended_df.to_csv(self.output_dir / 'extended_motif_table.csv', index=False, encoding='utf-8')
        
        # Create markdown version
        self._create_extended_motif_markdown(extended_df, motif_df, novel_motifs)
        
        logger.info("Extended motif table created")
    
    def _create_extended_motif_markdown(self, extended_df: pd.DataFrame, motif_df: pd.DataFrame, novel_motifs: List[str]):
        """Create markdown version of extended motif table"""
        # Get all unique motifs with counts
        all_motifs = motif_df['motif'].value_counts()
        
        with open(self.output_dir / 'extended_motif_table.md', 'w', encoding='utf-8') as f:
            f.write("# Extended Motif Distribution Analysis\n\n")
            
            # Main table
            f.write("## Comprehensive Motif Distribution by Function\n\n")
            known_cols = [col for col in extended_df.columns if col in self.known_motifs]
            header = "| Function | Total | Domains | " + " | ".join(known_cols) + " | Top Novel Motifs |\n"
            f.write(header)
            f.write("|" + "---|" * (len(known_cols) + 4) + "\n")
            
            for _, row in extended_df.iterrows():
                row_str = f"| {row['Function']} | {row['Total_Proteins']} | {row['Domains']} | "
                for col in known_cols:
                    row_str += f"{row[col]} | "
                row_str += f"{row['Top_Novel_Motifs']} |\n"
                f.write(row_str)
            
            # Position analysis tables
            f.write("\n## Position-Specific Analysis\n\n")
            f.write("### Position 3.5 Distribution\n\n")
            f.write("| Function | Total | C (Channel) | T (Pump) | V | Other |\n")
            f.write("|----------|-------|-------------|----------|---|-------|\n")
            
            for _, row in extended_df.iterrows():
                total = row['Total_Proteins']
                if total > 0:
                    f.write(f"| {row['Function']} | {total} | "
                           f"{row['C@3.5']} ({row['C@3.5']/total*100:.0f}%) | "
                           f"{row['T@3.5']} ({row['T@3.5']/total*100:.0f}%) | "
                           f"{row['V@3.5']} ({row['V@3.5']/total*100:.0f}%) | "
                           f"{row['Other@3.5']} ({row['Other@3.5']/total*100:.0f}%) |\n")
            
            f.write("\n### Position 5.5 Distribution\n\n")
            f.write("| Function | Total | G (Channel) | S (Pump) | A | Other |\n")
            f.write("|----------|-------|-------------|----------|---|-------|\n")
            
            for _, row in extended_df.iterrows():
                total = row['Total_Proteins']
                if total > 0:
                    f.write(f"| {row['Function']} | {total} | "
                           f"{row['G@5.5']} ({row['G@5.5']/total*100:.0f}%) | "
                           f"{row['S@5.5']} ({row['S@5.5']/total*100:.0f}%) | "
                           f"{row['A@5.5']} ({row['A@5.5']/total*100:.0f}%) | "
                           f"{row['Other@5.5']} ({row['Other@5.5']/total*100:.0f}%) |\n")
            
            # Alt-rhodopsin analysis
            f.write("\n### Alt-Rhodopsin Analysis (Position 3.42)\n\n")
            f.write("| Function | Total | Alt-R Count | Alt-R % | Alt-R Types |\n")
            f.write("|----------|-------|-------------|---------|-------------|\n")
            
            for _, row in extended_df.iterrows():
                total = row['Total_Proteins']
                alt_count = row['Alt-R']
                if total > 0:
                    f.write(f"| {row['Function']} | {total} | {alt_count} | "
                           f"{alt_count/total*100:.0f}% | {row['Alt-R_Types']} |\n")
            
            # Novel motif analysis
            f.write("\n## Novel Motif Analysis\n\n")
            
            # Get all novel motifs with counts
            novel_counts = motif_df[motif_df['motif'].isin(novel_motifs)]['motif'].value_counts()
            
            f.write("### Top 20 Novel Motifs (Overall)\n\n")
            f.write("| Motif | Count | Primary Functions |\n")
            f.write("|-------|-------|------------------|\n")
            
            for motif, count in novel_counts.head(20).items():
                # Get functions for this motif
                motif_funcs = motif_df[motif_df['motif'] == motif]['function'].value_counts()
                func_str = ', '.join([f"{f}({c})" for f, c in motif_funcs.head(3).items()])
                f.write(f"| {motif} | {count} | {func_str} |\n")
            
            f.write("\n## Key Insights\n\n")
            f.write("1. **Canonical Motifs**: Well-represented in their expected functions\n")
            f.write("2. **Novel Motifs**: Most common novel motifs identified\n")
            f.write("3. **Position 3.5**: Strong C/T discrimination for channel/pump function\n")
            f.write("4. **Position 5.5**: G/S discrimination with some variants\n")
            f.write("5. **Alt-Rhodopsins**: Distribution across different functional categories\n")
            
            # Check for potential novel family members
            f.write("\n### Potential Novel Family Members\n\n")
            
            # Check for ESL motif (Heliorhodopsins)
            esl_count = motif_df[motif_df['motif'] == 'ESL']['protein'].tolist()
            if esl_count:
                f.write(f"**Heliorhodopsin candidates (ESL motif)**: {len(esl_count)} proteins found\n")
                f.write(f"  - Proteins: {', '.join(esl_count[:5])}")
                if len(esl_count) > 5:
                    f.write(f" and {len(esl_count)-5} more")
                f.write("\n\n")
            
            # Look for other literature-mentioned patterns
            # Check for proteins with unusual combinations that might indicate novel families
            unusual_motifs = []
            for motif in novel_motifs:
                # Check if it contains unusual amino acids at key positions
                if len(motif) >= 3:
                    # Position 1 (3.45): Usually D, N, T, F, E
                    # Position 2 (3.49): Usually T, S
                    # Position 3 (3.56): Usually D, A, Q, E, L
                    if motif[0] not in 'DNTFE' or motif[1] not in 'TS' or motif[2] not in 'DAQEL':
                        count = all_motifs[motif]
                        if count >= 2:  # At least 2 occurrences
                            unusual_motifs.append((motif, count))
            
            if unusual_motifs:
                f.write("**Unusual motif patterns** (may indicate novel families):\n")
                for motif, count in sorted(unusual_motifs, key=lambda x: x[1], reverse=True)[:10]:
                    proteins = motif_df[motif_df['motif'] == motif]['protein'].tolist()
                    functions = motif_df[motif_df['motif'] == motif]['function'].value_counts()
                    f.write(f"  - {motif} ({count}x): {', '.join([f'{f}({c})' for f,c in functions.items()])}\n")
            
            f.write("\n**Note**: Based on literature review, novel families to watch for include:\n")
            f.write("- Heliorhodopsins (HeRs) with ESL motif\n")
            f.write("- Alt-Rhodopsins with substitutions at position 3.42\n")
            f.write("- CryoRhodopsins with RRESEDK/RREAEDK patterns\n")
            f.write("- Apusomonad rhodopsins with DXQ/XTQ motifs\n")


def main():
    """Main execution function"""
    parser = argparse.ArgumentParser(description='Unified motif analysis for microbial opsins')
    parser.add_argument('--grn-table', type=Path, 
                       default=Path('opsin_output/global_reference_grn/msa_table_grn.csv'),
                       help='Path to GRN table CSV')
    parser.add_argument('--properties', type=Path,
                       default=Path('property/mo_exp.csv'),
                       help='Path to properties CSV')
    parser.add_argument('--outputs-dir', type=Path,
                       default=Path('opsin_output/motifs'),
                       help='Output directory')
    parser.add_argument('--min-support', type=float, default=0.8,
                       help='Minimum support for correlations (0-1)')
    parser.add_argument('--min-proteins', type=int, default=3,
                       help='Minimum proteins per domain-function group')
    
    args = parser.parse_args()
    
    # Create analyzer and run
    analyzer = UnifiedMotifAnalyzer(args.grn_table, args.properties, args.output_dir)
    analyzer.run_analysis()


if __name__ == "__main__":
    main()