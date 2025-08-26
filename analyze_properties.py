
#!/usr/bin/env python3
"""
Property Distribution Analysis Script
Calculates all distributions and diversity metrics from raw property data
Generates statistics mentioned in diversity.md
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
from collections import defaultdict, Counter
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Set, Optional
import logging
import argparse
import sys
from scipy.stats import entropy
from datetime import datetime

# Add src to path for imports
sys.path.append('src')
from property_mapping import PropertyMapper

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

class PropertyAnalyzer:
    """Analyzes property distributions and diversity metrics"""
    
    def __init__(self, property_path: Path, grn_table_path: Optional[Path] = None, output_dir: Path = Path("property_analysis")):
        self.property_path = property_path
        self.grn_table_path = grn_table_path
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Initialize results storage
        self.results = {
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'dataset_overview': {},
            'function_distribution': {},
            'domain_distribution': {},
            'function_domain_combinations': {},
            'diversity_metrics': {},
            'missing_combinations': {}
        }
        
    def load_data(self):
        """Load property data using PropertyMapper"""
        logger.info("Loading property data...")
        
        # Initialize property mapper
        self.property_mapper = PropertyMapper(self.property_path)
        
        # Get all structure IDs
        if self.grn_table_path and self.grn_table_path.exists():
            # If GRN table provided, use its structure IDs
            grn_df = pd.read_csv(self.grn_table_path, index_col=0)
            structure_ids = grn_df.index.tolist()
            logger.info(f"Using {len(structure_ids)} structure IDs from GRN table")
        else:
            # Otherwise, use all available structures from property mapper
            # This will require getting all possible structure IDs
            structure_ids = []
            # Read property file to get all protein names
            property_df = pd.read_csv(self.property_path)
            for _, row in property_df.iterrows():
                # Generate all possible structure ID formats
                opsin_name = row.get('opsin_name', '')
                short_name = row.get('short_name', '')
                pdb_id = row.get('pdb_id', '')
                
                # Add various ID formats
                if opsin_name:
                    structure_ids.append(opsin_name)
                if short_name and short_name != opsin_name:
                    structure_ids.append(short_name)
                if pdb_id:
                    structure_ids.append(pdb_id)
                    structure_ids.append(f"{pdb_id}_exp")
                    structure_ids.append(f"{pdb_id}_pred")
            
            structure_ids = list(set(structure_ids))
            logger.info(f"Generated {len(structure_ids)} possible structure IDs from property file")
        
        # Map properties to proteins
        self.proteins = {}
        self.experimental_count = 0
        self.predicted_count = 0
        
        for structure_id in structure_ids:
            properties = self.property_mapper.get_properties(structure_id)
            if properties:
                # Get molecular function
                function = properties.get('molecular_function', 'Unknown')
                if not function or function == '?' or pd.isna(function):
                    function = 'Unknown'
                
                # Get domain
                domain = properties.get('domain', 'Unknown')
                if not domain or domain == '?' or pd.isna(domain):
                    domain = 'Unknown'
                
                # Determine if experimental or predicted
                is_experimental = '_exp' in structure_id or properties.get('is_experimental', False)
                if is_experimental:
                    self.experimental_count += 1
                else:
                    self.predicted_count += 1
                
                self.proteins[structure_id] = {
                    'function': function,
                    'domain': domain,
                    'is_experimental': is_experimental,
                    'opsin_name': properties.get('opsin_name', structure_id),
                    'short_name': properties.get('short_name', structure_id),
                    'pdb_id': properties.get('pdb_id', None)
                }
        
        logger.info(f"Loaded properties for {len(self.proteins)} structures")
        logger.info(f"Experimental: {self.experimental_count}, Predicted: {self.predicted_count}")
        
    def calculate_distributions(self):
        """Calculate function and domain distributions"""
        logger.info("Calculating distributions...")
        
        # Dataset overview
        total_structures = len(self.proteins)
        self.results['dataset_overview'] = {
            'total_structures': total_structures,
            'experimental_structures': self.experimental_count,
            'predicted_structures': self.predicted_count,
            'structures_with_annotations': sum(1 for p in self.proteins.values() 
                                             if p['function'] != 'Unknown' and p['domain'] != 'Unknown')
        }
        
        # Function distribution
        function_counts = Counter(p['function'] for p in self.proteins.values())
        self.results['function_distribution'] = {}
        for function, count in sorted(function_counts.items(), key=lambda x: x[1], reverse=True):
            self.results['function_distribution'][function] = {
                'count': count,
                'percentage': round(count / total_structures * 100, 1)
            }
        
        # Domain distribution
        domain_counts = Counter(p['domain'] for p in self.proteins.values())
        self.results['domain_distribution'] = {}
        for domain, count in sorted(domain_counts.items(), key=lambda x: x[1], reverse=True):
            self.results['domain_distribution'][domain] = {
                'count': count,
                'percentage': round(count / total_structures * 100, 1)
            }
        
        # Function-Domain combinations
        combination_counts = Counter((p['function'], p['domain']) for p in self.proteins.values())
        self.results['function_domain_combinations'] = {
            'combinations': {},
            'matrix': {},
            'top_20': []
        }
        
        # Store all combinations
        for (function, domain), count in combination_counts.items():
            key = f"{function}|{domain}"
            self.results['function_domain_combinations']['combinations'][key] = {
                'function': function,
                'domain': domain,
                'count': count,
                'percentage': round(count / total_structures * 100, 1)
            }
        
        # Create combination matrix
        all_functions = sorted(set(p['function'] for p in self.proteins.values()))
        all_domains = sorted(set(p['domain'] for p in self.proteins.values()))
        
        for function in all_functions:
            self.results['function_domain_combinations']['matrix'][function] = {}
            for domain in all_domains:
                count = combination_counts.get((function, domain), 0)
                self.results['function_domain_combinations']['matrix'][function][domain] = count
        
        # Top 20 combinations
        top_combinations = sorted(combination_counts.items(), key=lambda x: x[1], reverse=True)[:20]
        for (function, domain), count in top_combinations:
            self.results['function_domain_combinations']['top_20'].append({
                'function': function,
                'domain': domain,
                'count': count,
                'percentage': round(count / total_structures * 100, 1)
            })
        
    def calculate_diversity_metrics(self):
        """Calculate Shannon diversity indices and other metrics"""
        logger.info("Calculating diversity metrics...")
        
        # Shannon diversity for functions
        function_counts = [p['count'] for p in self.results['function_distribution'].values()]
        function_probs = np.array(function_counts) / sum(function_counts)
        shannon_function = entropy(function_probs, base=np.e)
        
        # Shannon diversity for domains
        domain_counts = [p['count'] for p in self.results['domain_distribution'].values()]
        domain_probs = np.array(domain_counts) / sum(domain_counts)
        shannon_domain = entropy(domain_probs, base=np.e)
        
        # Combination metrics
        unique_combinations = len(self.results['function_domain_combinations']['combinations'])
        n_functions = len(self.results['function_distribution'])
        n_domains = len(self.results['domain_distribution'])
        theoretical_max = n_functions * n_domains
        coverage = unique_combinations / theoretical_max * 100
        
        self.results['diversity_metrics'] = {
            'shannon_diversity_function': round(shannon_function, 3),
            'shannon_diversity_domain': round(shannon_domain, 3),
            'unique_combinations': unique_combinations,
            'theoretical_max_combinations': theoretical_max,
            'combination_coverage_percent': round(coverage, 1)
        }
        
        # Natural combinations (excluding Synthetic and Unknown)
        natural_functions = [f for f in self.results['function_distribution'].keys() if f != 'Unknown']
        natural_domains = [d for d in self.results['domain_distribution'].keys() if d != 'Synthetic']
        
        natural_combinations = sum(1 for (f, d), count in 
                                 Counter((p['function'], p['domain']) for p in self.proteins.values()).items()
                                 if f in natural_functions and d in natural_domains and count > 0)
        
        theoretical_natural = len(natural_functions) * len(natural_domains)
        natural_coverage = natural_combinations / theoretical_natural * 100 if theoretical_natural > 0 else 0
        
        self.results['diversity_metrics']['natural'] = {
            'functions': len(natural_functions),
            'domains': len(natural_domains),
            'combinations_observed': natural_combinations,
            'theoretical_combinations': theoretical_natural,
            'coverage_percent': round(natural_coverage, 1)
        }
        
    def analyze_missing_combinations(self):
        """Analyze which function-domain combinations are missing"""
        logger.info("Analyzing missing combinations...")
        
        # Get all observed combinations
        observed = set((p['function'], p['domain']) for p in self.proteins.values())
        
        # Get all possible combinations
        all_functions = sorted(set(p['function'] for p in self.proteins.values()))
        all_domains = sorted(set(p['domain'] for p in self.proteins.values()))
        all_possible = set((f, d) for f in all_functions for d in all_domains)
        
        # Find missing combinations
        missing = all_possible - observed
        
        self.results['missing_combinations'] = {
            'total_missing': len(missing),
            'missing_list': sorted(list(missing)),
            'by_function': {},
            'by_domain': {}
        }
        
        # Group missing by function
        for function in all_functions:
            missing_domains = [d for (f, d) in missing if f == function]
            if missing_domains:
                self.results['missing_combinations']['by_function'][function] = sorted(missing_domains)
        
        # Group missing by domain
        for domain in all_domains:
            missing_functions = [f for (f, d) in missing if d == domain]
            if missing_functions:
                self.results['missing_combinations']['by_domain'][domain] = sorted(missing_functions)
        
    def generate_visualizations(self):
        """Generate visualization plots"""
        logger.info("Generating visualizations...")
        
        # Create figure directory
        fig_dir = self.output_dir / 'figures'
        fig_dir.mkdir(exist_ok=True)
        
        # 1. Function distribution pie chart
        plt.figure(figsize=(10, 8))
        functions = list(self.results['function_distribution'].keys())
        counts = [self.results['function_distribution'][f]['count'] for f in functions]
        colors = plt.cm.Set3(np.linspace(0, 1, len(functions)))
        
        plt.pie(counts, labels=functions, autopct='%1.1f%%', startangle=90, colors=colors)
        plt.title('Molecular Function Distribution', fontsize=16, fontweight='bold')
        plt.tight_layout()
        plt.savefig(fig_dir / 'function_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 2. Domain distribution bar chart
        plt.figure(figsize=(10, 6))
        domains = list(self.results['domain_distribution'].keys())
        domain_counts = [self.results['domain_distribution'][d]['count'] for d in domains]
        
        bars = plt.bar(domains, domain_counts, color='steelblue', alpha=0.8)
        plt.xlabel('Domain', fontsize=12)
        plt.ylabel('Number of Structures', fontsize=12)
        plt.title('Domain Distribution', fontsize=16, fontweight='bold')
        
        # Add count labels on bars
        for bar, count in zip(bars, domain_counts):
            plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                    str(count), ha='center', va='bottom')
        
        plt.tight_layout()
        plt.savefig(fig_dir / 'domain_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 3. Function-Domain heatmap
        plt.figure(figsize=(12, 8))
        
        # Create matrix for heatmap
        functions = sorted(self.results['function_distribution'].keys())
        domains = sorted(self.results['domain_distribution'].keys())
        
        matrix = np.zeros((len(functions), len(domains)))
        for i, func in enumerate(functions):
            for j, dom in enumerate(domains):
                matrix[i, j] = self.results['function_domain_combinations']['matrix'].get(func, {}).get(dom, 0)
        
        # Create heatmap
        sns.heatmap(matrix, annot=True, fmt='g', cmap='YlOrRd',
                   xticklabels=domains, yticklabels=functions,
                   cbar_kws={'label': 'Number of Structures'})
        
        plt.title('Function-Domain Combination Matrix', fontsize=16, fontweight='bold')
        plt.xlabel('Domain', fontsize=12)
        plt.ylabel('Molecular Function', fontsize=12)
        plt.tight_layout()
        plt.savefig(fig_dir / 'function_domain_heatmap.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        # 4. Top 20 combinations bar chart
        plt.figure(figsize=(12, 8))
        
        top_20 = self.results['function_domain_combinations']['top_20']
        labels = [f"{item['function']}\n{item['domain']}" for item in top_20]
        counts = [item['count'] for item in top_20]
        
        bars = plt.barh(range(len(labels)), counts, color='darkgreen', alpha=0.7)
        plt.yticks(range(len(labels)), labels)
        plt.xlabel('Number of Structures', fontsize=12)
        plt.title('Top 20 Function-Domain Combinations', fontsize=16, fontweight='bold')
        plt.gca().invert_yaxis()
        
        # Add count labels
        for bar, count in zip(bars, counts):
            plt.text(bar.get_width() + 0.3, bar.get_y() + bar.get_height()/2,
                    str(count), ha='left', va='center')
        
        plt.tight_layout()
        plt.savefig(fig_dir / 'top_20_combinations.png', dpi=300, bbox_inches='tight')
        plt.close()
        
        logger.info(f"Saved visualizations to {fig_dir}")
        
    def generate_report(self):
        """Generate a comprehensive report"""
        logger.info("Generating report...")
        
        report_path = self.output_dir / 'property_analysis_report.md'
        
        with open(report_path, 'w', encoding='utf-8') as f:
            f.write("# Microbial Opsin Property Analysis Report\n\n")
            f.write(f"Generated on: {self.results['timestamp']}\n\n")
            
            # Dataset overview
            f.write("## Dataset Overview\n\n")
            overview = self.results['dataset_overview']
            f.write(f"- Total structures analyzed: {overview['total_structures']}\n")
            f.write(f"- Experimental structures: {overview['experimental_structures']}\n")
            f.write(f"- Predicted structures: {overview['predicted_structures']}\n")
            f.write(f"- Structures with both function and domain annotations: {overview['structures_with_annotations']}\n\n")
            
            # Function distribution
            f.write("## Molecular Function Distribution\n\n")
            for function, data in self.results['function_distribution'].items():
                f.write(f"- **{function}**: {data['count']} structures ({data['percentage']}%)\n")
            f.write("\n")
            
            # Domain distribution
            f.write("## Domain Distribution\n\n")
            for domain, data in self.results['domain_distribution'].items():
                f.write(f"- **{domain}**: {data['count']} structures ({data['percentage']}%)\n")
            f.write("\n")
            
            # Top combinations
            f.write("## Function-Domain Combinations\n\n")
            f.write("### Top 20 Most Common Combinations\n\n")
            f.write("| Rank | Molecular Function | Domain | Count | Percentage |\n")
            f.write("|------|-------------------|---------|--------|------------|\n")
            for i, combo in enumerate(self.results['function_domain_combinations']['top_20'], 1):
                f.write(f"| {i} | {combo['function']} | {combo['domain']} | {combo['count']} | {combo['percentage']}% |\n")
            f.write("\n")
            
            # Combination matrix
            f.write("### Complete Combination Matrix\n\n")
            f.write("Rows: Molecular Function, Columns: Domain\n\n")
            
            domains = sorted(self.results['domain_distribution'].keys())
            f.write("| Function / Domain |")
            for domain in domains:
                f.write(f" {domain} |")
            f.write("\n|-------------------|")
            for _ in domains:
                f.write("--------|")
            f.write("\n")
            
            for function in sorted(self.results['function_distribution'].keys()):
                f.write(f"| {function} |")
                for domain in domains:
                    count = self.results['function_domain_combinations']['matrix'].get(function, {}).get(domain, 0)
                    f.write(f" {count if count > 0 else '-'} |")
                f.write("\n")
            f.write("\n")
            
            # Diversity metrics
            f.write("## Diversity Metrics\n\n")
            metrics = self.results['diversity_metrics']
            f.write(f"- **Shannon Diversity Index (Molecular Function)**: {metrics['shannon_diversity_function']}\n")
            f.write(f"- **Shannon Diversity Index (Domain)**: {metrics['shannon_diversity_domain']}\n")
            f.write(f"- **Number of unique function-domain combinations**: {metrics['unique_combinations']}\n")
            f.write(f"- **Theoretical maximum combinations**: {metrics['theoretical_max_combinations']}\n")
            f.write(f"- **Combination coverage**: {metrics['combination_coverage_percent']}%\n\n")
            
            f.write("### Natural Combinations Only (Excluding Synthetic Domain and Unknown Function)\n\n")
            natural = metrics['natural']
            f.write(f"- **Natural functions**: {natural['functions']} (excluding Unknown)\n")
            f.write(f"- **Natural domains**: {natural['domains']} (excluding Synthetic)\n")
            f.write(f"- **Natural function-domain combinations observed**: {natural['combinations_observed']}\n")
            f.write(f"- **Theoretical natural combinations**: {natural['theoretical_combinations']} ")
            f.write(f"({natural['functions']} functions × {natural['domains']} domains)\n")
            f.write(f"- **Natural combination coverage**: {natural['coverage_percent']}%\n\n")
            
            # Missing combinations
            f.write("## Missing Function-Domain Combinations\n\n")
            missing = self.results['missing_combinations']
            f.write(f"Total missing combinations: {missing['total_missing']}\n\n")
            
            if missing['by_function']:
                f.write("### Missing by Function:\n")
                for function, domains in missing['by_function'].items():
                    f.write(f"- **{function}**: Missing in {', '.join(domains)}\n")
                f.write("\n")
            
            if missing['by_domain']:
                f.write("### Missing by Domain:\n")
                for domain, functions in missing['by_domain'].items():
                    f.write(f"- **{domain}**: Missing {', '.join(functions)}\n")
            
        # Save JSON results
        json_path = self.output_dir / 'property_analysis_results.json'
        with open(json_path, 'w', encoding='utf-8') as f:
            json.dump(self.results, f, indent=2, ensure_ascii=False)
        
        logger.info(f"Report saved to {report_path}")
        logger.info(f"JSON results saved to {json_path}")
        
    def run_analysis(self):
        """Run complete property analysis"""
        logger.info("Starting property analysis...")
        
        # Load data
        self.load_data()
        
        # Calculate distributions
        self.calculate_distributions()
        
        # Calculate diversity metrics
        self.calculate_diversity_metrics()
        
        # Analyze missing combinations
        self.analyze_missing_combinations()
        
        # Generate visualizations
        self.generate_visualizations()
        
        # Generate report
        self.generate_report()
        
        logger.info("Property analysis complete!")
        
        return self.results


def main():
    parser = argparse.ArgumentParser(description='Analyze property distributions in microbial opsin dataset')
    parser.add_argument('--property-file', type=Path, default=Path('property/mo_exp.csv'),
                       help='Path to property CSV file')
    parser.add_argument('--grn-table', type=Path, default=Path('opsin_output/global_reference_grn/msa_table_grn.csv'),
                       help='Path to GRN table (optional, for filtering structures)')
    parser.add_argument('--output-dir', type=Path, default=Path('property_analysis'),
                       help='Output directory for results')
    parser.add_argument('--no-grn', action='store_true',
                       help='Analyze all properties without GRN table filtering')
    
    args = parser.parse_args()
    
    # Check if files exist
    if not args.property_file.exists():
        logger.error(f"Property file not found: {args.property_file}")
        sys.exit(1)
    
    grn_path = None if args.no_grn else args.grn_table
    if grn_path and not grn_path.exists():
        logger.warning(f"GRN table not found: {grn_path}. Proceeding without it.")
        grn_path = None
    
    # Run analysis
    analyzer = PropertyAnalyzer(
        property_path=args.property_file,
        grn_table_path=grn_path,
        output_dir=args.output_dir
    )
    
    results = analyzer.run_analysis()
    
    # Print summary
    print("\n" + "="*60)
    print("PROPERTY ANALYSIS SUMMARY")
    print("="*60)
    print(f"Total structures: {results['dataset_overview']['total_structures']}")
    print(f"Unique functions: {len(results['function_distribution'])}")
    print(f"Unique domains: {len(results['domain_distribution'])}")
    print(f"Function-domain combinations: {results['diversity_metrics']['unique_combinations']}")
    print(f"Shannon diversity (function): {results['diversity_metrics']['shannon_diversity_function']}")
    print(f"Shannon diversity (domain): {results['diversity_metrics']['shannon_diversity_domain']}")
    print(f"\nResults saved to: {args.output_dir}")
    print("="*60)


if __name__ == "__main__":
    main()