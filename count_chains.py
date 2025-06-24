"""
Count the number of non-ligand chains in experimental opsin structures.
Uses the same loading method as opsin_analysis_workflow.py
"""

import os
import pandas as pd
from pathlib import Path
from collections import Counter
import json

# Import the same data loading function used in the workflow
from src.data_processing import load_experimental_dataset

def count_non_ligand_chains(cp_mo_exp):
    """
    Count the number of non-ligand chains for each structure.
    
    Args:
        cp_mo_exp: CifBaseProcessor with loaded experimental structures
        
    Returns:
        dict: Dictionary mapping PDB ID to chain count and chain list
    """
    chain_counts = {}
    
    # Common ligand residue names to exclude
    ligand_residues = {'RET', 'LIG', 'LYR', 'HOH', 'SO4', 'PO4', 'GOL', 'EDO', 
                      'PEG', 'CL', 'NA', 'K', 'CA', 'MG', 'ZN', 'FE', 'MN',
                      'ACT', 'FMT', 'ACE', 'NH2', 'TRS', 'MES', 'EPE', 'MPD',
                      'DMS', 'BME', 'DTT', 'CIT', 'MAL', 'SUC', 'TAR', 'MLI',
                      'FLC', 'ADP', 'ATP', 'GTP', 'GDP', 'NAD', 'FAD', 'FMN',
                      'HEM', 'BCL', 'CHL', 'PHO', 'U10', 'CDL', 'LMT', 'CLR',
                      'OLC', 'OLA', 'PLM', 'STE', 'PGE', 'PG4', 'PE4', 'P6G',
                      '12P', '1PE', 'DGA', 'DAG', 'MAG', 'POP', 'PCW', 'POV'}
    
    print(f"\nAnalyzing {len(cp_mo_exp.pdb_ids)} experimental structures...")
    
    for pdb_id in cp_mo_exp.pdb_ids:
        # Get all atoms for this PDB ID
        df_pdb = cp_mo_exp.data[cp_mo_exp.data['pdb_id'] == pdb_id]
        
        # Get unique chain IDs
        all_chains = df_pdb['auth_chain_id'].unique()
        
        # Filter out chains that only contain ligands/water
        protein_chains = []
        
        for chain_id in all_chains:
            # Get all residues in this chain
            chain_df = df_pdb[df_pdb['auth_chain_id'] == chain_id]
            residue_names = chain_df['res_name3l'].unique()
            
            # Check if chain contains any non-ligand residues (amino acids)
            amino_acids = {'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 
                          'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 
                          'THR', 'TRP', 'TYR', 'VAL', 'MSE', 'SEP', 'TPO', 'PTR',
                          'HYP', 'MLY', 'CME', 'CSO', 'CSD', 'LYR'}  # LYR is lysine-retinal
            
            # If chain contains amino acids, it's a protein chain
            if any(res in amino_acids for res in residue_names):
                protein_chains.append(chain_id)
        
        chain_counts[pdb_id] = {
            'total_chains': len(all_chains),
            'protein_chains': len(protein_chains),
            'protein_chain_ids': sorted(protein_chains),
            'all_chain_ids': sorted(all_chains)
        }
    
    return chain_counts

def analyze_chain_distribution(chain_counts):
    """
    Analyze the distribution of chain counts across structures.
    
    Args:
        chain_counts: Dictionary from count_non_ligand_chains
        
    Returns:
        dict: Statistics about chain distribution
    """
    # Count distribution
    protein_chain_distribution = Counter()
    total_chain_distribution = Counter()
    
    for pdb_id, info in chain_counts.items():
        protein_chain_distribution[info['protein_chains']] += 1
        total_chain_distribution[info['total_chains']] += 1
    
    # Find structures with multiple protein chains
    multi_chain_structures = {
        pdb_id: info for pdb_id, info in chain_counts.items() 
        if info['protein_chains'] > 1
    }
    
    # Statistics
    stats = {
        'total_structures': len(chain_counts),
        'protein_chain_distribution': dict(protein_chain_distribution),
        'total_chain_distribution': dict(total_chain_distribution),
        'multi_chain_count': len(multi_chain_structures),
        'multi_chain_pdb_ids': list(multi_chain_structures.keys()),
        'average_protein_chains': sum(info['protein_chains'] for info in chain_counts.values()) / len(chain_counts),
        'average_total_chains': sum(info['total_chains'] for info in chain_counts.values()) / len(chain_counts)
    }
    
    return stats, multi_chain_structures

def main():
    """Main function to run the chain counting analysis."""
    
    print("=" * 80)
    print("CHAIN COUNTING ANALYSIS FOR EXPERIMENTAL OPSIN STRUCTURES")
    print("=" * 80)
    
    # Load experimental structures using the same method as the workflow
    print("\nLoading experimental structures (mo_exp dataset)...")
    cp_mo_exp = load_experimental_dataset('mo_exp')
    
    # Count chains
    print("\nCounting non-ligand chains in each structure...")
    chain_counts = count_non_ligand_chains(cp_mo_exp)
    
    # Analyze distribution
    print("\nAnalyzing chain distribution...")
    stats, multi_chain_structures = analyze_chain_distribution(chain_counts)
    
    # Print results
    print("\n" + "=" * 80)
    print("RESULTS")
    print("=" * 80)
    
    print(f"\nTotal structures analyzed: {stats['total_structures']}")
    print(f"Average protein chains per structure: {stats['average_protein_chains']:.2f}")
    print(f"Average total chains per structure: {stats['average_total_chains']:.2f}")
    
    print("\nProtein chain distribution:")
    for n_chains, count in sorted(stats['protein_chain_distribution'].items()):
        print(f"  {n_chains} protein chain(s): {count} structures ({count/stats['total_structures']*100:.1f}%)")
    
    print(f"\nStructures with multiple protein chains: {stats['multi_chain_count']}")
    if stats['multi_chain_count'] > 0:
        print("Multi-chain structures:")
        for pdb_id, info in sorted(multi_chain_structures.items()):
            print(f"  {pdb_id}: {info['protein_chains']} chains {info['protein_chain_ids']}")
    
    # Save results
    output_dir = Path('opsin_output')
    output_dir.mkdir(exist_ok=True)
    
    # Save detailed chain counts
    with open(output_dir / 'chain_counts.json', 'w') as f:
        json.dump(chain_counts, f, indent=2)
    
    # Save statistics
    with open(output_dir / 'chain_statistics.json', 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"\nResults saved to {output_dir}/")
    print("  - chain_counts.json: Detailed chain information for each structure")
    print("  - chain_statistics.json: Summary statistics")
    
    # Create a summary CSV for easy viewing
    summary_data = []
    for pdb_id, info in sorted(chain_counts.items()):
        summary_data.append({
            'PDB_ID': pdb_id,
            'Protein_Chains': info['protein_chains'],
            'Total_Chains': info['total_chains'],
            'Protein_Chain_IDs': ', '.join(info['protein_chain_ids']),
            'All_Chain_IDs': ', '.join(info['all_chain_ids'])
        })
    
    summary_df = pd.DataFrame(summary_data)
    summary_df.to_csv(output_dir / 'chain_summary.csv', index=False)
    print("  - chain_summary.csv: Summary table of all structures")

if __name__ == '__main__':
    main()