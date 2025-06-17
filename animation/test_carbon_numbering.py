#!/usr/bin/env python3
"""
Test and validate carbon numbering inference for retinal.
Compare inferred carbon positions against experimental structures with canonical names.
"""

import os
import sys
import pandas as pd
import numpy as np
from typing import Dict, List, Tuple
from collections import defaultdict

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

from src.data_processing import load_opsin_structures
from src.retinal_carbon_mapping_correct import assign_retinal_carbons_correct
from scipy.spatial.distance import cdist


def get_retinal_atoms(structure_df: pd.DataFrame) -> pd.DataFrame:
    """Extract retinal atoms from structure DataFrame."""
    retinal_mask = structure_df['res_name3l'].isin(['RET', 'LYR'])
    retinal_df = structure_df[retinal_mask].copy()
    
    if 'LYR' in retinal_df['res_name3l'].values:
        lys_atoms = {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'}
        retinal_df = retinal_df[~retinal_df['res_atom_name'].isin(lys_atoms)]
    
    return retinal_df


def build_connectivity_graph(atoms_df: pd.DataFrame, bond_cutoff: float = 1.7) -> Dict[int, List[int]]:
    """Build connectivity graph based on distances between atoms."""
    coords = atoms_df[['x', 'y', 'z']].values
    distances = cdist(coords, coords)
    
    connectivity = defaultdict(list)
    
    for i in range(len(atoms_df)):
        for j in range(i + 1, len(atoms_df)):
            if distances[i, j] < bond_cutoff:
                connectivity[i].append(j)
                connectivity[j].append(i)
                
    return dict(connectivity)


def get_canonical_carbons(retinal_df: pd.DataFrame) -> Dict[str, pd.Series]:
    """Extract carbons with canonical names from retinal."""
    canonical_carbons = {}
    
    # Look for C1 through C20 (some retinals have extended numbering)
    for i in range(1, 21):
        carbon_name = f'C{i}'
        carbon_atoms = retinal_df[retinal_df['res_atom_name'] == carbon_name]
        if not carbon_atoms.empty:
            canonical_carbons[carbon_name] = carbon_atoms.iloc[0]
    
    return canonical_carbons


def find_ring_system(atoms_df: pd.DataFrame, connectivity: Dict[int, List[int]]) -> List[int]:
    """Find the beta-ionone ring system (6-membered ring)."""
    print(atoms_df.columns)
    carbon_mask = atoms_df['atom_name'] == 'C'
    carbon_indices = atoms_df.index[carbon_mask].tolist()
    
    # Look for 6-membered rings
    rings = []
    
    for start_idx in carbon_indices:
        # DFS to find rings starting from this carbon
        visited = set()
        path = []
        
        def dfs(current, target, depth):
            if depth > 6:  # Don't search too deep
                return False
            
            if current == target and depth == 6:
                rings.append(path[:])
                return True
            
            visited.add(current)
            path.append(current)
            
            for neighbor in connectivity.get(current, []):
                if neighbor in carbon_indices and (neighbor not in visited or (neighbor == target and depth >= 5)):
                    if dfs(neighbor, target, depth + 1):
                        break
            
            path.pop()
            return False
        
        dfs(start_idx, start_idx, 0)
    
    # Return the first ring found (should be the ionone ring)
    return rings[0] if rings else []


def infer_carbon_numbering_old(retinal_df: pd.DataFrame) -> Dict[str, pd.Series]:
    """
    Infer carbon numbering for retinal based on connectivity and structure.
    
    Returns:
        Dictionary mapping carbon names (C1-C20) to atom rows
    """
    # Reset index for easier handling
    retinal_df = retinal_df.reset_index(drop=True)

    print(retinal_df.columns)
    # Build connectivity
    connectivity = build_connectivity_graph(retinal_df)
    
    # Get all carbon atoms
    carbon_mask = retinal_df['atom_name'] == 'C'
    carbon_df = retinal_df[carbon_mask]
    carbon_indices = carbon_df.index.tolist()
    
    print(f"Found {len(carbon_indices)} carbon atoms in retinal")
    
    # Find the beta-ionone ring
    ring_indices = find_ring_system(retinal_df, connectivity)
    print(f"Found ring with {len(ring_indices)} carbons")
    
    # Identify different carbon types based on connectivity
    carbon_types = {}
    for idx in carbon_indices:
        neighbors = connectivity.get(idx, [])
        carbon_neighbors = [n for n in neighbors if n in carbon_indices]
        carbon_types[idx] = {
            'total_neighbors': len(neighbors),
            'carbon_neighbors': len(carbon_neighbors),
            'in_ring': idx in ring_indices
        }
    
    # Find methyl carbons (connected to only one carbon)
    methyl_carbons = [idx for idx, info in carbon_types.items() 
                      if info['carbon_neighbors'] == 1]
    print(f"Found {len(methyl_carbons)} methyl carbons")
    
    # The retinal numbering system:
    # C1, C5: methyl groups on the ring
    # C2, C3, C4, C6: ring carbons
    # C7-C15: polyene chain
    # C16, C17: methyl groups on C1
    # C18: methyl group on C5
    # C19, C20: methyl groups on C9 and C13
    
    inferred_carbons = {}
    
    # Start by finding C6 - the ring carbon connected to the polyene chain
    # It should be in the ring and connected to a carbon outside the ring
    c6_candidates = []
    for idx in ring_indices:
        neighbors = connectivity.get(idx, [])
        for n in neighbors:
            if n in carbon_indices and n not in ring_indices:
                c6_candidates.append(idx)
                break
    
    if c6_candidates:
        # C6 is typically the one with specific geometry
        c6_idx = c6_candidates[0]  # For now, take the first
        inferred_carbons['C6'] = retinal_df.loc[c6_idx]
        
        # Find C5 - connected to C6 in the ring
        for n in connectivity.get(c6_idx, []):
            if n in ring_indices:
                # Check if this could be C5 (should have a methyl group)
                for nn in connectivity.get(n, []):
                    if nn in methyl_carbons:
                        inferred_carbons['C5'] = retinal_df.loc[n]
                        inferred_carbons['C18'] = retinal_df.loc[nn]  # Methyl on C5
                        break
        
        # Trace the polyene chain starting from C6
        current = c6_idx
        visited = {c6_idx}
        visited.update(ring_indices)  # Don't go back into the ring
        
        chain_carbons = []
        
        # Find the carbon connected to C6 that's not in the ring (this is C7)
        for n in connectivity.get(c6_idx, []):
            if n in carbon_indices and n not in ring_indices:
                current = n
                break
        
        # Trace the polyene chain
        carbon_number = 7
        while current is not None and carbon_number <= 15:
            inferred_carbons[f'C{carbon_number}'] = retinal_df.loc[current]
            visited.add(current)
            
            # Find methyl groups attached to this carbon
            for n in connectivity.get(current, []):
                if n in methyl_carbons and n not in visited:
                    if carbon_number == 9:
                        inferred_carbons['C19'] = retinal_df.loc[n]
                    elif carbon_number == 13:
                        inferred_carbons['C20'] = retinal_df.loc[n]
            
            # Find next carbon in chain
            next_carbon = None
            for n in connectivity.get(current, []):
                if n in carbon_indices and n not in visited:
                    # Make sure it's part of the main chain, not a methyl
                    if carbon_types[n]['carbon_neighbors'] >= 2:
                        next_carbon = n
                        break
            
            current = next_carbon
            carbon_number += 1
    
    # Now work on the ring carbons
    # Find C1 - should have two methyl groups (C16, C17)
    for idx in ring_indices:
        if idx not in [v.name for v in inferred_carbons.values()]:
            methyl_count = 0
            methyl_indices = []
            for n in connectivity.get(idx, []):
                if n in methyl_carbons:
                    methyl_count += 1
                    methyl_indices.append(n)
            
            if methyl_count == 2:
                inferred_carbons['C1'] = retinal_df.loc[idx]
                inferred_carbons['C16'] = retinal_df.loc[methyl_indices[0]]
                inferred_carbons['C17'] = retinal_df.loc[methyl_indices[1]]
                break
    
    # Fill in remaining ring carbons (C2, C3, C4)
    remaining_ring = [idx for idx in ring_indices 
                      if idx not in [v.name for v in inferred_carbons.values()]]
    
    # Order them by connectivity
    if 'C1' in inferred_carbons and 'C6' in inferred_carbons:
        c1_idx = inferred_carbons['C1'].name
        c6_idx = inferred_carbons['C6'].name
        
        # Trace from C1 to C6
        current = c1_idx
        visited_ring = {c1_idx}
        carbon_number = 2
        
        while current != c6_idx and carbon_number <= 4:
            for n in connectivity.get(current, []):
                if n in remaining_ring and n not in visited_ring:
                    inferred_carbons[f'C{carbon_number}'] = retinal_df.loc[n]
                    visited_ring.add(n)
                    current = n
                    carbon_number += 1
                    break
    
    return inferred_carbons


def compare_carbon_assignments(canonical: Dict[str, pd.Series], 
                              inferred: Dict[str, pd.Series]) -> Dict[str, bool]:
    """
    Compare canonical and inferred carbon assignments.
    
    Returns:
        Dictionary mapping carbon names to whether they match
    """
    results = {}
    
    print("\nComparing carbon assignments:")
    print("-" * 60)
    print(f"{'Carbon':<10} {'Canonical Atom':<15} {'Inferred Atom':<15} {'Match':<10}")
    print("-" * 60)
    
    all_carbons = sorted(set(canonical.keys()) | set(inferred.keys()), 
                        key=lambda x: int(x[1:]))
    
    for carbon in all_carbons:
        if carbon in canonical and carbon in inferred:
            # Compare by coordinates (they should be the same atom)
            canonical_coords = canonical[carbon][['x', 'y', 'z']].values
            inferred_coords = inferred[carbon][['x', 'y', 'z']].values
            
            distance = np.linalg.norm(canonical_coords - inferred_coords)
            match = distance < 0.01  # Should be exactly the same atom
            
            results[carbon] = match
            
            canonical_name = canonical[carbon]['res_atom_name']
            inferred_name = inferred[carbon]['res_atom_name']
            
            print(f"{carbon:<10} {canonical_name:<15} {inferred_name:<15} {'✓' if match else '✗':<10}")
        elif carbon in canonical:
            results[carbon] = False
            print(f"{carbon:<10} {canonical[carbon]['res_atom_name']:<15} {'Not found':<15} {'✗':<10}")
        elif carbon in inferred:
            results[carbon] = False
            print(f"{carbon:<10} {'Not found':<15} {inferred[carbon]['res_atom_name']:<15} {'✗':<10}")
    
    print("-" * 60)
    
    # Summary
    matched = sum(1 for v in results.values() if v)
    total = len(results)
    print(f"\nMatched: {matched}/{total} ({matched/total*100:.1f}%)")
    
    return results


def test_carbon_numbering():
    """Test carbon numbering inference on structures with known canonical names."""
    
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
    
    # Find experimental structures (more likely to have canonical names)
    experimental_structures = {
        sid: sdata for sid, sdata in processed_structures.items()
        if '_pred' not in sid and '_model' not in sid
    }
    
    print(f"\nFound {len(experimental_structures)} experimental structures")
    
    # Test on each experimental structure
    test_results = {}
    
    for struct_id, struct_data in experimental_structures.items():
        print(f"\n{'='*80}")
        print(f"Testing structure: {struct_id}")
        print('='*80)
        
        # Get retinal atoms
        retinal_df = get_retinal_atoms(struct_data['df'])
        
        if retinal_df.empty:
            print(f"No retinal found in {struct_id}")
            continue
        
        # Get canonical carbons
        canonical = get_canonical_carbons(retinal_df)
        
        if not canonical:
            print(f"No canonical carbon names found in {struct_id}")
            continue
        
        print(f"Found {len(canonical)} carbons with canonical names")
        
        # Since the structures already have canonical names, 
        # the test should compare those names after removing them
        # and seeing if we can infer them back correctly
        
        # Create a copy without the canonical names
        retinal_df_copy = retinal_df.copy()
        retinal_df_copy['res_atom_name_original'] = retinal_df_copy['res_atom_name']
        retinal_df_copy['res_atom_name'] = 'C'  # Remove the canonical names
        
        # Infer carbon numbering using the corrected algorithm
        from src.retinal_carbon_mapping_fixed import assign_retinal_carbons_robust
        inferred_indices = assign_retinal_carbons_robust(retinal_df_copy)
        
        # Convert indices to Series for comparison
        inferred = {}
        for carbon_name, idx in inferred_indices.items():
            if idx < len(retinal_df_copy):
                inferred[carbon_name] = retinal_df_copy.iloc[idx]
        
        print(f"Inferred {len(inferred)} carbon positions")
        
        # Compare
        results = compare_carbon_assignments(canonical, inferred)
        test_results[struct_id] = results
    
    # Overall summary
    print(f"\n{'='*80}")
    print("OVERALL TEST RESULTS")
    print('='*80)
    
    all_carbons = set()
    for results in test_results.values():
        all_carbons.update(results.keys())
    
    carbon_accuracy = {}
    for carbon in sorted(all_carbons, key=lambda x: int(x[1:])):
        matches = 0
        total = 0
        for results in test_results.values():
            if carbon in results:
                total += 1
                if results[carbon]:
                    matches += 1
        
        if total > 0:
            carbon_accuracy[carbon] = matches / total * 100
            print(f"{carbon}: {matches}/{total} correct ({carbon_accuracy[carbon]:.1f}%)")
    
    # Save detailed results
    results_file = os.path.join(output_dir, 'carbon_numbering_test_results.txt')
    with open(results_file, 'w') as f:
        f.write("Carbon Numbering Test Results\n")
        f.write("="*80 + "\n\n")
        
        for struct_id, results in test_results.items():
            f.write(f"Structure: {struct_id}\n")
            f.write("-"*40 + "\n")
            for carbon, match in sorted(results.items(), key=lambda x: int(x[0][1:])):
                f.write(f"  {carbon}: {'MATCH' if match else 'MISMATCH'}\n")
            f.write("\n")
        
        f.write("\nOverall Accuracy by Carbon:\n")
        f.write("-"*40 + "\n")
        for carbon, accuracy in sorted(carbon_accuracy.items(), key=lambda x: int(x[0][1:])):
            f.write(f"{carbon}: {accuracy:.1f}%\n")
    
    print(f"\nDetailed results saved to: {results_file}")


if __name__ == "__main__":
    test_carbon_numbering()