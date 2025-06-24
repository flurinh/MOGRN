#!/usr/bin/env python3
"""
Check alignment between GRN table entries and helix boundaries
"""

import pandas as pd
import json
from pathlib import Path
from typing import Dict, List, Tuple

def load_grn_entry(grn_table_path: str, protein_id: str) -> Dict[str, str]:
    """Load a single GRN table entry and extract sequence info"""
    df = pd.read_csv(grn_table_path, index_col=0)
    
    if protein_id not in df.index:
        raise ValueError(f"Protein {protein_id} not found in GRN table")
    
    # Get the row as a dict
    row = df.loc[protein_id]
    sequence_dict = {}
    
    for col, value in row.items():
        if pd.notna(value) and value != '-':
            sequence_dict[col] = value
    
    return sequence_dict

def load_helix_boundaries(helices_path: str, protein_id: str) -> Dict[str, List[int]]:
    """Load helix boundaries for a specific protein"""
    with open(helices_path, 'r') as f:
        helices_data = json.load(f)
    
    if protein_id not in helices_data:
        raise ValueError(f"Protein {protein_id} not found in helices file")
    
    return helices_data[protein_id]

def extract_residue_info(grn_value: str) -> Tuple[str, int]:
    """Extract residue type and position from GRN value like 'K296'"""
    if not grn_value or grn_value == '-':
        return None, None
    
    residue = grn_value[0]
    position = int(grn_value[1:])
    return residue, position

def check_helix_alignment(sequence_dict: Dict[str, str], helix_boundaries: Dict[str, List[int]]) -> Dict:
    """Check if GRN positions align with helix boundaries"""
    
    results = {
        'helix_positions': {},
        'misalignments': [],
        'summary': {}
    }
    
    # For each helix
    for helix_num in ['1', '2', '3', '4', '5', '6', '7']:
        if helix_num not in helix_boundaries:
            continue
            
        start, end = helix_boundaries[helix_num]
        helix_positions = []
        
        # Find all GRN positions for this helix
        for grn_pos, grn_value in sequence_dict.items():
            if grn_pos.startswith(f"{helix_num}."):
                residue, position = extract_residue_info(grn_value)
                if position is not None:
                    helix_positions.append({
                        'grn': grn_pos,
                        'residue': residue,
                        'position': position,
                        'in_range': start <= position <= end
                    })
        
        # Sort by position
        helix_positions.sort(key=lambda x: x['position'])
        
        # Store results
        results['helix_positions'][helix_num] = {
            'boundaries': [start, end],
            'positions': helix_positions,
            'length': end - start + 1
        }
        
        # Check for misalignments
        if helix_positions:
            actual_start = helix_positions[0]['position']
            actual_end = helix_positions[-1]['position']
            
            if actual_start < start or actual_end > end:
                results['misalignments'].append({
                    'helix': helix_num,
                    'defined_range': [start, end],
                    'actual_range': [actual_start, actual_end],
                    'start_diff': actual_start - start,
                    'end_diff': actual_end - end
                })
    
    return results

def print_analysis(protein_id: str, results: Dict):
    """Print detailed analysis of helix alignment"""
    print(f"\n{'='*80}")
    print(f"Helix Alignment Analysis for: {protein_id}")
    print(f"{'='*80}\n")
    
    for helix_num in ['1', '2', '3', '4', '5', '6', '7']:
        if helix_num not in results['helix_positions']:
            continue
            
        helix_data = results['helix_positions'][helix_num]
        print(f"\nHelix {helix_num}:")
        print(f"  Defined boundaries: {helix_data['boundaries'][0]}-{helix_data['boundaries'][1]} (length: {helix_data['length']})")
        
        positions = helix_data['positions']
        if positions:
            print(f"  Actual range: {positions[0]['position']}-{positions[-1]['position']}")
            print(f"  Number of GRN positions: {len(positions)}")
            
            # Check .50 position
            fifty_pos = next((p for p in positions if p['grn'] == f"{helix_num}.50"), None)
            if fifty_pos:
                print(f"  Position {helix_num}.50: {fifty_pos['residue']}{fifty_pos['position']}")
            
            # Show any out-of-range positions
            out_of_range = [p for p in positions if not p['in_range']]
            if out_of_range:
                print(f"  ⚠️  Out-of-range positions:")
                for pos in out_of_range:
                    print(f"      {pos['grn']}: {pos['residue']}{pos['position']}")
    
    if results['misalignments']:
        print(f"\n{'⚠️  MISALIGNMENTS DETECTED ':-^80}")
        for mis in results['misalignments']:
            print(f"\nHelix {mis['helix']}:")
            print(f"  Defined: {mis['defined_range'][0]}-{mis['defined_range'][1]}")
            print(f"  Actual:  {mis['actual_range'][0]}-{mis['actual_range'][1]}")
            print(f"  Start difference: {mis['start_diff']:+d}")
            print(f"  End difference:   {mis['end_diff']:+d}")

def main():
    # Paths
    grn_table_path = "opsin_output/opsin_grn_tables/residue_table_grn.csv"
    helices_path = "property/helices_curated.json"
    
    # Test with first protein
    protein_id = "MerMAID1_model_0"
    
    try:
        # Load data
        print(f"Loading GRN entry for {protein_id}...")
        sequence_dict = load_grn_entry(grn_table_path, protein_id)
        
        print(f"Loading helix boundaries...")
        helix_boundaries = load_helix_boundaries(helices_path, protein_id)
        
        # Analyze alignment
        results = check_helix_alignment(sequence_dict, helix_boundaries)
        
        # Print results
        print_analysis(protein_id, results)
        
        # Also save detailed results
        output_file = f"helix_alignment_check_{protein_id}.json"
        with open(output_file, 'w') as f:
            json.dump({
                'protein_id': protein_id,
                'helix_boundaries': helix_boundaries,
                'analysis': results
            }, f, indent=2)
        
        print(f"\nDetailed results saved to: {output_file}")
        
    except Exception as e:
        print(f"Error: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main()