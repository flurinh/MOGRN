#!/usr/bin/env python3
"""
Debug script to check how transformations are stored and should be applied
"""

import pickle
import numpy as np

def load_rmsd_cache():
    """Load and examine the RMSD cache"""
    cache_file = "opsin_output/cache/rmsd_cache_a49caa6bee4e02eca6f1239bd483437e.pkl"
    
    with open(cache_file, 'rb') as f:
        cache_data = pickle.load(f)
    
    return cache_data

def main():
    cache_data = load_rmsd_cache()
    alignment_paths = cache_data.get('alignment_paths', {})
    
    print(f"Total alignment paths: {len(alignment_paths)}")
    
    # Find some example transformations involving MerMAID1_model_0
    mermaid_examples = []
    for key, value in alignment_paths.items():
        if 'MerMAID1_model_0' in key:
            mermaid_examples.append((key, value))
            if len(mermaid_examples) >= 5:
                break
    
    print(f"\nFound {len(mermaid_examples)} examples with MerMAID1_model_0:")
    
    for (struct1, struct2), alignment_info in mermaid_examples:
        print(f"\n--- {struct1} -> {struct2} ---")
        print(f"RMSD: {alignment_info['rmsd']:.3f}")
        
        R = np.array(alignment_info['rotation'])
        t = np.array(alignment_info['translation'])
        
        print(f"Rotation matrix shape: {R.shape}")
        print(f"Translation vector shape: {t.shape}")
        print(f"Rotation matrix (first row): {R[0]}")
        print(f"Translation vector: {t}")
        
        # Check if it's orthogonal (should be close to identity)
        should_be_identity = np.dot(R, R.T)
        print(f"R @ R.T should be identity - diagonal: {np.diag(should_be_identity)}")
        print(f"R @ R.T should be identity - off-diagonal max: {np.max(np.abs(should_be_identity - np.eye(3))):.6f}")
        
        # Check if we have both directions
        reverse_key = (struct2, struct1)
        if reverse_key in alignment_paths:
            reverse_info = alignment_paths[reverse_key]
            R_rev = np.array(reverse_info['rotation'])
            t_rev = np.array(reverse_info['translation'])
            
            print(f"Reverse rotation (first row): {R_rev[0]}")
            print(f"Reverse translation: {t_rev}")
            
            # Check if reverse is actually the inverse
            print(f"R_rev should equal R.T - max diff: {np.max(np.abs(R_rev - R.T)):.6f}")
            expected_t_rev = -np.dot(R.T, t)
            print(f"Expected reverse translation: {expected_t_rev}")
            print(f"Actual vs expected t_rev diff: {np.max(np.abs(t_rev - expected_t_rev)):.6f}")

if __name__ == "__main__":
    main()