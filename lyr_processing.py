"""
Helper module for processing LYR residues (lysine-retinal Schiff base complexes)
"""
import pandas as pd
import numpy as np
from typing import Dict, Any, Optional, Union, List

def convert_lyr_to_lys_ret(df, retinal_name='RET'):
    """
    Converts LYR residues (lysine-retinal Schiff base) to separate LYS and RET residues.
    
    Args:
        df: DataFrame containing atom data with LYR residues
        retinal_name: Name to use for retinal residues (default: 'RET')
        
    Returns:
        DataFrame with LYR residues split into LYS and RET
    """
    # If no LYR residues, return original dataframe
    if 'res_name3l' not in df.columns or not (df['res_name3l'] == 'LYR').any():
        return df
    
    # Make a copy of the dataframe to avoid modifying the original
    result_df = df.copy()
    
    # Find all LYR residues and get their unique residue IDs
    lyr_mask = result_df['res_name3l'] == 'LYR'
    lyr_data = result_df[lyr_mask].copy()
    
    # Process each LYR residue separately
    if 'auth_seq_id' in lyr_data.columns:
        lyr_residue_ids = lyr_data['auth_seq_id'].unique()
        print(f"Found {len(lyr_residue_ids)} LYR residues to convert")
        
        all_lys_atoms = []
        all_ret_atoms = []
        
        for lyr_id in lyr_residue_ids:
            # Get atoms for this specific LYR residue
            residue_mask = (lyr_data['auth_seq_id'] == lyr_id)
            residue_atoms = lyr_data[residue_mask].copy()
            
            # Create LYS atoms (protein backbone and sidechain atoms)
            lys_atoms = residue_atoms[residue_atoms['res_atom_name'].str.contains('N|CA|C|O|CB|CG|CD|CE|NZ')].copy()
            lys_atoms['res_name3l'] = 'LYS'
            # Keep group as ATOM
            if 'group' in lys_atoms.columns:
                lys_atoms['group'] = 'ATOM' 
            
            # Create RET atoms (retinal part - everything else)
            ret_atoms = residue_atoms[~residue_atoms['res_atom_name'].str.contains('N|CA|C|O|CB|CG|CD|CE|NZ')].copy()
            ret_atoms['res_name3l'] = retinal_name
            # Set group to HETATM for retinal
            if 'group' in ret_atoms.columns:
                ret_atoms['group'] = 'HETATM'
            
            # Assign a new residue number for retinal
            # Find maximum residue ID and add 1000 to ensure it doesn't conflict
            if 'auth_seq_id' in result_df.columns:
                max_resid = result_df['auth_seq_id'].astype(str).str.extract('(\d+)', expand=False).astype(float).max()
                new_ret_id = int(max_resid + 1000) if not pd.isna(max_resid) else 1000
                
                # Assign the new residue ID to retinal atoms
                ret_atoms['auth_seq_id'] = new_ret_id
                if 'label_seq_id' in ret_atoms.columns:
                    ret_atoms['label_seq_id'] = new_ret_id
                
                print(f"Converted LYR residue {lyr_id} to LYS + {retinal_name} (new ID: {new_ret_id})")
            
            all_lys_atoms.append(lys_atoms)
            all_ret_atoms.append(ret_atoms)
        
        # Combine all converted residues
        if all_lys_atoms and all_ret_atoms:
            lys_atoms_combined = pd.concat(all_lys_atoms, ignore_index=True)
            ret_atoms_combined = pd.concat(all_ret_atoms, ignore_index=True)
            
            # Remove original LYR residues from the result dataframe
            result_df = result_df[~lyr_mask].copy()
            
            # Add the converted atoms back
            result_df = pd.concat([result_df, lys_atoms_combined, ret_atoms_combined], ignore_index=True)
            
            print(f"Total: Converted {len(lyr_data)} LYR atoms to {len(lys_atoms_combined)} LYS + {len(ret_atoms_combined)} {retinal_name} atoms")
    else:
        # Simpler approach if we don't have auth_seq_id
        # Create LYS atoms (protein backbone and sidechain atoms)
        lys_atoms = lyr_data[lyr_data['res_atom_name'].str.contains('N|CA|C|O|CB|CG|CD|CE|NZ')].copy()
        lys_atoms['res_name3l'] = 'LYS'
        if 'group' in lys_atoms.columns:
            lys_atoms['group'] = 'ATOM'
        
        # Create RET atoms (retinal part - everything else)
        ret_atoms = lyr_data[~lyr_data['res_atom_name'].str.contains('N|CA|C|O|CB|CG|CD|CE|NZ')].copy()
        ret_atoms['res_name3l'] = retinal_name
        if 'group' in ret_atoms.columns:
            ret_atoms['group'] = 'HETATM'
        
        # Remove original LYR residues from the result dataframe
        result_df = result_df[~lyr_mask].copy()
        
        # Add the converted atoms back
        result_df = pd.concat([result_df, lys_atoms, ret_atoms], ignore_index=True)
        
        print(f"Converted {len(lyr_data)} LYR atoms to {len(lys_atoms)} LYS + {len(ret_atoms)} {retinal_name} atoms")
    
    return result_df

def process_lyr_in_structures(structures_dict, retinal_name='RET'):
    """
    Process all structures in the dictionary to convert LYR residues to LYS+RET.
    
    Args:
        structures_dict: Dictionary of structure data with DataFrames
        retinal_name: Name to use for retinal residues (default: 'RET')
        
    Returns:
        Dictionary with properly processed structure data
    """
    result = {}
    
    for pdb_id, data in structures_dict.items():
        # Create a copy to avoid modifying original
        result[pdb_id] = dict(data)
        
        # Process the main structure DataFrame
        if 'df' in data and isinstance(data['df'], pd.DataFrame) and not data['df'].empty:
            result[pdb_id]['df'] = convert_lyr_to_lys_ret(data['df'], retinal_name=retinal_name)
        
        # Process normalized DataFrames if present
        if 'df_norm' in data and isinstance(data['df_norm'], pd.DataFrame) and not data['df_norm'].empty:
            result[pdb_id]['df_norm'] = convert_lyr_to_lys_ret(data['df_norm'], retinal_name=retinal_name)
        
        if 'df_ca_norm' in data and isinstance(data['df_ca_norm'], pd.DataFrame) and not data['df_ca_norm'].empty:
            # For CA-only dataframe, we just need to rename any LYR to LYS (no splitting needed)
            df_ca = data['df_ca_norm'].copy()
            if 'res_name3l' in df_ca.columns:
                df_ca.loc[df_ca['res_name3l'] == 'LYR', 'res_name3l'] = 'LYS'
            result[pdb_id]['df_ca_norm'] = df_ca
            
    return result

def process_lyr_in_processor(processor, retinal_name='RET'):
    """
    Process LYR residues in a CifProcessor's data.
    
    Args:
        processor: CifProcessor or CifBaseProcessor instance
        retinal_name: Name to use for retinal residues (default: 'RET')
        
    Returns:
        None (modifies processor.data in place)
    """
    if not hasattr(processor, 'data') or processor.data is None or processor.data.empty:
        print("Warning: Processor has no data to process")
        return
    
    # Convert LYR residues in the data DataFrame
    processor.data = convert_lyr_to_lys_ret(processor.data, retinal_name=retinal_name)
    
    print(f"Processed LYR residues in {processor.name if hasattr(processor, 'name') else 'processor'}")