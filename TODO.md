# TODO: LYR (Lysine-Retinal Schiff Base) Processing

## Background
Some opsin structures have residues labeled as "LYR" which represent a lysine residue covalently bound to a retinal molecule (Schiff base). For proper analysis, these need to be split into separate LYS (protein) and RET (retinal) components.

## Implementation Tasks

1. **Create LYR Processing Module**
   - Create `lyr_processing.py` with functions to handle LYR residues:
     - `convert_lyr_to_lys_ret(df, retinal_name='RET')`: Separate LYR into LYS+RET atoms at the DataFrame level
     - `process_lyr_in_structures(structures_dict, retinal_name='RET')`: Process all structures in a dictionary
     - `process_lyr_in_processor(processor, retinal_name='RET')`: Process all data in a CifBaseProcessor

2. **Integration Points**
   - Update `data_processing.py`:
     - Import LYR processing utilities
     - Add `process_lyr_in_loaded_structures` function to apply after loading
   
   - Update `common_utils.py`:
     - Modify `find_retinal_within_cutoff` to also recognize LYR residues
     - Add LYR detection in both input and output

3. **Structure Processing Flow**
   - When structure is loaded, identify LYR residues
   - Split using atom names:
     - LYS atoms: N, CA, C, O, CB, CG, CD, CE, NZ (protein backbone/sidechain)
     - RET atoms: All other atoms (retinal part)
   - For LYS atoms:
     - Set `res_name3l = 'LYS'`
     - Keep original residue number
     - Mark as `ATOM` type
   - For RET atoms:
     - Set `res_name3l = 'RET'`
     - Assign a new unique residue number (e.g., max_resid + 1000)
     - Mark as `HETATM` type

4. **Testing**
   - Test with structures containing LYR residues
   - Verify that resulting LYS residues have correct atom counts
   - Verify that RET residues are properly identified in later processing

## Implementation Design Notes

1. **Identification**: Filter rows where `df['res_name3l'] == 'LYR'`

2. **Atom Classification**:
   ```python
   # Lysine atoms (backbone + sidechain)
   lys_atoms = lyr_data[lyr_data['res_atom_name'].str.contains('N|CA|C|O|CB|CG|CD|CE|NZ')].copy()
   # Retinal atoms (everything else)
   ret_atoms = lyr_data[~lyr_data['res_atom_name'].str.contains('N|CA|C|O|CB|CG|CD|CE|NZ')].copy()
   ```

3. **Residue Numbering**:
   ```python
   # Find maximum residue ID and add 1000 for the new retinal
   max_resid = df['auth_seq_id'].astype(str).str.extract('(\d+)', expand=False).astype(float).max()
   new_ret_id = int(max_resid + 1000) if not pd.isna(max_resid) else 1000
   ret_atoms['auth_seq_id'] = new_ret_id
   ```

4. **Integration Strategy**:
   - Apply LYR processing right after loading structures but before chain filtering
   - This ensures proper chain separation and retinal distance calculations