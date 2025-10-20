# MOGRN Project - TODO and Status

## Project Overview
MOGRN (Microbial Opsin Generic Residue Numbering) is a comprehensive framework for analyzing, comparing, and visualizing experimental and predicted microbial opsin structures using Generic Residue Numbering (GRN) to standardize structural comparisons.

## Recent Achievements (2025-09-03)

### 1. Clustered MO Structure Processing ✓
We successfully processed a new dataset of 21 clustered microbial opsin structures:

- **Created**: `prepare_clustered_mo_data.py` 
  - Loads CIF files from `structures/clustered_mo/`
  - Uses PROTOS framework to create a dataset
  - Extracts sequences from structures
  - **Key Finding**: The atom naming convention differs - CA atoms are identified by `res_atom_name == 'CA'` not `atom_name == 'CA'`

- **Created**: `clustered_mo_intermediate.json`
  - Contains sequence data and mapping information for all 21 structures
  - Maps MO_XXX identifiers to full organism names from `property/mo_small_name_mapping.txt`

### 2. Structure Alignment Pipeline ✓
We implemented a comprehensive alignment workflow to compare clustered structures against all processed structures:

- **Created**: `align_clustered_mo_step_by_step.py`
  - Performs N×M alignment matrix (21 clustered × 199 processed structures)
  - Uses multiple alignment strategies:
    1. Direct RMSD calculation for same-length sequences
    2. Truncated alignment for different lengths
    3. Sequence-based alignment as fallback
  - Successfully aligned all 21 structures with valid RMSD values

- **Results**: 
  - `clustered_mo_alignment_results_v2.json` - Contains best matches and all pairwise RMSDs
  - `clustered_mo_rmsd_matrix_v2.csv` - Full RMSD matrix (21×199)
  - Most clustered structures best match with `S13_Bin138_Proteo_SR_model_0`
  - RMSD values range from ~10-14 Å (reasonable for these structures)

### 3. Visualization Pipeline (In Progress)
Started implementing visualization of alignments using PROTOS visualization capabilities:

- **Created**: `visualize_clustered_alignments.py`
  - Loads alignment results and structures
  - Implements proper GRN filtering for processed structures
  - Applies rotation/translation transformations for alignment visualization
  - Creates both individual alignment visualizations and overview plots

**Current Issues**:
- Dataset loading path issues between different environments
- Need to ensure proper coordinate transformation is applied

## Key Technical Details

### Structure Data Format
- Clustered structures use different atom naming:
  - `atom_name` contains element symbols (C, N, O, S)
  - `res_atom_name` contains atom names (CA, CB, CG, etc.)
  - Must use `res_atom_name == 'CA'` to identify C-alpha atoms

### Alignment Approach
1. Filter structures by chain A
2. Extract CA atoms only
3. For processed structures: filter by GRN (remove residues without valid GRN)
4. Use QCPSuperimposer for optimal superposition
5. Apply transformation: `coords_transformed = coords @ rotation + translation`

### Property Table Requirements
The goal is to create `property/mo_clustered.csv` with same format as `mo_exp.csv`:
- Sequence information
- Organism/species (extracted from mapping file)
- Molecular function (inferred from best match)
- Structure metadata
- Alignment quality metrics

## Next Steps

### Immediate Tasks
1. [ ] Fix visualization script path issues
2. [ ] Complete property table generation with alignment-based function inference
3. [ ] Generate publication-quality alignment visualizations
4. [ ] Create summary statistics and analysis

### Analysis Goals
1. [ ] Determine functional classification of clustered structures based on best matches
2. [ ] Identify novel structural features or variations
3. [ ] Analyze distribution of RMSD values to understand structural diversity
4. [ ] Create phylogenetic/similarity tree based on structural alignment

### Documentation
1. [ ] Document the complete workflow for reproducibility
2. [ ] Create figure legends for visualizations
3. [ ] Summarize biological insights from structural comparisons

## File Structure Summary

```
clustered_mo_analysis/
├── Input Data
│   ├── structures/clustered_mo/*.cif (21 structures)
│   └── property/mo_small_name_mapping.txt
│
├── Processing Scripts
│   ├── prepare_clustered_mo_data.py
│   ├── align_clustered_mo_step_by_step.py
│   └── visualize_clustered_alignments.py
│
├── Results
│   ├── clustered_mo_intermediate.json
│   ├── clustered_mo_alignment_results_v2.json
│   ├── clustered_mo_rmsd_matrix_v2.csv
│   └── property/mo_clustered.csv (to be generated)
│
└── Visualizations
    └── opsin_output/clustered/ (to be populated)
```

## Success Metrics
- ✓ All 21 structures successfully loaded and processed
- ✓ 100% alignment success rate (21/21 structures aligned)
- ✓ All structures have identified best matches
- [ ] Property table generated with complete annotations
- [ ] Visualizations created for top alignments
- [ ] Biological insights documented