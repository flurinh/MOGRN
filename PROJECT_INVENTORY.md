# MOGRN Project Inventory

## Overview
This document provides an exhaustive inventory of all scripts and their functionality in the MOGRN (Microbial Opsin Generic Residue Numbering) project.

---

## CORE WORKFLOW FILES (Keep)

### `opsin_analysis_workflow.py`
**Purpose**: Main analysis pipeline orchestrating the entire workflow
**Functionality**:
- Loads opsin structures from protos datasets
- Runs structure prediction validation (Boltz-1 vs experimental)
- Performs structural conservation analysis
- Generates GRN (Generic Residue Numbering) tables
- Creates RMSD matrices and clustering
**Dependencies**: src/*, protos
**Status**: CORE - Keep

### `plot.py`
**Purpose**: Visualization script for publication figures
**Functionality**:
- Generates Figure 2: Prediction validation plots (RMSD distributions)
- Generates Figure 3: Structural clustering heatmaps
- Creates GRN distance profile plots
- Produces alignment visualizations
**Dependencies**: opsin_output/curated_grn.csv, protos
**Status**: CORE - Keep

### `prepare_data.py`
**Purpose**: Data preparation and protos dataset initialization
**Functionality**:
- Sets up ProtosPaths configuration
- Loads opsin property data from CSV/Excel
- Creates and registers structure datasets in protos
- Handles legacy dataset migration
**Dependencies**: protos, property/mo_exp.csv
**Status**: CORE - Needs refactoring for new protos API

---

## SOURCE MODULES (`src/`)

### `src/data_processing.py`
**Purpose**: Data loading and preprocessing utilities
**Functions**:
- `load_experimental_dataset()`: Load dataset via CifProcessor
- `filter_structures_by_chain()`: Extract specific chain
- `filter_structures_by_chain_and_retinal()`: Filter by chain + retinal proximity
- `load_opsin_property_data()`: Load property CSV
- `load_opsin_structures()`: Main structure loading function
**Status**: CORE

### `src/error_analysis.py`
**Purpose**: Error/RMSD calculation between structures
**Functions**:
- `calculate_atom_level_errors()`: Per-atom error statistics
- `calculate_structure_errors()`: Overall structure RMSD
- `compute_retinal_mean_closest_distance()`: Retinal L-RMSD
**Status**: CORE

### `src/structure_comparison.py`
**Purpose**: Structure alignment and comparison
**Functions**:
- `compute_all_vs_all_rmsd_improved()`: Pairwise RMSD matrix
- `compare_structures()`: Compare two structures
- `create_unified_structure_mapping()`: Create residue mapping
- `get_structure_alignment()`: CEalign wrapper
- `calculate_binding_pocket_rmsd_for_pairs()`: Pocket RMSD
**Constants**: `HELIX_RET_BANDS` - Distance bands for helix alignment
**Status**: CORE

### `src/helix_analysis.py`
**Purpose**: Helix identification and annotation
**Functions**:
- `define_reference_helices()`: Load helix boundaries from JSON
- `align_to_reference_and_annotate_helices()`: Propagate helix annotations
- `annotate_helices_from_reference()`: Apply helix labels
**Status**: CORE

### `src/assign_grns.py`
**Purpose**: GRN assignment workflow
**Functions**:
- `align_and_assign_grn()`: Main GRN assignment function
- Uses tree-based and reference-based alignment strategies
**Imports**: reference_alignment, tree_based_alignment, msa_grn
**Status**: CORE

### `src/reference_alignment.py`
**Purpose**: Reference-based structure alignment
**Functions**:
- `find_type_references()`: Find reference structures per functional type
- `find_global_reference()`: Select global reference (centroid)
- `create_seq_alignment_dicts_from_paths()`: Build alignment dictionaries
**Status**: CORE

### `src/tree_based_alignment.py`
**Purpose**: Guide tree-based alignment for GRN propagation
**Functions**:
- `build_similarity_tree()`: Build UPGMA tree from RMSD matrix
- `create_guide_tree()`: Create alignment guide tree
- `generate_transitive_alignment_paths()`: Path from reference to all structures
- `align_and_assign_grn_tree_based()`: Tree-based GRN assignment
**Status**: CORE

### `src/msa_grn.py`
**Purpose**: MSA-based GRN utilities
**Functions**:
- `analyze_residue_composition()`: Analyze residue conservation
- `generate_grn_msa_tables()`: Generate MSA from GRN positions
**Status**: CORE

### `src/property_mapping.py`
**Purpose**: Map structure IDs to property metadata
**Class**: `PropertyMapper`
- Handles 4 structure types: exp, pred, hideaki_exp, hideaki_pred
- Maps PDB IDs to functional classifications
**Status**: CORE

### `src/common_utils.py`
**Purpose**: Shared utility functions (avoid circular imports)
**Functions**:
- `compute_retinal_mean_closest_distance()`: Retinal RMSD
- `find_retinal_within_cutoff()`: Find retinal-proximal residues
**Status**: CORE

### `src/lyr_processing.py`
**Purpose**: Handle LYR (Lysine-Retinal) residue splitting
**Functions**:
- `convert_single_lyr_entry()`: Split LYR into LYS + RET
- `process_lyr_in_processor_data()`: Process all LYR residues
**Status**: CORE

### `src/visualization_functions.py`
**Purpose**: Visualization utilities
**Functions**:
- `create_and_visualize_similarity_tree()`: Dendrogram visualization
- `visualize_rmsd_heatmap()`: RMSD matrix heatmap
**Status**: CORE

### `src/visualize_alignment_grn.py`
**Purpose**: Interactive GRN alignment visualization
**Functions**:
- Creates HTML visualizations of GRN alignments
- Plotly-based interactive views
**Status**: CORE

### `src/opsin_color_scheme.py`
**Purpose**: Color schemes for opsin visualization
**Status**: CORE

### `src/retinal_carbon_mapping.py`
**Purpose**: Retinal atom naming/mapping utilities
**Status**: CORE

---

## ANALYSIS SCRIPTS (Review for consolidation)

### `analyze_grns.py` (84KB)
**Purpose**: Comprehensive GRN analysis tool
**Functionality**:
- GRN conservation analysis
- GRN pattern analysis
- Summary figure generation
- Consolidates multiple older scripts
**Status**: REVIEW - Large, may need modularization

### `analyze_motifs.py` (60KB)
**Purpose**: Motif analysis for functional residues
**Functionality**:
- Single position motif analysis (DTD, DTE, etc.)
- Literature-defined functional motifs
- 2-3 residue combination search
- Function-specific motif discovery
**Status**: REVIEW - Large, may need modularization

### `analyze_properties.py` (25KB)
**Purpose**: Property distribution analysis
**Functionality**:
- Calculates diversity metrics
- Generates statistics for manuscript
- Analyzes function/domain distributions
**Status**: REVIEW

### `analyze_validation_errors.py`
**Purpose**: Validation set error analysis for manuscript
**Functionality**:
- Loads Set A/B error CSVs
- Reports statistics (mean, std, min, max)
- Creates Figure 2a visualization
**Status**: REVIEW

---

## UTILITY SCRIPTS

### `scripts/generate_mo_exp_contact_yaml.py`
**Purpose**: Generate YAML configs for Boltz structure prediction
**Functionality**:
- Identifies Schiff base lysine
- Writes pocket constraints for Boltz
**Status**: UTILITY - Keep if using Boltz

### `scripts/update_opsin_reference_years.py`
**Purpose**: Augment property data with RCSB release dates
**Functionality**:
- Queries RCSB API for release dates
- Adds reference_year, dataset_split columns
**Status**: UTILITY - Keep

### `prepare_yaml.py`
**Purpose**: Generate YAML configs from Excel property file
**Functionality**:
- Reads mo_exp.xlsx
- Cleans sequences, handles duplicates
- Writes YAML configs for folding
**Status**: UTILITY - May overlap with scripts/

---

## CLUSTER ANALYSIS (`mo_clusters/`)

### `mo_clusters/prepare_clustered_mo_data.py`
**Purpose**: Process clustered MO structures
**Functionality**:
- Load structures using protos
- Parse mapping and extract sequences
**Status**: REVIEW

### `mo_clusters/align_clustered_mo_step_by_step.py`
**Purpose**: Step-by-step alignment of clustered structures
**Status**: REVIEW

### `mo_clusters/inspect_clustered_structures.py`
**Purpose**: Inspect clustered structure quality
**Status**: REVIEW

---

## VISUALIZATION SCRIPTS

### `visualize_clustered_alignments.py`
**Purpose**: Visualize alignments between clustered MO structures
**Functionality**:
- Uses protos visualization capabilities
- Creates plotly interactive visualizations
**Status**: REVIEW

---

## DEPRECATED/SCRATCH FILES (Move to deprecated/)

### `x.py`
**Content**: Single line printing character from hardcoded sequence
**Status**: DEPRECATED - Scratch file

### `xx.py`
**Content**: Boltz-2 ligand pickle generator from CIF
**Purpose**: Create Boltz-compatible ligand files
**Status**: DEPRECATED - One-off utility, move to scripts/ or deprecated/

### `flatten.py`
**Purpose**: Flatten Boltz output directories, select best model
**Functionality**:
- Scans new_opsins_outputs/ for CIF files
- Selects model with minimum ligand-protein distance
- Copies to flat directory
**Status**: DEPRECATED - One-off data processing

### `src/test_imports.py`
**Status**: DEPRECATED - Test file

### `src/test_protos.py`
**Status**: DEPRECATED - Test file

### `src/structure_alignment_subset.py`
**Status**: REVIEW - May be unused

---

## GENERATED FILES (Add to .gitignore)

### Large HTML files
- `overview_30.html` (8.9MB)
- `overview_45.html` (10.4MB)

### Output directories
- `flat_outputs/`
- `flat_new_opsins_outputs/`
- `new_opsins_outputs/`
- `outputs/`
- `opsin_output/cache/`
- `__pycache__/`

---

## DOCUMENTATION FILES

### Keep
- `README.md` - Project documentation
- `GUIDE.md` - Usage guide (English)
- `GUIDE_JP.md` - Usage guide (Japanese)
- `citations.md` - Literature citations
- `draft.md` / `draft.txt` - Methods manuscript

### Review
- `AGENTS.md` - Agent configuration
- `TODO.md` - Task tracking
- `diversity.md` - Diversity analysis notes
- `literature_review.txt` - Literature notes
- `REPLY` - Unknown

---

## RECOMMENDED ACTIONS

### 1. Move to `deprecated/`
```
x.py
xx.py
flatten.py
src/test_imports.py
src/test_protos.py
```

### 2. Add to `.gitignore`
```
overview_*.html
flat_outputs/
flat_new_opsins_outputs/
new_opsins_outputs/
outputs/
__pycache__/
opsin_output/cache/
*.pyc
```

### 3. Consolidate
- Consider merging `analyze_*.py` scripts or creating an `analysis/` module
- Move `prepare_yaml.py` to `scripts/`

### 4. Refactor
- Update `prepare_data.py` for new protos API
- Update all files importing from protos to use explicit `user_data_root`
