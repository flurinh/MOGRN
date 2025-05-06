# Opsin Analysis Modules Documentation

This document provides an overview of the different modules used in the opsin analysis workflow. Each module contains specialized functions for handling different aspects of the analysis pipeline.

## Core Workflow

The main workflow is orchestrated by **opsin_analysis_workflow.py**, which imports functionality from several specialized modules to perform a comprehensive analysis of opsin structures.

## Module Functions

### data_processing.py
Contains functions for loading, filtering, and preprocessing opsin structure data.

- **load_opsin_structures**: Loads opsin structures from prepared datasets with two-stage caching support.
- **filter_structures_by_chain_and_retinal**: Filters structures by chain and only includes the retinal closest to that chain.
- **ensure_structure_dtypes**: Ensures all structure DataFrames have proper data types for coordinates and other fields.
- **load_opsin_property_data**: Loads property data from CSV file and creates structure mappings.

### structure_comparison.py
Functions for comparing and aligning structures, calculating RMSD between structure pairs.

- **compute_all_vs_all_rmsd_improved**: Calculates RMSD between all pairs of structures using C-alpha atoms.
- **calculate_binding_pocket_rmsd_for_pairs**: Calculates RMSD for binding pocket residues between pairs of structures.
- **compare_structures**: Creates similarity matrices and visualizations based on structure RMSD.
- **create_unified_structure_mapping**: Creates a single mapping between experimental and predicted structures.

### helix_analysis.py
Functions for identifying, annotating, and analyzing helices in protein structures.

- **define_reference_helices**: Defines helix boundaries for reference structures.
- **annotate_helices_from_alignments**: Annotates helices in all structures based on reference alignments.
- **align_to_reference_and_annotate_helices**: Aligns structures to reference and annotates helices.
- **orient_and_annotate_structures**: Performs orientation and annotation in a single step.

### error_analysis.py
Functions for analyzing errors between experimental and predicted structures.

- **calculate_structure_errors**: Calculates RMSD errors between paired structures.
- **make_rmsd_table**: Builds a summary DataFrame of RMSD values.
- **compute_retinal_mean_closest_distance**: Calculates the mean distance between retinal atoms.
- **calculate_atom_level_errors**: Provides detailed error statistics for individual atoms.

### reference_alignment.py
Functions for selecting reference structures and aligning other structures to them.

- **find_type_references**: Identifies reference structures for each group/type.
- **find_global_reference**: Selects a global reference from type references.
- **create_seq_alignment_dicts_from_paths**: Converts alignment paths to sequence alignment dictionaries.

### msa_grn.py
Functions for multiple structure alignment and Generic Residue Numbering (GRN).

- **generate_grn_msa_tables**: Creates MSA tables with proper GRN column names.
- **create_msa_table**: Builds MSA-like tables from alignment dictionaries.
- **create_msa_distance_table**: Creates distance tables showing closest distances to retinal.
- **calculate_helix_distances**: Calculates mean distances to retinal for each helix.
- **analyze_residue_composition**: Analyzes amino acid composition at specific positions.

### assign_grns.py
Functions for assigning Generic Residue Numbers to structures.

- **align_and_assign_grn**: Handles structure alignment and GRN assignment.

### visualization_functions.py
Functions for visualizing structure alignments, distances, and conservation.

- **visualize_single_7tm_bundle**: Creates 3D visualization of a 7TM bundle.
- **visualize_msa_distances**: Visualizes multiple sequence alignment distances.
- **plot_average_distances_by_helix**: Creates plots of average distances by helix.
- **plot_distance_heatmap**: Generates heatmaps of distances.
- **create_residue_conservation_plot**: Creates plots of residue conservation.
- **create_sequence_logo**: Generates sequence logo plots of aligned positions.

### foldmason_helpers.py
Functions for interfacing with FoldMason structure alignment tools.

- **run_easy_msa**: Runs FoldMason's easy-msa command.
- **run_structuremsa**: Runs FoldMason's structuremsa command.
- **run_refinemsa**: Runs FoldMason's refinemsa command for iterative refinement.
- **align_with_foldmason**: Uses FoldMason for multiple structure alignment.

## Workflow Steps

The workflow in opsin_analysis_workflow.py follows these steps:

1. **Load Structures**: Loads opsin structures from prepared datasets.
2. **Calculate Errors**: Calculates RMSD errors between experimental and predicted structures.
3. **Annotate Helices**: Orients structures and annotates helices using alignments.
4. **Structure Comparison**: Computes RMSD matrices and visualizes structure similarities.
5. **Assign GRNs**: Assigns Generic Residue Numbers to structures using multiple sequence alignment.

## Data Flow

1. Raw structures are loaded from datasets using `load_opsin_structures`
2. Structures are filtered using `filter_structures_by_chain_and_retinal`
3. Structure errors are calculated with `calculate_structure_errors`
4. Helices are annotated using `align_to_reference_and_annotate_helices`
5. Structure comparisons are performed with `compare_structures`
6. Generic Residue Numbers are assigned using `align_and_assign_grn`

## Circular Import Resolution

To resolve circular imports between modules, ensure that:

1. Each module only imports functions from other modules that it directly depends on
2. If two modules need to import functions from each other, move the shared functions to a common module
3. Use delayed imports within functions when necessary to break import cycles
4. Consider reorganizing functionality to minimize cross-module dependencies

## Module Dependencies

- **data_processing.py**: Minimal dependencies, imported by most other modules
- **error_analysis.py**: Depends on structure_comparison.py
- **helix_analysis.py**: Depends on data_processing.py, visualization_functions.py
- **structure_comparison.py**: Depends on error_analysis.py, visualization_functions.py
- **msa_grn.py**: Depends on visualization_functions.py
- **reference_alignment.py**: No direct module dependencies
- **assign_grns.py**: Depends on reference_alignment.py, msa_grn.py, visualization_functions.py
- **visualization_functions.py**: Minimal dependencies
- **foldmason_helpers.py**: No direct module dependencies