# Opsin Analysis Project Documentation

## Overview

The Opsin Analysis project provides a comprehensive framework for analyzing, comparing, and visualizing experimental and predicted opsin structures with Generic Residue Numbering (GRN). The project focuses on the structural analysis of G-protein coupled receptors (GPCRs), particularly opsins, which are light-sensitive proteins involved in vision and other light-sensing functions.

## Structure Dataset and FoldMason Integration

This project utilizes the CifProcessor's dataset capabilities and FoldMason for structure alignment. The workflow includes:

1. **Find Available Datasets**
   - Examine datasets.json in data/structure_dataset/ to identify available structure datasets
   - Opsin-specific datasets include: mo_ref, mo_exp, mo_exp_refs, and hideaki_exp

2. **Load a Dataset**
   - Use CifProcessor.load_dataset(dataset_name) to load a specific dataset
   - This loads all structures in the dataset from CIF files

3. **Filter the Dataset**
   - Filter structures to select only chain A where appropriate
   - Apply additional filters as needed (by resolution, ligand presence, etc.)

4. **Write Temporary Structures**
   - Create temporary CIF files for the filtered structures
   - These temporary files will be used as input for structure alignment

5. **Create Structure Database**
   - Build a FoldMason database from temporary structure files
   - This database is required for structure-based MSA

6. **Run Multiple Structure Alignment with FoldMason**
   - Use the FoldMason.structuremsa method to perform structure alignment
   - Generate LDDT scores and visualization reports

7. **Extract Aligned Structures and Transformations**
   - Export aligned structures with updated coordinates using the superpose command
   - Extract rotation matrices and translations for all aligned proteins
   - Save transformations in JSON format for further analysis

## Project Structure

- **opsin_analysis_functions.py**: Core computational and data processing functions
- **visualization_functions.py**: Separate module for visualization capabilities
- **opsin_color_scheme.py**: Color schemes for consistent visualization
- **foldmason_helpers.py**: Helper functions for integrating with FoldMason models
- **opsin_analysis_workflow.py**: Complete end-to-end analysis workflow
- **dataset_analysis.py**: Script for analyzing structure datasets using CifProcessor and FoldMason

## Workflow Overview

The `opsin_analysis_workflow.py` implements a complete, modular pipeline with the following steps:

1. **Data Loading & Processing**: Loads experimental and predicted structures, filters and prepares them for analysis
2. **Structure Error Calculation**: Computes errors between experimental and predicted structures
3. **Structure Orientation & Normalization**: Orients structures consistently with the membrane normal along the Z-axis
4. **Helix Annotation**: Identifies and annotates helical segments in the structures
5. **Structure Comparison**: Calculates RMSD between all structure pairs and visualizes the relationships
6. **Generic Residue Numbering (GRN)**: Aligns structures and assigns standardized GRN position labels

## Key Components

### 1. Atom-Level Error Analysis

Functions to calculate atom-level errors and distances between experimental and predicted structures:

- `calculate_atom_level_errors`: Computes detailed error statistics for each atom in a residue
- `summarize_atom_errors`: Summarizes atom-level errors across all proteins and residues
- `compute_retinal_mean_closest_distance`: Finds the best 1:1 pairing between retinal atoms using the Hungarian method
- `find_retinal_within_cutoff`: Extracts retinal atoms within a specific distance of a protein chain
- `make_rmsd_table`: Builds a DataFrame summarizing backbone, pocket, and retinal RMSD values

### 2. Data Assembly & Preprocessing

Functions for loading, filtering, and rebuilding datasets:

- `load_experimental_dataset`: Loads experimental structural data using CifProcessor
- `filter_structures_by_chain`: Extracts structures for a specific chain from a dataset
- `rebuild_hideaki_predicted_dataset`: Loads predicted structures and creates a new dataset
- `extract_pdb_id`: Extracts PDB ID from a model file path
- `collect_cif_files`: Recursively searches for CIF files in a folder

### 3. Structure Orientation & Normalization

Functions to orient and normalize protein structures for consistent analysis:

- `orient_structure`: Orients a structure based on PCA of CA atom coordinates
- `orient_all_structures`: Applies orientation to all structures in a dictionary
- `compute_orientation_vector_pca`: Computes the principal component to determine orientation
- `align_vector_to_z`: Computes a rotation matrix to align a vector to the z-axis
- `orient_structures_n_terminus_up`: Orients structures so the N-terminus is up

### 4. Helix Annotation

Functions for identifying and annotating helices in protein structures:

- `calc_angle_between_vectors`: Calculates the angle between two vectors
- `compute_direction_vector`: Computes a direction vector from the first to last residue
- `identify_helical_like_residues`: Flags residues as helical based on phi/psi angles
- `assign_helix_numbers_with_direction_and_global`: Assigns helix numbers based on direction changes
- `merge_close_helices`: Merges helix segments separated by small gaps
- `annotate_helices_directional_merge`: Detects and merges helical segments
- `assign_tm_helix_flag`: Determines whether helices are transmembrane
- `finalize_helix_assignments`: Finalizes helix assignments based on membrane normal

### 5. Structure Alignment & RMSD Calculation

Functions for aligning structures and calculating RMSD values:

- `compute_all_vs_all_rmsd_improved`: Calculates RMSD between all structure pairs
- `calculate_binding_pocket_rmsd_for_pairs`: Calculates binding pocket RMSD for structure pairs

### 6. 7TM Helix Identification

Functions for identifying and visualizing the 7-transmembrane (7TM) helices:

- `select_correct_ret_cluster`: Selects the correct retinal cluster in a chain
- `select_correct_ret_cluster_full`: Selects the correct retinal cluster in a full structure
- `identify_7tm_helices`: Identifies which helices correspond to the 7TM bundle

### 7. Reference Structure Selection and Alignment

Functions for selecting reference structures and aligning others to them:

- `find_type_references`: Finds a reference structure for each group/type
- `find_global_reference`: Finds a global reference structure from type references
- `apply_alignment`: Applies rotation and translation to coordinates
- `create_sequence_alignment_dict`: Creates a mapping of sequence IDs based on alignment
- `align_structure_to_reference`: Aligns a target structure to a reference
- `update_structures_alignment`: Aligns all structures to appropriate references

### 8. Multiple Structure Alignment & GRN Assignment

Functions for creating and analyzing multiple sequence alignments with GRN:

- `create_msa_table`: Creates an MSA-like table from alignment dictionaries
- `process_alignment_df`: Converts alignment DataFrame to a position weight matrix
- `analyze_residue_composition`: Analyzes residue composition at specific positions
- `sort_grn_columns`: Sorts columns of a GRN-labeled DataFrame in the correct order
- `generate_grn_msa_tables`: Generates all MSA tables with GRN labeling

### 9. Visualization Functions

Separate visualization module with functions for:

- `plot_rmsd_heatmap`: Visualizes RMSD matrix as a heatmap
- `plot_similarity_tree`: Generates a dendrogram based on RMSD values
- `visualize_single_7tm_bundle`: Visualizes the identified 7TM helices for a structure
- `visualize_rmsd_matrix_improved`: Enhanced RMSD matrix visualization with metadata
- `create_and_visualize_similarity_tree`: Creates a similarity tree from RMSD matrix
- `create_sequence_logo`: Creates a sequence logo from an MSA DataFrame
- `print_residue_composition`: Prints the composition of residues at specific positions
- `plot_average_distances_by_helix`: Plots average distances to retinal by helix
- `plot_distance_heatmap`: Creates a heatmap of distances to retinal
- `visualize_binding_pocket`: Visualizes the binding pocket of a structure

## Running the Analysis

The complete analysis can be run using the `opsin_analysis_workflow.py` script:

```python
# Import the workflow
from projects.opsin_analysis.opsin_analysis_workflow import run_opsin_analysis_workflow

# Run the complete analysis with default parameters
results = run_opsin_analysis_workflow(
    output_dir='opsin_analysis_results',
    visualize=True
)
```

The workflow will:
1. Load and process all structures
2. Calculate errors between experimental and predicted structures
3. Orient structures and annotate helices
4. Compare structures and generate RMSD matrices
5. Align structures and assign GRN positions
6. Generate visualizations and statistical analyses
7. Save all results to the specified output directory

## Output Files

The analysis produces several key output files:

- **`rmsd_matrix.csv`**: Pairwise RMSD values between all structures
- **`rmsd_matrix.png`**: Visualization of the RMSD matrix as a heatmap
- **`msa_table_grn.csv`**: Multiple sequence alignment with GRN positions
- **`distance_table_grn.csv`**: Distances from each residue to retinal with GRN positions
- **`distances_by_helix.png`**: Plot of average distances to retinal by helix
- **`distance_heatmap.png`**: Heatmap showing residue distances to retinal
- **`*_7tm_bundle.html`**: Interactive 3D visualizations of the 7TM bundle for each structure

## Integration with Other Components

This project integrates with several other components in the PROTOS package:

- Uses `struct_processor.py` for loading and manipulating PDB/CIF files
- Leverages `struct_alignment.py` for structural alignment algorithms
- Connects with GRN assignment modules for residue numbering
- Employs visualization tools for structural analysis display
- Integrates with PropertyProcessor for metadata handling