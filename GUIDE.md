# MOGRN Project Guide

This guide explains how to run the microbial opsin analysis workflow and what each script does.

## Project Overview

The Opsin Analysis project provides a comprehensive framework for analyzing, comparing, and visualizing experimental and predicted opsin structures using Generic Residue Numbering (GRN). The project focuses on the structural analysis of microbial opsins, which are light-sensitive proteins involved in various sensing functions.

The analysis pipeline integrates multiple components:
- Structure data loading and preprocessing from various sources
- Structural alignment and error calculation
- Transmembrane helix identification and annotation
- Comparative analysis across multiple structures
- Generic Residue Number assignment for standardized comparison
- Comprehensive visualization tools for structural insights

## Module Architecture

The workflow is divided into specialized modules:

- **data_processing.py**: Loads, filters, and preprocesses opsin structures
- **structure_comparison.py**: Compares and aligns structures, calculates RMSD
- **helix_analysis.py**: Identifies and annotates transmembrane helices
- **error_analysis.py**: Analyzes errors between experimental and predicted structures
- **msa_grn.py**: Handles multiple sequence alignment and GRN assignment
- **visualization_functions.py**: Creates visualizations of alignments and analysis results
- **foldmason_helpers.py**: Interfaces with FoldMason structure alignment tools
- **opsin_analysis_workflow.py**: Orchestrates the complete analysis pipeline

## Step 1: Prepare Data Infrastructure and Initialize Datasets

**Command to run:** 
```bash
# Set up data structure and process datasets
python prepare_data_fixed.py
```

**What this step does:**
- Reviews and organizes property data for microbial opsins stored in CSV files
- Sets up the necessary directory structure for analysis
- Initializes the ProtosPaths system from the Protos framework to manage file paths
- Creates data directories for user-writable and reference data
- Validates the directory structure and reports any issues
- Sets up the CifBaseProcessor for managing protein structures
- Creates registry files and initializes dataset management
- Processes four types of datasets:
  1. Experimental microbial opsin structures (downloaded from PDB)
  2. Predicted microbial opsin structures (from local files)
  3. Experimental structures from Hideaki dataset
  4. Predicted structures from Hideaki dataset
- Validates datasets to ensure all structures are accessible
- Sets up caching for efficient data retrieval

## Step 2: Create YAML Configuration Files for Protein Sequences

**Command to run:**
```bash
# Generate YAML files from protein sequences
python prepare_yaml.py
```

**What this step does:**
- Processes protein sequence data from CSV files
- Creates YAML files for each protein with proper sequence formatting
- Extracts sequence segments using start/end positions when available
- Reports on sequences longer than 400 amino acids without range information
- Creates special configurations for proteins with multiple retinals
- Validates amino acid sequences and cleans invalid characters
- Generates statistics and visualizations of sequence length distributions

## Step 3: Run Analysis Workflow

**Command to run:**
```bash
# Run the complete analysis pipeline
python opsin_analysis_workflow.py
```

**What this step does:**
- Loads all structures from the datasets created in Step 1
- Performs structural alignment and standardization
- Calculates structural errors between experimental and predicted structures
- Identifies and annotates transmembrane helical regions
- Compares structures to analyze similarities and differences
- Creates a unified structure mapping for cross-comparison
- Assigns Generic Residue Numbers (GRNs) to standardize residue positions
- Caches results at each stage for faster reuse

The script uses default parameters that are already configured appropriately for the standard workflow.

The workflow follows these key processing steps:
1. **Structure Loading**: Loads opsin structures from prepared datasets with caching support
2. **Error Calculation**: Computes RMSD errors between paired experimental and predicted structures
3. **Helix Annotation**: Identifies and labels helical segments with transmembrane designation
4. **Structure Comparison**: Creates RMSD matrices and visualizes structural relationships
5. **GRN Assignment**: Aligns structures and assigns standard numbering based on evolutionary conservation

## Step 4: Generate Visualizations

**Command to run:**
```bash
# Create standardized visualizations and plots
python plot_figures.py
```

**What this step does:**
- Loads precomputed data from the analysis workflow
- Generates multiple visualization types:
  - RMSD heatmaps showing structural similarities
  - Similarity trees clustering structures by similarity
  - Distance plots showing proximity to retinal
  - Conservation plots showing amino acid conservation
  - Helix logo plots displaying sequence patterns
  - Overview plots summarizing key findings
- Filters out low-quality structures
- Creates protein summary CSV with key metrics
- Saves all visualizations to a 'figures' subdirectory

The script automatically detects and uses the data created by the analysis workflow in the default output directory.

The visualizations include:
- Structure similarity heatmaps and dendrograms
- Residue distance plots showing proximity to retinal
- Helix annotation visualizations
- Residue conservation analysis
- Binding pocket comparisons
- Statistical summaries of structural features

## Output Files

The analysis produces several key output files:

- **`rmsd_matrix.csv`**: Pairwise RMSD values between all structures
- **`msa_table_grn.csv`**: Multiple sequence alignment with GRN positions
- **`ca_distance_table_grn.csv`**: C-alpha distances from each residue to retinal
- **`distance_table_grn.csv`**: Sidechain distances to retinal with GRN positions
- **`protein_summary.csv`**: Summary of analyzed proteins with key metrics

## Optional Advanced Usage

If you need to customize the workflow beyond the default settings, the scripts support various command-line options. Here are some examples:

### Custom Directories
```bash
# Specify a custom output directory
python opsin_analysis_workflow.py --output-dir custom_output

# Specify input and output directories for visualizations
python plot_figures.py --input-dir custom_output --output-dir custom_figures
```

### Caching Control
```bash
# Run without using cached results (recompute everything)
python opsin_analysis_workflow.py --no-cache
```

### Visualization Control
```bash
# Run without generating visualizations during analysis
python opsin_analysis_workflow.py --no-visualize

# Generate higher quality figures (takes longer)
python plot_figures.py --quality high
```

These advanced options should only be used when you need to override the default behavior.