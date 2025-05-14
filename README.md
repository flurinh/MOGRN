# MOGRN: Microbial Opsin Generic Residue Numbering

## Project Overview

MOGRN is a comprehensive framework for analyzing, comparing, and visualizing experimental and predicted microbial opsin structures using Generic Residue Numbering (GRN). Microbial opsins are light-sensitive proteins involved in various sensing functions, and this framework provides tools for standardized structural analysis.

The analysis pipeline integrates multiple components:
- Structure data loading and preprocessing from various sources
- Structural alignment and error calculation
- Transmembrane helix identification and annotation
- Comparative analysis across multiple structures
- Generic Residue Number assignment for standardized comparison
- Comprehensive visualization tools for structural insights

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/MOGRN.git
cd MOGRN

# Install dependencies
pip install -r requirements.txt
```

### Installing Protos Framework (Required)

MOGRN depends on the Protos framework for protein structure analysis. To install it:

```bash
# Clone the Protos repository
git clone https://github.com/flurinh/protos.git
cd protos

# Install Protos in development mode
pip install -e .

# Return to the MOGRN directory
cd ..
```

## Required Data Folders (Not Included)

**IMPORTANT**: The repository does not include data files. You must create and populate these folders:

1. **`property/`** - Contains property data for microbial opsins:
   - `mo_exp.csv`: Main experimental data properties
   - `helices.json`: Helix definitions
   - Other supporting files

2. **`structures/`** - Contains structure files organized in subdirectories:
   - `hideaki_exp/`: Experimental structures from Hideaki dataset (CIF format)
   - `hideaki_pred/`: Predicted structures from Hideaki dataset (CIF format)
   - `mo_pred/`: Predicted microbial opsin structures (CIF format)

Directory structure example:
```
MOGRN/
├── property/
│   ├── mo_exp.csv
│   ├── helices.json
│   └── ...
├── structures/
│   ├── hideaki_exp/
│   │   └── [structure files (.cif)]
│   ├── hideaki_pred/
│   │   └── [structure files (.cif)]
│   └── mo_pred/
│       └── [structure files (.cif)]
└── ...
```

## Project Structure

The project is organized as follows:
- **Main workflow scripts** are in the root directory:
  - `prepare_data_fixed.py`: Data preparation and initialization
  - `prepare_yaml.py`: Configuration generation for sequences
  - `opsin_analysis_workflow.py`: Main analysis pipeline
  - `plot.py`: Visualization generation
- **Helper modules** are in the `src/` directory:
  - Processing utilities
  - Analysis tools
  - Visualization functions

## Workflow Steps

### Step 1: Prepare Data Infrastructure and Initialize Datasets

```bash
python prepare_data_fixed.py
```

This step:
- Reviews and organizes property data for microbial opsins
- Sets up the necessary directory structure for analysis
- Initializes the ProtosPaths system for file path management
- Creates data directories for user-writable and reference data
- Sets up the CifBaseProcessor for managing protein structures
- Processes experimental and predicted microbial opsin structures
- Validates datasets to ensure all structures are accessible
- Sets up caching for efficient data retrieval

### Step 2: Create YAML Configuration Files for Protein Sequences

```bash
python prepare_yaml.py
```

This step:
- Processes protein sequence data from CSV files
- Creates YAML files for each protein with proper sequence formatting
- Extracts sequence segments using start/end positions when available
- Creates special configurations for proteins with multiple retinals
- Validates amino acid sequences and cleans invalid characters

### Step 3: Run Analysis Workflow

```bash
python opsin_analysis_workflow.py
```

This step performs the complete analysis pipeline:
1. **Structure Loading**: Loads opsin structures from prepared datasets with caching
2. **Error Calculation**: Computes RMSD errors between paired structures
3. **Helix Annotation**: Identifies and labels helical segments 
4. **Structure Comparison**: Creates RMSD matrices and visualizes relationships
5. **GRN Assignment**: Aligns structures and assigns standard numbering

### Step 4: Generate Visualizations

```bash
python plot.py
```

This step creates visualizations from the analysis results:
- RMSD heatmaps showing structural similarities
- Similarity trees clustering structures by similarity
- Distance plots showing proximity to retinal
- Conservation plots showing amino acid conservation
- Helix logo plots displaying sequence patterns
- Overview plots summarizing key findings

## Output Files

The analysis produces several key output files:
- **rmsd_matrix.csv**: Pairwise RMSD values between all structures
- **msa_table_grn.csv**: Multiple sequence alignment with GRN positions
- **ca_distance_table_grn.csv**: C-alpha distances from each residue to retinal
- **distance_table_grn.csv**: Sidechain distances to retinal with GRN positions
- **protein_summary.csv**: Summary of analyzed proteins with key metrics

## Module Architecture

The workflow is divided into specialized modules:
- **data_processing.py**: Loads, filters, and preprocesses opsin structures
- **structure_comparison.py**: Compares and aligns structures, calculates RMSD
- **helix_analysis.py**: Identifies and annotates transmembrane helices
- **error_analysis.py**: Analyzes errors between experimental and predicted structures
- **msa_grn.py**: Handles multiple sequence alignment and GRN assignment
- **assign_grns.py**: Assigns Generic Residue Numbers to standardize residue positions
- **visualization_functions.py**: Creates visualizations of alignments and analysis results
- **foldmason_helpers.py**: Interfaces with FoldMason structure alignment tools
- **opsin_analysis_workflow.py**: Orchestrates the complete analysis pipeline

## Advanced Usage

For customized workflow options:

```bash
# Run with custom output directory
python opsin_analysis_workflow.py --output-dir custom_output

# Run without using cached results
python opsin_analysis_workflow.py --no-cache

# Generate visualizations with specific input and output directories
python plot.py --input-dir custom_output --output-dir custom_figures
```

## Dependencies

This project depends on several key packages:
- **Protos**: Framework for protein structure analysis (https://github.com/flurinh/protos.git)
- **BioPython**: For sequence and structure manipulation
- **NumPy/Pandas**: For data analysis and manipulation
- **Matplotlib/Seaborn**: For visualization
- **FoldMason**: For structure alignment (included in Protos)

## Detailed Documentation

For a complete explanation of the project's workflow and how to run each step, see the `GUIDE.md` file.