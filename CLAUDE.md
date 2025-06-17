# CLAUDE.md - MOGRN Project Reference

## Project Overview

**MOGRN (Microbial Opsin Generic Residue Numbering)** is a comprehensive framework for analyzing, comparing, and visualizing experimental and predicted microbial opsin structures using Generic Residue Numbering (GRN). This project provides standardized structural analysis tools for understanding light-sensitive proteins involved in various sensing functions.

### Key Features

- **Structure Processing**: Loads and processes opsin structures from multiple sources (experimental/predicted)
- **Generic Residue Numbering**: Assigns standardized residue positions across different structures
- **Structural Alignment**: Uses advanced alignment algorithms (CEalign, structure-based) for comparison
- **Transmembrane Helix Analysis**: Identifies and annotates 7TM helical regions
- **Distance Analysis**: Calculates distances to retinal chromophore for binding pocket analysis
- **Visualization**: Generates comprehensive plots including RMSD matrices, sequence logos, and distance plots
- **Caching System**: Multi-level caching for efficient data reuse and faster processing

### Project Architecture

The project follows a modular pipeline design with 6 main processing steps:

1. **Data Loading & Structure Processing** (`data_processing.py`)
2. **Error Analysis** (`error_analysis.py`) 
3. **Helix Annotation** (`helix_analysis.py`)
4. **Structure Comparison** (`structure_comparison.py`)
5. **GRN Assignment** (`assign_grns.py`, `msa_grn.py`)
6. **Visualization** (`plot.py`, `visualization_functions.py`)

## Workflow Execution

### Main Scripts

- **`prepare_data_fixed.py`**: Initialize datasets and infrastructure
- **`prepare_yaml.py`**: Generate configuration files for protein sequences  
- **`opsin_analysis_workflow.py`**: Main analysis pipeline orchestrator
- **`plot.py`**: Generate all visualizations from analysis results

### Quick Start

```bash
# 1. Set up data structure and datasets
python prepare_data_fixed.py

# 2. Generate sequence configurations
python prepare_yaml.py

# 3. Run complete analysis pipeline
python opsin_analysis_workflow.py

# 4. Generate visualizations
python plot.py
```

## Implementation Details

### Core Data Structures

**Processed Structures Dictionary**:
```python
processed_structures = {
    'pdb_id': {
        'df': pd.DataFrame,           # Full structure data
        'df_norm': pd.DataFrame,      # Normalized coordinates
        'df_ca_norm': pd.DataFrame,   # CA atoms only
        'df_ret': pd.DataFrame,       # Retinal atoms
        'chain_id': str,              # Chain identifier
        'structure_type': str,        # 'experimental' or 'predicted'
        'helix_definitions': dict,    # Helix boundary definitions
        'helix_assignments': dict,    # Residue to helix mapping
        'properties': dict            # Metadata and properties
    }
}
```

**Structure Mapping**:
```python
structure_mapping = {
    'experimental_id': 'predicted_id',  # Maps exp to pred structures
    # Used throughout workflow for paired analysis
}
```

### Key Processing Functions

#### 1. Data Loading (`data_processing.py`)

**`load_opsin_structures()`**:
- Two-stage caching system (raw + processed)
- Handles 4 datasets: mo_exp, mo_pred, hideaki_exp, hideaki_pred
- Chain filtering and retinal selection
- LYR residue processing (splits into LYS + RET)

**`filter_structures_by_chain_and_retinal()`**:
- Extracts specific chain and closest retinal
- Handles distance-based retinal selection
- Manages duplicate atom removal

#### 2. Structure Comparison (`structure_comparison.py`)

**`compute_all_vs_all_rmsd_improved()`**:
- Uses CEalign for structural alignment
- CA atom-based RMSD calculation
- Helix-only alignment option (TM regions 1-7)
- Comprehensive caching with hash-based keys
- Returns alignment paths for GRN mapping

**`calculate_binding_pocket_rmsd_for_pairs()`**:
- Compares experimental vs predicted binding pockets
- Uses alignment paths for residue mapping
- Separate retinal RMSD calculation
- Per-residue error analysis

#### 3. Helix Analysis (`helix_analysis.py`)

**`align_to_reference_and_annotate_helices()`**:
- Uses reference structure with predefined helix boundaries
- Transfers helix annotations via structural alignment
- Caches helix definitions in `property/helices.json`
- Validates and repairs corrupted cache files

**Helix Definition Format**:
```python
helix_definitions = {
    'structure_id': {
        '1': [start_pos, end_pos],  # Helix 1 boundaries
        '2': [start_pos, end_pos],  # Helix 2 boundaries
        # ... up to helix 7
    }
}
```

#### 4. GRN Assignment (`assign_grns.py`, `msa_grn.py`)

**`align_and_assign_grn()`**:
- Uses cached alignment paths from structure comparison
- Implements both reference-based and tree-based alignment
- Filters structures by RMSD quality (removes <0.1Å outliers)
- Generates MSA tables with GRN positions

**`generate_grn_msa_tables()`**:
- Creates comprehensive MSA with residue and distance information
- Supports both CA-only and all-atom analysis
- Filters by RMSD threshold and prioritizes experimental structures
- Outputs GRN-labeled position tables

#### 5. LYR Processing (`lyr_processing.py`)

**Special Handling for Lysine-Retinal Schiff Base**:
- Automatically detects LYR residues in structures
- Splits into separate LYS (protein) and RET (retinal) components
- Maintains all original atom properties and coordinates
- Processes at multiple levels: DataFrame, structures dict, processor

### Visualization System (`visualization_functions.py`, `plot.py`)

#### Generated Plots

1. **Opsin Overview Plot** (`01_opsin_overview.png`)
   - Structure counts by type and domain
   - Experimental vs predicted distributions

2. **RMSD Clustermap** (`02_rmsd_clustermap.png`)
   - Hierarchical clustering of structures
   - Color-coded by functional groups

3. **Distance to Retinal Plots** (`03_all_atom_distance_std.png`, `04_ca_atom_distance_std.png`)
   - Average distances with standard deviations
   - GRN-labeled positions
   - Identifies binding pocket residues

4. **Helix Logo Plots** (`05_helix_logos_x50.png`)
   - Sequence conservation within each TM helix
   - Scaled by conservation score

#### Visualization Architecture

**Color Schemes**:
- Consistent opsin family color mapping
- Functional group classification (ion pumps, channelrhodopsins, etc.)
- Helix-specific color coding

**Plot Styling**:
- Standardized figure sizes and DPI settings
- Consistent font families and sizes
- Professional publication-ready formatting

### File Organization

#### Directory Structure
```
MOGRN/
├── src/                          # Core implementation modules
│   ├── data_processing.py        # Structure loading and filtering
│   ├── structure_comparison.py   # RMSD calculation and alignment
│   ├── helix_analysis.py         # Helix annotation and transfer
│   ├── assign_grns.py            # GRN assignment logic
│   ├── msa_grn.py               # MSA table generation
│   ├── error_analysis.py         # Error calculation utilities
│   ├── lyr_processing.py         # LYR residue handling
│   ├── visualization_functions.py # Plot generation functions
│   ├── common_utils.py           # Shared utility functions
│   └── opsin_color_scheme.py     # Color mapping definitions
├── property/                     # Reference data and configurations
│   ├── mo_exp.csv               # Microbial opsin properties
│   ├── helices.json             # Cached helix definitions
│   └── helix_ref_*.json         # Reference helix boundaries
├── structures/                   # Structure files (not in repo)
│   ├── hideaki_exp/             # Experimental structures
│   ├── hideaki_pred/            # Predicted structures
│   └── mo_pred/                 # MO predicted structures
├── opsin_output/                 # Analysis results and cache
│   ├── cache/                   # Cached intermediate results
│   ├── figures/                 # Generated visualizations
│   └── *.csv                    # Output tables and matrices
└── protos/                       # Protos framework (dependency)
```

#### Output Files

**Key Result Files**:
- `rmsd_matrix.csv`: Pairwise RMSD values between structures
- `msa_table_grn.csv`: Multiple sequence alignment with GRN positions
- `distance_table_grn.csv`: Distances from each residue to retinal
- `ca_distance_table_grn.csv`: CA-specific distance measurements
- `analysis_summary.json`: Workflow execution summary

### Advanced Features

#### Caching System

**Multi-Level Caching**:
1. **Raw Structure Cache**: Stores unfiltered structure data (chain-independent)
2. **Processed Structure Cache**: Stores chain-specific filtered data
3. **Component Cache**: Individual workflow step results
4. **RMSD Cache**: Hash-based caching of expensive alignment calculations

**Cache Management**:
- Automatic cache validation and corruption detection
- Backup creation for critical cache files
- Atomic file operations to prevent corruption
- Hash-based keys for parameter-dependent caching

#### Error Handling

**Robust Error Recovery**:
- Graceful degradation when structures are missing
- Automatic fallback for failed alignments
- Comprehensive logging of processing issues
- Validation of input data integrity

**Quality Control**:
- RMSD outlier detection and filtering
- Structure completeness validation
- Coordinate data type enforcement
- Missing data imputation strategies

## Protos Framework Integration

### ProtosPaths System

ProtosPaths handles path management throughout the Protos framework:

```python
# Key classes in protos.io.paths
from protos.io.paths.path_config import ProtosPaths, DataSource

# Initialization
paths = ProtosPaths(
    user_data_root=None,  # Uses env var or defaults to 'data'
    ref_data_root=None,   # Uses package resources
    create_dirs=True,     # Creates directories if missing
    validate=True         # Validates directory structure
)

# Core functions
structure_path = paths.get_processor_path("structure")
abs_path = paths.resolve_path("dataset.json", source=DataSource.USER)
exists, source = paths.exists("structure/1abc.cif")
```

### DatasetManager

DatasetManager provides standardized dataset operations:

```python
# In protos.core.dataset_manager
manager = DatasetManager(
    processor_type="structure",  # Processor type
    paths=paths                  # Path resolver instance
)

# Dataset operations
dataset = manager.create_dataset(
    dataset_id="my_dataset",
    name="My Structure Dataset", 
    description="A collection of structure IDs",
    content=["1abc", "2xyz", "3def"],
    metadata={"source": "PDB"}
)
```

### CifBaseProcessor

Foundation for structure processing:

```python
# Initialization pattern for processors
processor = CifBaseProcessor(
    name="my_processor",     # Identifier
    data_root="~/data",      # Custom data location
    processor_data_dir="structure"  # Directory name
)

# Common methods
processor.load_dataset("dataset_name")
processor.load_structures(pdb_ids, apply_dtypes=True)
processor.save_dataset("output_dataset")
```

## Important Notes and Maintenance Requirements

### Critical Dependencies

**This project requires the Protos framework to be installed in development mode**:
```bash
cd protos
pip install -e .
```

**Required data directories** (not included in repository):
- `property/mo_exp.csv`: Experimental properties and mappings
- `structures/`: CIF structure files organized by type
- `property/helices.json`: Cached helix definitions (auto-generated)

### Code Maintenance Guidelines

#### When making substantial code changes or additions, this CLAUDE.md file MUST be updated to reflect:

1. **New Functions or Modules**: Add descriptions and usage examples
2. **Modified Data Structures**: Update structure definitions and schemas
3. **Changed File Formats**: Update input/output file documentation
4. **New Dependencies**: Add installation and configuration instructions
5. **Modified Workflow Steps**: Update process descriptions and sequence
6. **New Visualization Types**: Document plot types and interpretation
7. **Changed Cache Structures**: Update caching mechanism descriptions
8. **Modified Configuration Options**: Update parameter documentation

#### Update Requirements:

- **Function Signatures**: Document all new public functions with parameters and return values
- **Data Flow**: Update architecture diagrams and data structure documentation  
- **File Locations**: Update directory structure and file organization descriptions
- **Error Handling**: Document new error conditions and recovery mechanisms
- **Performance Considerations**: Note computational complexity and memory requirements

#### Documentation Standards:

- Use clear, descriptive section headers
- Include code examples for complex operations
- Maintain consistent formatting and style
- Cross-reference related sections
- Keep implementation details current with code

**Failure to maintain this documentation will result in knowledge loss and reduced project maintainability.**

### Performance Notes

**Memory Usage**:
- Large structures may require 8-16GB RAM for full workflow
- RMSD calculations are memory-intensive for >100 structures
- Use caching to avoid recomputation

**Computational Complexity**:
- RMSD calculation: O(n²) for n structures
- Structure alignment: O(n*m) for sequences of length n,m
- Helix annotation: O(n) per structure

**Optimization Strategies**:
- Use helix-only alignment for faster computation
- Enable all caching mechanisms
- Process structures in batches for large datasets
- Use chain filtering to reduce data volume

### Common Issues and Solutions

**Structure Loading Problems**:
- Ensure CIF files are properly formatted
- Check chain identifiers match expected values
- Verify retinal naming conventions (RET vs LIG vs LYR)

**Alignment Failures**:
- Structures may be too dissimilar for automatic alignment
- Missing CA atoms will cause alignment errors
- Very short sequences may fail alignment criteria

**Visualization Issues**:
- Large datasets may require plot parameter adjustment
- Memory limitations may require downsampling
- Color schemes may need adjustment for new structure types

### Development Guidelines

**Code Style**:
- Follow existing naming conventions
- Add comprehensive docstrings for all functions
- Include error handling for all file operations
- Use type hints where possible

**Testing**:
- Test with known good structure sets
- Validate output file formats
- Check visualization rendering on different systems
- Verify cache integrity and recovery mechanisms

## Project Documentation References

- **GUIDE.md**: Complete workflow explanation and step-by-step execution instructions
- **README.md**: Installation instructions and project overview
- **TODO.md**: Current development tasks and known issues