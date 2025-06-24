# MOGRN: Microbial Opsin Generic Residue Numbering

## Project Overview

MOGRN is a comprehensive framework for analyzing, comparing, and visualizing experimental and predicted microbial opsin structures using Generic Residue Numbering (GRN). Microbial opsins are light-sensitive proteins that function as ion pumps, channels, and sensors across diverse organisms. This framework provides advanced tools for understanding their structural conservation, functional determinants, and evolutionary relationships.

### Key Features

- **Generic Residue Numbering (GRN)**: Standardized numbering system enabling consistent comparison across diverse opsin structures
- **Conservation Analysis**: Identifies functionally critical positions and domain-specific patterns
- **Motif Detection**: Validates literature-known motifs and discovers novel functional patterns
- **Functional Discrimination**: Distinguishes structural features that determine pump vs channel function
- **Co-evolution Analysis**: Detects correlated residue changes indicating functional coupling

### Scientific Background

Microbial opsins contain seven transmembrane helices with a retinal chromophore. Key positions (labeled with .50 suffix in GRN) serve as structural anchors:
- **Position 3.50**: Critical functional switch (T in pumps → C in channels)
- **Position 6.50**: Conserved tryptophan forming the retinal binding pocket
- **Position 7.50**: Lysine forming Schiff base with retinal
- Other .50 positions show function-specific conservation patterns

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/MOGRN.git
cd MOGRN

# Create conda environment (recommended)
conda create -n mogrn python=3.10
conda activate mogrn

# Install dependencies
pip install -r requirements.txt
```

### Installing Protos Framework (Required)

MOGRN depends on the Protos framework for protein structure analysis:

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

**IMPORTANT**: Create and populate these folders before running the analysis:

1. **`property/`** - Property data files:
   - `mo_exp.csv`: Experimental opsin properties (required)
   - `helices_curated.json`: Helix boundary definitions (required)

2. **`structures/`** - Structure files in CIF format:
   - `hideaki_exp/`: Experimental structures from Hideaki dataset
   - `hideaki_pred/`: AlphaFold predictions for Hideaki dataset
   - `mo_pred/`: AlphaFold predictions for microbial opsins

## Analysis Pipeline

The complete analysis involves six main steps:

### Step 1: Data Preparation
```bash
python prepare_data.py
```
Initializes the data infrastructure, validates structure files, and sets up caching.

### Step 2: Sequence Configuration
```bash
python prepare_yaml.py
```
Generates YAML configuration files for protein sequences with proper formatting.

### Step 3: Core Analysis Workflow
```bash
python opsin_analysis_workflow.py
```
Performs structural alignment, RMSD calculations, helix identification, and GRN assignment.

### Step 4: Manuscript Visualizations
```bash
python plot.py
```
Generates RMSD heatmaps, similarity trees, distance plots, and conservation visualizations.

### Step 5: GRN Conservation Analysis
```bash
python analyze_grns.py
```
**New comprehensive analysis** including:
- Conservation patterns at functionally critical positions
- Functional group discrimination (pumps vs channels)
- Co-evolution network analysis
- Domain-specific conservation patterns
- Multi-panel summary figures

### Step 6: Motif Pattern Analysis
```bash
python analyze_motifs.py
```
**Advanced motif analysis** featuring:
- Single position conservation analysis
- Literature motif validation (DTD, etc.)
- Cross-functional pattern detection
- Novel motif family discovery
- Correlation analysis within functional groups

## Output Structure

```
MOGRN/
├── output/                          # Core analysis results
│   ├── rmsd_matrix.csv             # Pairwise structural similarities
│   ├── msa_table_grn.csv           # Aligned sequences with GRN
│   ├── distance_table_grn.csv      # Residue-retinal distances
│   └── protein_summary.csv         # Structure quality metrics
├── opsin_output/
│   ├── grn_analysis/               # GRN conservation analysis
│   │   ├── conservation/           # Conservation data files
│   │   ├── figures/                # Multi-panel visualizations
│   │   └── reports/                # Detailed text reports
│   └── motifs/                     # Motif analysis results
│       ├── single_position_analysis/
│       ├── literature_motifs/
│       └── correlation_analysis/
└── figures/                        # Standard visualization plots
```

## Key Results and Interpretation

### Conservation Analysis Results

The GRN conservation analysis reveals:
- **Universal Conservation**: K7.50 (Schiff base) and W6.50 (retinal pocket)
- **Functional Switches**: Position 3.50 discriminates pumps (T) from channels (C)
- **Domain Patterns**: Archaea show highest conservation; Eukaryotes most diverse

### Motif Analysis Results

Key functional motifs identified:

### Functional Discrimination

The analysis identifies positions that determine functional properties:
```
Pumps:    T3.50 + S5.50 + M4.50 → Tight structure, controlled transport
Channels: C3.50 + G5.50 + T4.50 → Flexible gating, ion selectivity
```

### Visualization Options

```bash
# Custom figure directories
python plot.py --input-dir custom_output --output-dir custom_figures

# Specific plot types
python plot.py --plots "rmsd,conservation,distance"

# High DPI output
python plot.py --dpi 600
```

## Module Architecture

### Core Analysis Modules
- **data_processing.py**: Structure loading and preprocessing
- **structure_comparison.py**: RMSD calculations and alignment
- **helix_analysis.py**: Transmembrane helix identification
- **msa_grn.py**: Multiple sequence alignment with GRN
- **assign_grns.py**: GRN assignment orchestration

### New Analysis Modules
- **analyze_grns.py**: Conservation and co-evolution analysis
- **analyze_motifs.py**: Motif detection and validation
- **property_mapping.py**: Unified property-structure mapping

### Visualization Modules
- **visualization_functions.py**: Core plotting functions
- **plot.py**: Main visualization script with extended capabilities
- **opsin_color_scheme.py**: Consistent coloring for functional groups

## Dependencies

- **Python 3.8+**
- **Protos**: Protein structure analysis framework
- **BioPython**: Sequence and structure manipulation
- **NumPy/Pandas**: Data analysis
- **Matplotlib/Seaborn**: Visualization
- **scikit-learn**: Clustering and statistics

## Citation

If you use MOGRN in your research, please cite:
```
[Citation information to be added]
```

## Detailed Documentation

For comprehensive step-by-step instructions, see:
- **GUIDE.md**: Detailed English guide
- **GUIDE_JP.md**: Japanese translation (日本語ガイド)

## Troubleshooting

Common issues and solutions:
1. **Missing structures**: Ensure all CIF files are in correct directories
2. **Memory errors**: Use `--no-cache` flag or process in batches
3. **Import errors**: Verify Protos installation with `python -c "import protos"`
