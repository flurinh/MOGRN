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
git clone https://github.com/flurinh/MOGRN.git
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

## Data Download

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18147121.svg)](https://doi.org/10.5281/zenodo.18147121)

The analysis data is hosted on Zenodo as a single archive (~650 MB compressed, ~2 GB uncompressed). Download and extract using:

```bash
# Download and extract all data
python download_data.py 18147121

# List available files without downloading
python download_data.py 18147121 --list

# Force overwrite existing data folders
python download_data.py 18147121 --overwrite

# Keep the archive after extraction
python download_data.py 18147121 --keep-archive
```

### Data Folders Downloaded from Zenodo

| Folder | Size | Description |
|--------|------|-------------|
| `opsin_output/` | ~1.4 GB | Main analysis output (figures, cache, GRN analysis) |
| `outputs/` | ~380 MB | Boltz prediction results |
| `new_opsins_outputs/` | ~140 MB | Boltz results for new opsins |
| `flat_outputs/` | ~66 MB | Flattened output files |
| `structures/` | ~32 MB | Structure files (hideaki_exp/pred, mo_pred) |
| `flat_new_opsins_outputs/` | ~6 MB | Flattened new opsin outputs |
| `property/` | ~1.5 MB | Property data files |
| `yaml_configs/` | ~500 KB | YAML configuration files |

**Note**: The `data/` folder (protos working directory) is regenerated locally when you run the analysis and is not included in the Zenodo archive.

### Required Data Contents

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
└── opsin_output/                    # All analysis outputs
    ├── rmsd_matrix.csv             # Pairwise structural similarities
    ├── curated_msa.csv             # Curated multiple sequence alignment
    ├── analysis_summary.json       # Summary of analysis parameters
    ├── hideaki_errors.csv          # Error analysis for Hideaki dataset
    ├── mo_exp_errors.csv           # Error analysis for experimental structures
    ├── cache/                      # Cached intermediate results
    │   ├── processed_structures_*.pkl
    │   ├── helix_annotations_*.pkl
    │   └── grn_assignment_*.pkl
    ├── global_reference_grn/       # Global GRN alignment results
    │   ├── msa_table_grn.csv      # Full MSA with GRN positions
    │   ├── ca_msa_table_grn.csv   # CA-only MSA
    │   ├── distance_table_grn.csv # Sidechain-retinal distances
    │   └── ca_distance_table_grn.csv # CA-retinal distances
    ├── tree_based_grn/             # Tree-based GRN alignment
    │   └── [same files as global_reference_grn]
    ├── grn_analysis/               # GRN conservation analysis
    │   ├── comprehensive_analysis_results.json
    │   ├── conservation/           # Conservation data by group
    │   │   ├── conservation_by_function.json
    │   │   ├── conservation_by_domain.json
    │   │   ├── motifs_by_function.json
    │   │   └── neighborhood_by_*.json
    │   ├── figures/                # Analysis visualizations
    │   │   ├── grn_conservation_comprehensive_summary.png
    │   │   ├── functional_discrimination_heatmap.png
    │   │   ├── coevolution_patterns.png
    │   │   └── *_conservation_heatmap.png
    │   └── reports/                # Detailed text reports
    │       ├── conservation_summary_report.txt
    │       └── grn_pattern_interpretation.txt
    ├── motifs/                     # Motif analysis results
    │   ├── comprehensive_motif_report.md
    │   ├── extended_motif_table.csv
    │   ├── literature_motifs.csv
    │   ├── single_position_discriminators.csv
    │   ├── correlations_*residue.csv
    │   ├── cross_function_patterns.csv
    │   └── comprehensive_motif_analysis.png
    └── paper_figures/              # Publication-ready figures
        ├── 01_opsin_overview.png
        ├── 02_rmsd_clustermap.png
        ├── 03_all_atom_distance_std.png
        ├── 04_ca_atom_distance_std.png
        ├── 05_helix_logos_x50.png
        ├── 06_property_analysis_natural_domains.png
        ├── 07_prediction_contribution_bars.png
        └── interactive_grn_alignment.html
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
python plot.py --input-dir custom_output --outputs-dir custom_figures

# Specific plot types
python plot.py --plots "rmsd,conservation,distance"

# High DPI outputs
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

- **Python 3.10+**
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

For comprehensive step-by-step instructions, see the GUIDE.md

## Troubleshooting

Common issues and solutions:
1. **Missing structures**: Ensure all CIF files are in correct directories
2. **Memory errors**: Use `--no-cache` flag or process in batches
3. **Import errors**: Verify Protos installation with `python -c "import protos"`

---

## Development Status (2026-01-05)

### Current Goals

We are working on improving the GRN alignment quality and visualization:

1. **GRN Alignment Quality Control**
   - Ensuring all structures are properly aligned to a global reference
   - Detecting misaligned residues by comparing to consensus positions
   - Validating lysine anchor at position 7.50 across all structures

2. **Global Reference Selection**
   - Changed global reference from `MerMAID1_model_0` to **`1jgj`** (Sensory Rhodopsin II)
   - 1jgj chosen for having best alignment quality (lowest mean distance from consensus)

3. **Interactive Visualization**
   - Centering coordinates on overall center for better rotation
   - Property mapping from `mo_exp_ST1.csv` for domain and molecular function
   - Handling Hideaki structure naming patterns (e.g., `TsChR_J132_refine3`)

### Recent Changes

- **Global reference**: Set to `1jgj` in all visualization and alignment modules
- **Property file**: Updated to use `mo_exp_ST1.csv`
- **PDB ID parsing**: Fixed scientific notation issue where `1e12` was read as `1.00E+12`
- **Coordinate centering**: Changed from retinal-center to overall-center for interactive plots
- **Misalignment detection**: Created `scripts/detect_misalignment.py` to identify outliers

### Completed Tasks

- [x] Set global reference to `1jgj` in source files
- [x] Delete all cache files and regenerate
- [x] Rerun structure comparison workflow with new reference
- [x] Run misalignment detection analysis

### Pending Tasks / Known Issues

- [ ] Investigate **7pl9** - severely misaligned (22 Å mean distance, 99% outliers)
- [ ] Investigate **OcHeR_model_0** - systematic issues (6 Å mean distance, 88% outliers)
- [ ] Investigate **7aky** and **AbHeR_model_0** - moderate issues (~5 Å mean distance)
- [ ] Review structures with localized issues at helix termini
- [ ] Validate postprocessor K-anchor shifts

### Misalignment Analysis Results

With `1jgj` as global reference:

| Structure | Mean Distance (Å) | Outlier % | Status |
|-----------|-------------------|-----------|--------|
| 7pl9 | 22.05 | 99.4% | Needs investigation |
| OcHeR_model_0 | 6.07 | 87.8% | Needs investigation |
| 7aky | 5.11 | 40% | Moderate issues |
| AbHeR_model_0 | 4.96 | 37% | Moderate issues |
| VbACR2_model_0 | 3.48 | 16% | Localized issues |
| Most structures | < 2 Å | < 5% | Good alignment |

### Key Files Modified

- `src/visualize_alignment_grn.py` - Reference defaults changed to `1jgj`
- `src/structure_alignment_subset.py` - Reference defaults changed to `1jgj`
- `scripts/detect_misalignment.py` - New script for alignment quality control
- `scripts/postprocess_grn_table.py` - GRN postprocessing with K-anchor validation
- `plot.py` - Property mapping and Hideaki structure handling

### Running the Analysis

```bash
# Full workflow with 1jgj reference
python opsin_analysis_workflow.py --global-ref 1jgj --skip-prepare

# Misalignment detection
python scripts/detect_misalignment.py

# Results saved to:
# - opsin_output/misalignment_analysis/structure_outliers.csv
# - opsin_output/misalignment_analysis/grn_statistics.csv
# - opsin_output/misalignment_analysis/distance_matrix.csv
```
