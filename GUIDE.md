# MOGRN Project Guide

This guide provides comprehensive step-by-step instructions for running the complete microbial opsin analysis workflow, including the new conservation and motif analysis capabilities.

## Table of Contents
1. [Project Overview](#project-overview)
2. [Prerequisites](#prerequisites)
3. [Environment Setup](#environment-setup)
4. [Step-by-Step Analysis](#step-by-step-analysis)
5. [Understanding the Results](#understanding-the-results)
6. [Advanced Options](#advanced-options)
7. [Troubleshooting](#troubleshooting)

## Project Overview

The MOGRN (Microbial Opsin Generic Residue Numbering) project provides a comprehensive framework for analyzing microbial opsin structures. The workflow now includes six main steps:

1. **Data Preparation**: Initialize infrastructure and validate datasets
2. **Sequence Configuration**: Generate YAML files for protein sequences
3. **Core Analysis**: Perform structural alignment and GRN assignment
4. **Visualization**: Create standard plots and heatmaps
5. **Conservation Analysis**: Analyze GRN conservation patterns
6. **Motif Analysis**: Detect and validate functional motifs

### What is Generic Residue Numbering (GRN)?

GRN provides a standardized way to compare residue positions across different opsin structures. Key positions end with ".50" (e.g., 3.50, 7.50) and serve as structural anchors. This allows consistent comparison even when proteins have different lengths or insertions/deletions.

## Prerequisites

### System Requirements
- Python 3.10+
- Ideally 16GB RAM or higher

### Installing Dependencies

1. **Create and activate conda environment** (recommended):
```bash
conda create -n mogrn python=3.10
conda activate mogrn
```

2. **Install Python packages**:
```bash
pip install -r requirements.txt
```

3. **Install Protos Framework**:
```bash
git clone https://github.com/flurinh/protos.git
cd protos
pip install -e .
cd ..
```

### Required Data Files

Create these directories and add your data files:

```
MOGRN/
├── property/
│   ├── mo_exp.csv        # Required: Experimental properties
│   └── helices.json      # Required: Helix definitions
└── structures/
    ├── hideaki_exp/      # Experimental structures (CIF format)
    ├── hideaki_pred/     # Predicted structures (CIF format)
    └── mo_pred/          # Additional predictions (CIF format)
```

## Environment Setup

Before running any analysis, activate the conda environment:

```bash
source /home/[username]/miniconda/etc/profile.d/conda.sh
conda activate mogrn
```

Verify the setup:
```bash
python -c "import protos; print('Protos installed successfully')"
python -c "import pandas, numpy, matplotlib; print('Core packages ready')"
```

## Step-by-Step Analysis

### Step 1: Data Preparation

**Command:**
```bash
python prepare_data.py
```

**What happens:**
- Initializes the ProtosPaths system for file management
- Creates necessary directory structure
- Sets up CifBaseProcessor for structure handling
- Validates all structure files are accessible
- Creates cache directories for efficient processing

**Expected output:**
```
Setting up data directories...
Initializing ProtosPaths...
Processing datasets:
  - mo_exp: X structures
  - mo_pred: Y structures
  - hideaki_exp: Z structures
  - hideaki_pred: W structures
Data preparation complete!
```

**Troubleshooting:**
- If you see "FileNotFoundError", check that property/ and structures/ directories exist
- If structures aren't found, verify CIF files are in the correct subdirectories

### Step 2: Sequence Configuration

**Command:**
```bash
python prepare_yaml.py
```

**What happens:**
- Reads protein sequences from property/mo_exp.csv
- Creates individual YAML files for each protein
- Handles special cases (multiple retinals, long sequences)
- Validates sequence integrity

**Expected output:**
```
Processing sequences from property/mo_exp.csv...
Created YAML files for X proteins
Special cases:
  - Proteins with multiple retinals: Y
  - Long sequences (>400 aa): Z
YAML preparation complete!
```

### Step 3: Core Analysis Workflow

**Command:**
```bash
python opsin_analysis_workflow.py
```

**What happens:**
1. **Structure Loading**: Loads all structures with caching
2. **RMSD Calculation**: Computes pairwise structural similarities
3. **Helix Assignment**: Identifies transmembrane regions
4. **GRN Assignment**: Performs multiple sequence alignment with standardized numbering
5. **Distance Analysis**: Calculates residue distances to retinal

### Step 4: Standard Visualizations

**Command:**
```bash
python plot.py
```

**What happens:**
- Loads analysis results from output directory
- Creates various visualization types
- Applies quality filters
- Generates publication-ready figures

**Key visualizations:**
- `rmsd_heatmap.png`: Structure similarity matrix
- `similarity_tree.png`: Hierarchical clustering
- `conservation_plots/`: Position-specific conservation
- `helix_logos/`: Sequence logos for each helix

### Step 5: GRN Conservation Analysis (NEW)

**Command:**
```bash
python analyze_grns.py
```

**What happens:**
1. **Conservation Analysis**: Examines conservation at .50 positions by functional group
2. **Pattern Detection**: Identifies discriminating positions between pumps and channels
3. **Co-evolution Analysis**: Detects correlated residue changes
4. **Report Generation**: Creates detailed interpretation reports

**Key results to examine:**
1. **Comprehensive summary figure**: Multi-panel visualization showing all key findings
2. **Discrimination heatmap**: Shows which positions distinguish functional groups
3. **Co-evolution patterns**: Networks of correlated positions

### Step 6: Motif Pattern Analysis (NEW)

**Command:**
```bash
python analyze_motifs.py
```

**What happens:**
1. **Single Position Analysis**: Conservation at each GRN position
2. **Literature Motif Validation**: Checks known motifs (DTD etc.)
3. **Correlation Analysis**: Finds new correlated patterns
4. **Novel Motif Detection**: Identifies previously unknown patterns

**Key outputs:**
1. **Extended motif table**: Comprehensive table of all motifs with statistics
2. **Correlation matrices**: Shows which positions co-vary
3. **Novel family detection**: New motif patterns specific to functional groups

## Understanding the Results

### Conservation Scores
- **>90%**: Universally conserved, critical for function
- **70-90%**: Highly conserved within functional groups
- **50-70%**: Moderately conserved, may determine specificity
- **<50%**: Variable positions, possibly involved in regulation

### Functional Discrimination
The analysis identifies "signature" positions that distinguish functional groups:
- **Position 3.50**: T (threonine) in pumps → C (cysteine) in channels
- **Position 5.50**: S (serine) in pumps → G (glycine) in channels
- **Position 4.50**: M (methionine) in pumps → T/R in channels

### Motif Interpretation
- **Conserved motifs**: Essential for basic opsin function
- **Modified motifs**: Determine specific functional properties
- **Novel motifs**: May represent undiscovered functional elements

## Troubleshooting

### Common Issues and Solutions

1. **ModuleNotFoundError: No module named 'protos'**
   - Solution: Ensure Protos is installed: `cd protos && pip install -e .`

2. **MemoryError during RMSD calculation**
   - Solution: Use `--no-cache` flag

3. **Empty figures or missing data**
   - Check that previous steps completed successfully
   - Verify output files exist in expected directories
   - Check protein property mappings in mo_exp.csv

4. **Slow performance**
   - Enable caching (default)
   - Get a better PC ;)


### Getting Help

Open an issue on GitHub with:
   - Error message
   - Steps to reproduce
   - System information