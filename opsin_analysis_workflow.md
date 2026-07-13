# MOGRN Opsin Analysis Workflow - Technical Documentation

## Overview

The MOGRN (Microbial Opsin Generic Residue Numbering) workflow is a comprehensive pipeline for analyzing microbial opsin structures. It compares experimental structures from the PDB with Boltz-1 predicted structures, calculates various RMSD metrics, and assigns Generic Residue Numbers (GRNs) for standardized position referencing.

---

## Step 0: Data Preparation and Structure Registration

**Script**: `prepare_data.py`

### Dataset Organization

The workflow organizes structures into four datasets based on experimental status and temporal splits:

| Dataset | Description | Source |
|---------|-------------|--------|
| `mo_exp_A` | Experimental structures pre-Sept 2021 (within Boltz training window) | RCSB PDB |
| `mo_exp_B` | Experimental structures post-Sept 2021 (true test set) + Hideaki experimental | RCSB PDB + Hideaki |
| `mo_pred_exp` | Boltz-1 predictions of experimental structures | Boltz-1 predictions |
| `mo_pred_novel` | Boltz-1 predictions of novel opsins (no experimental counterpart) | Boltz-1 predictions |

### Structure Sources

- **RCSB structures**: Downloaded mmCIF files stored in `data/structure/mmcif/`
- **Predicted structures**: Boltz-1 outputs stored in `structures/mo_pred/`
- **Hideaki structures**: Experimental and predicted stored in `structures/hideaki_exp/` and `structures/hideaki_pred/`

### Property Metadata

**File**: `property/mo_exp_ST1.csv`

Contains:
- PDB IDs and short names
- `experimentally_determined` flag (1 = experimental, 0 = predicted only)
- `dataset_split` (A = pre-Sept 2021, B = post-Sept 2021)
- Opsin family classifications

### Structure Registration

The protos library (`StructureProcessor`) is used to:
1. Parse mmCIF files
2. Register structures with unique IDs
3. Create datasets grouping related structures

**Output**: `opsin_output/structure_mapping.json` - Maps experimental PDB IDs to their predicted counterparts (e.g., `"7bmh" -> "7bmh_model_0"`)

---

## Step 1: Structure Loading and Retinal Standardization

**Module**: `src/lyr_processing.py`

### The Retinal Naming Problem

Retinal (the chromophore in opsins) appears in structures under different residue names:

| Residue Name | Description | Source |
|--------------|-------------|--------|
| `RET` | Free retinal | Standard naming |
| `LYR` | Lysine-Retinal Schiff base (covalently bound) | Some experimental structures |
| `LIG` | Generic ligand | Boltz-1 predicted structures |

### LYR Processing: Splitting Covalent Schiff Base

When retinal forms a covalent Schiff base with lysine, the entire complex is named `LYR`. The workflow splits this into separate components:

**LYR → LYS + RET**

```
LYR residue atoms:
├── Lysine backbone atoms: N, CA, C, O, CB, CG, CD, CE, NZ → becomes LYS (ATOM)
└── Retinal atoms: C1-C20, etc. → becomes RET (HETATM)
```

**Implementation**:
```python
LYS_ATOM_NAMES_IN_LYR = {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'}
# Any atom NOT in this set is assigned to RET
```

### LIG Renaming

Boltz-1 predictions use generic `LIG` naming for the retinal ligand:

**LIG → RET**

Simple residue name replacement: `df.loc[df['res_name3l'] == 'LIG', 'res_name3l'] = 'RET'`

### Processing Order

The standardization is applied **early** in the workflow (before caching) to ensure consistency:

1. LYR → LYS + RET (split covalent complex)
2. LIG → RET (rename generic ligand)

---

## Step 2: Helix Annotation

**File**: `property/helices_grn.json`

### Helix Boundary Definitions

Microbial opsins have 7 transmembrane helices (TM1-TM7). Helix boundaries are defined based on GRN positions:

**Rule**: Each helix spans GRN positions X.41 to X.59 (±9 from the conserved X.50 position)

```json
{
  "7bmh": {
    "1": [53, 71],    // Helix 1: residues 53-71
    "2": [82, 99],    // Helix 2: residues 82-99
    "3": [135, 153],  // Helix 3: residues 135-153
    "4": [163, 181],  // Helix 4: residues 163-181
    "5": [187, 205],  // Helix 5: residues 187-205
    "6": [228, 246],  // Helix 6: residues 228-246
    "7": [261, 279]   // Helix 7: residues 261-279
  }
}
```

### Helix Number Assignment

Each CA atom in a structure is assigned a `helix_num` (1-7) based on whether its `auth_seq_id` falls within the defined boundaries. Residues outside helix boundaries are assigned `helix_num = 0` (loops/tails).

---

## Step 3: All-vs-All RMSD Matrix Calculation

**Module**: `src/structure_comparison.py`
**Function**: `compute_all_vs_all_rmsd_improved()`

### Purpose

Creates a symmetric RMSD matrix comparing all structures pairwise to identify:
- Structural clusters
- The global reference structure (lowest average RMSD to all others)
- Outliers with unusual conformations

### Atoms Used

- **CA atoms only** (C-alpha backbone)
- **Helix residues only** (helix_num 1-7, excluding loops and tails)
- Chain A only

### Alignment Method

Uses **CEalign** (Combinatorial Extension) via protos:

```python
R, t, alignment_path_indices, rmsd = get_structure_alignment(
    coords1, coords2,
    window_size=8,  # AFP fragment size
    max_gap=30      # Maximum gap in alignment path
)
```

CEalign returns:
- `R`: 3x3 rotation matrix
- `t`: Translation vector
- `alignment_path_indices`: Matched residue pairs
- `rmsd`: RMSD of the aligned path

### Output

- `opsin_output/rmsd_matrix.csv`: N×N symmetric matrix
- Used to select global reference: structure with lowest mean RMSD to all others

---

## Step 4: Two-Pass CEalign for Exp-Pred Comparison

**Module**: `src/structure_comparison.py`
**Function**: `calculate_binding_pocket_rmsd_for_pairs()`

### Why Two-Pass Alignment?

Single-pass CEalign can be misled by:
- Flexible loop regions with high variability
- Terminal regions with poor predictions
- Local misalignments in peripheral regions

The two-pass approach **filters outliers** to focus alignment on well-predicted core regions.

### Two-Pass Algorithm

#### Pass 1: Initial Alignment

```python
# CEalign on ALL CA atoms
R1, t1, path1, rmsd1 = get_structure_alignment(
    exp_ca_coords, pred_ca_coords,
    window_size=8, max_gap=30
)
```

#### Inlier Filtering

After Pass 1, identify residues with good alignment:

```python
INLIER_THRESH = 5.0  # Angstroms

# Apply Pass-1 transformation
pred_aligned = pred_coords @ R1 + t1

# Filter to inliers (residual ≤ threshold)
inlier_pairs = []
for exp_idx, pred_idx in zip(path1[0], path1[1]):
    dist = ||exp_coords[exp_idx] - pred_aligned[pred_idx]||
    if dist <= INLIER_THRESH:
        inlier_pairs.append((exp_idx, pred_idx))
```

#### Pass 2: Refined Alignment

```python
# CEalign on INLIERS ONLY
exp_inlier_coords = exp_coords[inlier_indices]
pred_inlier_coords = pred_coords[inlier_indices]

R2, t2, path2, rmsd2 = get_structure_alignment(
    exp_inlier_coords, pred_inlier_coords,
    window_size=8, max_gap=30
)
```

### Why This Works

1. **Robustness**: Outliers (flexible loops, mispredicted regions) don't bias the final transformation
2. **Focus**: The refined R2, t2 optimizes alignment for the structural core
3. **Honest Reporting**: The reported RMSD reflects alignment quality on the consensus region, not inflated by outliers
4. **Similar to RANSAC**: This iterative outlier removal is analogous to RANSAC in computer vision - fit on all data, identify inliers, refit on inliers only

### Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `INLIER_THRESH` | 5.0 Å | Maximum residual for inlier classification |
| `window_size` | 8 | CEalign AFP fragment size |
| `max_gap` | 30 | Maximum gap in alignment path |

---

## Step 5: Binding Pocket RMSD Calculation

### Binding Pocket Definition

Residues within **6.0 Å** of any retinal atom (based on CA-to-retinal minimum distance):

```python
cutoff = 6.0  # Angstroms

# Calculate minimum distance from each CA to any RET atom
dist_to_ret = calculate_min_distances(ca_coords, ret_coords)

# Pocket residues are those within cutoff
pocket_residue_ids = ca_df[dist_to_ret <= cutoff]['auth_seq_id'].unique()
```

### Residue Mapping via CEalign Path

The CEalign alignment path provides structural correspondence between experimental and predicted residues:

```python
# For each experimental pocket residue, find corresponding predicted residue
for exp_pocket_id in experimental_pocket_ids:
    for exp_path_idx, pred_path_idx in zip(exp_path, pred_path):
        if exp_ca_df.iloc[exp_path_idx]['auth_seq_id'] == exp_pocket_id:
            pred_pocket_id = pred_ca_df.iloc[pred_path_idx]['auth_seq_id']
            pocket_pairs.append((exp_pocket_id, pred_pocket_id))
            break
```

### RMSD Calculation

**All-atom RMSD** for pocket residues (not just CA):

```python
pocket_rmsd_sum_sq = 0
pocket_atom_count = 0

for exp_res_id, pred_res_id in pocket_pairs:
    exp_atoms = exp_df[exp_df['auth_seq_id'] == exp_res_id]
    pred_atoms = pred_df_aligned[pred_df_aligned['auth_seq_id'] == pred_res_id]

    # Match atoms by name (N, CA, C, O, CB, etc.)
    for atom_name in exp_atoms['res_atom_name'].unique():
        exp_coord = exp_atoms[exp_atoms['res_atom_name'] == atom_name][['x','y','z']].values[0]
        pred_coord = pred_atoms[pred_atoms['res_atom_name'] == atom_name][['x','y','z']].values[0]

        dist_sq = sum((exp_coord - pred_coord)**2)
        pocket_rmsd_sum_sq += dist_sq
        pocket_atom_count += 1

pocket_rmsd = sqrt(pocket_rmsd_sum_sq / pocket_atom_count)
```

---

## Step 6: Retinal RMSD Calculation

**Module**: `src/common_utils.py`
**Function**: `compute_retinal_mean_closest_distance()`

### Challenge

Retinal atom naming may differ between experimental and predicted structures, making direct atom-by-atom matching unreliable.

### Solution: Nearest-Neighbor Matching

For each predicted retinal atom, find the closest experimental atom and compute RMSD:

```python
def compute_retinal_mean_closest_distance(exp_ret_coords, pred_ret_coords):
    # Build distance matrix: pred (M) vs exp (N)
    dist_mat = cdist(pred_ret_coords, exp_ret_coords)  # shape (M, N)

    # For each predicted atom, find closest experimental atom
    min_dists = np.min(dist_mat, axis=1)  # shape (M,)

    # Compute RMSD
    return np.sqrt(np.mean(min_dists ** 2))
```

### Transformation Pipeline

The predicted retinal is transformed using the **backbone alignment R, t** before computing RMSD:

```python
# Transform predicted retinal using backbone alignment
aligned_pred_ret_coords = pred_ret_coords @ R_use + t_use

# Compute RMSD via nearest-neighbor matching
retinal_rmsd = compute_retinal_mean_closest_distance(
    exp_ret_coords,
    aligned_pred_ret_coords
)
```

No additional alignment refinement is performed on the retinal - the backbone alignment is used directly.

---

## Step 7: GRN Assignment

**Module**: `src/assign_grns.py`

### Generic Residue Numbering (GRN)

GRN provides a standardized numbering scheme for opsins:
- Format: `H.PP` where H = helix number (1-7), PP = position within helix
- **X.50** is the most conserved position in each helix
- Positions range from X.41 to X.59 (±9 from X.50)

### Tree-Based Alignment Method

1. **Build similarity tree** from RMSD matrix
2. **Select global reference** (structure with lowest mean RMSD: `7bmh`)
3. **Generate transitive alignment paths** through the tree
4. **Create MSA-like table** mapping residues to GRN positions

### Output Files

| File | Description |
|------|-------------|
| `msa_table_grn.csv` | Residue identities at each GRN position |
| `distance_table_grn.csv` | Retinal distances for each position |
| `ca_msa_table_grn.csv` | CA atom coordinates at GRN positions |
| `ca_distance_table_grn.csv` | CA-retinal distances |

---

## Summary of RMSD Metrics

| Metric | Atoms | Alignment | Purpose |
|--------|-------|-----------|---------|
| **Backbone RMSD** | CA only (inliers) | Two-pass CEalign | Overall structural similarity |
| **Pocket RMSD** | All atoms in pocket residues | Via CEalign path mapping | Binding site accuracy |
| **Retinal RMSD** | All retinal atoms | Hungarian matching after backbone alignment | Ligand placement accuracy |

---

## Configuration Files

| File | Purpose |
|------|---------|
| `property/mo_exp_ST1.csv` | Opsin metadata (PDB IDs, families, splits) |
| `property/helices_grn.json` | Helix boundaries per structure |
| `opsin_output/structure_mapping.json` | Exp→Pred structure ID mapping |

---

## Key Parameters

| Parameter | Value | Location |
|-----------|-------|----------|
| Retinal cutoff | 6.0 Å | Binding pocket definition |
| Inlier threshold | 5.0 Å | Two-pass CEalign filtering |
| CEalign window | 8 | Fragment size |
| CEalign max_gap | 30 | Maximum path gap |
| GRN helix range | X.41-X.59 | ±9 from X.50 |
| Global reference | 7bmh | Lowest mean RMSD |
