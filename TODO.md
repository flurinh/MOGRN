# MOGRN Refactoring Progress

## Completed

### 1. Data Preparation (`prepare_data.py`)
- [x] Switched to `mo_exp_ST1.csv` as primary property file
- [x] Fixed Excel scientific notation parsing (e.g., `1E12` PDB ID)
- [x] Implemented proper dataset creation with 4 datasets:
  - **Dataset 1a (mo_exp_A)**: 42 experimental structures (Set A, pre-Sept 2021)
  - **Dataset 1b (mo_exp_B)**: 27 experimental structures (Set B + 8 Hideaki)
  - **Dataset 2a (mo_pred_exp)**: 71 Boltz predictions of experimental structures (63 + 8 Hideaki)
  - **Dataset 2b (mo_pred_novel)**: 58 Boltz predictions of novel opsins
- [x] Structure mapping (PDB ID -> prediction): 70 pairs (62 standard + 8 Hideaki)
- [x] Hideaki structure handling with special naming convention
- [x] Dataset validation with property mapping verification
- [x] All datasets: 100% mapped to predictions and properties

### 2. File Renames
- [x] Renamed `ChR2_model_0.cif` -> `CrChR2_model_0.cif` (to match ST1 naming)

### 3. Core Files Created
- [x] `prepare_data.py` - Data preparation with protos API
- [x] `opsin_analysis_workflow.py` - Analysis workflow (needs update)
- [x] `plot.py` - Visualization script (needs update)

## In Progress

### 4. Opsin Analysis Workflow (`opsin_analysis_workflow.py`)
- [ ] Update to use correct protos dataset API
- [ ] Use structure mapping from `prepare_data.py` output
- [ ] Integrate with new dataset structure (1a, 1b, 2a, 2b)
- [ ] Fix data loading to work with registered datasets
- [ ] Update error calculation for exp vs pred comparison
- [ ] Update helix annotation pipeline
- [ ] Update RMSD matrix calculation
- [ ] Update GRN assignment

## Pending

### 5. Plot Script (`plot.py`)
- [ ] Update to use new workflow outputs
- [ ] Verify all visualizations work with new data structure

### 6. Testing & Validation
- [ ] Run full workflow end-to-end
- [ ] Verify RMSD calculations
- [ ] Verify GRN assignments
- [ ] Generate paper figures

## Dataset Summary

| Dataset | Name | Count | Description |
|---------|------|-------|-------------|
| 1a | mo_exp_A | 42 | Experimental Set A (pre-Sept 2021, in Boltz training) |
| 1b | mo_exp_B | 27 | Experimental Set B (post-Sept 2021) + Hideaki exp |
| 2a | mo_pred_exp | 71 | Boltz predictions of experimental (63 + 8 Hideaki) |
| 2b | mo_pred_novel | 58 | Boltz predictions of novel opsins |

## Key Files

- `property/mo_exp_ST1.csv` - Primary property file (134 entries)
- `opsin_output/structure_mapping.json` - PDB ID -> prediction mapping (70 pairs)
- `structures/mo_pred/` - Prediction CIF files (121)
- `structures/hideaki_exp/` - Hideaki experimental (8)
- `structures/hideaki_pred/` - Hideaki predictions (8)
