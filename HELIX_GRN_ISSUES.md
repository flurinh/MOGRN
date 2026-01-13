# MOGRN Helix Boundaries and GRN Table Issues

## Overview

This document outlines critical issues identified in the GRN (Generic Residue Numbering) assignment workflow and proposes solutions.

---

## Issue 1: Truncated Helix Boundaries

### Problem

The `property/helices_curated.json` file contains helix boundary definitions that are **too strict**, causing:

1. **Missing residues at helix ends** - The GRN tables are missing residues at the N-terminal and C-terminal ends of each helix
2. **Over-confident alignments** - The alignment appears more accurate than it actually is because problematic edge residues are excluded
3. **Incomplete coverage** - The final GRN table only covers 171 TM positions when more should be included

### Current Format

```json
{
  "CnChR2_J230_refine9": {
    "1": [84, 109],   // Helix 1: residues 84-109
    "2": [120, 145],  // Helix 2: residues 120-145
    ...
  }
}
```

### Evidence

From the workflow output:
```
[INFO] TM residue positions: 171
[INFO] Helix 1: 25 positions
[INFO] Helix 2: 25 positions
...
```

The 25 positions per helix suggest tight boundaries. Actual TM helices typically have 20-30 residues, but the **interface regions** (where helices connect to loops) are being excluded.

### Solutions

#### Option A: Expand GRN Table Post-Creation (Recommended)

After GRN assignment, extend the table by:
1. For each structure, identify residues adjacent to assigned helix boundaries
2. Assign extended GRN positions (e.g., 1.36, 1.35 for N-terminal extension of helix 1)
3. Use structural alignment to validate extended positions

**Pros**: Doesn't require re-curating helix definitions; can be done programmatically
**Cons**: Extended positions may have lower confidence

#### Option B: Expand Helix Boundaries in JSON

Manually expand each helix boundary by 3-5 residues on each end:

```json
{
  "CnChR2_J230_refine9": {
    "1": [79, 114],   // Extended from [84, 109]
    "2": [115, 150],  // Extended from [120, 145]
    ...
  }
}
```

**Pros**: More accurate helix definitions
**Cons**: Requires manual curation for each structure; may introduce alignment errors at true loop regions

---

## Issue 2: Global Reference Selection

### Problem

The current global reference selection logic has issues:

1. ~~**Hard-coded reference** in `src/reference_alignment.py:295-296`~~ **FIXED**

2. **Tree-based method uses MerMAID1** which is a predicted structure without retinal coordinates

3. **Inconsistent selection** - Different methods (global vs tree-based) use different references

### FIXED: Removed Hardcoded Reference

The hardcoded reference `CnChR2_J230_refine9` was removed from `src/reference_alignment.py`.

Now the selection follows proper logic:
- **DataFrame input (RMSD matrix)**: Selects structure with lowest average RMSD to all others
- **Dict input**: Selects structure with best resolution from type references
- **Current default**: 4fbz is automatically selected as it has lowest avg RMSD

### Current Logic (`find_global_reference`)

```python
# From type references, find the one with lowest average RMSD to all structures
for struct_type, struct_id in type_references.items():
    avg_rmsd = all_structures.loc[struct_id].mean()
    if avg_rmsd < best_avg_rmsd:
        best_global_ref = struct_id
```

### Recommended Solution

**Two-pass approach:**

1. **Pass 1 (Bootstrap)**: Run initial analysis with approximate helix boundaries
   - Calculate full RMSD matrix
   - Identify the structure with **lowest average RMSD to all others** (most central)
   - This should be an **experimental structure** with retinal

2. **Pass 2 (Production)**: Use the identified reference
   - Re-run GRN assignment with the optimal reference
   - Ensure reference has complete retinal coordinates

### Selection Criteria for Global Reference

The ideal global reference should:
- Be an **experimental structure** (not predicted)
- Have **complete retinal (RET)** coordinates
- Have **low average RMSD** to all other structures
- Be from a **well-characterized opsin** (e.g., bacteriorhodopsin, channelrhodopsin)

### Suggested References by Quality

| Structure | Type | Retinal | Avg RMSD | Notes |
|-----------|------|---------|----------|-------|
| 4fbz | Experimental | Yes | Low | Current global ref (proton pump) |
| CnChR2_J230_refine9 | Hideaki | Yes | Low | Good channelrhodopsin reference |
| 1c3w | Experimental | Yes | Low | Classic bacteriorhodopsin |

---

## Issue 2b: Helix Length Inconsistencies in helices_curated.json

### Problem

Analysis of `property/helices_curated.json` reveals inconsistent helix definitions:

| Helix | Min Length | Max Length | Mean | Outliers |
|-------|------------|------------|------|----------|
| 1 | 1 | 37 | 24.1 | VbACR2_model_0 (1!), OtHKR_model_0 (37) |
| 2 | 21 | 41 | 25.5 | 8JH0 (41), B1ChR2_model_0 (39) |
| 3 | 19 | 32 | 22.2 | - |
| 4 | 19 | 32 | 22.5 | - |
| 5 | 18 | 35 | 24.4 | - |
| 6 | 17 | 33 | 25.3 | - |
| 7 | 23 | 29 | 25.8 | - |

### Critical Issues

1. **VbACR2_model_0 Helix 1**: `[46, 46]` - single residue! This is clearly incorrect.
2. Large variance in helix lengths suggests inconsistent annotation methods

### Known Non-7TM Structures

Some structures legitimately have fewer than 7 TM helices:

| Structure | Helices | Notes |
|-----------|---------|-------|
| S13_Bin138_Proteo_SR | 6 | Bacterial signal activator |

### New Approach: Phi/Psi-Based Helix Detection

Created `scripts/analyze_helices_phipsi.py` which:
1. Uses phi/psi dihedral angles to detect helical residues (α-helix: φ≈-60°, ψ≈-45°)
2. Applies sliding window smoothing to find continuous helical segments
3. Assigns helix numbers by N-to-C order (sequence alignment for mapping planned)

### SOLUTION: scripts/generate_helices_json.py

Created a new script that:
1. Uses phi/psi dihedral angles to detect helical residues (α-helix: φ≈-60°, ψ≈-45°)
2. Calculates phi/psi from coordinates when protos values are missing (fixes Hideaki structures)
3. Aligns sequences to CnChR2 reference to identify TM helix regions (1-7)
4. Finds continuous helical segments within each TM region

**Output**: `property/helices_phipsi_aligned.json`

**Results**:
- 198 structures processed (all datasets)
- 57 structures with all 7 helices detected
- Hideaki experimental structures now properly processed

**Limitations**:
- Sequence alignment doesn't work well for all opsins (variable N-terminal length)
- Structures with very different sequences may miss helix 1 and/or 2
- Some opsins have <7 TM helices biologically (e.g., S13_Bin138_Proteo_SR has 6)

**Recommended Next Steps**:
- Use profile-based alignment (PFAM rhodopsin HMM) instead of pairwise
- Or use structural superposition instead of sequence alignment
- Consider merging with curated helices for problematic structures

---

## Issue 2c: GRN Table Density (Sparse Alignment)

### Problem

The GRN table has 34.4% gaps (8,967 out of 26,051 cells).

Gap distribution by region:
| Region | Columns | Gap % | Notes |
|--------|---------|-------|-------|
| N-terminal (n.) | 5 | 99.1% | Only reference has these |
| Helix 1 | 25 | 33.8% | Variable N-terminal extensions |
| Helix 2 | 25 | 4.5% | Good coverage |
| Helix 3 | 22 | 2.1% | Good coverage |
| Helix 4 | 23 | 5.0% | Good coverage |
| Helix 5 | 25 | 4.6% | Good coverage |
| Helix 6 | 25 | 2.8% | Good coverage |
| Helix 7 | 26 | 7.2% | Some gaps |
| C-terminal (c.) | 13 | 99.1% | Only reference has these |
| Loops | ~50 | ~100% | Not propagated from alignment |

### Root Cause

The alignment paths from structure comparison don't cover:
1. N-terminal/C-terminal residues (not in TM helices)
2. Loop residues (alignment focuses on helices)
3. Helix edge residues (alignment may fail at boundaries)

### Recommended Solutions

1. **Post-processing for TM helices**: Fill gaps using:
   - Sequence alignment for missing helix residues
   - Nearest neighbor interpolation for edge cases

2. **Loop handling**: Mark loop positions as "unassigned" rather than gaps
   - Loops vary significantly between structures
   - Current approach is actually correct for loops

3. **N/C-terminal handling**: Only include if structure has the residues
   - These are structure-specific
   - Consider removing from final table or marking optional

---

## Issue 3: LYR/LIG/RET Naming Inconsistencies

### Problem

Retinal is represented differently across structures:
- **RET** - Standard retinal ligand
- **LIG** - Generic ligand (used in some predicted structures)
- **LYR** - Lysine-retinal Schiff base (covalently bound)

### FIXED: Chain ID Issue for Predicted Structures

**Root cause**: Boltz predictions place protein on chain A but LIG (retinal) on chain B. When filtering by chain A, retinal was lost.

**Fix applied** in `opsin_analysis_workflow.py` (lines 218-222):
```python
# For predicted structures (_model_0), set all chains to A
# (Boltz puts protein on chain A but LIG on chain B)
if "_model_0" in struct_id:
    df = df.copy()
    df["auth_chain_id"] = chain_id
```

**Result**: "No retinal found" messages reduced from 60+ to 0.

### Current Handling

The codebase has multiple conversion points:

#### 1. LYR Processing (`src/lyr_processing.py`)

Splits LYR into:
- **LYS** (lysine backbone atoms: N, CA, C, O, CB, CG, CD, CE, NZ)
- **RET** (remaining atoms = retinal moiety)

```python
LYS_ATOM_NAMES_IN_LYR = {'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'CE', 'NZ'}

def convert_single_lyr_entry(lyr_df, retinal_res_name='RET'):
    for atom_row in lyr_df.iterrows():
        if atom_row['res_atom_name'] in LYS_ATOM_NAMES_IN_LYR:
            atom_row['res_name3l'] = 'LYS'
        else:
            atom_row['res_name3l'] = retinal_res_name  # 'RET'
```

#### 2. LIG → RET Conversion (`src/msa_grn.py`, `src/data_processing.py`)

```python
# Multiple locations do this conversion:
df_norm.loc[df_norm['res_name3l'] == 'LIG', 'res_name3l'] = 'RET'
```

### Issues with Current Approach

1. **Scattered conversions** - LIG→RET conversion happens in multiple files
2. **Order dependency** - LYR must be processed before other operations
3. **Missing retinal warnings** - Many structures show "No retinal (RET or LIG) found"

### Recommended Solution

**Centralized preprocessing pipeline:**

```python
def standardize_retinal_naming(processor_data):
    """
    Standardize all retinal naming to 'RET'.
    Should be called early in the pipeline.

    Order of operations:
    1. Process LYR → LYS + RET (if LYR present)
    2. Convert LIG → RET (if LIG present)
    3. Validate retinal presence
    """
    # Step 1: LYR processing
    if 'LYR' in processor_data['res_name3l'].values:
        processor_data = process_lyr_in_dataframe(processor_data)

    # Step 2: LIG → RET
    processor_data.loc[
        processor_data['res_name3l'] == 'LIG', 'res_name3l'
    ] = 'RET'

    # Step 3: Validate
    ret_count = (processor_data['res_name3l'] == 'RET').sum()
    if ret_count == 0:
        logging.warning("No retinal found after standardization")

    return processor_data
```

---

## Issue 4: Structures Without Retinal

### Problem

Many predicted structures lack retinal coordinates:
```
[DEBUG] No retinal (RET or LIG) found for structure MerMAID1_model_0, setting all distances to NaN
[DEBUG] No retinal (RET or LIG) found for structure GtACR2_model_0, setting all distances to NaN
...
```

This affects:
- Distance calculations (all become NaN)
- Binding pocket analysis
- Selection of reference structures

### Solutions

1. **Exclude from distance calculations** but keep in MSA (current behavior)
2. **Inherit retinal** from closest experimental structure with known retinal
3. **Flag in output** which structures lack retinal data

---

## Implementation Priorities

### High Priority

1. **Fix global reference selection**
   - Remove hard-coded reference
   - Implement two-pass selection based on RMSD centrality
   - Ensure selected reference has retinal

2. **Centralize retinal preprocessing**
   - Single function for LYR/LIG/RET standardization
   - Call early in pipeline (in `prepare_data.py` or `opsin_analysis_workflow.py`)

### Medium Priority

3. **Expand helix boundaries**
   - Option A: Post-hoc GRN table expansion
   - Option B: Curate expanded boundaries in JSON

### Low Priority

4. **Document structures without retinal**
   - Add metadata field indicating retinal presence
   - Consider separate analysis tracks for structures with/without retinal

---

## Files to Modify

| File | Changes |
|------|---------|
| `src/reference_alignment.py` | Remove hard-coded reference (line 295-296), improve `find_global_reference()` |
| `src/assign_grns.py` | Add helix boundary expansion logic |
| `src/msa_grn.py` | Add GRN table post-expansion function |
| `opsin_analysis_workflow.py` | Add centralized retinal preprocessing step |
| `property/helices_curated.json` | (Optional) Expand helix boundaries |

---

## Testing Checklist

- [x] Global reference is selected automatically (not hard-coded) - **FIXED**
- [x] Selected reference is an experimental structure with retinal (4fbz auto-selected)
- [x] LYR residues are properly split into LYS + RET
- [x] LIG residues are renamed to RET
- [x] Predicted structures have chain IDs fixed (all set to A) - **FIXED**
- [ ] GRN table includes helix edge residues
- [ ] Distance tables have valid values for structures with retinal
- [x] Structures without retinal are flagged but included in MSA
- [x] Phi/psi angles calculated for structures missing them (Hideaki) - **NEW FIX**
- [x] New helix JSON generated with phi/psi-based detection - **NEW**

---

## References

- GRN notation: `N.XX` where N = helix number (1-7), XX = position relative to pivot (50)
- Pivot positions (X.50) are closest to retinal in each helix
- Loop notation: `NM.XXX` where N,M = flanking helices, XXX = distance from N, N is the nearer flanking helix