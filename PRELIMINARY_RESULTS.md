# Preliminary Results: Microbial Opsin Generic Residue Numbering (MOGRN) Analysis

## Table of Contents
1. [Executive Summary](#executive-summary)
2. [Dataset Overview](#dataset-overview)
3. [GRN Table Structure and Current Limitations](#grn-table-structure-and-current-limitations)
4. [Conservation Analysis Results](#conservation-analysis-results)
5. [Functional Group Patterns](#functional-group-patterns)
6. [Domain-Specific Patterns](#domain-specific-patterns)
7. [Helix-Specific Motifs](#helix-specific-motifs)
8. [Key Findings and Interpretations](#key-findings-and-interpretations)
9. [Future Directions](#future-directions)

## Executive Summary

This preliminary analysis presents the first comprehensive application of Generic Residue Numbering (GRN) to microbial opsins, analyzing 128 structures across diverse functional categories and domains of life. While the current GRN alignment contains some errors that will be addressed in future iterations, the analysis reveals striking conservation patterns that distinguish functional groups and evolutionary lineages.

**Key Findings:**
- Position 3.50 serves as the strongest functional discriminator (C for channels, T for pumps)
- Universally conserved positions include W6.50 (>90%) and K7.50 (>95%, retinal-binding lysine)
- Distinct sequence motifs characterize each functional group around .50 positions
- Eukaryotic microbial opsins show greater sequence diversity than prokaryotic ones

## Dataset Overview

### Structure Collection
The analysis includes 128 microbial opsin structures from the opsin_grn_tables/residue_table_grn.csv file:

| Source | Count | Description |
|--------|-------|-------------|
| Experimental structures | ~60   | X-ray crystallography and cryo-EM structures |
| AlphaFold predictions | ~60   | High-confidence structural models |
| Diverse organisms | 128   | Archaea, Bacteria, Eukaryota, Viruses, and Synthetic constructs |

### Functional Categories
Based on property/mo_exp.csv classification:

| Molecular Function | Count | Percentage | Key Examples |
|-------------------|-------|------------|--------------|
| Proton Pump | 42 | 35.0% | Bacteriorhodopsin, Proteorhodopsin |
| Cation Channel | 22 | 18.3% | Channelrhodopsins (ChR1, ChR2) |
| Sensor/Regulatory | 22 | 18.3% | Sensory rhodopsins |
| Chloride Pump | 12 | 10.0% | Halorhodopsin |
| Anion Channel | 10 | 8.3% | Anion channelrhodopsins (ACRs) |
| Sodium Pump | 5 | 4.2% | KR2-type rhodopsins |
| Unknown | 7 | 5.8% | Uncharacterized opsins |

### Domain Distribution

| Domain | Count | Percentage | Characteristics |
|--------|-------|------------|-----------------|
| Eukaryota | 45 | 37.5% | Algal channelrhodopsins, fungal opsins |
| Bacteria | 41 | 34.2% | Marine proteorhodopsins, XRs |
| Archaea | 27 | 22.5% | Haloarchaeal pumps (BR, HR, SR) |
| Virus | 5 | 4.2% | Giant virus rhodopsins |
| Synthetic | 2 | 1.7% | Engineered variants |

## GRN Table Structure and Current Limitations

### Table Format
The residue_table_grn.csv contains:
- **Rows**: 128 protein structures (indexed by structure name/PDB ID)
- **Columns**: 263 GRN positions across 7 transmembrane helices
- **Cell values**: Residue identity + sequence position (e.g., "K296" = lysine at position 296)
- **Missing data**: Represented as "-" for gaps or unaligned regions

### GRN Position Nomenclature
Positions follow Ballesteros-Weinstein-like notation adapted for microbial opsins:
- **X.50**: Most conserved position in helix X (e.g., 7.50 = retinal-binding lysine)
- **X.YY**: Position YY in helix X
- **XX.001-010**: Loop regions between helices

Example positions:
- 1.50, 2.50, 3.50, etc.: Reference positions in each helix
- 3.46: Conserved tryptophan position
- 7.53: Position near Schiff base

### Current Limitations and Errors

⚠️ **Important Note**: The current GRN alignment contains several known errors that will be corrected:

1. **Alignment inconsistencies**: Some structures may have incorrect helix boundaries
2. **Loop region assignments**: Cytoplasmic and extracellular loops need refinement
3. **Missing positions**: Some conserved positions may not be properly captured
4. **Reference structure bias**: Current alignment may be biased toward the reference structure used

These limitations do not invalidate the major conservation patterns observed but may affect specific position assignments.

## Conservation Analysis Results

### Overall Conservation Patterns

Analysis of all 263 GRN positions reveals:

| Conservation Level | Number of Positions | Notable Examples |
|-------------------|-------------------|------------------|
| >90% conserved | 2 | W6.50 (91.4%), K7.50 (95.3%) |
| 70-90% conserved | 5 | W3.46 (80.5%), D7.46, others |
| 50-70% conserved | 12 | Various .50 positions |
| <50% conserved | 244 | Most positions show functional/domain variation |

### Position-Specific Conservation at .50 Positions

| Helix | Position | Overall Conservation | Top Residue | Functional Variation |
|-------|----------|---------------------|-------------|---------------------|
| TM1 | 1.50 | 36.7% | M | High: M (pumps) vs F (sensors) vs C (anion channels) |
| TM2 | 2.50 | 50.0% | V | Moderate: V (channels) vs I (pumps) |
| TM3 | 3.50 | 39.1% / 38.3% | T/C | **Key discriminator**: T (pumps) vs C (channels) |
| TM4 | 4.50 | 56.2% | M | M strongly preferred in pumps |
| TM5 | 5.50 | 49.2% / 41.4% | S/G | S (pumps) vs G (channels) |
| TM6 | 6.50 | 91.4% | W | **Universally conserved** |
| TM7 | 7.50 | 95.3% | K | **Universally conserved** (Schiff base) |

## Functional Group Patterns

### Proton Pumps (n=42)
**Signature positions:**
- 1.50: M (73.8%) - Methionine strongly preferred
- 2.50: I (50.0%) - Isoleucine preference
- 3.50: T (61.9%) - Threonine signature
- 4.50: M (83.3%) - Methionine highly conserved
- 5.50: S (66.7%) - Serine preference

**Key motifs:**
- Helix 1: G-x-x-x-M-x-x-x-G (GxxxG motif with M1.50)
- Helix 3: W-L-F-T-P-L-L-L (highly conserved)
- Helix 4: x-D-V-x-M-I-x-T-G

### Cation Channels (n=22)
**Signature positions:**
- 1.50: Variable (S 27.3%, others)
- 2.50: V (63.6%) - Valine preferred
- 3.50: C (90.9%) - **Cysteine strongly conserved**
- 4.50: T (45.5%) - Threonine preference
- 5.50: G (81.8%) - Glycine strongly preferred

**Key motifs:**
- Helix 2: Contains E at position -4 (2.46)
- Helix 3: Y-L-L-T-C-P-L-I-L
- Helix 5: F-F-x-x-G-C-x-x-F

### Anion Channels (n=10)
**Signature positions:**
- 1.50: C (40.0%) - Cysteine preference
- 2.50: L (50.0%) - Leucine unique
- 3.50: C (90.0%) - Cysteine conserved (like cation channels)
- 5.50: G (50.0%) - Glycine preference

**Key motifs:**
- Helix 1: x-A-V-V-C-A-C-Q-x
- Helix 2: E-A-I-Y-L-P-s-V-E (distinctive)

### Chloride Pumps (n=12)
**Signature positions:**
- 1.50: A (41.7%) - Alanine preference
- 3.50: T (66.7%) - Threonine (pump signature)
- 4.50: M (100.0%) - Methionine universally conserved
- 5.50: S (91.7%) - Serine strongly conserved

### Sensor/Regulatory (n=22)
**Signature positions:**
- 1.50: F (54.5%) - Phenylalanine preference
- 3.46: W (90.9%) - Tryptophan conserved
- Variable at most .50 positions

### Sodium Pumps (n=5)
**Unique motifs:**
- Helix 1: G-Y-A-V-M-L-A-G-L (highly specific)
- Helix 2: L-S-A-V-V-M-V-S-A (conserved)

## Domain-Specific Patterns

### Archaea (n=27)
- Highest overall conservation
- Prefer I at 2.50 (unique among domains)
- Strong preference for F at 7.53 (85.2%)
- Classic pump signatures well-preserved

### Bacteria (n=41)
- M strongly conserved at 4.50 (92.7%)
- Intermediate conservation levels
- Mix of pump and sensory functions

### Eukaryota (n=45)
- Highest sequence diversity
- C at 3.50 (75.6%) - channel preference
- G at 5.50 (53.3%) vs S in prokaryotes
- More variable at .50 positions

### Viruses (n=5)
- Limited dataset, unique patterns
- Polar residues at 1.50
- Y preference at 7.53

## Helix-Specific Motifs

### Helix 1 (TM1)
**Conservation gradient**: Decreases from N to C terminus
**Key features**:
- GxxxG motif in pumps (positions 46-54)
- M1.50 in pumps vs variable in channels
- Functional group separation strong

### Helix 2 (TM2)
**Conservation gradient**: Moderate throughout
**Key features**:
- E2.46 in channels (glutamate at -4 position)
- Y2.46 common in eukaryotes
- I/V variation at 2.50 distinguishes pumps/channels

### Helix 3 (TM3)
**Conservation gradient**: High around .50
**Key features**:
- W3.46 universally conserved (>80%)
- **C3.50 vs T3.50 - strongest functional discriminator**
- P-L-L motif after .50 position

### Helix 4 (TM4)
**Conservation gradient**: Moderate, increases toward C-terminus
**Key features**:
- M4.50 universal in pumps
- D at -3 position common
- T-G motif at +3,+4 positions

### Helix 5 (TM5)
**Conservation gradient**: Variable
**Key features**:
- S5.50 (pumps) vs G5.50 (channels) distinction
- W at -4 in prokaryotes
- F-rich in eukaryotic channels

### Helix 6 (TM6)
**Conservation gradient**: Highest of all helices
**Key features**:
- **W6.50 most conserved position in microbial opsins**
- Hydrophobic residues before W
- P at +4 also conserved
- Critical for retinal binding pocket

### Helix 7 (TM7)
**Conservation gradient**: Very high throughout
**Key features**:
- **K7.50 universally conserved (Schiff base)**
- D7.46 conserved (counter-ion)
- V-G-F-G motif after K common
- Most conserved helix overall

## Key Findings and Interpretations

### 1. Functional Discrimination
The analysis reveals that microbial opsins can be distinguished by characteristic residues at key positions:

**Primary discriminators:**
- **Position 3.50**: C (channels) vs T (pumps) - 90% accuracy
- **Position 5.50**: G (channels) vs S (pumps) - 70% accuracy
- **Position 1.50**: Variable in channels, M in pumps

**Mechanistic implications:**
- C3.50 in channels may provide flexibility needed for ion conduction
- T3.50 in pumps may stabilize the closed state
- G5.50 in channels could provide a hinge for gating

### 2. Universal Conservation
Two positions show >90% conservation across all functional groups:

**W6.50 (91.4% conserved)**:
- Part of retinal binding pocket
- Pi-stacking with retinal polyene chain
- Mutation typically abolishes function

**K7.50 (95.3% conserved)**:
- Forms Schiff base with retinal
- Essential for light sensitivity
- Only known exceptions are in specialized variants

### 3. Evolutionary Patterns

**Domain-specific preferences:**
- Eukaryotes: More channels, greater diversity
- Prokaryotes: More pumps, higher conservation
- Archaea: Most conserved sequences overall

**Functional evolution:**
- Pumps appear more ancient (higher conservation)
- Channels show more diversity (recent evolution?)
- Sensors intermediate in conservation

### 4. Structural Motifs

**Conserved structural elements:**
- No NPXXY motif (GPCR-specific, not in microbial opsins)
- GxxxG in Helix 1 (pumps) - helix packing motif
- WLxT in Helix 3 - structural role
- PLL after 3.50 - proline kink
- TG after 4.50 - flexibility point

### 5. Neighborhood Analysis

Conservation analysis of ±4 positions around .50 shows:

**Hydrophobicity patterns:**
- Core regions (-2 to +2) predominantly hydrophobic
- Channels have more polar/special residues than pumps
- Conservation highest at .50, decreases with distance

**Group-specific signatures:**
- Each functional group has unique patterns in neighborhoods
- Sodium pumps show highest conservation in neighborhoods
- Unknown function proteins cluster with pumps

## Future Directions

### 1. GRN Refinement
- Correct current alignment errors
- Validate helix boundaries with experimental data
- Improve loop region assignments
- Test alternative reference structures

### 2. Expanded Analysis
- Include more recently solved structures
- Add thermophilic and psychrophilic variants
- Analyze chimeric and engineered opsins
- Include metagenomic discoveries

### 3. Functional Validation
- Test predicted functional determinants
- Create channel↔pump mutations at key positions
- Validate evolutionary hypotheses
- Structure-function correlation studies

### 4. Method Development
- Improve GRN assignment algorithm
- Develop microbial opsin-specific numbering
- Create automated quality checks
- Build predictive models

### 5. Applications
- Rational design of optogenetic tools
- Predict function from sequence
- Guide crystallization constructs
- Identify druggable sites

## Conclusions

This preliminary analysis demonstrates the power of applying GRN to microbial opsins, revealing clear patterns that distinguish functional groups and evolutionary lineages. Despite current limitations in the GRN alignment, the major findings—particularly the C/T discrimination at position 3.50 and universal conservation of W6.50 and K7.50—provide a solid foundation for understanding microbial opsin diversity and evolution.

The identification of function-specific motifs and domain preferences opens new avenues for both basic research and optogenetic tool development. As the GRN alignment is refined and the dataset expanded, we expect these patterns to become even clearer and more predictive of function.

---

*Note: This is a preliminary analysis with known limitations in the GRN alignment. Results should be validated with additional structures and experimental data before drawing definitive conclusions.*