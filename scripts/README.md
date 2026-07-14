# Maintenance and audit scripts

Scripts in this directory are command-line entry points, so most are not imported by the
core package. They are retained only when they reproduce a current table, input, audit,
or figure.

## Canonical reference maintenance

- `apply_visual_grn_curations.py` applies the approved PsChR2 and HulaCCR1 register
  corrections with strict before/after checks.
- `build_tara_reference_rows.py` builds TARA_A/TARA_B rows, validates the Protos
  round-trip, and synchronizes MOGRN runtime copies.
- `audit_tara_tm5_consensus.py` evaluates candidate TARA_A TM5 registers.
- `audit_type_i_nearest_neighbors.py` audits curated GRNs against structurally aligned
  neighbours.

## Tandem-rhodopsin preprocessing

- `preprocess_tandem_rhodopsins.py` exposes the sequence-level tandem adapter without
  importing Protos.
- `detect_annotate_dual_rhodopsins.py` runs tandem detection and annotation diagnostics.

The production structure split is integrated into `prepare_data.py` through
`src/tandem_structure_preprocessing.py` and configured in
`src/resources/tandem_structure_domains.json`.

## Analysis and quality control

- `analyze_grn_distance.py`: retinal-distance tables for GRN residues.
- `analyze_rmsd.py`: per-structure summaries from the RMSD matrix.
- `detect_misalignment.py`: coordinate outlier diagnostics for the interactive alignment.
- `generate_helices_grn.py`: helix boundaries from curated GRN intervals.
- `analyze_helices_phipsi.py` and `generate_helices_json.py`: optional structural
  helix-boundary diagnostics.
- `postprocess_grn_table_v2.py`: historical full-residue table reconstruction retained
  for reproducibility; it is not run by the standard pipeline.
- `visualize_aligned_helices.py`: optional helix-boundary inspection figure.

## Prediction-input helpers

- `generate_mo_exp_contact_yaml.py` creates experimental-structure Boltz inputs with
  ligand-pocket constraints.
- `new_opsin_yaml/generate_new_opsin_yaml.py` creates sequence/prediction inputs for new
  opsins.

All outputs belong under ignored runtime directories such as `opsin_output/`, `property/`,
or `yaml_configs/`; scripts must not write generated artifacts into the tracked source
tree.
