<h1 align="center">MOGRN</h1>

<p align="center">
  <b>Microbial Opsin Generic Residue Numbering</b><br>
  A structure-guided coordinate system for comparing type-I opsins residue by residue.
</p>

<p align="center"><img src="docs/grn-positions.jpg" alt="Generic residue positions on bacteriorhodopsin" width="380"></p>

## Overview

MOGRN aligns experimental and predicted microbial-opsin structures, assigns Generic
Residue Numbers (GRNs), compares transmembrane bundles, and produces conservation,
retinal-pocket, RMSD, and interactive structural analyses. The repository contains the
analysis code and the hand-curated type-I reference table; large structures, metadata,
caches, and generated figures are distributed separately.

The canonical reference input is [`type_I.csv`](type_I.csv). It currently contains 130
single-domain entities. Tandem TARA/bestrhodopsin parents are split during preprocessing
into ordinary A/B entities before the standard pipeline runs.

## Installation

MOGRN targets Python 3.10 or newer and uses
[ProtOS](https://github.com/flurinh/protos) for structure and GRN storage.

```bash
git clone https://github.com/flurinh/protos.git ../protos
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ../protos
python -m pip install -r requirements.txt
```

For tests, also install `requirements-dev.txt`.

## Data

Runtime data are intentionally not committed. Restore the associated Zenodo archive
before running the full workflow:

```bash
python download_data.py <ZENODO_RECORD_ID>
```

The archive must provide at least:

- `property/mo_exp_ST1.csv` and `property/helices_grn.json`;
- source structures and predictions under `structures/`;
- any Protos runtime data required by the selected release.

`type_I.csv` is the tracked, hand-curated GRN reference. Generated copies under
`data/grn/reference/` and `opsin_output/` must be treated as derived artifacts.

Release maintainers can build a deterministic archive with an internal file manifest
and an external SHA-256 checksum:

```bash
python create_zenodo_archive.py
```

## Reproduce the analysis

```bash
# Register structures, split configured tandem parents, and build datasets
python prepare_data.py --rebuild

# Run structure comparison and GRN assignment
python opsin_analysis_workflow.py

# Validate and synchronize TARA/manual reference rows
python scripts/build_tara_reference_rows.py

# Generate static and interactive figures
python plot.py --input-dir opsin_output --output-dir opsin_output/paper_figures
```

The standard interactive outputs are:

- `opsin_output/paper_figures/interactive_grn_alignment.html`;
- `opsin_output/paper_figures/interactive_grn_alignment_b.html`.

These files are large, reproducible outputs and are not version-controlled.

## Curated tandem and register behavior

- Structures shorter than 500 residues are never considered configured tandem inputs.
- Known TARA parents are split before dataset construction, with local residue numbering
  and parent-coordinate provenance retained.
- Retinal is assigned to the nearest split domain.
- TARA_A retains the real `5.451` insertion while its internal `5.43` gap is closed.
- HulaCCR1 H5 is continuous and has no `5.451` insertion.
- PsChR2 H1 is continuous from `T2=1.32` through `W33=1.63`.

Manual corrections are encoded by `scripts/apply_visual_grn_curations.py` and guarded by
regression tests; they are not hidden changes to ProtOS.

## Tests

```bash
python -m pip install -r requirements-dev.txt
python -m pytest -q
```

The focused scientific-integrity suite covers tandem sequence/structure splitting,
TARA reference round-trips, and the curated GRN registers.

## Repository layout

```text
src/                         Core analysis and preprocessing modules
src/resources/               Curated tandem-domain configuration
scripts/                     Reproducible audits and maintenance tools
type_I.csv                   Canonical 130-entity type-I GRN table
prepare_data.py              Structure registration and dataset preparation
opsin_analysis_workflow.py   Main analysis pipeline
plot.py                      Static and interactive visualization entry point
docs/                        Small publication-facing documentation assets
```

Large inputs and all generated outputs are excluded by `.gitignore`. The historical
clustered-MO experiment and stale generated figures were removed from the source tree;
they remain recoverable from Git history and should be deposited as data only if they are
part of a cited analysis.

## Citation and release status

Citation metadata are provided in [`CITATION.cff`](CITATION.cff). The manuscript
*A Generic Residue-Numbering System for Microbial Rhodopsins* is under review.

Before an archival release, complete the remaining items in
[`docs/PUBLICATION_CHECKLIST.md`](docs/PUBLICATION_CHECKLIST.md), notably the software
license, immutable Zenodo record, metadata-table freeze, and tested ProtOS revision.

## Related projects

- [ProtOS](https://github.com/flurinh/protos): structure and GRN framework
- [Lambda](https://github.com/flurinh/lambda): downstream opsin-colour modelling
