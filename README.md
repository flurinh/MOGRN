<h1 align="center">MOGRN</h1>

<p align="center">
  <b>Microbial Opsin Generic Residue Numbering</b><br>
  A structure-guided coordinate system for comparing type-I opsins residue by residue.
</p>

<p align="center"><img src="docs/grn-positions.jpg" alt="Generic residue positions on bacteriorhodopsin" width="380"></p>

## Overview

MOGRN has two connected roles. First, it explains how the microbial-opsin GRN
system is anchored and runs a reproducible structure-alignment workflow that
produces a raw, uncurated baseline. Manual scientific curation happens between
that baseline and the clean ProtOS reference table. Second, it demonstrates the
finished GRN system by applying that curated table to structures for analysis and
archival release.

The canonical 130-entity reference input is
`protos/src/protos/reference_data/grn/reference/type_I_opsins.csv` from the pulled
ProtOS `master` branch. Tandem TARA/bestrhodopsin parents are split during
preprocessing into ordinary A/B entities before the standard pipeline runs.

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

Runtime data are intentionally not committed. Restore an older complete runtime
snapshot, when needed, with:

```bash
python download_data.py <ZENODO_RECORD_ID>
```

That runtime snapshot must provide at least:

- `property/mo_exp_ST1.csv` and `property/helices_grn.json`;
- source structures and predictions under `structures/`;
- any Protos runtime data required by the selected release.

Runtime copies under `data/grn/reference/` and `opsin_output/grn_reference.csv`
are derived byte-for-byte from the authoritative ProtOS table.

Release maintainers build the canonical ground-truth deposit only after the clean
workflow has completed:

```bash
python build_groundtruth_bundle.py
```

The builder normalizes `mo_exp_ST5_HEK1.xlsx`, copies the 130-row ProtOS reference
and GRN config, verifies the function-card CSV against its V2 workbook, rewrites all
130 exported structures' `grn` columns from the ProtOS table, and validates every
identifier and residue join. It writes
`zenodo_upload/grn_opsins_groundtruth.zip` only when all release checks pass.
Cell-level validation reports are written beside the staging directory when a
conflict blocks release.

## Reproduce the analysis

```bash
# Register structures, split configured tandem parents, and build datasets
python prepare_data.py --rebuild

# Run structure comparison and raw, provisional GRN assignment without caches.
python opsin_analysis_workflow.py --skip-prepare --no-cache

# After manual curation has been incorporated into the pulled ProtOS table,
# clear provisional labels and persist that table as the final annotation.
python apply_curated_grns.py

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

Final manual corrections live in the ProtOS `type_I_opsins.csv` table. MOGRN's
maintenance scripts audit the TARA/register behavior but cannot replace that source.

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
protos/                      ProtOS checkout and canonical type-I GRN table
prepare_data.py              Structure registration and dataset preparation
opsin_analysis_workflow.py   Main analysis pipeline
apply_curated_grns.py        Separate curated ProtOS annotation application
build_groundtruth_bundle.py  Strict canonical Zenodo bundle builder
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
