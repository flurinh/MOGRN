# Repository Guidelines

## Project Structure & Module Organization
Core Python modules live in `src/`, organized by workflow stage (e.g., `data_processing.py` for IO prep, `structure_comparison.py` for RMSD metrics, `helix_analysis.py` for TM labeling). Top-level scripts such as `prepare_data.py`, `opsin_analysis_workflow.py`, `analyze_grns.py`, and `analyze_motifs.py` orchestrate standard pipeline steps. Domain resources sit under `property/` and `structures/`, while intermediate configs land in `yaml_configs/` and results collect under `opsin_output/`. Keep experimental notebooks and papers inside `paper/` and stash large cached artefacts under `data_clustered_mo/`.

## Build, Test, and Development Commands
Create a Python 3.10 environment and install dependencies with `pip install -r requirements.txt` plus `pip install -e protos`. Run `python prepare_data.py` to seed directory scaffolding, `python prepare_yaml.py` to emit per-protein configs, and `python opsin_analysis_workflow.py` for the full analysis. Visualization and post-processing live in `python plot.py`, `python analyze_grns.py`, and `python analyze_motifs.py`. Sanity-check module imports with `python -m pytest src/test_imports.py -s`; adjust the `PROTOS_DATA_ROOT` env var per `src/test_protos.py` before exercising Protos-dependent flows.

## Coding Style & Naming Conventions
The codebase follows PEP 8 with four-space indentation, module-level docstrings, and snake_case for functions, variables, and filenames (`assign_grns.py`). Prefer descriptive helper names and keep public APIs concentrated in `src/` modules. Use type hints where practical, mirror existing docstring format, and keep plotting color dictionaries in `opsin_color_scheme.py`. Avoid committing notebook output or large data; treat anything under `opsin_output/` as generated artefacts.

## Testing Guidelines
Lightweight checks rely on pytest-compatible scripts in `src/`. Add focused tests beside the module under test (e.g., `src/test_structure_comparison.py`) and run with `python -m pytest src -k <module>`. When adding analysis steps, include command-line smoke tests that read from `property/mo_exp.csv` and record expected figures in `opsin_output/figures`. Document any data prerequisites in README-adjacent files so collaborators can recreate fixtures.

## Commit & Pull Request Guidelines
Recent history uses short imperative subjects (`update`, `xO`); continue the imperative mood but expand summaries to describe scope (e.g., `Refine helix alignment cache`). Squash noisy intermediate commits before opening a PR. Each PR should state the pipeline stage touched, data assumptions, and include before/after artifact paths or screenshots when altering figures. Link related issues and flag any downstream scripts that must be rerun after the change.
