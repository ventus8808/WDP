<!-- Guidance for AI coding agents working on the WDP repository -->
# WDP — Copilot instructions (concise)

Be brief and make minimal, well-scoped changes. Prefer edits that follow existing project conventions (use `config.yaml` for paths, write outputs to `Data/Processed/` or `Result/`), and preserve reproducibility.

Key concepts — high level
- Project root contains `config.yaml` (single source of path configuration). All scripts read paths from it. Do not hard-code absolute paths; use the keys under `data_sources`, `data_directories`, and `result_directories`.
- Primary workflows:
  - Data ingest/download: `Code/Download/*` (Python)
  - Data cleaning: `Code/Clean/*` (Python). Example: `Code/Clean/County_Adjacency.py` reads `config.yaml` and writes CSVs into `Data/Processed`.
  - Covariate generation: `Code/PCA/*` (Python) produces `Data/Processed/PCA/Master_Covariates.csv`.
  - Modeling/analysis: `Code/INLA/*` (R + shell) and `Code/Test/*` (Python legacy). INLA workflows use SLURM submit scripts and assume R packages (INLA, dplyr) are installed.

Developer workflows & commands (explicit)
- Run a data-cleaning script from project root:
  - python Code/Clean/County_Adjacency.py
- Run PCA workflow:
  - python Code/PCA/PCA_Systematic.py
- Run Bayesian Python runner (legacy):
  - python Code/Test/BSTM_Run.py --aggregated
- INLA single-compound tests (preferred on cluster):
  - sbatch Code/INLA/submit_single_compound_test.sh
  - For local debugging: bash Code/INLA/run_single_compound_test.sh "5" "Atrazine" "C81-C96"

Project-specific conventions and patterns
- Centralized configuration: always read `config.yaml` from project root (scripts locate it using relative parents). Example: `project_root = Path(__file__).resolve().parents[2]` then `yaml.safe_load(config_path)`.
- Output locations: follow `config.yaml` `processed` and `result_directories.filter` values (e.g., `Data/Processed/Pesticide`, `Result/Filter`).
- Filenames: many scripts emit timestamped CSVs (Result/Filter/Results_*.csv). When adding code that consumes these outputs, match existing filename patterns.
- R/SLURM in INLA: INLA scripts expect system R with INLA package; SLURM submit scripts (e.g., `submit_single_compound_test.sh`) build sbatch array jobs and call `run_single_compound_test.sh`.

Integration points & external dependencies
- External data sources referenced in `config.yaml`: CDC WONDER, USGS PNSP, NLDAS, GEE, CACES LUR. Download scripts in `Code/Download/` follow these URLs.
- R packages required for INLA runs: at minimum `INLA` and common tidy tools. Python dependencies include `pandas`, `geopandas`, `pyyaml`.

When editing or adding code
- Follow existing style: small, single-purpose scripts. Prefer adding CLI args or reading `config.yaml` rather than changing global path resolution.
- Tests: there are no formal unit tests. If you add logic, include a small smoke-check (e.g., a short `if __name__ == '__main__':` demo) and document expected input/output files in the docstring.
- Logging: scripts use simple prints; maintain this pattern unless adding a broader logging facility.

Files and places to inspect for context (examples)
- `config.yaml` — canonical paths and data source URLs
- `Code/Clean/County_Adjacency.py` — example pattern reading config and writing CSV
- `Code/INLA/README_Single_Compound_Test.md` — cluster usage, parameters, expected outputs
- `Code/INLA/run_single_compound_test.sh` & `submit_single_compound_test.sh` — how R scripts are invoked under SLURM
- `Code/PCA/PCA_Systematic.py` — covariate generation pipeline

If unsure, prefer non-disruptive changes and ask the repo owner for permission to modify major analysis scripts. When adding or changing outputs, mirror the directory structure under `Data/Processed/` and `Result/`.

If you update this file: keep it short (20–50 lines), concrete, and reference the exact filenames above.
