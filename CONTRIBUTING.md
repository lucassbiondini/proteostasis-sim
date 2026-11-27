# Contributing to `proteostasis-sim`

Thank you for your interest in contributing to `proteostasis-sim`! Contributions are welcome across model extensions, empirical transcriptomic/proteomic calibration pipelines, numerical solvers, documentation, and test coverage.

## Reporting Issues

Before opening a new issue, please search existing [GitHub Issues](https://github.com/lucassbiondini/proteostasis-sim/issues) to avoid duplicates. Use the structured issue templates:
- **Bug Report:** Include Python/R versions, operating system, minimal reproducible snippet, and full traceback.
- **Feature Request / Model Extension:** Describe the biological or statistical motivation, mathematical formulation (if applicable), and expected API behavior.

## Development Setup

1. Fork the repository and clone your fork locally:
   ```bash
   git clone https://github.com/<YOUR_USER>/proteostasis-sim.git
   cd proteostasis-sim
   ```
2. Create the Conda environment (recommended if running R cross-validation) or a Python virtual environment:
   ```bash
   conda env create -f environment.yml
   conda activate proteostasis-sim
   # OR:
   pip install -e ".[omics,dev]"
   ```
3. Create a feature or bugfix branch from `main`:
   ```bash
   git checkout -b feat/my-extension
   ```

## Coding & Verification Standards

Before submitting a Pull Request, please ensure all checks pass locally:

1. **Linting & formatting (`ruff`):**
   ```bash
   ruff check src tests scripts
   ```
2. **Unit tests & coverage (`pytest`):**
   ```bash
   pytest --cov=proteostasis_sim --cov-report=term-missing
   ```
3. **R biostatistical cross-validation (when touching `wp3_biostats.py`):**
   ```bash
   Rscript scripts/validate_python_vs_r.R data/processed/crossval_cohort.csv data/processed/python_reference_metrics.json
   ```
4. **Figure reproduction smoke check (when touching plotting or default parameters):**
   ```bash
   python scripts/reproduce_figures.py --quick --out-dir /tmp/proteostasis_figs_check
   ```

## Pull Request Process

- Keep PRs focused on a single logical change and update docstrings, `README.md`, and `CHANGELOG.md` (under an `[Unreleased]` section) when adding user-facing functionality.
- Ensure GitHub Actions CI passes on your PR branch before requesting review.
- By contributing, you agree that your contributions will be licensed under the repository's [MIT License](LICENSE).
