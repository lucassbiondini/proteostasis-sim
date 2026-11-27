# `proteostasis-sim`

[![CI — Python Unit Tests & R Cross-Validation](https://github.com/lucassbiondini/proteostasis-sim/actions/workflows/ci.yml/badge.svg)](https://github.com/lucassbiondini/proteostasis-sim/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-3776ab.svg)](pyproject.toml)
[![R 4.3+ Cross-Validated](https://img.shields.io/badge/R-4.3%2B%20survival%2Fcmprsk-276DC3.svg)](scripts/validate_python_vs_r.R)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24575804.svg)](https://doi.org/10.5281/zenodo.24575804)

**Multiscale *in silico* framework for separating genuine physiological antagonism (adaptive cellular hormesis) from competing-risks statistical masking in genetic proteinopathies**, empirically calibrated from human hepatic transcriptomics (*SERPINA1* PiZZ archetype).

---

## Overview

`proteostasis-sim` couples three mechanistic scales strictly bottom-up without hardcoded inter-scale shortcuts:

1. **Empirical Omics Anchoring & WP1 (`src/proteostasis_sim/omics_geo.py`, `wp1_molecular.py`):**
   - Queries and scores transcription factor regulons (`HSF1`, `HIF1A`, `NFE2L2`, `ATF4`) across four public human hepatic NCBI GEO datasets (`GSE66076`, `GSE41666`, `GSE18494`, `GSE31286`), anchoring stress-response turnover to $\delta_H = 0.3825\,\mathrm{h}^{-1}$ ($t_{1/2} = 1.81\,\mathrm{h}$).
   - Solves a stiff 5-ODE/SDE intracellular system tracking misfolded monomers ($M$), toxic aggregates ($A$), adaptive proteostasis capacity ($H$), oxidative damage ($D$), and mitochondrial ATP reserve ($B$), capturing saddle-node bistability ($E \in [0.15, 1.55]$), intermittent pulsatile chronobiology (45% duty cycle), global LHS-PRCC sensitivity, and Bliss Independence drug synergy.
2. **WP2 Spatial Tissue Remodeling (`src/proteostasis_sim/wp2_tissue.py`):**
   - Embeds log-normally heterogeneous hepatocytes onto a 2D lattice coupled to an explicit reaction-diffusion profibrotic cytokine PDE, bystander ROS injury, MMP-mediated ECM remodeling, and replicative senescence, fusing single-cell and tissue failure times via right-censored Weibull MLE.
3. **WP3 Shared-Frailty Competing-Risks Causal Audit (`src/proteostasis_sim/wp3_biostats.py`):**
   - Simulates virtual cohorts under **Hormesis** and **Null Masking** counterfactuals, computing exact Aalen–Johansen cumulative incidence (CIF), Cox cause-specific hazards (CSH), exact $O(N\log N)$ Fine–Gray subdistribution hazards (IPCW), and the counterfactual genuine-antagonism index ($\Phi$).

---

## Installation

### Option 1: Full Conda Environment (Python 3.11 + R 4.3 `survival`/`cmprsk`)

```bash
git clone https://github.com/lucassbiondini/proteostasis-sim.git
cd proteostasis-sim
conda env create -f environment.yml
conda activate proteostasis-sim
```

### Option 2: `pip` Installation (Python 3.10+)

```bash
git clone https://github.com/lucassbiondini/proteostasis-sim.git
cd proteostasis-sim
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -e ".[omics,dev]"
```

*(Optional)* For R cross-validation, ensure R $\ge 4.3$ is installed with `survival`, `cmprsk`, and `jsonlite`:

```bash
Rscript -e "install.packages(c('survival', 'cmprsk', 'jsonlite'), repos='https://cloud.r-project.org')"
```

---

## Reproducing Manuscript Results & Figures 1–4

### 1. Run Automated Unit Tests & Python-vs.-R Cross-Validation

```bash
# Run full test suite (WP1, WP2, WP3, GEO regulon parser, and R cross-validation)
pytest --cov=proteostasis_sim --cov-report=term-missing

# Or run the standalone R cross-validation benchmark directly:
Rscript scripts/validate_python_vs_r.R data/processed/crossval_cohort.csv data/processed/python_reference_metrics.json
```

### 2. Fetch & Calibrate from Real NCBI GEO Datasets (`GSE66076`, `GSE41666`, `GSE18494`, `GSE31286`)

```bash
python scripts/run_geo_calibration.py
```

Outputs are written to `data/processed/`:
- `GSE66076_SERPINA1_regulons.csv`
- `GSE41666_HepG2_regulons.csv`
- `GSE18494_HepG2_kinetics.csv` and `GSE18494_kinetic_calibration.json`

### 3. Reproduce Figures 1–4 (Publication PDF + High-DPI PNG)

```bash
# Full publication resolution (exports figures/fig1_*.pdf..fig4_*.pdf and 600 DPI PNGs)
python scripts/reproduce_figures.py --out-dir figures --dpi 600

# Download fresh GEO series before plotting (requires internet access to NCBI GEO)
python scripts/reproduce_figures.py --fetch-geo --out-dir figures --dpi 600

# Fast smoke verification (~15 seconds)
python scripts/reproduce_figures.py --quick --out-dir figures_preview --dpi 150
```

| Manuscript Figure | Generated Files | Underlying Scale / Content |
| :--- | :--- | :--- |
| **Figure 1** | `figures/fig1_empirical_geo.{pdf,png}` | Empirical human hepatic GEO regulon scoring (`GSE66076`, `GSE41666`, `GSE31286`) and $\delta_H = 0.3825\,\mathrm{h}^{-1}$ turnover fit (`GSE18494`) |
| **Figure 2** | `figures/fig2_multiscale_validation.{pdf,png}` | WP1 saddle-node bistability, chronic vs. 45% pulsatile survival ($7.93\times$), `PiMZ`/`PiZZ`/`PRS_H` stratification, and Python-vs.-R verification |
| **Figure 3** | `figures/fig3_wp2_wp3_spatial_causal.{pdf,png}` | WP2 censored Weibull MLE bridge and WP3 shared-frailty causal separation ($HR_{\mathrm{CSH}}$ vs. $sHR_{\mathrm{FG}}$ and $\Phi$) |
| **Figure 4** | `figures/fig4_prcc_bliss_screening.{pdf,png}` | Global LHS-PRCC sensitivity ($n=160$) and genotype-specific Bliss Independence synergy (`HSF1 + Mito`) |

---

## Repository Structure

```text
proteostasis-sim/
├── .github/
│   ├── ISSUE_TEMPLATE/          # Structured bug report & feature request templates
│   └── workflows/ci.yml         # Automated Python + R continuous integration
├── data/processed/              # Benchmark cohorts & calibrated GEO summary tables
├── paper/                       # JOSS software paper (paper.md & paper.bib)
├── scripts/
│   ├── reproduce_figures.py     # End-to-end reproduction of Figures 1-4 (PDF + PNG)
│   ├── run_geo_calibration.py   # NCBI GEO fetch, MyGene.info mapping & kinetic fit
│   └── validate_python_vs_r.R   # Independent R (survival/cmprsk) numerical audit
├── src/proteostasis_sim/
│   ├── omics_geo.py             # GEOparse ingestion, TF regulon z-scores & curve fitting
│   ├── wp1_molecular.py         # 5-ODE/SDE proteostasis-bioenergetics, PRCC & Bliss
│   ├── wp2_tissue.py            # 2D hybrid ABM-PDE tissue lattice & Weibull MLE fusion
│   ├── wp3_biostats.py          # Shared-frailty cohorts, CIF, Cox CSH, Fine-Gray IPCW & Phi
│   └── platform.py              # Unified InSilicoMultiscalePlatform facade
├── tests/                       # Pytest suite covering WP1, WP2, WP3 & GEO parsing
├── CHANGELOG.md
├── CITATION.cff
├── CONTRIBUTING.md
├── environment.yml
├── LICENSE
└── pyproject.toml
```

---

## Contributing & Support

- **Bug reports & feature requests:** Please open an issue via [GitHub Issues](https://github.com/lucassbiondini/proteostasis-sim/issues) using the provided templates.
- **Pull requests:** See [`CONTRIBUTING.md`](CONTRIBUTING.md) for development setup, code formatting (`ruff`), and testing guidelines.

---

## Citation

If you use `proteostasis-sim` in your research, please cite the software archive (see [`CITATION.cff`](CITATION.cff)) and accompanying manuscript:

```bibtex
@software{biondini2026proteostasis_sim,
  author       = {Biondini, Lucas Silveira},
  title        = {{proteostasis-sim: Multiscale in silico framework for physiological antagonism and competing risks in genetic proteinopathies}},
  year         = {2026},
  version      = {0.2.0},
  publisher    = {Zenodo},
  doi          = {10.5281/zenodo.24575804},
  url          = {https://github.com/lucassbiondini/proteostasis-sim}
}
```

## License

Released under the [MIT License](LICENSE).
