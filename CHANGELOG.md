# Changelog

All notable changes to `proteostasis-sim` will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] - 2026-11-26

### Added
- **Empirical NCBI GEO calibration pipeline (`proteostasis_sim.omics_geo`):** Automated ingestion of `GSE66076`, `GSE41666`, `GSE18494`, and `GSE31286` via `GEOparse`, median probe-to-HGNC symbol collapse, MyGene.info v3 reverse mapping for custom array `GPL9419`, weighted transcription factor regulon scoring (`HSF1`, `HIF1A`, `NFE2L2`, `ATF4`), and nonlinear least-squares turnover fitting ($\delta_H = 0.3825\,\mathrm{h}^{-1}$, $t_{1/2} = 1.81\,\mathrm{h}$).
- **End-to-end figure reproduction CLI (`scripts/reproduce_figures.py`):** Generates publication-ready vector PDF and high-DPI PNG files for Figures 1–4 (`--quick` and `--fetch-geo` modes).
- **Exact Python-vs.-R cross-validation (`scripts/validate_python_vs_r.R`):** Automated numerical verification of Aalen–Johansen CIF ($\le 2.00 \times 10^{-15}$), Cox cause-specific hazards ($\le 1.37 \times 10^{-7}$), and vectorized $O(N\log N)$ Fine–Gray IPCW subdistribution hazards ($\le 1.79 \times 10^{-9}$) against R 4.3 `survival` and `cmprsk`.
- **Repository governance & JOSS readiness:** Added `CITATION.cff`, `CHANGELOG.md`, `CONTRIBUTING.md`, structured `.github/ISSUE_TEMPLATE/` forms, and JOSS paper draft (`paper/paper.md`, `paper/paper.bib`).

### Changed
- Anchored default WP1 stress-response turnover parameter `delta_H` from initial exploratory value (`0.42`) to empirically calibrated human hepatic turnover (`0.3825 h^-1`).
- Unified pulsatile chronobiology benchmark (45% duty cycle, period $T = 8.0$) reporting $7.93\times$ survival gain at $E = 1.75$.

### Fixed
- Resolved GPL9419 reporter-to-HGNC mapping gaps in `GSE18494` via batch REST resolution.
- Enforced explicit CFL stability checks ($D_C \Delta t / \Delta x^2 = 0.04 \ll 0.25$) in the 2D WP2 cytokine diffusion solver.

## [0.1.0] - 2026-11-18

### Added
- Initial release of the three-scale `proteostasis-sim` architecture:
  - WP1: 5-ODE/SDE nonlinear intracellular proteostasis-bioenergetics ($M, A, H, D, B$), saddle-node hysteresis solver, LHS-PRCC sensitivity, and Bliss synergy screening.
  - WP2: 2D hybrid ABM–PDE hepatocyte lattice with log-normal parameter heterogeneity, bystander cytokine signaling, MMP fibrosis remodeling, and censored Weibull MLE bottom-up surrogate coupling.
  - WP3: Shared-frailty competing-risks cohort simulator and counterfactual genuine-antagonism index ($\Phi$).
