---
title: 'proteostasis-sim: A Multiscale In Silico Framework for Physiological Antagonism, Hormesis, and Competing Risks in Proteinopathies'
tags:
  - Python
  - R
  - systems biology
  - proteostasis
  - hormesis
  - competing risks
  - SERPINA1
authors:
  - name: Lucas Silveira Biondini
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Instituto de Ciências Biológicas, Universidade Federal de Minas Gerais (ICB/UFMG), Belo Horizonte, Brazil
    index: 1
date: 25 November 2026
bibliography: paper.bib
---

# Summary

Late-onset genetic proteinopathies—such as Alpha-1 Antitrypsin Deficiency (*SERPINA1* PiZZ) and neurodegenerative aggregation disorders—exhibit pronounced inter-individual heterogeneity in age of onset and fibrotic progression despite shared causal mutations [@mann2025deep; @wilson2015emergence]. Concurrently, epidemiological studies frequently observe paradoxical reductions in the cumulative incidence of genetic diseases under chronic environmental stressors, raising a fundamental causal dilemma: does an environmental stressor induce genuine biological protection via cellular hormesis (physiological antagonism), or does it merely mask the genetic phenotype through competing-risk attrition?

`proteostasis-sim` is an open-source, modular Python package with automated R (`survival`/`cmprsk`) cross-validation that couples three biological scales bottom-up: (1) a 5-variable nonlinear molecular ODE/SDE proteostasis-bioenergetics model with saddle-node hysteresis and intermittent vs. chronic chronobiology (`wp1_molecular`); (2) a vectorized 2D hybrid agent-based and reaction-diffusion PDE tissue model capturing spatial heterogeneity, bystander cytokine signaling, fibrotic remodeling, and replicative senescence (`wp2_tissue`); and (3) a shared-frailty virtual cohort engine implementing exact Aalen-Johansen cumulative incidence, cause-specific Cox proportional hazards (CSH), exact vectorized Fine-Gray subdistribution hazard ratios via IPCW [@fine1999proportional], and a counterfactual antagonism index ($\Phi$) (`wp3_biostats`). An integrated NCBI GEO calibration module (`omics_geo`) anchors molecular activation parameters to real human hepatic transcriptomic datasets (`GSE66076`, `GSE41666`, `GSE18494`, `GSE31286`).

# Statement of Need

Existing software ecosystems typically treat molecular proteostasis networks, spatial tissue remodeling, and competing-risks survival analysis in isolated silos. Mechanistic ODE/ABM simulators rarely propagate failure-time distributions into competing-risks causal diagnostics, whereas epidemiological tools (`lifelines`, `cmprsk`, `mstate`) lack mechanistic generators capable of separating etiological hazard reduction from competitive masking and survivor-selection bias under shared frailty. `proteostasis-sim` fills this gap by providing a reproducible *in silico* clinical trial and hypothesis-testing platform for systems biologists, biostatisticians, and translational researchers studying protein misfolding disorders, chronopharmacology, and combinatorial drug synergy (Bliss independence).

# Architecture and Validation

The package is organized into four core modules under `src/proteostasis_sim/` exposed through the `InSilicoMultiscalePlatform` facade:

1. **Molecular Scale (`wp1_molecular.py`):** Models misfolded monomers ($M$), toxic aggregates ($A$), proteostasis capacity ($H$), oxidative/environmental damage ($D$), and mitochondrial ATP reserve ($B$) with age-dependent basal decline, oxidative cross-acceleration of nucleation, and chaperone sequestration. Supports stiff ODE integration (`Radau`), stochastic ensemble simulations (SDE), hysteresis continuation, genotype/PRS stratification, LHS-PRCC sensitivity, and Bliss combination screening.
2. **Tissue Scale (`wp2_tissue.py`):** Embeds heterogeneous intracellular dynamics into a 2D hepatocyte lattice coupled to a diffusing pro-fibrotic cytokine field solved via discrete Laplacian convolution, linking tissue failure times to WP3 through right-censored Weibull maximum-likelihood estimation (MLE).
3. **Population Scale (`wp3_biostats.py`):** Simulates virtual cohorts with shared and cause-specific frailty, contrasting hormetic vs. null-masking counterfactual worlds via Cox CSH and exact $O(N \log N)$ Fine-Gray IPCW estimation.
4. **Empirical Omics Calibration (`omics_geo.py`):** Downloads and maps NCBI GEO hepatocyte datasets to HGNC symbols, computes transcription factor regulon scores (`HSF1`, `HIF1A`, `NFE2L2`, `ATF4`), and fits bootstrap-validated Hill and kinetic relaxation curves.

Numerical integrity is enforced via automated continuous integration (`pytest`) covering dynamical positivity, saddle-node bistability, tissue regime transitions, and direct Python-versus-R cross-validation (`scripts/validate_python_vs_r.R`) verifying agreement between Python's Fine-Gray/Aalen-Johansen/Cox CSH estimators and R's gold-standard `cmprsk` and `survival` implementations.

# References
