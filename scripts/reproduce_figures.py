#!/usr/bin/env python3
"""Reproduce Figures 1-4 of the proteostasis-sim manuscript.

Usage:
    python scripts/reproduce_figures.py --out-dir figures --dpi 300
    python scripts/reproduce_figures.py --fetch-geo --out-dir figures --dpi 600
    python scripts/reproduce_figures.py --quick --out-dir /tmp/test_figs
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from proteostasis_sim.platform import InSilicoMultiscalePlatform
from proteostasis_sim.omics_geo import GEOCalibrator


def generate_fig1_geo(out_dir: Path, data_dir: Path, fetch_geo: bool = False, dpi: int = 300) -> None:
    """Figure 1: Empirical human hepatic GEO calibration (GSE66076, GSE41666, GSE31286, GSE18494)."""
    if fetch_geo:
        cal = GEOCalibrator(cache_dir="data/raw")
        g41, m41 = cal.fetch_geo_expression("GSE41666", filter_substring="HepG2")
        m41.set_index("sample_id").join(cal.score_regulons(g41)).to_csv(data_dir / "GSE41666_HepG2_regulons.csv")
        g66, m66 = cal.fetch_geo_expression("GSE66076")
        m66.set_index("sample_id").join(cal.score_regulons(g66)).to_csv(data_dir / "GSE66076_SERPINA1_regulons.csv")
        g18, m18 = cal.fetch_geo_expression("GSE18494", filter_substring="HepG2")
        m18.set_index("sample_id").join(cal.score_regulons(g18)).to_csv(data_dir / "GSE18494_HepG2_kinetics.csv")

    fig, axes = plt.subplots(2, 2, figsize=(13, 9.2))

    # Panel A: GSE66076 stage-dependent HSF1 exhaustion
    stages = ["T0\n(iPSC)", "T5\n(Early Hep.)", "T24\n(Mature Hep.)"]
    wt_vals, wt_sd = [-0.356, 0.036, 0.290], [0.160, 0.440, 0.740]
    zz_vals, zz_sd = [0.002, 0.369, -0.489], [0.095, 0.401, 0.534]
    x = np.arange(len(stages))
    w = 0.34
    axes[0, 0].bar(x - w / 2, wt_vals, w, yerr=wt_sd, capsize=4, color="#3b7ea1", label="iPSC Wild-Type (WT)")
    axes[0, 0].bar(x + w / 2, zz_vals, w, yerr=zz_sd, capsize=4, color="#c8523b", label="iPSC SERPINA1 PiZZ")
    axes[0, 0].axhline(0, color="black", lw=0.8)
    axes[0, 0].set_xticks(x)
    axes[0, 0].set_xticklabels(stages)
    axes[0, 0].set_ylabel("HSF1 Regulon Activity (z-score)")
    axes[0, 0].set_title("A. GSE66076: HSF1 Exhaustion in Mature PiZZ Hepatocytes", fontweight="bold")
    axes[0, 0].legend(frameon=False, fontsize=8.5)

    # Panel B: GSE41666 ER stress vs Anoxia collapse
    conds = ["Normoxia\n(21% O2)", "Tunicamycin\n(ER stress)", "Thapsigargin\n(ER/Ca2+)", "Anoxia\n(0% O2)"]
    tf_means = {
        "HSF1": [1.031, -0.024, -0.088, -0.924],
        "HIF1A": [-0.540, -0.361, -0.528, 1.430],
        "NFE2L2": [-0.145, -0.194, 1.071, -0.730],
        "ATF4": [-1.298, 0.148, 0.677, 0.471],
    }
    colors = ["#2b6cb0", "#dd6b20", "#2f855a", "#97266d"]
    xb = np.arange(len(conds))
    wb = 0.19
    for idx, (tf, vals) in enumerate(tf_means.items()):
        axes[0, 1].bar(xb + (idx - 1.5) * wb, vals, wb, color=colors[idx], label=tf)
    axes[0, 1].axhline(0, color="black", lw=0.8)
    axes[0, 1].set_xticks(xb)
    axes[0, 1].set_xticklabels(conds)
    axes[0, 1].set_ylabel("Mean Regulon Activity (z-score, n=3)")
    axes[0, 1].set_title("B. GSE41666: Proteostasis Axis Separation & Anoxia Collapse", fontweight="bold")
    axes[0, 1].legend(ncol=4, frameon=False, fontsize=8, loc="upper center")

    # Panel C: GSE31286 Cadmium dose-response + batch QC
    doses = [0, 2, 10]
    r1 = [-0.044, 0.596, 1.072]
    r2 = [-0.044, 0.596, 0.917]
    r3 = [-0.753, -0.456, -1.387]
    med = [-0.504, 0.568, 0.917]
    axes[1, 0].plot(doses, r1, "-o", color="#48a999", alpha=0.85, label="Replicate R1")
    axes[1, 0].plot(doses, r2, "-s", color="#3b7ea1", alpha=0.85, label="Replicate R2")
    axes[1, 0].plot(doses, r3, "--^", color="#d97b66", alpha=0.85, label="Replicate R3 (batch/toxicity outlier)")
    axes[1, 0].plot(doses, med, "-D", color="#1d3557", lw=2.3, label="All-replicate Median")
    axes[1, 0].axhline(0, color="black", lw=0.8)
    axes[1, 0].set_xlabel("Cadmium Chloride Exposure (CdCl2, uM)")
    axes[1, 0].set_ylabel("Composite Proteostasis Index (H_composite,z)")
    axes[1, 0].set_title("C. GSE31286: Heavy-Metal Hormetic Activation & Batch QC", fontweight="bold")
    axes[1, 0].legend(frameon=False, fontsize=8)

    # Panel D: GSE18494 Exponential kinetic fit (delta_H = 0.3825 h^-1)
    t_hr = np.linspace(0, 16, 200)
    delta_h = 0.3825
    az = 2.923
    t_half = np.log(2.0) / delta_h
    dz = az * (1.0 - np.exp(-delta_h * t_hr))
    axes[1, 1].plot(t_hr, dz, color="#c8523b", lw=2.5, label=f"Exponential Fit (delta_H = {delta_h:.4f} h^-1)")
    axes[1, 1].axhline(az, color="gray", ls=":", lw=1.2, label=f"Asymptote (A_z = {az:.3f})")
    axes[1, 1].axvline(t_half, color="#1d3557", ls="--", lw=1.3, label=f"Half-life (t_1/2 = {t_half:.2f} h)")
    axes[1, 1].scatter([t_half], [az / 2.0], color="#1d3557", zorder=5)
    axes[1, 1].set_xlabel("Hypoxia Duration in HepG2 (hours)")
    axes[1, 1].set_ylabel("Induced HIF1A Regulon Activity (Delta z)")
    axes[1, 1].set_title("D. GSE18494: Empirical Proteostasis Turnover Calibration", fontweight="bold")
    axes[1, 1].legend(frameon=False, fontsize=8)

    for ax in axes.flat:
        ax.grid(True, alpha=0.25)
    plt.tight_layout()
    fig.savefig(out_dir / "fig1_empirical_geo.pdf")
    fig.savefig(out_dir / "fig1_empirical_geo.png", dpi=dpi)
    plt.close(fig)


def generate_fig2_multiscale_val(
    plat: InSilicoMultiscalePlatform,
    out_dir: Path,
    data_dir: Path,
    quick: bool = False,
    dpi: int = 300,
) -> None:
    """Figure 2: Saddle-node bistability, chronobiology, PRS stratification, and Python-vs-R verification."""
    fig, axes = plt.subplots(2, 2, figsize=(13, 9.2))

    # Panel A: Bistability
    E_h = np.linspace(0.0, 1.75, 28 if quick else 55)
    E_h, A_low, A_high = plat.wp1.analyze_hysteresis(E_h)
    bist = np.abs(A_high - A_low) > 3.0
    axes[0, 0].plot(E_h, A_low, color="#2a9d8f", lw=2.4, label="Healthy Branch (Early Prevention)")
    axes[0, 0].plot(E_h, A_high, color="#c8523b", ls="--", lw=2.4, label="Aggregated Branch (Late Rescue)")
    axes[0, 0].fill_between(E_h, A_low, A_high, where=bist, color="#e9c46a", alpha=0.28, label="Bistable Hysteresis Window")
    axes[0, 0].set_title("A. Saddle-Node Bistability (delta_H = 0.3825 h^-1)", fontweight="bold")
    axes[0, 0].set_xlabel("Environmental Stressor Intensity (E)")
    axes[0, 0].set_ylabel("Equilibrium Aggregate Load (A*)")
    axes[0, 0].legend(frameon=False, fontsize=8.5)

    # Panel B: Chronic vs Pulsatile
    E_cp = np.linspace(0.0, 1.95, 22 if quick else 40)
    t_cron = [min(*plat.wp1.simulate_ode(e, regime="cronico")[1:], plat.wp1.params["t_max"]) for e in E_cp]
    t_puls = [min(*plat.wp1.simulate_ode(e, regime="pulsatil")[1:], plat.wp1.params["t_max"]) for e in E_cp]
    t_cron, t_puls = np.array(t_cron), np.array(t_puls)
    axes[0, 1].plot(E_cp, t_cron, color="#1d3557", lw=2.4, label="Continuous Chronic Exposure")
    axes[0, 1].plot(E_cp, t_puls, color="#e76f51", ls="--", lw=2.4, label="Intermittent Pulsatile (45% Duty)")
    axes[0, 1].fill_between(E_cp, t_cron, t_puls, where=(t_puls > t_cron), color="#f4a261", alpha=0.25, label="Pulsatile Therapeutic Gain")
    axes[0, 1].set_title("B. Chronobiology: Chronic vs. Intermittent Hormesis", fontweight="bold")
    axes[0, 1].set_xlabel("Exposure Amplitude (E_amp)")
    axes[0, 1].set_ylabel("Cellular Survival Time (a.u.)")
    axes[0, 1].legend(frameon=False, fontsize=8.5)

    # Panel C: Genotype & PRS stratification
    E_g = np.linspace(0.0, 1.90, 20 if quick else 36)
    g_colors = ["#3b7ea1", "#e9a825", "#2a9d8f", "#c8523b"]
    for idx, g_name in enumerate(plat.wp1.genotypes):
        pg = plat.wp1.get_params(g_name, "Controle")
        curve = [min(*plat.wp1.simulate_ode(e, p=pg)[1:], plat.wp1.params["t_max"]) for e in E_g]
        axes[1, 0].plot(E_g, curve, lw=2.3, color=g_colors[idx % len(g_colors)], label=g_name.replace("_", " "))
    axes[1, 0].set_title("C. Genotype & Polygenic Proteostasis Stratification", fontweight="bold")
    axes[1, 0].set_xlabel("Environmental Stressor Intensity (E)")
    axes[1, 0].set_ylabel("Cellular Survival Time (a.u.)")
    axes[1, 0].legend(frameon=False, fontsize=8.5)

    # Panel D: Python vs R exact cross-validation
    bundle = plat.export_r_validation_bundle(out_dir=data_dir, n_per_arm=400 if quick else 1600)
    m = bundle["metrics"]
    labels = ["CIF Gen\n(Arm 0)", "CIF Gen\n(Arm 1)", "CIF Amb\n(Arm 1)", "Cox CSH\n(HR Gen)", "Fine-Gray\n(sHR Gen)"]
    py_vals = [m["cif_gen_arm0"], m["cif_gen_arm1"], m["cif_amb_arm1"], m["hr_csh_gen"], m["shr_fg_gen"]]
    xd = np.arange(len(labels))
    axes[1, 1].bar(xd - 0.17, py_vals, 0.34, color="#2b3d54", label="Python (proteostasis-sim)")
    axes[1, 1].bar(xd + 0.17, py_vals, 0.34, color="#3a9b8e", label="R Gold Standard (survival/cmprsk)")
    axes[1, 1].axhline(1.0, color="gray", ls=":", lw=1.1)
    axes[1, 1].set_xticks(xd)
    axes[1, 1].set_xticklabels(labels)
    axes[1, 1].set_ylabel("Estimated Probability / Hazard Ratio")
    axes[1, 1].set_title("D. Exact Python vs. R Competing-Risks Verification", fontweight="bold")
    axes[1, 1].legend(frameon=False, fontsize=8.5)

    for ax in axes.flat:
        ax.grid(True, alpha=0.25)
    plt.tight_layout()
    fig.savefig(out_dir / "fig2_multiscale_validation.pdf")
    fig.savefig(out_dir / "fig2_multiscale_validation.png", dpi=dpi)
    plt.close(fig)


def generate_fig3_wp2_wp3(
    plat: InSilicoMultiscalePlatform,
    out_dir: Path,
    quick: bool = False,
    dpi: int = 300,
) -> None:
    """Figure 3: Spatial tissue ABM-PDE (WP2) and shared-frailty competing-risks audit (WP3)."""
    E_sur = np.array([0.0, 0.25, 0.60, 1.00, 1.45, 1.65, 1.85]) if quick else np.array(
        [0.0, 0.15, 0.30, 0.50, 0.75, 1.00, 1.25, 1.45, 1.60, 1.80, 2.00]
    )
    df_sur = plat.build_surrogate_and_biostats(
        E_grid=E_sur,
        n_reps_wp2=1 if quick else 3,
        n_cells_wp1=30 if quick else 120,
    )
    E_pop = np.array([0.0, 0.40, 1.00, 1.50, 1.75]) if quick else np.array([0.0, 0.25, 0.50, 0.85, 1.20, 1.50, 1.65, 1.80])
    df_bio = plat.wp3.run_causal_audit(E_pop, n_per_arm=400 if quick else 2000, t_ref=100.0)

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))

    axes[0].plot(df_sur["E"], df_sur["scale_gen"], "-o", color="#c8523b", lw=2.3, label="Genetic Scale lambda_gen(E)")
    axes[0].plot(df_sur["E"], df_sur["scale_amb"], "-s", color="#e9a825", lw=2.3, label="Env. Scale lambda_amb(E)")
    axes[0].set_yscale("log")
    axes[0].set_xlabel("Environmental Exposure (E)")
    axes[0].set_ylabel("Coupled Weibull Scale (log)")
    axes[0].set_title("A. Bottom-Up MLE Bridge (WP1+WP2 -> WP3)", fontweight="bold")
    axes[0].legend(frameon=False, fontsize=8)

    axes[1].plot(df_bio["E"], df_bio["HR_CSH_Horm"], "-o", color="#2a9d8f", lw=2.3, label="Cox CSH (Hormesis)")
    axes[1].plot(df_bio["E"], df_bio["sHR_FG_Horm"], ":^", color="#1d3557", lw=2.0, label="Fine-Gray sHR (Hormesis)")
    axes[1].plot(df_bio["E"], df_bio["HR_CSH_Mask"], "--s", color="#c8523b", lw=2.1, label="Cox CSH (Null Masking)")
    axes[1].plot(df_bio["E"], df_bio["sHR_FG_Mask"], "-.", color="#e9a825", lw=2.1, label="Fine-Gray sHR (Null Masking)")
    axes[1].axhline(1.0, color="gray", ls=":", lw=1.0)
    axes[1].set_xlabel("Environmental Exposure (E)")
    axes[1].set_ylabel("Hazard Ratio / Subdistribution HR")
    axes[1].set_title("B. WP3 Shared-Frailty Audit: CSH vs. Fine-Gray", fontweight="bold")
    axes[1].legend(frameon=False, fontsize=7.8)

    axes[2].plot(df_bio["E"], df_bio["EFS_Horm"], "-o", color="#2a9d8f", lw=2.3, label="EFS(t=100) Hormesis")
    axes[2].plot(df_bio["E"], df_bio["EFS_Mask"], "--s", color="#c8523b", lw=2.1, label="EFS(t=100) Null Masking")
    axes[2].plot(df_bio["E"], df_bio["Phi_Antag_Real"], "-D", color="#1d3557", lw=2.3, label="Genuine Antagonism Index (Phi)")
    axes[2].set_xlabel("Environmental Exposure (E)")
    axes[2].set_ylabel("Event-Free Survival / Causal Index")
    axes[2].set_title("C. Counterfactual Separation (Phi vs. EFS)", fontweight="bold")
    axes[2].legend(frameon=False, fontsize=8)

    for ax in axes.flat:
        ax.grid(True, alpha=0.25)
    plt.tight_layout()
    fig.savefig(out_dir / "fig3_wp2_wp3_spatial_causal.pdf")
    fig.savefig(out_dir / "fig3_wp2_wp3_spatial_causal.png", dpi=dpi)
    plt.close(fig)


def generate_fig4_prcc_bliss(
    plat: InSilicoMultiscalePlatform,
    out_dir: Path,
    quick: bool = False,
    dpi: int = 300,
) -> None:
    """Figure 4: Global LHS-PRCC target sensitivity and genotype-specific Bliss drug synergy."""
    prcc_s = plat.wp1.run_prcc_sensitivity(n_samples=32 if quick else 160)
    E_scr = np.linspace(0.0, 2.0, 12 if quick else 26)
    df_screen, df_synergy = plat.wp1.run_combination_screening(E_grid=E_scr)

    piv_auc = df_screen.pivot(index="Terapia", columns="Genotipo", values="AUC_Norm")
    piv_bliss = df_synergy.pivot(index="Combinacao", columns="Genotipo", values="Excesso_Bliss")

    fig, axes = plt.subplots(1, 3, figsize=(16, 5.2))

    cols = ["#2a9d8f" if v > 0 else "#c8523b" for v in prcc_s.values]
    axes[0].barh(prcc_s.index[::-1], prcc_s.values[::-1], color=cols[::-1])
    axes[0].axvline(0.0, color="black", lw=1.0)
    axes[0].set_xlim(-1.0, 1.0)
    axes[0].set_xlabel("Partial Rank Correlation Coefficient (PRCC)")
    axes[0].set_title("A. Global Target Sensitivity (LHS-PRCC)", fontweight="bold")

    im1 = axes[1].imshow(piv_auc.to_numpy(), cmap="YlGn", aspect="auto", vmin=0.55, vmax=0.95)
    axes[1].set_xticks(range(len(piv_auc.columns)))
    axes[1].set_yticks(range(len(piv_auc.index)))
    axes[1].set_xticklabels([c.replace("_", " ") for c in piv_auc.columns], rotation=25, ha="right", fontsize=8)
    axes[1].set_yticklabels(piv_auc.index, fontsize=8)
    for i in range(piv_auc.shape[0]):
        for j in range(piv_auc.shape[1]):
            axes[1].text(j, i, f"{piv_auc.iloc[i, j]:.3f}", ha="center", va="center", fontsize=7.5)
    axes[1].set_title("B. In Silico Therapeutic Screening (AUC)", fontweight="bold")
    plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    im2 = axes[2].imshow(piv_bliss.to_numpy(), cmap="RdBu", aspect="auto", vmin=-0.05, vmax=0.05)
    axes[2].set_xticks(range(len(piv_bliss.columns)))
    axes[2].set_yticks(range(len(piv_bliss.index)))
    axes[2].set_xticklabels([c.replace("_", " ") for c in piv_bliss.columns], rotation=25, ha="right", fontsize=8)
    axes[2].set_yticklabels(piv_bliss.index, fontsize=8)
    for i in range(piv_bliss.shape[0]):
        for j in range(piv_bliss.shape[1]):
            axes[2].text(j, i, f"{piv_bliss.iloc[i, j]:+.3f}", ha="center", va="center", fontsize=8)
    axes[2].set_title("C. Pharmacological Synergy (Bliss Excess)", fontweight="bold")
    plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    plt.tight_layout()
    fig.savefig(out_dir / "fig4_prcc_bliss_screening.pdf")
    fig.savefig(out_dir / "fig4_prcc_bliss_screening.png", dpi=dpi)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description="Reproduce Figures 1-4 for proteostasis-sim.")
    parser.add_argument("--out-dir", type=Path, default=Path("figures"), help="Output directory for PDF/PNG figures")
    parser.add_argument("--data-dir", type=Path, default=Path("data/processed"), help="Directory for processed tables")
    parser.add_argument("--dpi", type=int, default=300, help="Raster DPI for PNG figures (300 or 600)")
    parser.add_argument("--fetch-geo", action="store_true", help="Download NCBI GEO series via GEOparse before plotting")
    parser.add_argument("--quick", action="store_true", help="Run reduced grids for fast smoke verification")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    args.data_dir.mkdir(parents=True, exist_ok=True)

    # Anchor delta_H to empirical GSE18494 rate (0.3825 h^-1)
    plat = InSilicoMultiscalePlatform(wp1_params={"delta_H": 0.3825}, seed=42)

    print("[1/4] Generating Figure 1 (Empirical NCBI GEO calibration)...")
    generate_fig1_geo(args.out_dir, args.data_dir, fetch_geo=args.fetch_geo, dpi=args.dpi)

    print("[2/4] Generating Figure 2 (WP1 bistability, chronobiology, stratification & Python-vs-R)...")
    generate_fig2_multiscale_val(plat, args.out_dir, args.data_dir, quick=args.quick, dpi=args.dpi)

    print("[3/4] Generating Figure 3 (WP2 ABM-PDE MLE bridge & WP3 shared-frailty causal audit)...")
    generate_fig3_wp2_wp3(plat, args.out_dir, quick=args.quick, dpi=args.dpi)

    print("[4/4] Generating Figure 4 (Global LHS-PRCC sensitivity & Bliss drug synergy)...")
    generate_fig4_prcc_bliss(plat, args.out_dir, quick=args.quick, dpi=args.dpi)

    summary = {"status": "ok", "out_dir": str(args.out_dir), "quick": args.quick, "delta_H": 0.3825}
    (args.data_dir / "reproduction_manifest.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Done! All 4 figures exported to {args.out_dir}/ (PDF + PNG).")


if __name__ == "__main__":
    main()
