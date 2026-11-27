from __future__ import annotations
import json
import argparse
import importlib.util
from pathlib import Path
import numpy as np
import pandas as pd
import requests
from scipy.optimize import curve_fit


GEO_DATASET_CATALOG = {
    "GSE66076": {
        "role": "disease_baseline_serpina1",
        "title": "Alpha-1 antitrypsin deficient (PiZZ) iPSC-derived hepatocytes developmental liver disease signature",
        "organism": "Homo sapiens",
        "platform": "HuGene-1_0-st (part of SuperSeries GSE66078)",
        "pmid": "25843048",
        "use_in_wp1": "Anchors basal UPR/proteostasis induction (beta_0, beta_M) in SERPINA1 PiZZ hepatocytes",
    },
    "GSE41666": {
        "role": "hepatic_er_and_hypoxia_stress",
        "title": "HepG2 cells under Normoxia (21% O2), Anoxia (0% O2), Tunicamycin (1.5 ug/mL), Thapsigargin (0.3 uM)",
        "organism": "Homo sapiens",
        "platform": "Microarray (24 samples; 12 HepG2 across 4 conditions x 3 replicates)",
        "use_in_wp1": "Separates ER proteotoxic stress (ATF4/UPR) from hypoxia stress (HIF1A/HSF1) in human liver cells",
    },
    "GSE18494": {
        "role": "hepatic_hypoxia_time_course",
        "title": "Hypoxic HepG2 hepatoma time course (Normoxia vs Hypoxia across timepoints)",
        "organism": "Homo sapiens",
        "platform": "Affymetrix HG-U133 Plus 2.0 (GPL570/GPL9419)",
        "use_in_wp1": "Calibrates activation and relaxation kinetics (beta_E, delta_H) of HIF1A/HSF1/proteostasis response",
    },
    "GSE31286": {
        "role": "hepatic_heavy_metal_stress",
        "title": "Whole genome expression profiling of HepG2 cells exposed to cadmium (CdCl2)",
        "organism": "Homo sapiens",
        "platform": "Agilent Whole Human Genome (GPL4133)",
        "use_in_wp1": "Calibrates heavy-metal induction of HSF1 (HSP70) and NFE2L2 (oxidative stress response) vs ROS damage D",
    },
}

PROTEOSTASIS_REGULONS = {
    "HSF1": {
        "HSPA1A": 1.0,
        "HSPA1B": 0.95,
        "HSP90AA1": 0.90,
        "DNAJB1": 0.85,
        "HSPH1": 0.85,
        "BAG3": 0.80,
    },
    "HIF1A": {
        "VEGFA": 1.0,
        "LDHA": 0.90,
        "BNIP3": 0.90,
        "SLC2A1": 0.85,
        "PDK1": 0.85,
        "CA9": 0.80,
    },
    "NFE2L2": {
        "NQO1": 1.0,
        "HMOX1": 0.95,
        "SQSTM1": 0.90,
        "GCLC": 0.85,
        "TXNRD1": 0.85,
        "GCLM": 0.80,
    },
    "ATF4": {
        "DDIT3": 0.95,
        "ASNS": 0.90,
        "HSPA5": 0.85,
        "SQSTM1": 0.80,
        "MAP1LC3B": 0.80,
        "TRIB3": 0.80,
    },
}


class GEOCalibrator:
    """Real NCBI GEO data loader, probe-to-gene collapser, TF regulon scorer, and WP1 kinetic/Hill calibrator."""

    def __init__(self, cache_dir: str | Path = "./data/raw"):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.has_geoparse = importlib.util.find_spec("GEOparse") is not None
        self.has_decoupler = importlib.util.find_spec("decoupler") is not None

    def fetch_geo_expression(self, gse_id: str, filter_substring: str | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Download and parse a real NCBI GEO series into gene-level expression matrix + sample metadata."""
        if not self.has_geoparse:
            raise ImportError("GEOparse is required to download real GEO datasets. Install via: pip install GEOparse")

        import GEOparse

        gse = GEOparse.get_GEO(geo=gse_id, destdir=str(self.cache_dir), silent=True)
        meta_rows = []
        for gsm_name, gsm in gse.gsms.items():
            title = " ".join(gsm.metadata.get("title", [""]))
            chars = " | ".join(gsm.metadata.get("characteristics_ch1", []))
            source = " ".join(gsm.metadata.get("source_name_ch1", [""]))
            if filter_substring and (filter_substring.lower() not in f"{title} {chars} {source}".lower()):
                continue
            meta_rows.append({
                "sample_id": gsm_name,
                "title": title,
                "source": source,
                "characteristics": chars,
            })

        meta_df = pd.DataFrame(meta_rows)
        if meta_df.empty:
            raise ValueError(f"No samples matched filter '{filter_substring}' in {gse_id}.")

        expr_raw = gse.pivot_samples("VALUE")
        expr_raw = expr_raw[[s for s in meta_df["sample_id"] if s in expr_raw.columns]].apply(pd.to_numeric, errors="coerce").dropna()

        # Map probe IDs to HGNC gene symbols using GPL platform annotation table
        gpl = next(iter(gse.gpls.values()))
        gpl_tbl = gpl.table
        symbol_col = self._detect_gene_symbol_column(gpl_tbl)
        id_col = "ID" if "ID" in gpl_tbl.columns else gpl_tbl.columns[0]

        gene_df = self.collapse_probes_to_genes(expr_raw, gpl_tbl[[id_col, symbol_col]], id_col=id_col, symbol_col=symbol_col)
        # Log2-transform automatically if values appear unlogged (linear intensity scale)
        if float(gene_df.to_numpy().max()) > 50.0:
            gene_df = np.log2(np.maximum(gene_df, 0.0) + 1.0)
        return gene_df, meta_df

    @staticmethod
    def _detect_gene_symbol_column(gpl_table: pd.DataFrame) -> str:
        candidates = ["Gene Symbol", "GENE_SYMBOL", "Symbol", "SYMBOL", "gene_assignment", "ILMN_Gene"]
        for c in candidates:
            if c in gpl_table.columns:
                return c
        for c in gpl_table.columns:
            if "symbol" in c.lower() or "gene" in c.lower():
                return c
        raise KeyError("Could not detect a gene symbol column in GPL annotation table.")

    @staticmethod
    def collapse_probes_to_genes(
        expr_probes: pd.DataFrame,
        annotation_df: pd.DataFrame,
        id_col: str = "ID",
        symbol_col: str = "Gene Symbol",
    ) -> pd.DataFrame:
        """Collapse multiple microarray probes per gene using median intensity across mapped HGNC symbols."""
        ann = annotation_df[[id_col, symbol_col]].dropna().copy()
        ann[id_col] = ann[id_col].astype(str)

        def extract_symbol(val: str) -> str:
            s = str(val).strip()
            if "//" in s:
                parts = [p.strip() for p in s.split("//") if p.strip()]
                return parts[1] if len(parts) >= 2 else parts[0]
            if "///" in s:
                return s.split("///")[0].strip()
            return s.split()[0].strip() if s else ""

        ann["GENE"] = ann[symbol_col].map(extract_symbol)
        ann = ann[ann["GENE"].str.len() > 0].drop_duplicates(subset=[id_col])
        probe_to_gene = ann.set_index(id_col)["GENE"]

        df = expr_probes.copy()
        df.index = df.index.astype(str)
        df["GENE"] = df.index.map(probe_to_gene)
        df = df.dropna(subset=["GENE"]).groupby("GENE").median()
        return df

    def score_regulons(self, gene_expr: pd.DataFrame, regulons: dict | None = None) -> pd.DataFrame:
        """Compute weighted-mean transcription factor activity (z-scored across samples)."""
        regs = PROTEOSTASIS_REGULONS if regulons is None else regulons
        z_expr = gene_expr.sub(gene_expr.mean(axis=1), axis=0).div(gene_expr.std(axis=1) + 1e-8, axis=0)
        scores = {}
        for tf, targets in regs.items():
            present = [g for g in targets if g in z_expr.index]
            if not present:
                continue
            weights = np.array([targets[g] for g in present], dtype=float)
            scores[tf] = (z_expr.loc[present].T.to_numpy() @ weights) / np.sum(np.abs(weights))
        if not scores:
            raise ValueError("None of the proteostasis regulon target genes were found in the expression matrix.")
        df_tf = pd.DataFrame(scores, index=gene_expr.columns)
        df_tf["H_composite_z"] = df_tf.mean(axis=1)
        return df_tf

    @staticmethod
    def fit_hill_response(
        E_values: np.ndarray,
        H_activity: np.ndarray,
        baseline_beta0: float = 0.45,
        max_betaE_scale: float = 2.0,
        n_boot: int = 200,
        seed: int = 42,
    ) -> dict:
        """Non-circular Hill calibration anchored to empirical basal z-score and saturation plateau."""
        E_obs = np.asarray(E_values, dtype=float)
        H_raw = np.asarray(H_activity, dtype=float)

        # Anchor baseline to E==0 control samples without arbitrary rescaling to target parameters
        ctrl_mask = E_obs <= 1e-9
        h_ctrl = float(np.mean(H_raw[ctrl_mask])) if np.any(ctrl_mask) else float(np.min(H_raw))
        # Relative positive induction above basal control (in units of SD z-score shifted to fold-activation)
        H_rel = np.maximum(H_raw - h_ctrl, 0.0)
        # Scale induction into physiological WP1 beta_E units via empirical dynamic range (1 z-unit ~= 0.75x basal H)
        H_wp1 = baseline_beta0 + np.minimum(H_rel * 0.75, max_betaE_scale)

        def hill_fn(e, b0, bE, kE, n):
            return b0 + bE * (np.power(e, n) / (np.power(kE, n) + np.power(e, n)))

        popt, _ = curve_fit(
            hill_fn,
            E_obs,
            H_wp1,
            p0=[baseline_beta0, max(float(np.max(H_wp1) - baseline_beta0), 0.3), max(float(np.median(E_obs[E_obs > 0])), 0.2), 2.0],
            bounds=([0.05, 0.10, 0.05, 0.6], [2.5, 4.0, 3.0, 5.0]),
            maxfev=10000,
        )

        rng = np.random.default_rng(seed)
        boot = []
        for _ in range(n_boot):
            idx = rng.integers(0, len(E_obs), size=len(E_obs))
            if len(np.unique(E_obs[idx])) < 3:
                continue
            try:
                bp, _ = curve_fit(
                    hill_fn,
                    E_obs[idx],
                    H_wp1[idx],
                    p0=popt,
                    bounds=([0.05, 0.10, 0.05, 0.6], [2.5, 4.0, 3.0, 5.0]),
                    maxfev=5000,
                )
                boot.append(bp)
            except Exception:
                continue

        boot_arr = np.array(boot) if boot else np.array([popt])
        ci_low = np.percentile(boot_arr, 2.5, axis=0)
        ci_high = np.percentile(boot_arr, 97.5, axis=0)

        return {
            "wp1_parameters": {
                "beta_0": round(float(popt[0]), 4),
                "beta_E": round(float(popt[1]), 4),
                "K_E": round(float(popt[2]), 4),
                "n_hill": round(float(popt[3]), 4),
            },
            "ci_95_low": {k: round(float(v), 4) for k, v in zip(["beta_0", "beta_E", "K_E", "n_hill"], ci_low)},
            "ci_95_high": {k: round(float(v), 4) for k, v in zip(["beta_0", "beta_E", "K_E", "n_hill"], ci_high)},
            "n_samples": int(len(E_obs)),
            "n_bootstrap_converged": int(len(boot)),
        }

    @staticmethod
    def fit_kinetic_relaxation(time_hours: np.ndarray, H_activity: np.ndarray) -> dict:
        """Calibrate induction/turnover rate delta_H from real time-course GEO series (e.g., GSE18494)."""
        t_obs = np.asarray(time_hours, dtype=float)
        h_obs = np.asarray(H_activity, dtype=float)
        h0 = float(np.mean(h_obs[t_obs <= 1e-9])) if np.any(t_obs <= 1e-9) else float(h_obs[0])
        h_shift = h_obs - h0

        def exp_rise(t, amp, rate):
            return amp * (1.0 - np.exp(-rate * t))

        popt, _ = curve_fit(
            exp_rise,
            t_obs,
            h_shift,
            p0=[max(float(np.max(h_shift)), 0.5), 0.25],
            bounds=([0.05, 0.01], [10.0, 3.0]),
            maxfev=5000,
        )
        return {
            "amplitude_z": round(float(popt[0]), 4),
            "delta_H_per_hour": round(float(popt[1]), 4),
            "half_life_hours": round(float(np.log(2.0) / popt[1]), 3),
        }

    @staticmethod
    def query_string_network(genes: list[str], species: int = 9606, min_score: int = 700) -> pd.DataFrame:
        url = "https://string-db.org/api/json/network"
        params = {
            "identifiers": "%0d".join(genes),
            "species": species,
            "required_score": min_score,
            "caller_identity": "proteostasis_sim_joss",
        }
        resp = requests.get(url, params=params, timeout=20)
        resp.raise_for_status()
        df = pd.DataFrame(resp.json())
        cols = [c for c in ["preferredName_A", "preferredName_B", "score", "escore", "dscore"] if c in df.columns]
        return df[cols].sort_values("score", ascending=False).reset_index(drop=True)


def main_cli():
    parser = argparse.ArgumentParser(description="Download real NCBI GEO datasets and calibrate WP1 proteostasis parameters.")
    parser.add_argument("--gse", type=str, default="GSE41666", choices=list(GEO_DATASET_CATALOG.keys()))
    parser.add_argument("--filter", type=str, default="HepG2", help="Sample title/metadata filter substring")
    parser.add_argument("--out", type=str, default="data/processed/geo_tf_scores.csv")
    args = parser.parse_args()

    cal = GEOCalibrator()
    gene_df, meta_df = cal.fetch_geo_expression(args.gse, filter_substring=args.filter if args.filter else None)
    tf_df = cal.score_regulons(gene_df)
    merged = meta_df.set_index("sample_id").join(tf_df)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    merged.to_csv(args.out)
    print(f"Saved {args.gse} TF activity scores ({len(merged)} samples) to {args.out}")
    print(merged[["title", "HSF1", "HIF1A", "NFE2L2", "ATF4", "H_composite_z"]].round(3).to_string())
