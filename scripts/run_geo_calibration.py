import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
from proteostasis_sim.omics_geo import GEOCalibrator, GEO_DATASET_CATALOG


def run_gse41666_hepg2_audit(cal: GEOCalibrator, out_dir: Path):
    gene_df, meta_df = cal.fetch_geo_expression("GSE41666", filter_substring="HepG2")
    tf_df = cal.score_regulons(gene_df)
    res = meta_df.set_index("sample_id").join(tf_df)
    res.to_csv(out_dir / "GSE41666_HepG2_regulons.csv")
    print("\n[GSE41666] HepG2 Normoxia vs Anoxia vs Tunicamycin vs Thapsigargin:")
    print(res[["title", "HSF1", "HIF1A", "NFE2L2", "ATF4", "H_composite_z"]].round(3).to_string())
    return res


def run_gse66076_serpina1_audit(cal: GEOCalibrator, out_dir: Path):
    gene_df, meta_df = cal.fetch_geo_expression("GSE66076")
    tf_df = cal.score_regulons(gene_df)
    res = meta_df.set_index("sample_id").join(tf_df)
    res.to_csv(out_dir / "GSE66076_SERPINA1_regulons.csv")
    print("\n[GSE66076] SERPINA1 PiZZ iPSC-Hepatocytes Proteostasis Signature:")
    print(res[["title", "HSF1", "ATF4", "NFE2L2", "H_composite_z"]].head(12).round(3).to_string())
    return res


def run_gse18494_kinetics(cal: GEOCalibrator, out_dir: Path):
    gene_df, meta_df = cal.fetch_geo_expression("GSE18494", filter_substring="HepG2")
    tf_df = cal.score_regulons(gene_df)
    res = meta_df.set_index("sample_id").join(tf_df)

    def parse_hours(txt: str) -> float:
        m = re.search(r"(\d+(?:\.\d+)?)\s*h", txt.lower())
        if m:
            return float(m.group(1))
        return 0.0 if ("normox" in txt.lower() or "control" in txt.lower() or "0h" in txt.lower()) else np.nan

    res["hours"] = [parse_hours(f"{t} {c}") for t, c in zip(res["title"], res["characteristics"])]
    valid = res.dropna(subset=["hours"])
    kin = cal.fit_kinetic_relaxation(valid["hours"].to_numpy(), valid["HIF1A"].to_numpy()) if len(valid) >= 4 else {}
    res.to_csv(out_dir / "GSE18494_HepG2_kinetics.csv")
    (out_dir / "GSE18494_kinetic_calibration.json").write_text(json.dumps(kin, indent=2), encoding="utf-8")
    print("\n[GSE18494] Calibracao Cinetica Temporal (HIF1A/Proteostase em HepG2):")
    print(json.dumps(kin, indent=2))
    return res, kin


if __name__ == "__main__":
    out_dir = Path("data/processed")
    out_dir.mkdir(parents=True, exist_ok=True)
    cal = GEOCalibrator(cache_dir="data/raw")
    print("Catalogo GEO verificado:", list(GEO_DATASET_CATALOG.keys()))
    run_gse41666_hepg2_audit(cal, out_dir)
    run_gse66076_serpina1_audit(cal, out_dir)
    run_gse18494_kinetics(cal, out_dir)
