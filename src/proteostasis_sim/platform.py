from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .wp1_molecular import MolecularProteostasisModel
from .wp2_tissue import TissueABMPDEModel
from .wp3_biostats import CompetingRisksBiostats


class InSilicoMultiscalePlatform:
    """Unified bottom-up facade integrating WP1 (Molecular), WP2 (Tissue ABM-PDE), and WP3 (Biostats)."""

    def __init__(self, wp1_params: dict | None = None, tissue_params: dict | None = None, seed: int = 42):
        self.seed = seed
        self.wp1 = MolecularProteostasisModel(params=wp1_params, seed=seed)
        self.wp2 = TissueABMPDEModel(self.wp1, tissue_params=tissue_params)
        self.wp3: CompetingRisksBiostats | None = None

    def calibrate_from_geo_json(self, json_path: str | Path) -> dict:
        return self.wp1.load_omics_calibration(json_path)

    def build_surrogate_and_biostats(
        self,
        E_grid: np.ndarray | None = None,
        n_reps_wp2: int = 4,
        n_cells_wp1: int = 120,
    ) -> pd.DataFrame:
        E_grid = np.array([0.0, 0.15, 0.30, 0.50, 0.75, 1.00, 1.25, 1.45, 1.60, 1.80, 2.00]) if E_grid is None else E_grid
        df_sur = self.wp2.build_bottom_up_surrogate(E_grid, n_reps_wp2=n_reps_wp2, n_cells_wp1=n_cells_wp1)
        self.wp3 = CompetingRisksBiostats(self.wp2.interp_gen, self.wp2.interp_amb)
        return df_sur

    def export_r_validation_bundle(self, out_dir: str | Path = "data/processed", n_per_arm: int = 1500) -> dict:
        """Generate cohort CSV + Python reference metrics JSON for independent R (survival/cmprsk) cross-validation."""
        if self.wp3 is None:
            self.build_surrogate_and_biostats(n_reps_wp2=2, n_cells_wp1=80)

        out_path = Path(out_dir)
        out_path.mkdir(parents=True, exist_ok=True)
        # Use moderate exposure arms where both genetic and competing environmental events occur with random censoring
        E_arms = np.array([0.0, 0.20, 1.55])
        df_coh = self.wp3.simulate_virtual_cohort(
            E_arms,
            n_per_arm=n_per_arm,
            mode="mascaramento",
            sigma_shared=0.0,
            sigma_ind=0.25,
            t_cens=150.0,
            random_censor_rate=0.004,
            seed=2026,
        )
        # Keep a binary comparison arm (E=0.0 vs E=1.55) where competing events are abundant and censoring is non-trivial
        df_bin = df_coh[df_coh["E"].isin([0.0, 1.55])].copy()
        df_bin["arm"] = (df_bin["E"] > 0.0).astype(int)
        csv_file = out_path / "crossval_cohort.csv"
        df_bin[["time", "event", "arm"]].to_csv(csv_file, index=False)

        df0 = df_bin[df_bin["arm"] == 0]
        df1 = df_bin[df_bin["arm"] == 1]

        metrics = {
            "t_ref": 100.0,
            "cif_gen_arm0": self.wp3.aalen_johansen_cif(df0["time"].to_numpy(), df0["event"].to_numpy(), 1, 100.0),
            "cif_amb_arm0": self.wp3.aalen_johansen_cif(df0["time"].to_numpy(), df0["event"].to_numpy(), 2, 100.0),
            "cif_gen_arm1": self.wp3.aalen_johansen_cif(df1["time"].to_numpy(), df1["event"].to_numpy(), 1, 100.0),
            "cif_amb_arm1": self.wp3.aalen_johansen_cif(df1["time"].to_numpy(), df1["event"].to_numpy(), 2, 100.0),
            "hr_csh_gen": self.wp3.fit_cox_csh(df0, df1, target_cause=1, l2_penalty=0.0),
            "shr_fg_gen": self.wp3.fit_fine_gray_shr(df0, df1, target_cause=1, l2_penalty=0.0),
        }
        json_file = out_path / "python_reference_metrics.json"
        json_file.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
        return {"csv": str(csv_file), "json": str(json_file), "metrics": metrics}
