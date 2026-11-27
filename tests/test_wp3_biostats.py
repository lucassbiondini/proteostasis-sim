import shutil
import subprocess
import unittest
from pathlib import Path
import numpy as np
from proteostasis_sim.platform import InSilicoMultiscalePlatform
from proteostasis_sim.wp3_biostats import CompetingRisksBiostats


class TestWP3BiostatsAndCrossVal(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.platform = InSilicoMultiscalePlatform(
            tissue_params={"grid_size": 18, "dt_abm": 0.35},
            seed=42,
        )
        cls.platform.build_surrogate_and_biostats(
            E_grid=np.array([0.0, 0.3, 0.7, 1.3, 1.6, 1.9]),
            n_reps_wp2=2,
            n_cells_wp1=60,
        )

    def test_aalen_johansen_probability_bounds(self):
        df = self.platform.wp3.simulate_virtual_cohort(np.array([0.0, 1.55]), n_per_arm=600, mode="mascaramento", seed=1)
        sub = df[df["E"] == 1.55]
        c1 = CompetingRisksBiostats.aalen_johansen_cif(sub["time"].to_numpy(), sub["event"].to_numpy(), 1, 100.0)
        c2 = CompetingRisksBiostats.aalen_johansen_cif(sub["time"].to_numpy(), sub["event"].to_numpy(), 2, 100.0)
        self.assertGreaterEqual(c1, 0.0)
        self.assertGreaterEqual(c2, 0.0)
        self.assertLessEqual(c1 + c2, 1.0 + 1e-9)

    def test_masking_vs_hormesis_causal_separation(self):
        df_audit = self.platform.wp3.run_causal_audit(np.array([0.0, 0.5, 1.6]), n_per_arm=800, t_ref=100.0)
        row_horm = df_audit[df_audit["E"] == 0.5].iloc[0]
        self.assertLess(row_horm["HR_CSH_Horm"], 0.15)
        self.assertGreater(row_horm["HR_CSH_Mask"], 0.85)
        self.assertGreater(row_horm["Phi_Antag_Real"], 0.85)

    def test_fine_gray_masking_divergence_under_competing_toxicity(self):
        bundle = self.platform.export_r_validation_bundle(out_dir="data/processed", n_per_arm=800)
        m = bundle["metrics"]
        # Under pure masking (no biological effect on T_gen), Cox CSH ~ 1.0 while Fine-Gray sHR < 1.0
        self.assertGreater(m["hr_csh_gen"], 0.80)
        self.assertLess(m["hr_csh_gen"], 1.25)
        self.assertLess(m["shr_fg_gen"], m["hr_csh_gen"])

    def test_optional_r_cmprsk_cross_validation(self):
        bundle = self.platform.export_r_validation_bundle(out_dir="data/processed", n_per_arm=800)
        rscript = shutil.which("Rscript")
        if rscript is None:
            self.skipTest("Rscript not available in current environment; run scripts/validate_python_vs_r.R locally.")
        check_pkgs = subprocess.run(
            [rscript, "-e", "quit(status = ifelse(all(c('survival','cmprsk','jsonlite') %in% rownames(installed.packages())), 0, 1))"],
            capture_output=True,
        )
        if check_pkgs.returncode != 0:
            self.skipTest("R packages survival/cmprsk/jsonlite not installed in active R environment.")
        res = subprocess.run(
            [rscript, "scripts/validate_python_vs_r.R", bundle["csv"], bundle["json"], "0.01"],
            capture_output=True,
            text=True,
        )
        self.assertEqual(res.returncode, 0, msg=f"R cross-validation failed:\n{res.stdout}\n{res.stderr}")
