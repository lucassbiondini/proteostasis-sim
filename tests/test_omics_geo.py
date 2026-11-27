import tempfile
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from proteostasis_sim.omics_geo import GEOCalibrator, GEO_DATASET_CATALOG
from proteostasis_sim.wp1_molecular import MolecularProteostasisModel


class TestOmicsGEO(unittest.TestCase):
    def test_catalog_contains_verified_accessions(self):
        for acc in ["GSE66076", "GSE41666", "GSE18494", "GSE31286"]:
            self.assertIn(acc, GEO_DATASET_CATALOG)

    def test_probe_collapse_and_regulon_scoring(self):
        expr_probes = pd.DataFrame(
            {
                "GSM1": [6.0, 6.2, 5.5, 5.0, 5.1, 5.0],
                "GSM2": [6.1, 6.3, 5.4, 5.2, 5.0, 5.1],
                "GSM3": [8.5, 8.7, 7.9, 7.6, 7.5, 7.8],
                "GSM4": [9.0, 9.2, 8.4, 8.1, 8.0, 8.3],
            },
            index=["P1", "P2", "P3", "P4", "P5", "P6"],
        )
        ann = pd.DataFrame({
            "ID": ["P1", "P2", "P3", "P4", "P5", "P6"],
            "Gene Symbol": ["NM_005345 // HSPA1A // heat shock", "HSPA1A", "HSP90AA1", "VEGFA", "NQO1", "DDIT3"],
        })
        cal = GEOCalibrator()
        gene_df = cal.collapse_probes_to_genes(expr_probes, ann)
        self.assertIn("HSPA1A", gene_df.index)
        self.assertAlmostEqual(float(gene_df.loc["HSPA1A", "GSM1"]), 6.1, places=5)

        tf_df = cal.score_regulons(gene_df)
        self.assertGreater(float(tf_df.loc["GSM4", "HSF1"]), float(tf_df.loc["GSM1", "HSF1"]))

    def test_non_circular_hill_calibration_and_wp1_injection(self):
        cal = GEOCalibrator()
        E_vals = np.array([0.0, 0.0, 0.2, 0.2, 0.5, 0.5, 1.0, 1.0, 1.6, 1.6])
        H_z = np.array([-1.2, -1.1, -0.7, -0.6, 0.1, 0.2, 0.8, 0.9, 1.1, 1.2])
        res = cal.fit_hill_response(E_vals, H_z, n_boot=40)
        p = res["wp1_parameters"]
        self.assertGreater(p["beta_E"], 0.5)
        self.assertGreater(p["K_E"], 0.1)

        with tempfile.TemporaryDirectory() as tmp:
            jpath = Path(tmp) / "cal.json"
            jpath.write_text(__import__("json").dumps(res), encoding="utf-8")
            wp1 = MolecularProteostasisModel()
            updated = wp1.load_omics_calibration(jpath)
            self.assertEqual(wp1.params["beta_E"], updated["beta_E"])
