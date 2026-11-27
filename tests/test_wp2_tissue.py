import unittest
import numpy as np
from proteostasis_sim.wp1_molecular import MolecularProteostasisModel
from proteostasis_sim.wp2_tissue import TissueABMPDEModel


class TestWP2Tissue(unittest.TestCase):
    def setUp(self):
        self.wp1 = MolecularProteostasisModel(seed=42)
        self.wp2 = TissueABMPDEModel(self.wp1, tissue_params={"grid_size": 20, "dt_abm": 0.3})

    def test_three_tissue_regimes(self):
        out_ctrl = self.wp2.simulate_tissue(0.0, seed=10)
        out_horm = self.wp2.simulate_tissue(0.6, seed=10)
        out_tox = self.wp2.simulate_tissue(1.7, seed=10)
        self.assertEqual(out_ctrl["cause"], "genetica")
        self.assertEqual(out_horm["cause"], "censurado")
        self.assertEqual(out_tox["cause"], "ambiental")
        self.assertLess(out_tox["t_fail"], out_ctrl["t_fail"])

    def test_bottom_up_surrogate_hormetic_peak(self):
        E_grid = np.array([0.0, 0.5, 1.7])
        df_sur = self.wp2.build_bottom_up_surrogate(E_grid, n_reps_wp2=2, n_cells_wp1=50)
        s0 = float(df_sur.loc[df_sur["E"] == 0.0, "scale_gen"].iloc[0])
        s_mid = float(df_sur.loc[df_sur["E"] == 0.5, "scale_gen"].iloc[0])
        s_hi = float(df_sur.loc[df_sur["E"] == 1.7, "scale_amb"].iloc[0])
        self.assertGreater(s_mid, s0)
        self.assertLess(s_hi, 100.0)
