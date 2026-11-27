import unittest
import numpy as np
from proteostasis_sim.wp1_molecular import MolecularProteostasisModel


class TestWP1Molecular(unittest.TestCase):
    def setUp(self):
        self.model = MolecularProteostasisModel(seed=42)

    def test_state_positivity_and_bounds(self):
        sol, tg, ta = self.model.simulate_ode(0.6)
        self.assertTrue(np.all(sol.y[:4] >= -1e-6))
        self.assertTrue(np.all((sol.y[4] >= 0.0) & (sol.y[4] <= 1.01)))
        self.assertTrue(np.isinf(tg) and np.isinf(ta))

    def test_bistability_hysteresis_exists(self):
        E_grid, A_low, A_high = self.model.analyze_hysteresis(np.linspace(0.0, 1.7, 25))
        bistable = np.abs(A_high - A_low) > 3.0
        self.assertTrue(np.any(bistable), "Model must exhibit saddle-node hysteresis.")
        self.assertLess(E_grid[bistable][0], 0.35)
        self.assertGreater(E_grid[bistable][-1], 1.20)

    def test_pulsatile_hormesis_outperforms_chronic_at_high_E(self):
        _, tg_c, ta_c = self.model.simulate_ode(1.60, regime="cronico")
        _, tg_p, ta_p = self.model.simulate_ode(1.60, regime="pulsatil")
        t_obs_c = min(tg_c, ta_c, 150.0)
        t_obs_p = min(tg_p, ta_p, 150.0)
        self.assertGreater(t_obs_p, t_obs_c)

    def test_genotype_stratification_ordering(self):
        p_mz = self.model.get_params("PiMZ_Moderado")
        p_zz = self.model.get_params("PiZZ_Padrao")
        p_low = self.model.get_params("PiZZ_PRS_BaixoH")
        t_mz = min(*self.model.simulate_ode(0.0, p=p_mz)[1:], 150.0)
        t_zz = min(*self.model.simulate_ode(0.0, p=p_zz)[1:], 150.0)
        t_low = min(*self.model.simulate_ode(0.0, p=p_low)[1:], 150.0)
        self.assertGreater(t_mz, t_zz)
        self.assertGreater(t_zz, t_low)

    def test_combination_screening_runs_and_bliss_valid(self):
        df_screen, df_syn = self.model.run_combination_screening(np.linspace(0.0, 2.0, 11))
        self.assertEqual(len(df_screen), len(self.model.genotypes) * len(self.model.therapies))
        self.assertTrue(np.all(np.isfinite(df_syn["Excesso_Bliss"])))
