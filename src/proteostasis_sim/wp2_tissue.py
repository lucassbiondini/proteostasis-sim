from __future__ import annotations
import numpy as np
import pandas as pd
from scipy.ndimage import laplace, binary_dilation
from scipy.interpolate import PchipInterpolator
from .wp1_molecular import MolecularProteostasisModel


DEFAULT_WP2_PARAMS = {
    "grid_size": 30,
    "dt_abm": 0.25,
    "D_C": 0.16,
    "s_A": 0.045,
    "s_D": 0.085,
    "lambda_C": 0.24,
    "bystander_D": 0.11,
    "C_fib_thresh": 0.55,
    "k_mmp": 0.003,
    "p_prolif": 0.065,
    "max_div": 5,
    "fail_thresh": 0.35,
    "chronic_fib_thresh": 0.18,
}


class TissueABMPDEModel:
    """WP2: 2D Hybrid Agent-Based + Reaction-Diffusion PDE Tissue Model."""

    def __init__(self, wp1_model: MolecularProteostasisModel, tissue_params: dict | None = None):
        self.wp1 = wp1_model
        self.tissue_params = DEFAULT_WP2_PARAMS.copy()
        if tissue_params:
            self.tissue_params.update(tissue_params)
        self.interp_gen = None
        self.interp_amb = None
        self.df_surrogate = None

    def simulate_tissue(
        self,
        E_amp: float,
        wp1_params: dict | None = None,
        regime: str = "cronico",
        seed: int = 42,
        save_snapshot: bool = False,
    ) -> dict:
        p1 = self.wp1.params if wp1_params is None else wp1_params
        p2 = self.tissue_params
        rng = np.random.default_rng(seed)
        N = int(p2["grid_size"])
        dt = float(p2["dt_abm"])
        n_steps = int(p1["t_max"] / dt)

        alive = np.ones((N, N), dtype=bool)
        fibrosis = np.zeros((N, N), dtype=bool)
        div_count = np.zeros((N, N), dtype=np.int16)

        M = np.full((N, N), self.wp1.y0[0], dtype=float)
        A = np.full((N, N), self.wp1.y0[1], dtype=float)
        H = np.full((N, N), self.wp1.y0[2], dtype=float)
        D = np.zeros((N, N), dtype=float)
        B = np.ones((N, N), dtype=float)
        C = np.zeros((N, N), dtype=float)

        def sample_lognorm(base, cv, shape):
            return rng.lognormal(mean=np.log(base) - 0.5 * cv**2, sigma=cv, size=shape)

        k_syn_g = sample_lognorm(p1["k_syn"], 0.22, (N, N))
        gamma_E_g = sample_lognorm(p1["gamma_E"], 0.28, (N, N))
        beta_E_g = sample_lognorm(p1["beta_E"], 0.22, (N, N))

        gen_deaths = 0
        amb_deaths = 0
        t_fail = np.inf
        cause_fail = "censurado"
        t_chronic_amb = np.inf
        snapshot = None
        kernel_4 = np.array([[0, 1, 0], [1, 0, 1], [0, 1, 0]], dtype=bool)

        for step in range(n_steps + 1):
            t = step * dt
            func_mask = alive & (A < 0.5 * p1["A_crit"]) & (B > 0.35)
            frac_func = float(func_mask.mean())
            frac_fib = float(fibrosis.mean())

            if np.isinf(t_chronic_amb) and (frac_fib >= p2["chronic_fib_thresh"]) and (amb_deaths > 0 or E_amp >= 0.4):
                t_chronic_amb = t

            if save_snapshot and step == int(60.0 / dt):
                st = np.zeros((N, N), dtype=np.int8)
                st[func_mask & (H < 2.0)] = 1
                st[func_mask & (H >= 2.0)] = 2
                st[alive & (~func_mask)] = 3
                st[fibrosis] = 4
                snapshot = (st.copy(), C.copy())

            if np.isinf(t_fail) and frac_func < p2["fail_thresh"]:
                t_fail = t
                cause_fail = "genetica" if gen_deaths >= amb_deaths else "ambiental"
                if not save_snapshot:
                    break

            if step == n_steps:
                break

            E_t = self.wp1.exposure_signal(t, E_amp, regime=regime)
            hill_E = (E_t**p1["n_hill"]) / (p1["K_E"]**p1["n_hill"] + E_t**p1["n_hill"]) if E_t > 0 else 0.0
            age_factor = np.exp(-p1["alpha_age"] * t)

            S_prod = np.where(
                alive,
                p2["s_A"] * np.maximum(0.0, A - 3.0) + p2["s_D"] * np.maximum(0.0, D - 2.0),
                0.0,
            )
            C = np.maximum(0.0, C + dt * (p2["D_C"] * laplace(C, mode="nearest") + S_prod - p2["lambda_C"] * C))

            if np.any(alive):
                M_a, A_a, H_a, D_a, B_a, C_a = M[alive], A[alive], H[alive], D[alive], B[alive], C[alive]
                ks_a, ge_a, be_a = k_syn_g[alive], gamma_E_g[alive], beta_E_g[alive]

                atp_eff = (1.25 * B_a) / (0.25 + B_a)
                ros_inhib = 1.0 / (1.0 + 0.07 * D_a)
                hill_M = (M_a**p1["m_hill"]) / (p1["K_PM"]**p1["m_hill"] + M_a**p1["m_hill"])

                dM = ks_a - p1["k_deg"] * ros_inhib * M_a - (p1["k_clr"] * atp_eff * H_a * M_a) / (p1["K_M"] + M_a) - 2.0 * p1["k_nuc"] * (1.0 + p1["k_ox"] * D_a) * (M_a**2) - p1["k_elong"] * M_a * A_a
                dA = p1["k_nuc"] * (1.0 + p1["k_ox"] * D_a) * (M_a**2) + p1["k_elong"] * M_a * A_a - (p1["k_auto"] * atp_eff * ros_inhib * H_a * A_a) / (p1["K_A"] + A_a) - p1["delta_A"] * A_a
                dH = (p1["beta_0"] * age_factor + be_a * hill_E + p1["beta_M"] * hill_M) * atp_eff - p1["delta_H"] * H_a - p1["k_seq"] * H_a * A_a
                dD = ge_a * (E_t**p1["p_tox"]) + p2["bystander_D"] * C_a - (p1["gamma_H"] * atp_eff * H_a * D_a) / (p1["K_D"] + D_a) - p1["delta_D"] * D_a
                dB = p1["r_B"] * (1.0 - B_a) - p1["c_H"] * H_a * B_a - p1["c_D"] * D_a * B_a - p1["c_A"] * A_a * B_a

                M[alive] = np.maximum(0.0, M_a + dt * dM)
                A[alive] = np.maximum(0.0, A_a + dt * dA)
                H[alive] = np.maximum(0.0, H_a + dt * dH)
                D[alive] = np.maximum(0.0, D_a + dt * dD)
                B[alive] = np.clip(B_a + dt * dB, 0.01, 1.0)

                morte_gen = alive & (A >= p1["A_crit"])
                morte_amb = alive & (~morte_gen) & ((D >= p1["D_crit"]) | (B <= p1["B_crit"]))
                morte_tot = morte_gen | morte_amb

                if np.any(morte_tot):
                    gen_deaths += int(morte_gen.sum())
                    amb_deaths += int(morte_amb.sum())
                    alive[morte_tot] = False
                    fibrosis[morte_tot & (C >= p2["C_fib_thresh"])] = True
                    M[morte_tot] = A[morte_tot] = H[morte_tot] = D[morte_tot] = B[morte_tot] = 0.0

            if np.any(fibrosis):
                resolve = fibrosis & (C < 0.35 * p2["C_fib_thresh"]) & (rng.random((N, N)) < p2["k_mmp"] * dt)
                fibrosis[resolve] = False

            free_space = (~alive) & (~fibrosis)
            if np.any(free_space):
                donors = alive & (A < 0.3 * p1["A_crit"]) & (D < 0.5 * p1["D_crit"]) & (B > 0.55) & (div_count < p2["max_div"])
                cand = free_space & binary_dilation(donors, structure=kernel_4)
                if np.any(cand):
                    born = cand & (rng.random((N, N)) < p2["p_prolif"] * dt)
                    n_born = int(born.sum())
                    if n_born > 0:
                        alive[born] = True
                        M[born] = 1.0
                        A[born] = 0.1
                        H[born] = 1.0 + 1.2 * hill_E
                        D[born] = 0.0
                        B[born] = 0.95
                        div_count[born] = np.minimum(div_count[born] + 1, p2["max_div"])
                        k_syn_g[born] = sample_lognorm(p1["k_syn"], 0.22, n_born)
                        gamma_E_g[born] = sample_lognorm(p1["gamma_E"], 0.28, n_born)
                        beta_E_g[born] = sample_lognorm(p1["beta_E"], 0.22, n_born)

        return {
            "E": float(E_amp),
            "t_fail": float(t_fail),
            "cause": cause_fail,
            "t_chronic_amb": float(t_chronic_amb),
            "gen_deaths": int(gen_deaths),
            "amb_deaths": int(amb_deaths),
            "frac_fib_final": float(fibrosis.mean()),
            "snapshot": snapshot,
        }

    @staticmethod
    def fit_censored_weibull_mle(times: np.ndarray, events: np.ndarray, fixed_shape: float = 3.0) -> tuple[float, float]:
        times = np.maximum(np.asarray(times, dtype=float), 1e-3)
        events = np.asarray(events, dtype=float)
        k = float(fixed_shape)
        denom = max(float(events.sum()), 0.25)
        lam = (np.sum(times**k) / denom) ** (1.0 / k)
        return float(np.clip(lam, 5.0, 2500.0)), k

    def build_bottom_up_surrogate(self, E_grid: np.ndarray, n_reps_wp2: int = 4, n_cells_wp1: int = 120) -> pd.DataFrame:
        rows = []
        t_max = float(self.wp1.params["t_max"])
        for i, E_val in enumerate(E_grid):
            tg_sde, ta_sde = self.wp1.simulate_sde_ensemble(E_val, n_cells=n_cells_wp1, seed=200 + i)
            ev_g_sde, ev_a_sde = np.isfinite(tg_sde), np.isfinite(ta_sde)
            t_g_sde = np.where(ev_g_sde, tg_sde, t_max)
            t_a_sde = np.where(ev_a_sde, ta_sde, t_max)
            lam_g_wp1, _ = self.fit_censored_weibull_mle(t_g_sde, ev_g_sde, fixed_shape=3.2)
            lam_a_wp1, _ = self.fit_censored_weibull_mle(t_a_sde, ev_a_sde, fixed_shape=2.5)

            t_wp2, c_wp2, fib_wp2 = [], [], []
            for r in range(n_reps_wp2):
                out = self.simulate_tissue(E_val, seed=500 + i * 17 + r)
                t_wp2.append(min(out["t_fail"], t_max))
                c_wp2.append(out["cause"])
                fib_wp2.append(out["frac_fib_final"])

            t_wp2 = np.array(t_wp2)
            ev_g_wp2 = np.array([c == "genetica" for c in c_wp2], dtype=float)
            ev_a_wp2 = np.array([c == "ambiental" for c in c_wp2], dtype=float)
            lam_g_wp2, _ = self.fit_censored_weibull_mle(t_wp2, ev_g_wp2, fixed_shape=3.2)
            lam_a_wp2, _ = self.fit_censored_weibull_mle(t_wp2, ev_a_wp2, fixed_shape=2.5)

            w_g = ev_g_wp2.sum() / (ev_g_wp2.sum() + 1.0)
            w_a = ev_a_wp2.sum() / (ev_a_wp2.sum() + 1.0)
            lam_g_coupled = w_g * lam_g_wp2 + (1.0 - w_g) * (1.35 * lam_g_wp1)
            lam_a_coupled = w_a * lam_a_wp2 + (1.0 - w_a) * (1.15 * lam_a_wp1)

            rows.append({
                "E": float(E_val),
                "wp2_t_mean": float(t_wp2.mean()),
                "wp2_t_std": float(t_wp2.std()),
                "wp2_fib_mean": float(np.mean(fib_wp2)),
                "scale_gen": float(lam_g_coupled),
                "scale_amb": float(lam_a_coupled),
            })

        df_sur = pd.DataFrame(rows)
        self.interp_gen = PchipInterpolator(df_sur["E"], df_sur["scale_gen"])
        self.interp_amb = PchipInterpolator(df_sur["E"], df_sur["scale_amb"])
        self.df_surrogate = df_sur
        return df_sur
