from __future__ import annotations
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.stats import qmc, rankdata


DEFAULT_WP1_PARAMS = {
    "k_syn": 2.35,
    "k_deg": 0.20,
    "k_clr": 1.05,
    "K_M": 2.0,
    "k_ox": 0.12,
    "k_nuc": 0.0042,
    "k_elong": 0.042,
    "k_auto": 0.90,
    "K_A": 3.0,
    "delta_A": 0.045,
    "beta_0": 0.45,
    "beta_E": 1.85,
    "K_E": 0.45,
    "n_hill": 2.5,
    "beta_M": 0.15,
    "K_PM": 5.0,
    "m_hill": 2.0,
    "delta_H": 0.42,
    "k_seq": 0.09,
    "alpha_age": 0.0015,
    "gamma_E": 0.95,
    "p_tox": 2.6,
    "gamma_H": 0.28,
    "K_D": 2.5,
    "delta_D": 0.18,
    "r_B": 0.45,
    "c_H": 0.025,
    "c_D": 0.065,
    "c_A": 0.015,
    "A_crit": 15.0,
    "D_crit": 15.0,
    "B_crit": 0.22,
    "t_max": 150.0,
}

DEFAULT_GENOTYPES = {
    "PiMZ_Moderado": {"k_syn": 0.86, "k_nuc": 0.80, "beta_E": 1.00},
    "PiZZ_Padrao": {"k_syn": 1.00, "k_nuc": 1.00, "beta_E": 1.00},
    "PiZZ_PRS_AltoH": {"k_syn": 1.00, "k_nuc": 1.00, "beta_E": 1.25},
    "PiZZ_PRS_BaixoH": {"k_syn": 1.08, "k_nuc": 1.20, "beta_E": 0.65},
}

DEFAULT_THERAPIES = {
    "Controle": {},
    "Indutor_Autofagia": {"k_auto": 1.45},
    "Antioxidante_Mito": {"c_D": 0.55, "k_ox": 0.50, "gamma_H": 1.35},
    "Ativador_HSF1": {"beta_E": 1.35, "beta_0": 1.25},
    "Anti_Agregante": {"k_nuc": 0.65, "k_seq": 0.70},
    "Combo_Auto_Mito": {"k_auto": 1.45, "c_D": 0.55, "k_ox": 0.50, "gamma_H": 1.35},
    "Combo_HSF1_Mito": {"beta_E": 1.35, "beta_0": 1.25, "c_D": 0.55, "k_ox": 0.50, "gamma_H": 1.35},
    "Combo_Auto_AntiAgg": {"k_auto": 1.45, "k_nuc": 0.65, "k_seq": 0.70},
}


class MolecularProteostasisModel:
    """WP1: 5-ODE nonlinear proteostasis, bioenergetics, aging, and SDE model."""

    def __init__(self, params: dict | None = None, seed: int = 42):
        self.seed = seed
        self.params = DEFAULT_WP1_PARAMS.copy()
        if params:
            self.params.update(params)
        self.y0 = np.array([1.0, 0.05, 1.0, 0.0, 1.0], dtype=float)
        self.genotypes = {k: v.copy() for k, v in DEFAULT_GENOTYPES.items()}
        self.therapies = {k: v.copy() for k, v in DEFAULT_THERAPIES.items()}

    def load_omics_calibration(self, json_path: str | Path) -> dict:
        """Inject real GEO-calibrated kinetic/Hill parameters into WP1."""
        payload = json.loads(Path(json_path).read_text(encoding="utf-8"))
        calibrated = payload.get("wp1_parameters", payload)
        allowed = {"beta_0", "beta_E", "K_E", "n_hill", "delta_H", "beta_M"}
        updated = {}
        for k, v in calibrated.items():
            if k in allowed and np.isfinite(v):
                self.params[k] = float(v)
                updated[k] = float(v)
        return updated

    def get_params(self, genotype_name: str = "PiZZ_Padrao", therapy_name: str = "Controle") -> dict:
        p = self.params.copy()
        for k, mult in self.genotypes[genotype_name].items():
            p[k] = p[k] * mult
        for k, mult in self.therapies[therapy_name].items():
            p[k] = p[k] * mult
        return p

    @staticmethod
    def exposure_signal(t: float, E_amp: float, regime: str = "cronico", duty: float = 0.45, period: float = 8.0) -> float:
        if E_amp <= 0.0:
            return 0.0
        if regime == "cronico":
            return float(E_amp)
        return float(E_amp) if ((t % period) / period) < duty else 0.0

    def odes(self, t: float, y: np.ndarray, E_amp: float, p: dict | None = None, regime: str = "cronico") -> list[float]:
        p = self.params if p is None else p
        M, A, H, D, B = np.maximum(y, 0.0)
        B = float(np.clip(B, 0.01, 1.0))
        E = self.exposure_signal(t, E_amp, regime=regime)

        hill_E = (E**p["n_hill"]) / (p["K_E"]**p["n_hill"] + E**p["n_hill"]) if E > 0 else 0.0
        hill_M = (M**p["m_hill"]) / (p["K_PM"]**p["m_hill"] + M**p["m_hill"])
        age_factor = np.exp(-p["alpha_age"] * t)
        atp_eff = (1.25 * B) / (0.25 + B)
        ros_inhib = 1.0 / (1.0 + 0.07 * D)

        dM = (
            p["k_syn"]
            - p["k_deg"] * ros_inhib * M
            - (p["k_clr"] * atp_eff * H * M) / (p["K_M"] + M)
            - 2.0 * p["k_nuc"] * (1.0 + p["k_ox"] * D) * (M**2)
            - p["k_elong"] * M * A
        )
        dA = (
            p["k_nuc"] * (1.0 + p["k_ox"] * D) * (M**2)
            + p["k_elong"] * M * A
            - (p["k_auto"] * atp_eff * ros_inhib * H * A) / (p["K_A"] + A)
            - p["delta_A"] * A
        )
        dH = (
            (p["beta_0"] * age_factor + p["beta_E"] * hill_E + p["beta_M"] * hill_M) * atp_eff
            - p["delta_H"] * H
            - p["k_seq"] * H * A
        )
        dD = (
            p["gamma_E"] * (E**p["p_tox"])
            - (p["gamma_H"] * atp_eff * H * D) / (p["K_D"] + D)
            - p["delta_D"] * D
        )
        dB = (
            p["r_B"] * (1.0 - B)
            - p["c_H"] * H * B
            - p["c_D"] * D * B
            - p["c_A"] * A * B
        )
        return [dM, dA, dH, dD, dB]

    def simulate_ode(
        self,
        E_amp: float,
        p: dict | None = None,
        regime: str = "cronico",
        t_max: float | None = None,
        y0: np.ndarray | None = None,
    ):
        p = self.params if p is None else p
        t_end = float(p["t_max"] if t_max is None else t_max)
        y_init = self.y0 if y0 is None else np.asarray(y0, dtype=float)

        def ev_gen(t, y):
            return y[1] - p["A_crit"]
        ev_gen.direction = 1

        def ev_amb(t, y):
            return y[3] - p["D_crit"]
        ev_amb.direction = 1

        def ev_bio(t, y):
            return y[4] - p["B_crit"]
        ev_bio.direction = -1

        method = "Radau" if regime == "cronico" else "RK45"
        max_step = np.inf if regime == "cronico" else 0.45
        sol = solve_ivp(
            lambda t, y: self.odes(t, y, E_amp, p, regime=regime),
            (0.0, t_end),
            y_init,
            method=method,
            events=[ev_gen, ev_amb, ev_bio],
            max_step=max_step,
            rtol=1e-5,
            atol=1e-8,
        )
        t_gen = float(sol.t_events[0][0]) if len(sol.t_events[0]) > 0 else np.inf
        t_amb = float(sol.t_events[1][0]) if len(sol.t_events[1]) > 0 else np.inf
        t_bio = float(sol.t_events[2][0]) if len(sol.t_events[2]) > 0 else np.inf
        return sol, t_gen, min(t_amb, t_bio)

    def analyze_hysteresis(self, E_grid: np.ndarray | None = None):
        E_grid = np.linspace(0.0, 1.7, 35) if E_grid is None else np.asarray(E_grid, dtype=float)
        p_auto = self.params.copy()
        p_auto["alpha_age"] = 0.0
        y_low_init = np.array([0.8, 0.05, 1.8, 0.0, 1.0])
        y_high_init = np.array([2.5, 30.0, 0.15, 0.0, 0.5])
        A_low, A_high = [], []
        for E in E_grid:
            s_l = solve_ivp(lambda t, y: self.odes(t, y, E, p_auto), (0, 260), y_low_init, method="Radau")
            s_h = solve_ivp(lambda t, y: self.odes(t, y, E, p_auto), (0, 260), y_high_init, method="Radau")
            A_low.append(max(float(s_l.y[1, -1]), 0.0))
            A_high.append(max(float(s_h.y[1, -1]), 0.0))
        return E_grid, np.array(A_low), np.array(A_high)

    def simulate_sde_ensemble(
        self,
        E_amp: float,
        n_cells: int = 120,
        omega: float = 45.0,
        dt: float = 0.1,
        t_max: float | None = None,
        p: dict | None = None,
        seed: int | None = None,
    ):
        p = self.params if p is None else p
        t_end = float(p["t_max"] if t_max is None else t_max)
        rng = np.random.default_rng(self.seed if seed is None else seed)
        n_steps = int(t_end / dt)
        sqrt_dt = np.sqrt(dt)
        inv_sqrt_omega = 1.0 / np.sqrt(omega)

        M = np.full(n_cells, self.y0[0], dtype=float)
        A = np.full(n_cells, self.y0[1], dtype=float)
        H = np.full(n_cells, self.y0[2], dtype=float)
        D = np.full(n_cells, self.y0[3], dtype=float)
        B = np.full(n_cells, self.y0[4], dtype=float)

        t_gen = np.full(n_cells, np.inf)
        t_amb = np.full(n_cells, np.inf)
        hill_E = (E_amp**p["n_hill"]) / (p["K_E"]**p["n_hill"] + E_amp**p["n_hill"]) if E_amp > 0 else 0.0

        for step in range(1, n_steps + 1):
            t = step * dt
            age_factor = np.exp(-p["alpha_age"] * t)
            atp_eff = (1.25 * B) / (0.25 + B)
            ros_inhib = 1.0 / (1.0 + 0.07 * D)
            hill_M = (M**p["m_hill"]) / (p["K_PM"]**p["m_hill"] + M**p["m_hill"])

            v_syn = p["k_syn"]
            v_degM = p["k_deg"] * ros_inhib * M + (p["k_clr"] * atp_eff * H * M) / (p["K_M"] + M)
            v_nuc = p["k_nuc"] * (1.0 + p["k_ox"] * D) * (M**2)
            v_elo = p["k_elong"] * M * A
            v_degA = (p["k_auto"] * atp_eff * ros_inhib * H * A) / (p["K_A"] + A) + p["delta_A"] * A
            v_synH = (p["beta_0"] * age_factor + p["beta_E"] * hill_E + p["beta_M"] * hill_M) * atp_eff
            v_degH = p["delta_H"] * H + p["k_seq"] * H * A
            v_genD = p["gamma_E"] * (E_amp**p["p_tox"])
            v_degD = (p["gamma_H"] * atp_eff * H * D) / (p["K_D"] + D) + p["delta_D"] * D
            v_recB = p["r_B"] * (1.0 - B)
            v_useB = (p["c_H"] * H + p["c_D"] * D + p["c_A"] * A) * B

            dW = rng.normal(0.0, sqrt_dt, size=(5, n_cells)) * inv_sqrt_omega
            M = np.maximum(0.0, M + (v_syn - v_degM - 2.0 * v_nuc - v_elo) * dt + np.sqrt(np.maximum(v_syn + v_degM + 4.0 * v_nuc + v_elo, 1e-8)) * dW[0])
            A = np.maximum(0.0, A + (v_nuc + v_elo - v_degA) * dt + np.sqrt(np.maximum(v_nuc + v_elo + v_degA, 1e-8)) * dW[1])
            H = np.maximum(0.0, H + (v_synH - v_degH) * dt + np.sqrt(np.maximum(v_synH + v_degH, 1e-8)) * dW[2])
            D = np.maximum(0.0, D + (v_genD - v_degD) * dt + np.sqrt(np.maximum(v_genD + v_degD, 1e-8)) * dW[3])
            B = np.clip(B + (v_recB - v_useB) * dt + 0.15 * np.sqrt(np.maximum(v_recB + v_useB, 1e-8)) * dW[4], 0.01, 1.0)

            hit_g = np.isinf(t_gen) & (A >= p["A_crit"])
            hit_a = np.isinf(t_amb) & ((D >= p["D_crit"]) | (B <= p["B_crit"]))
            t_gen[hit_g] = t
            t_amb[hit_a] = t

        return t_gen, t_amb

    def run_prcc_sensitivity(self, n_samples: int = 80) -> pd.Series:
        target_keys = ["k_auto", "k_clr", "k_seq", "k_ox", "beta_E", "gamma_H", "c_D", "alpha_age"]
        sampler = qmc.LatinHypercube(d=len(target_keys), seed=self.seed)
        u = sampler.random(n=n_samples)
        X = np.zeros_like(u)
        Y = np.zeros(n_samples)

        for j, k in enumerate(target_keys):
            X[:, j] = self.params[k] * (0.6 + 0.8 * u[:, j])

        for i in range(n_samples):
            p_i = self.params.copy()
            for j, k in enumerate(target_keys):
                p_i[k] = float(X[i, j])
            _, tg1, ta1 = self.simulate_ode(0.20, p=p_i)
            _, tg2, ta2 = self.simulate_ode(1.55, p=p_i)
            Y[i] = min(tg1, ta1, self.params["t_max"]) + min(tg2, ta2, self.params["t_max"])

        R = np.column_stack([rankdata(X[:, j]) for j in range(len(target_keys))] + [rankdata(Y)])
        C_mat = np.corrcoef(R, rowvar=False)
        P_mat = np.linalg.pinv(C_mat)
        prcc = {
            k: float(-P_mat[j, -1] / np.sqrt(abs(P_mat[j, j] * P_mat[-1, -1]) + 1e-12))
            for j, k in enumerate(target_keys)
        }
        return pd.Series(prcc).sort_values(key=abs, ascending=False)

    def run_combination_screening(self, E_grid: np.ndarray | None = None):
        E_grid = np.linspace(0.0, 2.0, 21) if E_grid is None else np.asarray(E_grid, dtype=float)
        t_max = float(self.params["t_max"])
        records = []
        for g_name in self.genotypes:
            for th_name in self.therapies:
                p_gt = self.get_params(g_name, th_name)
                surv = np.array([min(*self.simulate_ode(e, p=p_gt, t_max=t_max)[1:], t_max) for e in E_grid])
                auc_norm = float(np.trapezoid(surv, E_grid) / (t_max * (E_grid[-1] - E_grid[0]))) if hasattr(np, "trapezoid") else float(np.trapz(surv, E_grid) / (t_max * (E_grid[-1] - E_grid[0])))
                safe_mask = surv >= (0.98 * t_max)
                width = float(E_grid[safe_mask][-1] - E_grid[safe_mask][0]) if np.sum(safe_mask) >= 2 else 0.0
                records.append({
                    "Genotipo": g_name,
                    "Terapia": th_name,
                    "AUC_Norm": auc_norm,
                    "Janela_Segura_dE": width,
                    "T_obs_E0": float(surv[0]),
                    "T_obs_E16": float(surv[np.argmin(np.abs(E_grid - 1.60))]),
                })

        df_screen = pd.DataFrame(records)
        ctrl = df_screen[df_screen["Terapia"] == "Controle"].set_index("Genotipo")["AUC_Norm"].to_dict()
        df_screen["Ganho_AUC"] = df_screen.apply(lambda r: r["AUC_Norm"] - ctrl[r["Genotipo"]], axis=1)

        combos = [
            ("Combo_Auto_Mito", "Indutor_Autofagia", "Antioxidante_Mito"),
            ("Combo_HSF1_Mito", "Ativador_HSF1", "Antioxidante_Mito"),
            ("Combo_Auto_AntiAgg", "Indutor_Autofagia", "Anti_Agregante"),
        ]
        syn_rows = []
        for g_name in self.genotypes:
            sub = df_screen[df_screen["Genotipo"] == g_name].set_index("Terapia")
            for c_name, a_name, b_name in combos:
                ga = max(float(sub.loc[a_name, "Ganho_AUC"]), 0.0)
                gb = max(float(sub.loc[b_name, "Ganho_AUC"]), 0.0)
                gc = float(sub.loc[c_name, "Ganho_AUC"])
                bliss_exp = ga + gb - ga * gb
                syn_rows.append({
                    "Genotipo": g_name,
                    "Combinacao": c_name,
                    "Ganho_A": ga,
                    "Ganho_B": gb,
                    "Bliss_Esperado": bliss_exp,
                    "Ganho_Combo": gc,
                    "Excesso_Bliss": gc - bliss_exp,
                })
        return df_screen, pd.DataFrame(syn_rows)
